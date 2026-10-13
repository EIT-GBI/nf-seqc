//
// Subworkflow: ILLUMINA_QC
//
// NextSeq demultiplexes on-instrument, so FASTQ already exist inside the run
// folder. This branch therefore only reads: per-sample FASTQ quality (FastQC)
// and run-level metrics from the binary InterOp files (%PF, %occupied, Q30,
// phasing, PhiX error rate) — none of which reach the FASTQ at all. The
// on-instrument BCL Convert reports (index assignment, undetermined barcodes,
// per-sample Q30) go to MultiQC's bclconvert module, the Illumina counterpart
// of the bases2fastq section on AVITI.
//
// The run-level half mirrors EIT-GBI/Automate-Seq-Run-Metrics-Collection: same
// metrics, same column names, same %Occupied vs %PF plot. Two behaviours from
// that tool are handled differently here:
//
//   * Its `RTAComplete.txt` + `CopyComplete.txt` gate is kept, but lives in
//     `requireRunComplete()` in main.nf: raising it there, where run_dir is a
//     plain value, gives one clean error instead of an abort inside a channel
//     operator followed by a spurious "empty List" message.
//   * Its `monitoring_complete.txt` marker exists to stop cron reprocessing a
//     run. Nextflow's -resume covers that, so no marker is written; nothing in
//     the run directory is modified.
//

include { FASTQC_FASTQC        } from '../../modules/fastqc/fastqc/main.nf'
include { INTEROP_RUN_METRICS  } from '../../modules/interop/run_metrics/main.nf'

// BCL Convert's report directory, wherever the FASTQ came from: Analysis/<n>/Data/Reports/
// on-instrument (NextSeq 1000/2000), Reports/ in a bcl-convert output dir. On-instrument
// runs also leave a Demultiplex_Stats.csv under RunInstrumentAnalyticsMetrics/ with the
// sample names and indices anonymised (0001, ELIDED), so a directory called Reports wins.
// Several analyses can leave several; the first one sorted wins.
def bclConvertReports( root ) {
    def found = ( files("${root}/**/Demultiplex_Stats.csv") + files("${root}/Demultiplex_Stats.csv") )
        .findAll { f -> f.exists() }
        .unique()
        .sort()
    def reports = found.findAll { f -> f.parent.name == 'Reports' } ?: found
    if ( !reports ) {
        log.warn "ILLUMINA_QC: no Demultiplex_Stats.csv under ${root}; the report will have no demultiplexing sections."
        return null
    }
    if ( reports.size() > 1 ) {
        log.warn "ILLUMINA_QC: ${reports.size()} BCL Convert reports under ${root}; using ${reports[0].parent}."
    }
    return reports[0].parent
}

// Inputs for MultiQC's bclconvert module. It pairs Demultiplex_Stats.csv with a RunInfo.xml
// in the same directory (all of qc/ once staged), and skips both without either one, so the
// run's RunInfo.xml goes along only when there is a demux file for it. MultiQC 1.35 does not
// read Adapter_Metrics.csv.
def bclConvertFiles( reports, run_dir ) {
    def demux = reports.resolve('Demultiplex_Stats.csv')
    def extra = [ 'Quality_Metrics.csv', 'Top_Unknown_Barcodes.csv' ]
        .collect { name -> reports.resolve(name) }
        .findAll { f -> f.exists() }
    return [ demux, run_dir.resolve('RunInfo.xml') ] + extra
}

// Cell backgrounds for flag columns: the occupancy plot's zone colours
// (steelblue / green / firebrick) semi-transparent, so the text stays readable on light and dark themes.
// Header `bgcols` is how MultiQC's own modules flag categories; the same map
// is in nf-mod-interop's Run Metrics table.
def bgcols( Map colours ) {
    return [ "#         bgcols:" ] + colours.collect { value, colour -> "#             ${value}: '${colour}'" }
}
def loadingBgcols() { return [ under: 'rgba(70,130,180,0.4)', optimal: 'rgba(0,128,0,0.35)', over: 'rgba(178,34,34,0.4)' ] }

// Reads per SampleID from Demultiplex_Stats.csv, summed over lanes. Returns
// [ samples: [id: reads], undetermined: reads, total: reads, even: % ], or null
// when there is nothing to report on. `even` is one sample's share of all
// reads if the demultiplexed reads were split evenly.
def demuxReads( stats ) {
    def reads = [:]
    stats.splitCsv( header: true ).each { row ->
        reads[row.SampleID] = ( reads[row.SampleID] ?: 0L ) + ( row['# Reads'] as long )
    }
    def total    = reads.values().sum() ?: 0L
    def samples  = reads.findAll { id, _n -> id != 'Undetermined' }
    def assigned = samples.values().sum() ?: 0L
    if ( !total || !samples ) {
        return null
    }
    return [
        samples:      samples,
        undetermined: reads['Undetermined'] ?: 0L,
        total:        total,
        even:         100.0 * assigned / total / samples.size(),
    ]
}

// A sample's loading in the pool: under below half an even split, over above
// twice one, optimal in between.
def sampleLoading( demux, long reads ) {
    def pct = 100.0 * reads / demux.total
    return pct < demux.even / 2 ? 'under' : pct > demux.even * 2 ? 'over' : 'optimal'
}

// MultiQC General Statistics columns: each sample's share of the run's reads
// (Undetermined counted in the total but given no row), and its Loading. The
// cut-offs depend on how many samples share the run, so they are computed
// here and given in the column description.
def readShareTable( demux ) {
    def under = String.format( '%.2f', demux.even / 2 )
    def over  = String.format( '%.2f', demux.even * 2 )
    def header = [
        "# id: 'illumina_read_share'",
        "# plot_type: 'generalstats'",
        "# headers:",
        "#     pct_reads:",
        "#         title: '% Reads'",
        "#         namespace: 'Demultiplexing'",
        "#         description: 'Share of all reads in the run, Undetermined included'",
        "#         suffix: '%'",
        "#         min: 0",
        "#         format: '{:,.2f}'",
        "#     loading:",
        "#         title: 'Loading'",
        "#         namespace: 'Demultiplexing'",
        "#         description: 'Share of reads against an even split of the pool: under below ${under}%, over above ${over}%'",
    ] + bgcols( loadingBgcols() ) + [
        "Sample\tpct_reads\tloading",
    ]
    def rows = demux.samples.collect { id, n -> "${id}\t${100.0 * n / demux.total}\t${sampleLoading( demux, n )}" }
    return ( header + rows ).join('\n') + '\n'
}

// Run-level demultiplexing summary, one row, placed under Run Metrics. A high
// Undetermined share is invisible in General Statistics, which only compares
// samples with each other, so it gets its own flag above --max_undetermined,
// and a warning in the Nextflow log. No colour-scale bars: in a one-row
// table every bar is full width, so each number turns into a solid block.
def demuxSummaryTable( run_id, demux, max_undetermined ) {
    def pct_undetermined = 100.0 * demux.undetermined / demux.total
    def level = pct_undetermined > max_undetermined ? 'high' : 'ok'
    if ( level == 'high' ) {
        log.warn "ILLUMINA_QC: ${String.format('%.1f', pct_undetermined)}% of reads in ${run_id} are Undetermined (above --max_undetermined ${max_undetermined}%)."
    }
    def counts = [ under: 0, optimal: 0, over: 0 ]
    demux.samples.each { _id, n -> counts[ sampleLoading( demux, n ) ] += 1 }
    def header = [
        "# id: 'illumina_demux_summary'",
        "# section_name: 'Illumina Demultiplexing'",
        "# description: 'Run-level summary of BCL Convert Demultiplex_Stats.csv. Undetermined is high above ${max_undetermined}% (--max_undetermined). Sample loading is the General Statistics Loading column.'",
        "# plot_type: 'table'",
        "# pconfig:",
        "#     id: 'illumina_demux_summary_table'",
        "#     namespace: 'Demultiplexing'",
        "# headers:",
        "#     pct_undetermined:",
        "#         title: 'Undetermined'",
        "#         scale: false",
        "#         description: 'Share of all reads not assigned to any sample'",
        "#         suffix: '%'",
        "#         min: 0",
        "#         max: 100",
        "#         format: '{:,.2f}'",
        "#     undetermined_level:",
        "#         title: 'Undetermined level'",
        "#         scale: false",
        "#         description: 'high above ${max_undetermined}% (--max_undetermined)'",
    ] + bgcols( [ ok: loadingBgcols().optimal, high: loadingBgcols().over ] ) + [
        "#     samples:",
        "#         title: 'Samples'",
        "#         scale: false",
        "#         format: '{:,.0f}'",
        "#     under:",
        "#         title: 'Under'",
        "#         scale: false",
        "#         description: 'Samples under half an even split'",
        "#         format: '{:,.0f}'",
        "#     optimal:",
        "#         title: 'Optimal'",
        "#         scale: false",
        "#         format: '{:,.0f}'",
        "#     over:",
        "#         title: 'Over'",
        "#         scale: false",
        "#         description: 'Samples over twice an even split'",
        "#         format: '{:,.0f}'",
        "Run\tpct_undetermined\tundetermined_level\tsamples\tunder\toptimal\tover",
        "${run_id}\t${pct_undetermined}\t${level}\t${demux.samples.size()}\t${counts.under}\t${counts.optimal}\t${counts.over}",
    ]
    return header.join('\n') + '\n'
}

workflow ILLUMINA_QC {

    take:
    ch_run               // channel: [ val(meta), path(run_dir) ]
    fastq_dir            // file or null: --fastq_dir, where the FASTQ live if not in the run dir
    include_undetermined // boolean: --include_undetermined
    max_undetermined     // number: --max_undetermined, % of reads

    main:
    ch_versions = channel.empty()

    // FASTQ location varies by layout: Analysis/<n>/Data/... on DRAGEN-enabled
    // NextSeq 1000/2000, Alignment_*/.../Fastq on older setups. Glob instead of
    // hardcoding; --fastq_dir overrides when the FASTQ live outside the run dir.
    ch_pairs = ch_run.flatMap { meta, run_dir ->
        def root = fastq_dir ?: run_dir
        def r1s  = ( files("${root}/**/*_R1_001.fastq.gz") + files("${root}/*_R1_001.fastq.gz") )
            .findAll { f -> f.exists() }
            .unique()

        if ( !include_undetermined ) {
            r1s = r1s.findAll { f -> !f.name.startsWith('Undetermined') }
        }
        if ( !r1s ) {
            error "ILLUMINA_QC: no *_R1_001.fastq.gz found under ${root}. Has the run been demultiplexed? Point --fastq_dir at the FASTQ directory if they live elsewhere."
        }

        r1s.sort().collect { r1 ->
            def r2  = file( r1.toString().replace('_R1_001.fastq.gz', '_R2_001.fastq.gz') )
            def sid = r1.name
                .replaceAll(/_S\d+(_L\d{3})?_R1_001\.fastq\.gz$/, '')
                .replaceAll(/_R1_001\.fastq\.gz$/, '')
            tuple( meta + [ id: sid, single_end: !r2.exists() ], r1, r2.exists() ? r2 : [] )
        }
    }

    FASTQC_FASTQC( ch_pairs )
    INTEROP_RUN_METRICS( ch_run )

    // flatMap, not map + filter: ch_run is a value channel, and filtering one
    // down to nothing never closes it, which leaves MULTIQC waiting forever on
    // any run without a Demultiplex_Stats.csv.
    ch_reports = ch_run
        .flatMap { meta, run_dir ->
            def reports = bclConvertReports( fastq_dir ?: run_dir )
            reports ? [ tuple( meta, run_dir, reports ) ] : []
        }

    ch_demux = ch_reports
        .map { meta, _run_dir, reports -> tuple( meta, demuxReads( reports.resolve('Demultiplex_Stats.csv') ) ) }
        .filter { _meta, demux -> demux }

    // bclconvert's own tables report each sample's share of the *assigned* reads;
    // this column adds its share of all reads, which is what drops when many
    // reads go Undetermined.
    ch_read_share = ch_demux
        .map { _meta, demux -> readShareTable( demux ) }
        .collectFile( name: 'illumina_read_share_mqc.tsv' )

    ch_demux_summary = ch_demux
        .map { meta, demux -> demuxSummaryTable( meta.run, demux, max_undetermined ) }
        .collectFile( name: 'illumina_demux_summary_mqc.tsv' )

    ch_bclconvert = ch_reports
        .map { _meta, run_dir, reports -> bclConvertFiles( reports, run_dir ) }

    // The metrics table and occupancy plot go to MultiQC as custom content;
    // the CSV itself is the machine-readable log row and is published separately.
    ch_qc = channel.empty()
        .mix( FASTQC_FASTQC.out.zip.map { _meta, zip -> zip } )
        .mix( INTEROP_RUN_METRICS.out.mqc_table.map { _meta, tsv -> tsv } )
        .mix( INTEROP_RUN_METRICS.out.plot.map { _meta, png -> png } )
        .mix( ch_read_share )
        .mix( ch_demux_summary )
        .mix( ch_bclconvert )

    ch_versions = ch_versions
        .mix( FASTQC_FASTQC.out.versions_fastqc )
        .mix( INTEROP_RUN_METRICS.out.versions_interop )

    emit:
    reads    = ch_pairs.map { meta, r1, r2 -> tuple( meta, [ r1, r2 ].flatten().findAll { f -> f } ) }
    qc       = ch_qc                          // channel: path — flat files for MultiQC
    metrics  = INTEROP_RUN_METRICS.out.metrics // channel: [ meta, run_metrics.csv ]
    versions = ch_versions                     // channel: version records collected from modules
}
