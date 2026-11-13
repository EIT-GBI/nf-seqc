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
include { demuxSummary; barcodeName; readShareTable; demuxSummaryTable; unknownBarcodesTable } from './demux_report.nf'

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

// Reads per SampleID from Demultiplex_Stats.csv, summed over lanes, and the
// unknown barcodes from the Top_Unknown_Barcodes.csv beside it (if BCL Convert
// wrote one), as demuxSummary() takes them. null when there is nothing to
// report on.
def demuxReads( stats ) {
    def reads = [:]
    stats.splitCsv( header: true ).each { row ->
        reads[row.SampleID] = ( reads[row.SampleID] ?: 0L ) + ( row['# Reads'] as long )
    }
    def unknown  = stats.resolveSibling('Top_Unknown_Barcodes.csv')
    def barcodes = unknown.exists()
        ? unknown.splitCsv( header: true ).collect { row ->
            [ lane: row.Lane, index: barcodeName( row.index, row.index2 ), reads: row['# Reads'] as long ]
        }
        : []
    return demuxSummary( reads.findAll { id, _n -> id != 'Undetermined' }, reads['Undetermined'] ?: 0L, barcodes )
}

// Reads BCL Convert could not assign to a sample: always delivered, QC'd only
// with --include_undetermined.
def isUndetermined( f ) {
    return f.name.startsWith('Undetermined')
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

        if ( !r1s.findAll { f -> !isUndetermined( f ) } ) {
            error "ILLUMINA_QC: no *_R1_001.fastq.gz found under ${root}. Has the run been demultiplexed? Point --fastq_dir at the FASTQ directory if they live elsewhere."
        }

        r1s.sort().collect { r1 ->
            def r2  = file( r1.toString().replace('_R1_001.fastq.gz', '_R2_001.fastq.gz') )
            def sid = r1.name
                .replaceAll(/_S\d+(_L\d{3})?_R1_001\.fastq\.gz$/, '')
                .replaceAll(/_R1_001\.fastq\.gz$/, '')
            tuple( meta + [ id: sid, single_end: !r2.exists(), undetermined: isUndetermined( r1 ) ], r1, r2.exists() ? r2 : [] )
        }
    }

    FASTQC_FASTQC( ch_pairs.filter { meta, _r1, _r2 -> include_undetermined || !meta.undetermined } )
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
        .map { _meta, demux -> readShareTable( demux, 'Undetermined' ) }
        .collectFile( name: 'demux_read_share_mqc.tsv' )

    ch_demux_summary = ch_demux
        .map { meta, demux -> demuxSummaryTable( meta.run, demux, max_undetermined, 'BCL Convert Demultiplex_Stats.csv', 'Undetermined', 'ILLUMINA_QC' ) }
        .collectFile( name: 'demux_summary_mqc.tsv' )

    ch_unknown_barcodes = ch_demux
        .flatMap { _meta, demux ->
            def table = unknownBarcodesTable( demux, 'BCL Convert Top_Unknown_Barcodes.csv', 'Undetermined', 10, 'ILLUMINA_QC' )
            table ? [ table ] : []
        }
        .collectFile( name: 'demux_unknown_barcodes_mqc.tsv' )

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
        .mix( ch_unknown_barcodes )
        .mix( ch_bclconvert )

    ch_versions = ch_versions
        .mix( FASTQC_FASTQC.out.versions_fastqc )
        .mix( INTEROP_RUN_METRICS.out.versions_interop )

    emit:
    files    = FASTQC_FASTQC.out.html.map { _meta, html -> html } // channel: path — FastQC HTML, published in qc/
    reads    = ch_pairs.map { meta, r1, r2 -> tuple( meta, [ r1, r2 ].flatten().findAll { f -> f } ) }
    qc       = ch_qc                          // channel: path — flat files for MultiQC
    metrics  = INTEROP_RUN_METRICS.out.metrics // channel: [ meta, run_metrics.csv ]
    versions = ch_versions                     // channel: version records collected from modules
}
