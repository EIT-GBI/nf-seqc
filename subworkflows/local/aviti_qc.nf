//
// Subworkflow: AVITI_QC
//
// The AVITI run folder holds no FASTQ, so conversion is part of this pipeline's
// job rather than a precondition: bases2fastq produces the FASTQ scientists get
// plus the run-level metrics (RunStats.json / Metrics.csv) that MultiQC reads.
//

include { BASES2FASTQ   } from '../../modules/local/bases2fastq/main.nf'
include { FASTQC_FASTQC } from '../../modules/fastqc/fastqc/main.nf'
include { demuxSummary; withPhixReads; phixSample; barcodeName; readShareTable; demuxSummaryTable; unknownBarcodesTable } from './demux_report.nf'

// Element's own run QC report (<run>_QC.html): published next to the MultiQC
// report rather than with the QC files, because scientists open it directly.
def isElementReport( f ) {
    return f.name.endsWith('_QC.html')
}

// PhiX is the run's spike-in control, not a sample: no FASTQ delivered, no
// FastQC, no sample row. Its share of all reads is a column of the
// Demultiplexing table; the rest stays in bases2fastq's own reports.
// Same rule as phixSample() on the sample name: PhiX_R1.fastq.gz,
// PhiX_v3_R1.fastq.gz.
def isPhiX( f ) {
    return f.name.startsWith('PhiX_')
}

// Reads bases2fastq could not assign to a sample: always delivered, QC'd only
// with --include_undetermined.
def isUnassigned( f ) {
    return f.name.startsWith('Unassigned')
}

// A process output declared as path("*.x") is a single Path when one file
// matches and a List when several do. `[] + path` does NOT wrap a single Path:
// Path is Iterable over its directory components, so it yields every parent
// directory (home, runner, work, ...) as if it were a file. Wrap explicitly.
def asFileList( files ) {
    return files instanceof List ? files : [ files ]
}

// Polonies per sample from bases2fastq's RunStats.json, and each lane's
// UnassignedSequences, as demuxSummary() takes them; Unassigned is the run
// total less every assigned sample, PhiX included. PhiX is the spike-in, not
// a sample: it gets no % Reads or Loading row and does not count in the even
// split, and its share of all reads goes to the Demultiplexing table instead.
// null when there is nothing to report on (e.g. the empty file a stub run
// writes).
def bases2fastqReads( run_stats ) {
    if ( !run_stats.size() ) {
        return null
    }
    def stats    = new groovy.json.JsonSlurper().parse( run_stats.toFile() )
    def assigned = stats.SampleStats.collectEntries { s -> [ ( s.SampleName ): s.NumPolonies as long ] }
    def total    = stats.NumPolonies as long
    def phix     = assigned.findAll { name, _n -> phixSample( name ) }.values().sum() ?: 0L
    def samples  = assigned.findAll { name, _n -> !phixSample( name ) }
    def barcodes = ( stats.Lanes ?: [] ).collectMany { lane ->
        ( lane.UnassignedSequences ?: [] ).collect { u ->
            [ lane: lane.Lane, index: barcodeName( u.I1, u.I2 ), reads: u.Count as long ]
        }
    }
    def demux = demuxSummary( samples, total - ( assigned.values().sum() ?: 0L ), barcodes )
    return withPhixReads( demux, phix as long, 'assigned to the PhiX control by bases2fastq' )
}

// Polonies per sample, lane and index pair, from RunStats.json's per-sample
// Occurrences, as MultiQC custom content. It replaces the bases2fastq module's
// own Index Assignment Metrics table (removed in assets/multiqc_config.yml),
// which MultiQC (1.35) keys on the index pair alone: when the flow cell is
// split into two independent lanes that reuse the same indexes, the lane 2
// sample got no row and its polonies were added to the lane 1 sample's
// (SEQEBA000015: 10 of 53 samples had no row; lab feedback 2026-11-23).
// Occurrences give the lane and the concatenated I1+I2; RunManifest.json gives
// the split, matched on sample, lane and sequence. Without an exact match
// (e.g. an index stored reverse-complemented, or no manifest) the sequence is
// split at the sample's I1 length from the manifest, else the run's, and a
// warning names the rows; with no manifest at all it all goes in Index 1. An
// occurrence without a Lane counts as lane 1. Rows run lane by lane in
// manifest order, so a sample sequenced on both lanes (1+2) has a row in
// each. PhiX keeps its rows, one per index and lane: it is in Element's
// table too. null when RunStats.json has no Occurrences (e.g. a stub run's
// empty file).
def indexAssignmentTable( run_stats, run_manifest ) {
    if ( !run_stats.size() ) {
        return null
    }
    def stats    = new groovy.json.JsonSlurper().parse( run_stats.toFile() )
    def manifest = run_manifest?.size() ? new groovy.json.JsonSlurper().parse( run_manifest.toFile() ) : [:]
    def indexes  = [:]
    def i1_len   = [:]
    def order    = [:]
    ( manifest.Samples ?: [] ).eachWithIndex { s, i ->
        order[ s.SampleName ] = s.SampleNumber != null ? s.SampleNumber as int : i
        ( s.Indexes ?: [] ).each { ix ->
            def i1 = ix.Index1 ?: ''
            def i2 = ix.Index2 ?: ''
            indexes[ "${s.SampleName}\t${ix.Lane}\t${i1}${i2}".toString() ] = [ i1, i2 ]
            i1_len.putIfAbsent( s.SampleName, i1.length() )
        }
    }
    def run_i1_len = i1_len ? i1_len.values().countBy { n -> n }.max { e -> e.value }.key : null
    def lane_total = ( stats.Lanes ?: [] ).collectEntries { l -> [ ( ( l.Lane ?: 1 ) as int ): ( l.NumPoloniesBeforeTrimming ?: l.NumPolonies ?: 0 ) as long ] }
    def run_total  = ( stats.NumPoloniesBeforeTrimming ?: stats.NumPolonies ?: 0 ) as long
    def unmatched  = []
    def rows = ( stats.SampleStats ?: [] ).collectMany { s ->
        ( s.Occurrences ?: [] ).collect { o ->
            def seq  = o.ExpectedSequence ?: ''
            def lane = ( o.Lane ?: 1 ) as int
            def ix   = indexes[ "${s.SampleName}\t${lane}\t${seq}".toString() ]
            if ( ix == null ) {
                def n = i1_len.containsKey( s.SampleName ) ? i1_len[ s.SampleName ] : run_i1_len
                ix = n != null && n < seq.length() ? [ seq.substring( 0, n ), seq.substring( n ) ] : [ seq, '' ]
                if ( manifest ) {
                    unmatched << "${s.SampleName} L${lane} ${seq}"
                }
            }
            [ sample: s.SampleName, lane: lane, i1: ix[0], i2: ix[1],
              polonies: ( o.NumPoloniesBeforeTrimming ?: o.NumPolonies ?: 0 ) as long,
              order: order.containsKey( s.SampleName ) ? order[ s.SampleName ] : ( s.SampleNumber ?: 0 ) as int ]
        }
    }
    if ( unmatched ) {
        log.warn "AVITI_QC: ${unmatched.size()} index sequence(s) in RunStats.json match no index in RunManifest.json; Index 1/2 split at the manifest's I1 length: ${unmatched.take( 5 ).join( ', ' )}${unmatched.size() > 5 ? ', ...' : ''}"
    }
    if ( !rows ) {
        return null
    }
    rows = rows.sort { a, b -> a.lane <=> b.lane ?: a.order <=> b.order ?: a.i1 <=> b.i1 ?: a.i2 <=> b.i2 }
    // A sample's row name is unique per lane, except PhiX's (several indexes).
    def per_lane = rows.countBy { r -> "${r.sample}\t${r.lane}".toString() }
    def header = [
        "# id: 'bases2fastq_index_assignment'",
        "# section_name: 'Index Assignment Metrics'",
        "# description: 'Polonies assigned to each sample in each lane, by index pair, from bases2fastq RunStats.json. One row per sample per lane: a sample sequenced on both lanes (1+2) appears twice, and samples sharing indexes in different lanes each get their own row.'",
        "# plot_type: 'table'",
        "# pconfig:",
        "#     id: 'bases2fastq_index_assignment_table'",
        "#     namespace: 'bases2fastq'",
        "#     col1_header: 'Sample (lane)'",
        "#     sort_rows: false",
        "# headers:",
        "#     lane:",
        "#         title: 'Lane'",
        "#         scale: false",
        "#         format: '{:,.0f}'",
        "#     index_1:",
        "#         title: 'Index 1'",
        "#         description: 'Sample index 1 (I1)'",
        "#     index_2:",
        "#         title: 'Index 2'",
        "#         description: 'Sample index 2 (I2)'",
        "#     polonies:",
        "#         title: 'Assigned Polonies'",
        "#         description: 'Polonies assigned to this sample and index pair in this lane, before trimming'",
        "#         scale: 'RdYlGn'",
        "#         format: '{:,.0f}'",
        "#     pct_lane:",
        "#         title: '% of Lane'",
        "#         description: 'Share of all polonies in this lane'",
        "#         suffix: '%'",
        "#         min: 0",
        "#         scale: 'RdYlGn'",
        "#         format: '{:,.2f}'",
        "#     pct_run:",
        "#         title: '% of Run'",
        "#         description: 'Share of all polonies in the run, every lane'",
        "#         suffix: '%'",
        "#         min: 0",
        "#         scale: 'RdYlGn'",
        "#         format: '{:,.2f}'",
        "Sample\tlane\tindex_1\tindex_2\tpolonies\tpct_lane\tpct_run",
    ]
    def lines = rows.collect { r ->
        def name     = "${r.sample} (L${r.lane})" + ( per_lane[ "${r.sample}\t${r.lane}".toString() ] > 1 ? " ${barcodeName( r.i1, r.i2 )}" : '' )
        def pct_lane = lane_total[ r.lane ] ? 100.0 * r.polonies / lane_total[ r.lane ] : ''
        def pct_run  = run_total ? 100.0 * r.polonies / run_total : ''
        "${name}\t${r.lane}\t${r.i1}\t${r.i2}\t${r.polonies}\t${pct_lane}\t${pct_run}"
    }
    return ( header + lines ).join('\n') + '\n'
}

workflow AVITI_QC {

    take:
    ch_run               // channel: [ val(meta), path(run_dir) ]
    include_undetermined // boolean: --include_undetermined
    max_undetermined     // number: --max_undetermined, % of reads Unassigned

    main:
    ch_versions = channel.empty()

    BASES2FASTQ( ch_run )

    // bases2fastq writes Samples/<project>/<sample>/<sample>_R{1,2}.fastq.gz.
    ch_pairs = BASES2FASTQ.out.fastq.flatMap { meta, fastqs ->
        def files = asFileList( fastqs )

        files = files.findAll { f -> !isPhiX( f ) }

        def r1s = files.findAll { f -> f.name ==~ /.*_R1(_001)?\.fastq\.gz$/ }
        if ( !r1s.findAll { f -> !isUnassigned( f ) } ) {
            error "AVITI_QC: bases2fastq produced no R1 FASTQ. Check the run manifest in ${meta.id}."
        }

        r1s.sort().collect { r1 ->
            def r2  = files.find { f -> f.name == r1.name.replace('_R1', '_R2') }
            def sid = r1.name.replaceAll(/_R1(_001)?\.fastq\.gz$/, '')
            tuple( meta + [ id: sid, single_end: r2 == null, undetermined: isUnassigned( r1 ) ], r1, r2 ?: [] )
        }
    }

    FASTQC_FASTQC( ch_pairs.filter { meta, _r1, _r2 -> include_undetermined || !meta.undetermined } )

    // flatMap: a missing or empty RunStats.json gives no tables rather than an error
    ch_demux = BASES2FASTQ.out.stats.flatMap { meta, stats ->
        def run_stats = asFileList( stats ).find { f -> f.name == 'RunStats.json' }
        def demux     = run_stats ? bases2fastqReads( run_stats ) : null
        demux ? [ tuple( meta, demux ) ] : []
    }
    ch_read_share = ch_demux
        .map { _meta, demux -> readShareTable( demux, 'Unassigned' ) }
        .collectFile( name: 'demux_read_share_mqc.tsv' )
    ch_demux_summary = ch_demux
        .map { meta, demux -> demuxSummaryTable( meta.run, demux, max_undetermined, 'bases2fastq RunStats.json', 'Unassigned', 'AVITI_QC' ) }
        .collectFile( name: 'demux_summary_mqc.tsv' )
    ch_unknown_barcodes = ch_demux
        .flatMap { _meta, demux ->
            def table = unknownBarcodesTable( demux, 'bases2fastq RunStats.json', 'Unassigned', 10, 'AVITI_QC' )
            table ? [ table ] : []
        }
        .collectFile( name: 'demux_unknown_barcodes_mqc.tsv' )
    ch_index_assignment = BASES2FASTQ.out.stats
        .flatMap { _meta, stats ->
            def files     = asFileList( stats )
            def run_stats = files.find { f -> f.name == 'RunStats.json' }
            def table     = run_stats ? indexAssignmentTable( run_stats, files.find { f -> f.name == 'RunManifest.json' } ) : null
            table ? [ table ] : []
        }
        .collectFile( name: 'bases2fastq_index_assignment_mqc.tsv' )

    ch_qc = channel.empty()
        .mix( FASTQC_FASTQC.out.zip.map { _meta, zip -> zip } )
        .mix( BASES2FASTQ.out.stats.flatMap { _meta, stats -> asFileList( stats ).findAll { f -> !isElementReport( f ) } } )
        .mix( BASES2FASTQ.out.sample_stats.flatMap { _meta, stats ->
            asFileList( stats ).findAll { f -> !isPhiX( f ) && ( include_undetermined || !isUnassigned( f ) ) }
        } )
        .mix( ch_read_share )
        .mix( ch_demux_summary )
        .mix( ch_unknown_barcodes )
        .mix( ch_index_assignment )

    ch_versions = ch_versions
        .mix( FASTQC_FASTQC.out.versions_fastqc )
        .mix( BASES2FASTQ.out.versions_bases2fastq )

    emit:
    // channel: path — published in qc/ for scientists, not read by MultiQC
    files    = FASTQC_FASTQC.out.html.map { _meta, html -> html }
                   .mix( BASES2FASTQ.out.run_files.map { _meta, f -> f } )
    report   = BASES2FASTQ.out.stats.map { _meta, stats -> asFileList( stats ).find { f -> isElementReport( f ) } }
                   .filter { f -> f } // channel: path — Element's <run>_QC.html
    reads    = ch_pairs.map { meta, r1, r2 -> tuple( meta, [ r1, r2 ].flatten().findAll { f -> f } ) }
    qc       = ch_qc
    versions = ch_versions
}
