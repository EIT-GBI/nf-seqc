//
// Subworkflow: AVITI_QC
//
// The AVITI run folder holds no FASTQ, so conversion is part of this pipeline's
// job rather than a precondition: bases2fastq produces the FASTQ scientists get
// plus the run-level metrics (RunStats.json / Metrics.csv) that MultiQC reads.
//

include { BASES2FASTQ   } from '../../modules/local/bases2fastq/main.nf'
include { FASTQC_FASTQC } from '../../modules/fastqc/fastqc/main.nf'
include { demuxSummary; barcodeName; readShareTable; demuxSummaryTable; unknownBarcodesTable } from './demux_report.nf'

// Element's own run QC report (<run>_QC.html): published next to the MultiQC
// report rather than with the QC files, because scientists open it directly.
def isElementReport( f ) {
    return f.name.endsWith('_QC.html')
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
// total less the samples. null when there is nothing to report on (e.g. the
// empty file a stub run writes).
def bases2fastqReads( run_stats ) {
    if ( !run_stats.size() ) {
        return null
    }
    def stats    = new groovy.json.JsonSlurper().parse( run_stats.toFile() )
    def samples  = stats.SampleStats.collectEntries { s -> [ ( s.SampleName ): s.NumPolonies as long ] }
    def total    = stats.NumPolonies as long
    def barcodes = ( stats.Lanes ?: [] ).collectMany { lane ->
        ( lane.UnassignedSequences ?: [] ).collect { u ->
            [ lane: lane.Lane, index: barcodeName( u.I1, u.I2 ), reads: u.Count as long ]
        }
    }
    return demuxSummary( samples, total - ( samples.values().sum() ?: 0L ), barcodes )
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

        // Unassigned reads are always delivered; they are QC'd (FastQC, a sample
        // row) only with --include_undetermined.
        // PhiX is the run's spike-in control, not a sample: no FASTQ delivered,
        // no FastQC. Its numbers stay in bases2fastq's own reports.
        files = files.findAll { f -> !f.name.startsWith('PhiX_') }

        def r1s = files.findAll { f -> f.name ==~ /.*_R1(_001)?\.fastq\.gz$/ }
        if ( !r1s.findAll { f -> !f.name.startsWith('Unassigned') } ) {
            error "AVITI_QC: bases2fastq produced no R1 FASTQ. Check the run manifest in ${meta.id}."
        }

        r1s.sort().collect { r1 ->
            def r2  = files.find { f -> f.name == r1.name.replace('_R1', '_R2') }
            def sid = r1.name.replaceAll(/_R1(_001)?\.fastq\.gz$/, '')
            tuple( meta + [ id: sid, single_end: r2 == null, undetermined: r1.name.startsWith('Unassigned') ], r1, r2 ?: [] )
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

    ch_qc = channel.empty()
        .mix( FASTQC_FASTQC.out.zip.map { _meta, zip -> zip } )
        .mix( BASES2FASTQ.out.stats.flatMap { _meta, stats -> asFileList( stats ).findAll { f -> !isElementReport( f ) } } )
        .mix( BASES2FASTQ.out.sample_stats.flatMap { _meta, stats -> asFileList( stats ) } )
        .mix( ch_read_share )
        .mix( ch_demux_summary )
        .mix( ch_unknown_barcodes )

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
