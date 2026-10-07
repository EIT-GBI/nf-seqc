//
// Subworkflow: AVITI_QC
//
// The AVITI run folder holds no FASTQ, so conversion is part of this pipeline's
// job rather than a precondition: bases2fastq produces the FASTQ scientists get
// plus the run-level metrics (RunStats.json / Metrics.csv) that MultiQC reads.
//

include { BASES2FASTQ   } from '../../modules/local/bases2fastq/main.nf'
include { FASTQC_FASTQC } from '../../modules/fastqc/fastqc/main.nf'

// A process output declared as path("*.x") is a single Path when one file
// matches and a List when several do. `[] + path` does NOT wrap a single Path:
// Path is Iterable over its directory components, so it yields every parent
// directory (home, runner, work, ...) as if it were a file. Wrap explicitly.
def asFileList( files ) {
    return files instanceof List ? files : [ files ]
}

workflow AVITI_QC {

    take:
    ch_run // channel: [ val(meta), path(run_dir) ]

    main:
    ch_versions = channel.empty()

    BASES2FASTQ( ch_run )

    // bases2fastq writes Samples/<project>/<sample>/<sample>_R{1,2}.fastq.gz.
    ch_pairs = BASES2FASTQ.out.fastq.flatMap { meta, fastqs ->
        def files = asFileList( fastqs )

        if ( !params.include_undetermined ) {
            files = files.findAll { f -> !f.name.startsWith('Unassigned') }
        }

        def r1s = files.findAll { f -> f.name ==~ /.*_R1(_001)?\.fastq\.gz$/ }
        if ( !r1s ) {
            error "AVITI_QC: bases2fastq produced no R1 FASTQ. Check the run manifest in ${meta.id}."
        }

        r1s.sort().collect { r1 ->
            def r2  = files.find { f -> f.name == r1.name.replace('_R1', '_R2') }
            def sid = r1.name.replaceAll(/_R1(_001)?\.fastq\.gz$/, '')
            tuple( meta + [ id: sid, single_end: r2 == null ], r1, r2 ?: [] )
        }
    }

    FASTQC_FASTQC( ch_pairs )

    ch_qc = channel.empty()
        .mix( FASTQC_FASTQC.out.zip.map { _meta, zip -> zip } )
        .mix( BASES2FASTQ.out.stats.flatMap { _meta, stats -> asFileList( stats ) } )
        .mix( BASES2FASTQ.out.sample_stats.flatMap { _meta, stats -> asFileList( stats ) } )

    ch_versions = ch_versions
        .mix( FASTQC_FASTQC.out.versions_fastqc )
        .mix( BASES2FASTQ.out.versions_bases2fastq )

    emit:
    reads    = ch_pairs.map { meta, r1, r2 -> tuple( meta, [ r1, r2 ].flatten().findAll { f -> f } ) }
    qc       = ch_qc
    versions = ch_versions
}
