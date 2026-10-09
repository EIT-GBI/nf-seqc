//
// Subworkflow: PACBIO_QC
//
// Revio basecalls on-instrument, so HiFi reads arrive as unaligned BAM and are
// delivered as-is (no FASTQ conversion requested). QC is reference-free read
// statistics plus the CCS reports the instrument already wrote.
//

include { SAMTOOLS_FLAGSTAT } from '../../modules/samtools/flagstat/main.nf'
include { NANOPLOT          } from '../../modules/local/nanoplot/main.nf'

workflow PACBIO_QC {

    take:
    ch_run // channel: [ val(meta), path(run_dir) ]

    main:
    ch_versions = channel.empty()

    // One BAM per SMRT Cell, or per barcode on a barcoded cell.
    ch_bam = ch_run.flatMap { meta, run_dir ->
        def root = params.bam_dir ? file(params.bam_dir, checkIfExists: true) : run_dir
        def bams = ( files("${root}/**/*.hifi_reads*.bam") + files("${root}/*.hifi_reads*.bam") )
            .findAll { f -> f.exists() }
            .unique()

        if ( !params.include_failed ) {
            bams = bams.findAll { f -> !f.name.contains('unassigned') && !f.name.contains('fail_reads') }
        }
        if ( !bams ) {
            error "PACBIO_QC: no *.hifi_reads*.bam found under ${root}. Point --bam_dir at the SMRT Cell output if it lives elsewhere."
        }

        bams.sort().collect { bam ->
            tuple( meta + [ id: bam.name.replaceAll(/\.bam$/, '') ], bam )
        }
    }

    // Instrument-written run reports: CCS yield/quality and the dataset manifest.
    ch_reports = ch_run.flatMap { _meta, run_dir ->
        ( files("${run_dir}/**/*.ccs_report*") + files("${run_dir}/**/*.consensusreadset.xml") )
            .findAll { f -> f.exists() }
            .unique()
    }

    SAMTOOLS_FLAGSTAT( ch_bam )
    NANOPLOT( ch_bam )

    ch_qc = channel.empty()
        .mix( NANOPLOT.out.stats.map { _meta, stats -> stats } )
        .mix( SAMTOOLS_FLAGSTAT.out.flagstat.map { _meta, txt -> txt } )
        .mix( ch_reports )

    ch_versions = ch_versions
        .mix( NANOPLOT.out.versions_nanoplot )

    emit:
    reads    = ch_bam // HiFi uBAM delivered unchanged
    qc       = ch_qc
    versions = ch_versions
}
