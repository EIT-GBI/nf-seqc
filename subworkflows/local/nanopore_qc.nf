//
// Subworkflow: NANOPORE_QC
//
// POD5 arrive already split per barcode, so no dorado demux step is needed —
// each directory holding *.pod5 is one sample. Dorado basecalls on GPU to an
// unaligned BAM (keeping MM/ML tags when params.basecall.modified_bases is set),
// which is then converted to FASTQ for delivery and QC'd reference-free.
//

include { DORADO_BASECALLER } from '../../modules/dorado/basecaller/main.nf'
include { SAMTOOLS_FASTQ    } from '../../modules/samtools/fastq/main.nf'
include { NANOPLOT          } from '../../modules/local/nanoplot/main.nf'

workflow NANOPORE_QC {

    take:
    ch_run // channel: [ val(meta), path(run_dir) ]

    main:
    ch_versions = channel.empty()

    // Group POD5 by containing directory. A barcode-named directory becomes the
    // sample id; a generic container (pod5/, pod5_pass/) means the run wasn't
    // barcoded, so the run id is used instead.
    ch_pod5 = ch_run.flatMap { meta, run_dir ->
        def root  = params.pod5_dir ? file(params.pod5_dir, checkIfExists: true) : run_dir
        def pod5s = ( [] + file("${root}/**/*.pod5") + file("${root}/*.pod5") )
            .findAll { f -> f.exists() }
            .unique()

        if ( !params.include_failed ) {
            pod5s = pod5s.findAll { f -> !f.toString().contains('_fail') }
        }
        if ( !pod5s ) {
            error "NANOPORE_QC: no *.pod5 found under ${root}. Point --pod5_dir at the POD5 directory if they live elsewhere."
        }

        def generic = [ 'pod5', 'pod5_pass', 'pod5_fail', 'pod5_skip' ]
        pod5s.groupBy { f -> f.parent }.collect { dir, _files ->
            tuple( meta + [ id: dir.name in generic ? meta.id : dir.name ], dir )
        }.sort { t -> t[0].id }
    }

    DORADO_BASECALLER( ch_pod5 )
    SAMTOOLS_FASTQ( DORADO_BASECALLER.out.bam )
    NANOPLOT( DORADO_BASECALLER.out.bam )

    ch_qc = channel.empty()
        .mix( NANOPLOT.out.stats.map { _meta, stats -> stats } )

    ch_versions = ch_versions
        .mix( NANOPLOT.out.versions_nanoplot )

    emit:
    reads    = SAMTOOLS_FASTQ.out.reads          // FASTQ for scientists
    bam      = DORADO_BASECALLER.out.bam         // uBAM, carries modified-base tags
    qc       = ch_qc
    versions = ch_versions
}
