//
// Runs COPY_READS for real (not -stub) on small files, so the cp command and
// the renaming are tested; the stub only touches the target names.
//
//     nextflow run tests/copy_reads.nf --src <dir> --outdir <dir>
//
// <src> must hold sample.bam, sample.bam.pbi and Undetermined_R1.fastq.gz.
// CI then checks the published names and contents.
//

include { COPY_READS } from '../modules/local/copy_reads/main.nf'

params.src    = null
params.outdir = 'results-copy-test'

workflow {
    main:
    def src = file( params.src, checkIfExists: true )
    ch_in = channel.of(
        // Revio: a BAM and its index, renamed after the sample.
        tuple( [ id: 'LibA' ], [ src.resolve('sample.bam'), src.resolve('sample.bam.pbi') ], [ 'LibA.hifi_reads.bam', 'LibA.hifi_reads.bam.pbi' ] ),
        // Illumina: a single file, keeping its name.
        tuple( [ id: 'Undetermined' ], [ src.resolve('Undetermined_R1.fastq.gz') ], [ 'Undetermined_R1.fastq.gz' ] ),
    )
    COPY_READS( ch_in )

    publish:
    reads = COPY_READS.out.reads
}

output {
    reads {
        path '.'
    }
}
