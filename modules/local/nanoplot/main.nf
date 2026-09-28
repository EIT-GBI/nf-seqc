// Long-read QC plots and summary stats, reference-free. Accepts either an
// unaligned BAM (ONT dorado output, PacBio HiFi reads) or FASTQ.
//
// Local module — promote to EIT-GBI/nf-mod-nanoplot when a second pipeline
// needs it. Publishing is handled by the workflow `output {}` block.

process NANOPLOT {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}/*.html"), emit: html, optional: true
    tuple val(meta), path("${meta.id}/*.png"), emit: png, optional: true
    tuple val(meta), path("${meta.id}/${meta.id}_NanoStats.txt"), emit: stats
    tuple val("${task.process}"), val('nanoplot'), eval('NanoPlot --version | sed "s/^NanoPlot //"'), emit: versions_nanoplot, topic: versions

    script:
    def args = task.ext.args ?: ''
    def input_flag = reads.name.endsWith('.bam') ? '--ubam' : '--fastq'
    """
    NanoPlot \\
        ${args} \\
        --threads ${task.cpus} \\
        --prefix ${meta.id}_ \\
        --outdir ${meta.id} \\
        ${input_flag} ${reads}
    """

    stub:
    """
    mkdir -p ${meta.id}
    touch ${meta.id}/${meta.id}_NanoStats.txt
    touch ${meta.id}/${meta.id}_LengthvsQualityScatterPlot_dot.html
    """
}
