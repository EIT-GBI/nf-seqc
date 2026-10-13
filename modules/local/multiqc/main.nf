// Single batch-level report aggregating every per-sample and run-level QC file.
//
// Local module — promote to EIT-GBI/nf-mod-multiqc when a second pipeline needs
// it. Publishing is handled by the workflow `output {}` block.

process MULTIQC {
    tag "${run_id}"
    label 'process_low'

    input:
    path qc_files, stageAs: 'qc/*'
    val  run_id
    // One or more configs, applied in order; staged under numbered directories
    // so two files with the same name cannot collide.
    path multiqc_configs, stageAs: 'config?/*'
    path software_versions

    output:
    path "${run_id}_multiqc_report.html", emit: report
    path "${run_id}_multiqc_report_data", emit: data
    tuple val("${task.process}"), val('multiqc'), eval('multiqc --version | sed "s/^multiqc, version //"'), emit: versions_multiqc, topic: versions

    script:
    def args   = task.ext.args ?: ''
    // A single staged file arrives as a Path, which iterates as its directory
    // components, so wrap it before building the arguments.
    def config_list = multiqc_configs instanceof List ? multiqc_configs : [ multiqc_configs ]
    def config = config_list.collect { c -> "--config ${c}" }.join(' ')
    """
    multiqc \\
        ${args} \\
        ${config} \\
        --title "${run_id}" \\
        --filename ${run_id}_multiqc_report.html \\
        .
    """

    stub:
    """
    touch ${run_id}_multiqc_report.html
    mkdir -p ${run_id}_multiqc_report_data
    """
}
