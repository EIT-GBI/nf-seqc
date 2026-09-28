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
    path multiqc_config
    path software_versions

    output:
    path "${run_id}_multiqc_report.html", emit: report
    path "${run_id}_multiqc_report_data", emit: data
    tuple val("${task.process}"), val('multiqc'), eval('multiqc --version | sed "s/^multiqc, version //"'), emit: versions_multiqc, topic: versions

    script:
    def args   = task.ext.args ?: ''
    def config = multiqc_config.name != 'NO_FILE' ? "--config ${multiqc_config}" : ''
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
