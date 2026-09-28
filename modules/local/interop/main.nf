// Illumina run-level metrics from the binary InterOp files, RunParameters.xml
// and PrimaryAnalysisMetrics.csv.
//
// Mirrors EIT-GBI/Automate-Seq-Run-Metrics-Collection: same metric set, same
// column names, same %Occupied vs %PF tile plot, so rows from either tool line
// up. The extraction itself lives in collect_illumina_metrics.py alongside this
// file (same pattern as nf-mod-utils/samplesheet).
//
// Local module — promote to EIT-GBI/nf-mod-interop when a second pipeline needs
// it. Publishing is handled by the workflow `output {}` block, so no publishDir.

process INTEROP_RUN_METRICS {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(run_dir)

    output:
    tuple val(meta), path("${meta.id}.run_metrics.csv"), emit: metrics
    tuple val(meta), path("${meta.id}_run_metrics_mqc.tsv"), emit: mqc_table
    tuple val(meta), path("${meta.id}_occupied_vs_pf_mqc.png"), emit: plot, optional: true
    // The `|| echo` matters: eval commands run even under -stub, and a non-zero
    // exit fails the whole task rather than just leaving the version blank.
    tuple val("${task.process}"), val('illumina-interop'), eval('python3 -c "import interop; print(getattr(interop, \'__version__\', \'unknown\'))" 2>/dev/null || echo unknown'), emit: versions_interop, topic: versions

    script:
    def args = task.ext.args ?: ''
    """
    python3 ${moduleDir}/collect_illumina_metrics.py \\
        ${args} \\
        --run-dir ${run_dir} \\
        --run-id ${meta.id} \\
        --out-prefix ${meta.id}
    """

    stub:
    """
    touch ${meta.id}.run_metrics.csv
    touch ${meta.id}_run_metrics_mqc.tsv
    touch ${meta.id}_occupied_vs_pf_mqc.png
    """
}
