// Element Biosciences bases2fastq: converts an AVITI run folder into
// demultiplexed FASTQ plus run-level metrics. CPU-only (no GPU support).
//
// Local module — promote to EIT-GBI/nf-mod-bases2fastq when a second pipeline
// needs it. Publishing is handled by the workflow `output {}` block.

process BASES2FASTQ {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(run_dir)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    // RunManifest.json gives MultiQC (>= 1.29) its lane-settings and per-sample
    // index-assignment sections; *_QC.html is Element's own run report.
    tuple val(meta), path("{RunStats.json,Metrics.csv,RunManifest.csv,RunManifest.json,*_QC.html}"), emit: stats, optional: true
    tuple val(meta), path("*_stats.json"), emit: sample_stats, optional: true
    // For scientists, not MultiQC: reads and index sequences per sample, the
    // most common unmatched index sequences, the run log, and the run-manifest
    // settings bases2fastq assumed.
    tuple val(meta), path("{IndexAssignment.csv,UnassignedSequences.csv,Bases2Fastq.log,RunManifestErrors.json}"), emit: run_files, optional: true
    // bases2fastq prints "bases2fastq version: 2.4.0.2499042781" (note the
    // colon); take the first version-like token after "version".
    tuple val("${task.process}"), val('bases2fastq'), eval('bases2fastq --version 2>&1 | sed -nE "s/.*version:?[[:space:]]*([0-9][^ ,]*).*/\\1/p" | head -n1'), emit: versions_bases2fastq, topic: versions

    script:
    def args = task.ext.args ?: ''
    // --skip-multi-qc: bases2fastq otherwise runs its own bundled MultiQC
    // (1.31), which duplicates BATCH_REPORT and crashed on a real run.
    """
    bases2fastq \\
        ${args} \\
        --skip-multi-qc \\
        --num-threads ${task.cpus} \\
        ${run_dir} \\
        b2f

    # Flatten before emitting. bases2fastq nests FASTQ under
    # Samples/<project>/<sample>/, and the output block preserves those relative
    # subpaths, so without this the published tree becomes
    # reads/b2f/Samples/<project>/<sample>/… instead of reads/<sample>_R1.fastq.gz.
    find b2f/Samples -name '*.fastq.gz' -exec mv {} . \\;
    find b2f/Samples -name '*_stats.json' -exec mv {} . \\;
    for m in RunStats.json Metrics.csv RunManifest.csv RunManifest.json; do
        if [ -f "b2f/\$m" ]; then mv "b2f/\$m" .; fi
    done
    # Run-level QC report only: the per-project ones under Samples/ stay behind.
    find b2f -maxdepth 1 -name '*_QC.html' -exec mv {} . \\;
    for m in IndexAssignment.csv UnassignedSequences.csv info/Bases2Fastq.log info/RunManifestErrors.json; do
        if [ -f "b2f/\$m" ]; then mv "b2f/\$m" .; fi
    done
    """

    stub:
    """
    touch SAMPLE1_R1.fastq.gz SAMPLE1_R2.fastq.gz PhiX_R1.fastq.gz PhiX_R2.fastq.gz Unassigned_R1.fastq.gz Unassigned_R2.fastq.gz
    touch SAMPLE1_stats.json PhiX_stats.json Unassigned_stats.json
    touch RunStats.json Metrics.csv RunManifest.csv RunManifest.json
    touch ${meta.id}_QC.html
    touch IndexAssignment.csv UnassignedSequences.csv Bases2Fastq.log RunManifestErrors.json
    """
}
