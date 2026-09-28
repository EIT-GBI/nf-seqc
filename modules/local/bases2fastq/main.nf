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
    tuple val(meta), path("{RunStats.json,Metrics.csv,RunManifest.csv}"), emit: stats, optional: true
    tuple val(meta), path("*_stats.json"), emit: sample_stats, optional: true
    tuple val("${task.process}"), val('bases2fastq'), eval('bases2fastq --version 2>&1 | sed "s/^.*version //;s/ .*//" | head -n1'), emit: versions_bases2fastq, topic: versions

    script:
    def args = task.ext.args ?: ''
    """
    bases2fastq \\
        ${args} \\
        --num-threads ${task.cpus} \\
        ${run_dir} \\
        b2f

    # Flatten before emitting. bases2fastq nests FASTQ under
    # Samples/<project>/<sample>/, and the output block preserves those relative
    # subpaths, so without this the published tree becomes
    # reads/b2f/Samples/<project>/<sample>/… instead of reads/<sample>_R1.fastq.gz.
    find b2f/Samples -name '*.fastq.gz' -exec mv {} . \\;
    find b2f/Samples -name '*_stats.json' -exec mv {} . \\;
    for m in RunStats.json Metrics.csv RunManifest.csv; do
        if [ -f "b2f/\$m" ]; then mv "b2f/\$m" .; fi
    done
    """

    stub:
    """
    touch SAMPLE1_R1.fastq.gz SAMPLE1_R2.fastq.gz
    touch SAMPLE1_stats.json
    touch RunStats.json Metrics.csv RunManifest.csv
    """
}
