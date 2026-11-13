// Per-barcode read statistics (reads, yield, mean length, N50, median Q) and
// the Nanopore QC plots of the MultiQC report, from one flow cell's MinKNOW
// sequencing_summary_*.txt and report_*.json. No BAM is read for QC. The
// summaries describe MinKNOW's live basecalls and barcoding, not the
// pipeline's re-basecalled reads.
//
// meta.flowcell names flow-cell-level series; meta.aliases ([ barcode: alias ]),
// meta.plain_name (an unbarcoded flow cell's sample) and meta.multi_flowcell
// make the plot series match the pipeline's sample names (ontSamples() in
// nanopore_qc.nf). Barcodes with fewer than min_reads reads are left out of
// the plots, as they are left out of the samples.
//
// Local module — nanopore-specific and pipeline-shaped; promote to an nf-mod
// repo only if another pipeline needs it. Publishing is handled by the
// workflow `output {}` block. The script ships as a module binary
// (resources/usr/bin/), so editing it does not invalidate -resume.

process ONT_SUMMARY_STATS {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(summaries, stageAs: 'summary?/*'), path(reports, stageAs: 'report?/*')
    val min_reads

    output:
    tuple val(meta), path("${meta.id}.summary_stats.tsv"), emit: stats
    path "mqc/*_mqc.json", emit: mqc, optional: true
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/^Python //"'), emit: versions_python, topic: versions

    script:
    // Single-quoted for the shell, so a ' in an alias is closed and escaped.
    def aliases = groovy.json.JsonOutput.toJson( meta.aliases ?: [:] ).replace("'", "'\\''")
    def report  = reports ? "--report ${reports}" : ''
    def plain   = meta.plain_name ? "--plain-name '${meta.plain_name}'" : ''
    def suffix  = meta.multi_flowcell ? '--suffix-flowcell' : ''
    """
    ont_summary_stats.py \\
        -o ${meta.id}.summary_stats.tsv \\
        --flowcell '${meta.flowcell}' \\
        --aliases '${aliases}' \\
        --min-reads ${min_reads} \\
        --mqc-dir mqc \\
        ${report} ${plain} ${suffix} \\
        ${summaries}
    """

    stub:
    """
    printf 'barcode\\tpass_reads\\tfail_reads\\tpass_bases\\tmean_length\\tn50\\tmedian_q\\n' > ${meta.id}.summary_stats.tsv
    printf 'barcode01\\t1000\\t50\\t5000000\\t5000\\t8000\\t18.5\\n' >> ${meta.id}.summary_stats.tsv
    printf 'barcode02\\t2000\\t80\\t8000000\\t4000\\t7000\\t19.0\\n' >> ${meta.id}.summary_stats.tsv
    printf 'unclassified\\t100\\t300\\t200000\\t2000\\t3000\\t12.0\\n' >> ${meta.id}.summary_stats.tsv
    printf 'all\\t3100\\t430\\t13200000\\t4258\\t7500\\t18.7\\n' >> ${meta.id}.summary_stats.tsv
    mkdir mqc
    echo '{"id": "ont_pass_fail", "plot_type": "bargraph", "data": {"barcode01": {"pass": 1000, "fail": 50}}}' > mqc/${meta.id}.ont_pass_fail_mqc.json
    """
}
