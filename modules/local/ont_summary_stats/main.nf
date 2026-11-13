// Per-barcode read statistics (reads, yield, mean length, N50, median Q) and
// the Nanopore QC plots of the MultiQC report, from one flow cell's MinKNOW
// sequencing_summary_*.txt and report_*.json. No BAM is read for QC. The
// summaries describe MinKNOW's live basecalls and barcoding, not the
// pipeline's re-basecalled reads.
//
// meta.flowcell names flow-cell-level series. meta.samples ([ barcode: id ],
// 'all' for an unbarcoded flow cell) is the pipeline's own sample ids, so the
// plot series match General Statistics and only reported samples are
// plotted. meta.offsets gives, per summary, how many seconds its acquisition
// started after the flow cell's first, so a restart follows it on the time
// axis.
//
// Local module — nanopore-specific and pipeline-shaped; promote to an nf-mod
// repo only if another pipeline needs it. Publishing is handled by the
// workflow `output {}` block. The script ships as a module binary
// (resources/usr/bin/), so editing it does not invalidate -resume.

// A value single-quoted for the shell, any ' in it closed and escaped.
def shellQuote( value ) {
    return "'" + value.toString().replace("'", "'\\''") + "'"
}

process ONT_SUMMARY_STATS {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(summaries, stageAs: 'summary?/*'), path(reports, stageAs: 'report?/*')

    output:
    tuple val(meta), path("${meta.id}.summary_stats.tsv"), emit: stats
    path "mqc/*_mqc.json", emit: mqc, optional: true
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/^Python //"'), emit: versions_python, topic: versions

    script:
    def samples = groovy.json.JsonOutput.toJson( meta.samples ?: [:] )
    def offsets = ( meta.offsets ?: [] ).collect { o -> "--offset ${o}" }.join(' ')
    def report  = ( reports instanceof List ? reports : [ reports ] ).findAll { r -> r }.collect { r -> "--report ${r}" }.join(' ')
    """
    ont_summary_stats.py \\
        -o ${meta.id}.summary_stats.tsv \\
        --flowcell ${shellQuote( meta.flowcell )} \\
        --samples ${shellQuote( samples )} \\
        --mqc-dir mqc \\
        ${offsets} ${report} \\
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
