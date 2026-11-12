// Per-barcode read statistics (reads, yield, mean length, N50, median Q) from
// a flow cell's MinKNOW sequencing_summary_*.txt files, so no BAM is read for
// QC. The summaries describe MinKNOW's live basecalls and barcoding, not the
// pipeline's re-basecalled reads.
//
// Local module — nanopore-specific and pipeline-shaped; promote to an nf-mod
// repo only if another pipeline needs it. Publishing is handled by the
// workflow `output {}` block. The script ships as a module binary
// (resources/usr/bin/), so editing it does not invalidate -resume.

process ONT_SUMMARY_STATS {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(summaries, stageAs: 'summary?/*')

    output:
    tuple val(meta), path("${meta.id}.summary_stats.tsv"), emit: stats
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/^Python //"'), emit: versions_python, topic: versions

    script:
    """
    ont_summary_stats.py -o ${meta.id}.summary_stats.tsv ${summaries}
    """

    stub:
    """
    printf 'barcode\\tpass_reads\\tfail_reads\\tpass_bases\\tmean_length\\tn50\\tmedian_q\\n' > ${meta.id}.summary_stats.tsv
    printf 'barcode01\\t1000\\t50\\t5000000\\t5000\\t8000\\t18.5\\n' >> ${meta.id}.summary_stats.tsv
    printf 'barcode02\\t2000\\t80\\t8000000\\t4000\\t7000\\t19.0\\n' >> ${meta.id}.summary_stats.tsv
    printf 'unclassified\\t100\\t300\\t200000\\t2000\\t3000\\t12.0\\n' >> ${meta.id}.summary_stats.tsv
    printf 'all\\t3100\\t430\\t13200000\\t4258\\t7500\\t18.7\\n' >> ${meta.id}.summary_stats.tsv
    """
}
