// A delivered BAM's integrity and modified-base tags, before scientists get it:
//
//   - samtools quickcheck -u: the header and the BGZF end-of-file block are
//     there, so the file is not truncated (-u: unaligned BAM, no @SQ needed).
//     A read error while scanning the first reads (a block damaged in the
//     middle, which quickcheck cannot see) also counts as FAIL.
//   - MM/ML tags in the first `max_reads` reads: which modifications the reads
//     carry (6mA, 5mC, ...) and the share of reads tagged. Modified bases are
//     only called when asked for (--basecall.modified_bases on Nanopore), and
//     an assay such as DiMeLo-seq is useless without them.
//
// Writes one tab-separated row: sample, quickcheck, reads_scanned,
// reads_with_mm, reads_with_ml, mods (comma-separated MM codes, e.g. A+a,C+m).
// Never fails the run on a bad BAM: the row says FAIL and the pipeline warns.
//
// Local module on the nf-mod-samtools image; promote to a quickcheck/modtags
// subcommand of nf-mod-samtools if another pipeline needs it.

process BAM_CHECK {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(bam, stageAs: 'in/*')
    val max_reads

    output:
    tuple val(meta), path("${meta.id}.bam_check.tsv"), emit: check
    tuple val("${task.process}"), val('samtools'), eval('samtools --version | head -n1 | sed "s/^samtools //"'), emit: versions_samtools, topic: versions

    script:
    // Only the BAM, not a .pbi next to it.
    def bams = ( bam instanceof List ? bam : [ bam ] ).findAll { f -> f.name.endsWith('.bam') }
    def b    = bams ? bams[0] : bam
    """
    if samtools quickcheck -u -v "${b}" > quickcheck.txt 2>&1; then quick=PASS; else quick=FAIL; fi
    # head closes the pipe early; samtools' resulting SIGPIPE is expected, and
    # so is the error it then writes about the closed pipe.
    { samtools view "${b}" 2> view.err || true; } | head -n ${max_reads} | awk -F '\\t' '
        {
            n++
            mm = ""; ml = 0
            for ( i = 12; i <= NF; i++ ) {
                if ( substr( \$i, 1, 5 ) == "MM:Z:" || substr( \$i, 1, 5 ) == "Mm:Z:" ) mm = substr( \$i, 6 )
                if ( substr( \$i, 1, 5 ) == "ML:B:" || substr( \$i, 1, 5 ) == "Ml:B:" ) ml = 1
            }
            if ( mm != "" ) {
                with_mm++
                k = split( mm, parts, ";" )
                for ( j = 1; j <= k; j++ ) {
                    code = parts[j]; sub( /,.*/, "", code ); sub( /[.?]\$/, "", code )
                    if ( code != "" ) codes[code] = 1
                }
            }
            if ( ml ) with_ml++
        }
        END {
            list = ""
            for ( c in codes ) list = list ( list == "" ? "" : "," ) c
            printf "%d\\t%d\\t%d\\t%s\\n", n + 0, with_mm + 0, with_ml + 0, list
        }' > tags.tsv
    # Any other samtools error is a read error in the reads scanned: a BAM
    # damaged in the middle, which quickcheck (header and end only) misses.
    if grep -viqE 'broken pipe|error closing standard output|writing to standard output|failed to write' view.err; then quick=FAIL; fi
    printf '%s\\t%s\\t' '${meta.id}' "\$quick" > ${meta.id}.bam_check.tsv
    cat tags.tsv >> ${meta.id}.bam_check.tsv
    """

    stub:
    """
    printf '%s\\tPASS\\t100\\t0\\t0\\t\\n' '${meta.id}' > ${meta.id}.bam_check.tsv
    """
}
