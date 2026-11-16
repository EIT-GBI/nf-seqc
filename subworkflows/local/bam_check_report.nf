//
// General Statistics columns from BAM_CHECK, for the long-read platforms,
// whose reads are delivered as BAM: whether each delivered BAM is intact, and
// which modified bases its reads carry.
//

include { bgcols; loadingBgcols } from './demux_report.nf'

// An MM tag code (base, strand, modification: C+m, A+a, T-a) as people say the
// modification. The strand and base are dropped: Revio calls 6mA on both
// strands (A+a and T-a), and that is one modification. Unknown codes are shown
// as written.
def modName( String code ) {
    def names = [ a: '6mA', m: '5mC', h: '5hmC', '21839': '4mC', '17596': 'inosine', '17802': 'pseudouridine', o: '8oxoG' ]
    def m = code =~ /^[ACGTUN][+-](.+)$/
    return m.matches() && names[ m.group(1) ] ? names[ m.group(1) ] : code
}

// One BAM_CHECK row: [ id:, quickcheck:, reads:, with_mm:, with_ml:, codes: [] ].
def bamCheckRow( String tsv_text ) {
    def f = ( tsv_text.trim() + '\t\t\t\t\t' ).split('\t', -1)
    return [
        id:         f[0],
        quickcheck: f[1],
        reads:      f[2] ? f[2] as long : 0L,
        with_mm:    f[3] ? f[3] as long : 0L,
        with_ml:    f[4] ? f[4] as long : 0L,
        codes:      f[5] ? f[5].tokenize(',').sort() : [],
    ]
}

// Warns in the Nextflow log about every truncated BAM, and, when modified
// bases were asked for (`expect_mods`), every BAM whose reads carry none.
def bamCheckWarnings( List rows, boolean expect_mods, String caller ) {
    rows.findAll { r -> r.quickcheck != 'PASS' }.each { r ->
        log.warn "${caller}: ${r.id}: samtools quickcheck failed. The delivered BAM is truncated or corrupt."
    }
    if ( expect_mods ) {
        rows.findAll { r -> r.quickcheck == 'PASS' && r.reads && !r.with_mm }.each { r ->
            log.warn "${caller}: ${r.id}: modified bases were requested but none of its first ${r.reads} reads carry an MM tag."
        }
    }
}

// MultiQC General Statistics columns, one row per sample: BAM check
// (PASS/FAIL), modified bases found, and the share of the scanned reads
// carrying an MM tag. null when there are no rows.
def bamCheckTable( List rows, int max_reads ) {
    if ( !rows ) {
        return null
    }
    def header = [
        "# id: 'bam_check'",
        "# plot_type: 'generalstats'",
        "# headers:",
        "#     bam_ok:",
        "#         title: 'BAM'",
        "#         namespace: 'BAM check'",
        "#         description: 'samtools quickcheck of the delivered BAM: FAIL means truncated or corrupt'",
        "#         scale: false",
    ] + bgcols( [ PASS: loadingBgcols().optimal, FAIL: loadingBgcols().over ] ) + [
        "#     mods:",
        "#         title: 'Mod. bases'",
        "#         namespace: 'BAM check'",
        "#         description: 'Modified bases in the reads'' MM tags (first ${max_reads} reads); none unless modified-base calling was requested'",
        "#         scale: false",
        "#     pct_mods:",
        "#         title: '% Mod. reads'",
        "#         namespace: 'BAM check'",
        "#         description: 'Share of the first ${max_reads} reads carrying an MM tag'",
        "#         suffix: '%'",
        "#         min: 0",
        "#         max: 100",
        "#         format: '{:,.1f}'",
        "#         hidden: true",
        "Sample\tbam_ok\tmods\tpct_mods",
    ]
    def lines = rows.collect { r ->
        def mods = r.codes ? r.codes.collect { c -> modName( c ) }.unique().sort().join(', ') : 'none'
        [ r.id, r.quickcheck, mods, r.reads ? 100.0 * r.with_mm / r.reads : '' ].join('\t')
    }
    return ( header + lines ).join('\n') + '\n'
}
