//
// Demultiplexing report, shared by ILLUMINA_QC and AVITI_QC: each platform
// reads its own demultiplexer's counts (BCL Convert, bases2fastq) into the
// map demuxSummary() returns, and gets the same MultiQC tables.
//
//   - General Statistics: % Reads and Loading per sample
//   - Demultiplexing: one run-level row, the unassigned share and its level
//   - Undetermined barcodes: the commonest index pairs that matched no sample
//

// Cell backgrounds for flag columns: the occupancy plot's zone colours
// (steelblue / green / firebrick), semi-transparent so the text stays readable
// in both the light and the dark theme. Header `bgcols` is how MultiQC's own
// modules flag categories; the same map is in nf-mod-interop's Run Metrics table.
def bgcols( Map colours ) {
    return [ "#         bgcols:" ] + colours.collect { value, colour -> "#             ${value}: '${colour}'" }
}
def loadingBgcols() { return [ under: 'rgba(70,130,180,0.4)', optimal: 'rgba(0,128,0,0.35)', over: 'rgba(178,34,34,0.4)' ] }

// samples: [ name: reads ] for the assigned samples; unassigned: reads that
// went to no sample; barcodes: the demultiplexer's top unassigned index pairs,
// as [ lane:, index:, reads: ] maps (empty when it lists none). Returns
// [ samples, unassigned, total, even, barcodes ], or null when there is nothing
// to report on. `even` is one sample's share of all reads if the assigned reads
// were split evenly.
def demuxSummary( Map samples, long unassigned, List barcodes ) {
    def assigned = samples.values().sum() ?: 0L
    def total    = assigned + unassigned
    if ( !total || !samples ) {
        return null
    }
    return [
        samples:    samples,
        unassigned: unassigned,
        total:      total,
        even:       100.0 * assigned / total / samples.size(),
        barcodes:   barcodes,
    ]
}

// An index pair as one string, i7+i5; a single-index run has no second index.
def barcodeName( index, index2 ) {
    return index2 ? "${index}+${index2}" : "${index}"
}

// A sample's loading in the pool: under below half an even split, over above
// twice one, optimal in between.
def sampleLoading( demux, long reads ) {
    def pct = 100.0 * reads / demux.total
    return pct < demux.even / 2 ? 'under' : pct > demux.even * 2 ? 'over' : 'optimal'
}

// MultiQC General Statistics columns: each sample's share of the run's reads
// (unassigned reads counted in the total but given no row), and its Loading.
// The cut-offs depend on how many samples share the run, so they are computed
// here and given in the column description.
def readShareTable( demux, String unassigned_label ) {
    def under = String.format( '%.2f', demux.even / 2 )
    def over  = String.format( '%.2f', demux.even * 2 )
    def header = [
        "# id: 'demux_read_share'",
        "# plot_type: 'generalstats'",
        "# headers:",
        "#     pct_reads:",
        "#         title: '% Reads'",
        "#         namespace: 'Demultiplexing'",
        "#         description: 'Share of all reads in the run, ${unassigned_label} included'",
        "#         suffix: '%'",
        "#         min: 0",
        "#         format: '{:,.2f}'",
        "#     loading:",
        "#         title: 'Loading'",
        "#         namespace: 'Demultiplexing'",
        "#         description: 'Share of reads against an even split of the pool: under below ${under}%, over above ${over}%'",
    ] + bgcols( loadingBgcols() ) + [
        "Sample\tpct_reads\tloading",
    ]
    def rows = demux.samples.collect { id, n -> "${id}\t${100.0 * n / demux.total}\t${sampleLoading( demux, n )}" }
    return ( header + rows ).join('\n') + '\n'
}

// Run-level demultiplexing summary, one row, placed under the run metrics. A
// high unassigned share is invisible in General Statistics, which only
// compares samples with each other, so it gets its own flag above
// --max_undetermined, and a warning in the Nextflow log (prefixed with
// `caller`). No colour-scale bars: in a one-row table every bar is full
// width, so each number turns into a solid block.
def demuxSummaryTable( run_id, demux, max_undetermined, String source, String unassigned_label, String caller ) {
    def pct_unassigned = 100.0 * demux.unassigned / demux.total
    def level = pct_unassigned > max_undetermined ? 'high' : 'ok'
    if ( level == 'high' ) {
        log.warn "${caller}: ${String.format('%.1f', pct_unassigned)}% of reads in ${run_id} are ${unassigned_label} (above --max_undetermined ${max_undetermined}%)."
    }
    def counts = [ under: 0, optimal: 0, over: 0 ]
    demux.samples.each { _id, n -> counts[ sampleLoading( demux, n ) ] += 1 }
    def header = [
        "# id: 'demux_summary'",
        "# section_name: 'Demultiplexing'",
        "# description: 'Run-level summary of ${source}. ${unassigned_label} is high above ${max_undetermined}% (--max_undetermined). Sample loading is the General Statistics Loading column.'",
        "# plot_type: 'table'",
        "# pconfig:",
        "#     id: 'demux_summary_table'",
        "#     namespace: 'Demultiplexing'",
        "# headers:",
        "#     pct_unassigned:",
        "#         title: '${unassigned_label}'",
        "#         scale: false",
        "#         description: 'Share of all reads not assigned to any sample'",
        "#         suffix: '%'",
        "#         min: 0",
        "#         max: 100",
        "#         format: '{:,.2f}'",
        "#     unassigned_level:",
        "#         title: '${unassigned_label} level'",
        "#         scale: false",
        "#         description: 'high above ${max_undetermined}% (--max_undetermined)'",
    ] + bgcols( [ ok: loadingBgcols().optimal, high: loadingBgcols().over ] ) + [
        "#     samples:",
        "#         title: 'Samples'",
        "#         scale: false",
        "#         format: '{:,.0f}'",
        "#     under:",
        "#         title: 'Under'",
        "#         scale: false",
        "#         description: 'Samples under half an even split'",
        "#         format: '{:,.0f}'",
        "#     optimal:",
        "#         title: 'Optimal'",
        "#         scale: false",
        "#         format: '{:,.0f}'",
        "#     over:",
        "#         title: 'Over'",
        "#         scale: false",
        "#         description: 'Samples over twice an even split'",
        "#         format: '{:,.0f}'",
        "Run\tpct_unassigned\tunassigned_level\tsamples\tunder\toptimal\tover",
        "${run_id}\t${pct_unassigned}\t${level}\t${demux.samples.size()}\t${counts.under}\t${counts.optimal}\t${counts.over}",
    ]
    return header.join('\n') + '\n'
}

// The commonest unassigned index pairs, top `per_lane` in each lane, as a
// table of its own. MultiQC cannot move a section out of its module, and the
// demultiplexers' own lists sit inside bclconvert / bases2fastq, so this one
// can be placed at the top of the report. Rows keep lane-then-reads order
// (sort_rows: false). null when the demultiplexer listed no barcodes.
//
// A barcode with at least as many reads as the Loading column's `under`
// cut-off (half an even split) is as big as a real sample, so it is flagged
// as a likely missing sample, and warned about in the Nextflow log (prefixed
// with `caller`).
def unknownBarcodesTable( demux, String source, String unassigned_label, int per_lane, String caller ) {
    if ( !demux.barcodes ) {
        return null
    }
    def rows = demux.barcodes
        .groupBy { b -> b.lane as int }
        .sort()
        .collectMany { _lane, bs -> bs.sort { b -> -b.reads }.take( per_lane ) }
    def sample_sized = demux.even / 2
    def missing = rows.findAll { b -> 100.0 * b.reads / demux.total >= sample_sized }
    if ( missing ) {
        log.warn "${caller}: ${missing.size()} ${unassigned_label} barcode(s) have as many reads as a sample (at least ${String.format('%.2f', sample_sized)}% of the run): ${missing.collect { b -> "L${b.lane} ${b.index}" }.join(', ')}. Check the sample sheet."
    }
    def header = [
        "# id: 'demux_unknown_barcodes'",
        "# section_name: '${unassigned_label} barcodes'",
        "# description: 'The ${per_lane} commonest index pairs per lane that matched no sample, from ${source}. A barcode with many reads here is often a sample missing from the sample sheet, or one with its index reverse-complemented. Flagged when a barcode has at least ${String.format('%.2f', sample_sized)}% of all reads, the share below which a sample counts as under-loaded.'",
        "# plot_type: 'table'",
        "# pconfig:",
        "#     id: 'demux_unknown_barcodes_table'",
        "#     namespace: 'Demultiplexing'",
        "#     col1_header: 'Lane and barcode'",
        "#     sort_rows: false",
        "# headers:",
        "#     reads:",
        "#         title: 'Reads'",
        "#         description: 'Reads with this index pair in this lane'",
        "#         format: '{:,.0f}'",
        "#     pct_total:",
        "#         title: '% of all reads'",
        "#         description: 'Share of all reads in the run'",
        "#         suffix: '%'",
        "#         min: 0",
        "#         format: '{:,.3f}'",
        "#     pct_unassigned:",
        "#         title: '% of ${unassigned_label}'",
        "#         description: 'Share of all ${unassigned_label} reads in the run'",
        "#         suffix: '%'",
        "#         min: 0",
        "#         format: '{:,.2f}'",
        "#     flag:",
        "#         title: 'Flag'",
        "#         scale: false",
        "#         description: 'likely missing sample: at least ${String.format('%.2f', sample_sized)}% of all reads, as many as an under-loaded sample'",
    ] + bgcols( [ 'likely missing sample': loadingBgcols().over ] ) + [
        "Barcode\treads\tpct_total\tpct_unassigned\tflag",
    ]
    def lines = rows.collect { b ->
        def pct_total      = 100.0 * b.reads / demux.total
        def pct_unassigned = demux.unassigned ? 100.0 * b.reads / demux.unassigned : 0
        def flag           = pct_total >= sample_sized ? 'likely missing sample' : ''
        "L${b.lane} ${b.index}\t${b.reads}\t${pct_total}\t${pct_unassigned}\t${flag}"
    }
    return ( header + lines ).join('\n') + '\n'
}
