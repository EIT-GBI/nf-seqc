//
// Demultiplexing report, shared by ILLUMINA_QC and AVITI_QC: each platform
// reads its own demultiplexer's counts (BCL Convert, bases2fastq) into the
// map demuxSummary() returns, and gets the same MultiQC tables.
//
//   - General Statistics: % Reads and Loading per sample
//   - Demultiplexing: one run-level row, the unassigned share and its level
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
// went to no sample. Returns [ samples, unassigned, total, even ], or null when
// there is nothing to report on. `even` is one sample's share of all reads if
// the assigned reads were split evenly.
def demuxSummary( Map samples, long unassigned ) {
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
    ]
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
