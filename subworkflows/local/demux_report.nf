//
// Demultiplexing report, shared by ILLUMINA_QC and AVITI_QC: each platform
// reads its own demultiplexer's counts (BCL Convert, bases2fastq) into the
// map demuxSummary() returns, and gets the same MultiQC tables.
//
//   - General Statistics: % Reads and Loading per sample
//   - Demultiplexing: one run-level row, the unassigned share and its level,
//     and the PhiX spike-in share where the platform reports it
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
// [ samples, unassigned, total, even, barcodes, pools ], or null when there is
// nothing to report on. `even` is one sample's share of all reads if the
// assigned reads were split evenly. A platform that knows the PhiX spike-in
// share adds phix_pct (% of all reads) and phix_source to the map (see
// withPhix()).
//
// pools: [ sample: [ pool: reads ] ], where a pool is what was loaded
// separately: an AVITI or Illumina lane (PacBio and Nanopore pass none, so
// their SMRT Cells and flow cells share the run-wide split). Two lanes
// can hold two different pools of different sizes, so a sample's even split
// is that of the samples sharing its pool(s) (expectedShare()); [:] or one
// pool gives the run-wide split.
def demuxSummary( Map samples, long unassigned, List barcodes, Map pools = [:] ) {
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
        pools:      pools.findAll { id, _p -> samples.containsKey( id ) },
    ]
}

// Each pool's even split, as a share of all reads: its assigned reads over the
// samples in it, [ pool: % ]. Empty when the run has one pool or none known.
def poolEven( demux ) {
    def pools = demux.pools ?: [:]
    def names = pools.values().collectMany { p -> p.keySet() as List }.unique()
    if ( names.size() < 2 ) {
        return [:]
    }
    return names.collectEntries { pool ->
        def in_pool = pools.findAll { _id, p -> p.containsKey( pool ) }
        def reads   = in_pool.values().sum { p -> p[ pool ] as long } ?: 0L
        [ ( pool ): 100.0 * reads / demux.total / in_pool.size() ]
    }
}

// The share of all reads a sample would have if every pool it was in were
// split evenly between that pool's samples; the run-wide even split for a
// run of one pool.
def expectedShare( demux, String id ) {
    def even = poolEven( demux )
    def mine = demux.pools?.get( id )
    return even && mine ? mine.keySet().sum { pool -> even[ pool ] } : demux.even
}

// The demux map with the PhiX spike-in's share of all reads, and where it came
// from, for the Demultiplexing table. pct null (unknown) leaves it unchanged.
// For a platform whose PhiX reads are already in the demux total (Illumina:
// they are part of Undetermined).
def withPhix( demux, pct, String source ) {
    return demux && pct != null ? demux + [ phix_pct: pct as double, phix_source: source ] : demux
}

// A demultiplexer sample that is the PhiX spike-in rather than a sample: PhiX,
// or PhiX_<anything> (AVITI's control, or PhiX given an index in an Illumina
// sample sheet). It is kept out of the samples (no % Reads or Loading row, not
// in the even split) and its reads are added back with withPhixSampleReads().
def phixSample( name ) {
    def n = name.toString()
    return n == 'PhiX' || n.startsWith('PhiX_')
}

// The demux map with reads demuxSummary() did not count, because they went to
// a PhiX sample, added to the run total, so % Reads, the unassigned share and
// the even split stay shares of all reads.
def withPhixSampleReads( demux, long reads ) {
    if ( !demux || !reads ) {
        return demux
    }
    def total    = demux.total + reads
    def assigned = demux.samples.values().sum() ?: 0L
    return demux + [ total: total, even: 100.0 * assigned / total / demux.samples.size() ]
}

// The same for a platform whose PhiX reads all went to a PhiX sample (AVITI's
// control), with phix_pct their share of that total.
def withPhixReads( demux, long reads, String source ) {
    def with_reads = withPhixSampleReads( demux, reads )
    return with_reads ? withPhix( with_reads, 100.0 * reads / with_reads.total, source ) : with_reads
}

// An index pair as one string, i7+i5; a single-index run has no second index.
def barcodeName( index, index2 ) {
    return index2 ? "${index}+${index2}" : "${index}"
}

// A sample's loading in its pool: under below half an even split, over above
// twice one, optimal in between.
def sampleLoading( demux, String id, long reads ) {
    def pct  = 100.0 * reads / demux.total
    def even = expectedShare( demux, id )
    // A sample with no reads is under whatever its pool's split: a lane whose
    // only sample got nothing has an even split of 0.
    return !reads || pct < even / 2 ? 'under' : pct > even * 2 ? 'over' : 'optimal'
}

// MultiQC General Statistics columns: each sample's share of the run's reads
// (unassigned reads counted in the total but given no row), and its Loading.
// The cut-offs depend on how many samples share the run, so they are computed
// here and given in the column description.
def readShareTable( demux, String unassigned_label ) {
    def cutoffs = poolEven( demux )
        ? 'of the samples in its lane (lanes holding pools of different sizes get different cut-offs): under below half of it, over above twice'
        : "of the pool: under below ${String.format( '%.2f', demux.even / 2 )}%, over above ${String.format( '%.2f', demux.even * 2 )}%"
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
        "#         description: 'Share of reads against an even split ${cutoffs}'",
    ] + bgcols( loadingBgcols() ) + [
        "Sample\tpct_reads\tloading",
    ]
    // Sorted, so the file (a MultiQC input) is the same on every run and
    // -resume does not rerun MULTIQC for nothing.
    def rows = demux.samples.sort().collect { id, n -> "${id}\t${100.0 * n / demux.total}\t${sampleLoading( demux, id, n )}" }
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
    def phix = demux.phix_pct
    if ( level == 'high' ) {
        def phix_note = phix != null ? "; PhiX is ${String.format('%.1f', phix)}% of all reads" : ''
        log.warn "${caller}: ${String.format('%.1f', pct_unassigned)}% of reads in ${run_id} are ${unassigned_label} (above --max_undetermined ${max_undetermined}%)${phix_note}."
    }
    def counts = [ under: 0, optimal: 0, over: 0 ]
    demux.samples.each { id, n -> counts[ sampleLoading( demux, id, n ) ] += 1 }
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
    ] + bgcols( [ ok: loadingBgcols().optimal, high: loadingBgcols().over ] ) + ( phix == null ? [] : [
        "#     pct_phix:",
        "#         title: 'PhiX'",
        "#         scale: false",
        "#         description: 'Share of all reads that are PhiX spike-in, ${demux.phix_source}. A high ${unassigned_label} share with PhiX close to it is mostly spike-in; with little PhiX it points to the sample sheet.'",
        "#         suffix: '%'",
        "#         min: 0",
        "#         max: 100",
        "#         format: '{:,.2f}'",
    ] ) + [
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
        "Run\tpct_unassigned\tunassigned_level" + ( phix == null ? '' : '\tpct_phix' ) + "\tsamples\tunder\toptimal\tover",
        "${run_id}\t${pct_unassigned}\t${level}" + ( phix == null ? '' : "\t${phix}" ) + "\t${demux.samples.size()}\t${counts.under}\t${counts.optimal}\t${counts.over}",
    ]
    return header.join('\n') + '\n'
}

// An index read of only G and N is no index at all: on two-colour chemistry a
// dark cycle reads as G (unindexed PhiX, empty clusters).
def noIndex( index ) {
    return index.toString() ==~ /[GN+]+/
}

// The commonest unassigned index pairs, top `per_lane` in each lane, as a
// table of its own. MultiQC cannot move a section out of its module, and the
// demultiplexers' own lists sit inside bclconvert / bases2fastq, so this one
// can be placed at the top of the report. Rows keep lane-then-reads order
// (sort_rows: false). null when the demultiplexer listed no barcodes.
//
// A barcode with at least as many reads as the Loading column's `under`
// cut-off (half an even split) is as big as a real sample, so it is flagged
// as missing from the sample sheet, and warned about in the Nextflow log (prefixed
// with `caller`).
def unknownBarcodesTable( demux, String source, String unassigned_label, int per_lane, String caller ) {
    if ( !demux.barcodes ) {
        return null
    }
    def rows = demux.barcodes
        .groupBy { b -> b.lane as int }
        .sort()
        .collectMany { _lane, bs -> bs.sort { b -> -b.reads }.take( per_lane ) }
    // As many reads as an under-loaded sample of its lane: half that lane's
    // even split (the run's, with one lane). A barcode of only G and N is no
    // index at all (unindexed PhiX, empty clusters read as dark cycles), so
    // it is never taken for a missing sample.
    def lane_even = poolEven( demux )
    def missing   = rows.findAll { b -> !noIndex( b.index ) && 100.0 * b.reads / demux.total >= ( lane_even[ b.lane.toString() ] ?: demux.even ) / 2 }
    def cutoff = lane_even
        ? 'half the even split of its lane'
        : "${String.format('%.2f', demux.even / 2)}% of all reads"
    if ( missing ) {
        log.warn "${caller}: ${missing.size()} ${unassigned_label} barcode(s) have as many reads as an under-loaded sample (${cutoff}): ${missing.collect { b -> "L${b.lane} ${b.index}" }.join(', ')}. Check the sample sheet."
    }
    def header = [
        "# id: 'demux_unknown_barcodes'",
        "# section_name: '${unassigned_label} barcodes'",
        "# description: 'The ${per_lane} commonest index pairs per lane that matched no sample, from ${source}. A barcode with many reads here is often a sample missing from the sample sheet, or one with its index reverse-complemented. Flagged when a barcode has at least ${cutoff}, the share below which a sample counts as under-loaded; never for an index of only G and N (no index read).'",
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
        "#         description: 'missing from samplesheet: at least ${cutoff}, as many as an under-loaded sample, so likely a sample left out of the sample sheet or entered with the wrong index'",
    ] + bgcols( [ 'missing from samplesheet': loadingBgcols().over ] ) + [
        "Barcode\treads\tpct_total\tpct_unassigned\tflag",
    ]
    def lines = rows.collect { b ->
        def pct_total      = 100.0 * b.reads / demux.total
        def pct_unassigned = demux.unassigned ? 100.0 * b.reads / demux.unassigned : 0
        def flag           = b in missing ? 'missing from samplesheet' : ''
        "L${b.lane} ${b.index}\t${b.reads}\t${pct_total}\t${pct_unassigned}\t${flag}"
    }
    return ( header + lines ).join('\n') + '\n'
}
