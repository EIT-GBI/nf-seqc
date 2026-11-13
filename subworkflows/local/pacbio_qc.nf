//
// Subworkflow: PACBIO_QC
//
// Revio basecalls on-instrument, so HiFi reads arrive as unaligned BAM and are
// delivered as-is (no FASTQ conversion requested). The instrument also writes
// every QC number this report needs, per SMRT Cell and per barcode, so no task
// reads the BAM: QC is the instrument's own reports, read here.
//
// A Revio run folder holds one directory per SMRT Cell:
//
//   r84334_<date>_<time>/<plate>_<well>/        e.g. 1_A01
//       hifi_reads/<movie>.hifi_reads[.<barcode>|.unassigned].bam
//       fail_reads/<movie>.fail_reads[...].bam
//       pb_formats/<movie>.hifi_reads.consensusreadset.xml   sample names
//       statistics/<movie>.summary.json, .ccs_report.json   per cell
//       statistics/<movie>.reports.zip                      loading, and
//                                     barcodes.report.json: per barcode
//       statistics/<movie>.{hifi,fail}_reads.lima_counts.txt
//       metadata/<movie>.metadata.xml                       software versions
//       metadata/<movie>.transferdone                       written last
//

include { demuxSummary; readShareTable; demuxSummaryTable } from './demux_report.nf'

// <movie>.<kind>[.<barcode>].bam, kind hifi_reads or fail_reads; barcode is
// e.g. bc2053, bc1001--bc1002, or unassigned. null for any other name.
def revioBamName( bam ) {
    def m = bam.name =~ /^(.+?)\.(hifi_reads|fail_reads)(?:\.(.+))?\.bam$/
    return m.matches() ? [ movie: m.group(1), kind: m.group(2), barcode: m.group(3) ] : null
}

// The SMRT Cell directory a BAM belongs to: <cell>/hifi_reads/x.bam. null for
// a BAM outside the Revio layout (e.g. a flat --bam_dir).
def revioCell( bam ) {
    return bam.parent.name in [ 'hifi_reads', 'fail_reads' ] ? bam.parent.parent : null
}

// The consensusreadset.xml names barcodes bc2053--bc2053; the BAM names a
// symmetric pair bc2053. Key both by the BAM's form.
def barcodeKey( String name ) {
    def parts = name.split('--')
    return parts.size() == 2 && parts[0] == parts[1] ? parts[0] : name
}

def localName( node ) {
    return node.name().tokenize(':')[-1]
}

// The strict syntax has no GPath `.@attr`, so read attributes from the map.
def xmlName( node ) {
    return node.attributes()['Name']?.toString()
}

// Sample names from a cell's <movie>.hifi_reads.consensusreadset.xml: each
// BioSample's name by its barcode, and the well sample name, which is all an
// unbarcoded cell has. Empty when the file is missing.
def revioSampleNames( cell, String movie ) {
    def xml = cell ? cell.resolve("pb_formats/${movie}.hifi_reads.consensusreadset.xml") : null
    if ( !xml?.exists() ) {
        return [ barcodes: [:], well_sample: null ]
    }
    def doc     = new groovy.xml.XmlSlurper( false, false ).parse( xml.toFile() )
    def samples = [:]
    doc.depthFirst().findAll { n -> localName( n ) == 'BioSample' }.each { bs ->
        bs.depthFirst().findAll { n -> localName( n ) == 'DNABarcode' }.each { bc ->
            samples[ barcodeKey( xmlName( bc ) ) ] = xmlName( bs )
        }
    }
    def well = doc.depthFirst().find { n -> localName( n ) == 'WellSample' }
    return [ barcodes: samples, well_sample: well ? xmlName( well ) : null ]
}

// Sample ids are MultiQC row names; keep them to one safe token.
def safeId( String name ) {
    return name.replaceAll(/[^A-Za-z0-9._-]+/, '_')
}

// [ meta, bam ] per HiFi BAM, meta.id the BioSample name (barcoded cell), the
// well sample name (unbarcoded cell) or the BAM name (no
// consensusreadset.xml). Only HiFi BAMs come here (failed reads are counted
// per barcode from lima, not delivered). A name used by more than one cell (a
// library re-sequenced) is suffixed with the well, so General Statistics keeps
// one row per cell.
def revioSamples( meta, List bams ) {
    def names = [:]
    def rows = bams.collect { bam ->
        def b    = revioBamName( bam )
        def cell = revioCell( bam )
        if ( !b ) {
            return [ meta: meta + [ id: safeId( bam.name - ~/\.bam$/ ), well: null, movie: null, barcode: null ], bam: bam ]
        }
        def key = "${cell}/${b.movie}".toString()
        if ( !names.containsKey( key ) ) {
            names[key] = revioSampleNames( cell, b.movie )
        }
        def n    = names[key]
        def name = b.barcode
            ? ( n.barcodes[ b.barcode ] ?: b.barcode )
            : ( n.well_sample ?: b.movie )
        def well = cell ? cell.name : b.movie
        [ meta: meta + [ id: safeId( name ), well: well, movie: b.movie, barcode: b.barcode ], bam: bam ]
    }
    def counts = rows.countBy { r -> r.meta.id }
    return rows.collect { r ->
        def id = counts[ r.meta.id ] > 1 ? "${r.meta.id}_${r.meta.well}" : r.meta.id
        tuple( r.meta + [ id: id ], r.bam )
    }
}

// Cells whose transfer has not finished: Revio writes
// metadata/<movie>.transferdone once a cell is fully copied off the
// instrument. Checked from the entry workflow, before any task runs.
def revioUnfinishedCells( root ) {
    def cells = files("${root}/*/metadata", type: 'dir') + files("${root}/metadata", type: 'dir')
    return cells
        .findAll { m -> m.exists() && !files("${m}/*.transferdone") }
        .collect { m -> m.parent.name }
        .unique()
}

// HiFi reads per sample and the Not Barcoded count, from each barcoded cell's
// lima_counts.txt, as demuxSummary() takes them. null when no cell was
// barcoded. `samples` is revioSamples()' output, which maps each cell's
// barcode to its sample id.
def limaReads( List samples ) {
    def by_cell = samples.findAll { meta, _bam -> meta.movie && meta.barcode && meta.barcode != 'unassigned' }
        .groupBy { meta, bam -> [ revioCell( bam ), meta.movie ] }
    def reads      = [:]
    def unassigned = 0L
    by_cell.each { cell_movie, rows ->
        def counts = cell_movie[0]?.resolve("statistics/${cell_movie[1]}.hifi_reads.lima_counts.txt")
        if ( !counts?.exists() ) {
            return
        }
        def ids = rows.collectEntries { meta, _bam -> [ ( meta.barcode ): meta.id ] }
        counts.splitCsv( header: true, sep: '\t' ).each { row ->
            def n = row.Counts as long
            if ( row.IdxCombinedNamed == 'Not Barcoded' ) {
                unassigned += n
            }
            else if ( ids[ barcodeKey( row.IdxCombinedNamed ) ] ) {
                def id = ids[ barcodeKey( row.IdxCombinedNamed ) ]
                reads[id] = ( reads[id] ?: 0L ) + n
            }
        }
    }
    return reads ? demuxSummary( reads, unassigned, [] ) : null
}

// pbreports attributes as [ id: value ], from a report JSON's text.
def pbAttributes( String json ) {
    def report = new groovy.json.JsonSlurper().parseText( json )
    return ( report.attributes ?: [] ).collectEntries { a -> [ ( a.id ): a.value ] }
}

// The text of one entry of a cell's reports.zip, or null. No try/finally in
// the strict syntax; withCloseable closes the zip.
def zipEntryText( zip, String name ) {
    if ( !zip?.exists() ) {
        return null
    }
    return new java.util.zip.ZipFile( zip.toFile() ).withCloseable { z ->
        def entry = z.getEntry( name )
        entry ? z.getInputStream( entry ).text : null
    }
}

// pbreports' per-barcode table (barcodes.report.json in reports.zip) as
// [ barcode: [ column: value ] ], keyed by the BAM's barcode form (bc2053).
// Columns are the ids' last part: number_of_reads, number_of_bases,
// mean_read_length, median_read_quality ('Q55'), mean_bcqual. Empty when the
// cell has no such report (unbarcoded, or an older SMRT Link).
def revioBarcodeTable( zip ) {
    def text = zipEntryText( zip, 'barcodes.report.json' )
    def table = text ? new groovy.json.JsonSlurper().parseText( text ).tables?.find { t -> t.id == 'barcode.barcode_table' } : null
    if ( !table ) {
        return [:]
    }
    def cols = table.columns.collectEntries { c -> [ ( c.id.tokenize('.')[-1] ): c.values ] }
    def rows = [:]
    cols.barcode?.eachWithIndex { bc, i ->
        rows[ barcodeKey( bc.toString() ) ] = cols.collectEntries { k, v -> [ ( k ): v[i] ] }
    }
    return rows
}

// Read counts per barcode from a lima_counts.txt, [ bc2053: n ]. Empty when
// the file is missing.
def limaCounts( counts ) {
    if ( !counts?.exists() ) {
        return [:]
    }
    return counts.splitCsv( header: true, sep: '\t' )
        .findAll { row -> row.IdxCombinedNamed != 'Not Barcoded' }
        .collectEntries { row -> [ ( barcodeKey( row.IdxCombinedNamed ) ): row.Counts as long ] }
}

// 'Q55' -> 55; numbers pass through; '' when absent.
def qvNumber( value ) {
    if ( value == null ) {
        return ''
    }
    def m = value.toString() =~ /^Q?(\d+(?:\.\d+)?)$/
    return m.matches() ? m.group(1) as BigDecimal : ''
}

// One sample's HiFi metrics, every value '' when its source is missing.
// Barcoded: the cell's barcodes.report.json row (falling back to lima_counts
// for the read count). Unbarcoded: the cell's ccs_report.json Q20 figures,
// since the whole cell is the sample. Failed reads per barcode come from
// fail_reads.lima_counts.txt; an unbarcoded cell has no such count.
def revioSampleRow( Map meta, Map cell_data ) {
    def bc = meta.barcode
    if ( bc ) {
        def r = cell_data.barcodes[ bc ] ?: [:]
        return [
            hifi_reads:  r.number_of_reads ?: cell_data.hifi_counts[ bc ] ?: '',
            hifi_yield:  r.number_of_bases ? r.number_of_bases / 1e6 : '',
            mean_length: r.mean_read_length ?: '',
            median_qv:   qvNumber( r.median_read_quality ),
            bc_quality:  r.mean_bcqual ?: '',
            fail_reads:  cell_data.fail_counts.containsKey( bc ) ? cell_data.fail_counts[ bc ] : '',
        ]
    }
    def c = cell_data.ccs
    def bases = c['ccs_processing.total_number_of_ccs_bases_q20']
    return [
        hifi_reads:  c['ccs_processing.number_of_ccs_reads_q20'] ?: '',
        hifi_yield:  bases ? bases / 1e6 : '',
        mean_length: c['ccs_processing.mean_ccs_readlength_q20'] ?: '',
        median_qv:   qvNumber( c['ccs_processing.median_qv_q20'] ),
        bc_quality:  '',
        fail_reads:  '',
    ]
}

// MultiQC General Statistics columns, one row per HiFi sample, from the
// instrument's reports: what NanoPlot used to compute by reading every BAM.
// null when there is nothing to show (no cell in the Revio layout).
def revioSampleTable( List samples ) {
    def cells = [:]
    def rows = samples.findAll { meta, bam -> meta.movie && revioCell( bam ) }.collect { meta, bam ->
        def cell = revioCell( bam )
        def key  = "${cell}/${meta.movie}".toString()
        if ( !cells.containsKey( key ) ) {
            def stats = cell.resolve('statistics')
            def ccs   = stats.resolve("${meta.movie}.ccs_report.json")
            cells[key] = [
                barcodes:    revioBarcodeTable( stats.resolve("${meta.movie}.reports.zip") ),
                hifi_counts: limaCounts( stats.resolve("${meta.movie}.hifi_reads.lima_counts.txt") ),
                fail_counts: limaCounts( stats.resolve("${meta.movie}.fail_reads.lima_counts.txt") ),
                ccs:         ccs.exists() ? pbAttributes( ccs.text ) : [:],
            ]
        }
        [ id: meta.id ] + revioSampleRow( meta, cells[key] )
    }
    if ( !rows ) {
        return null
    }
    def header = [
        "# id: 'pacbio_sample_metrics'",
        "# plot_type: 'generalstats'",
        "# headers:",
        "#     hifi_reads:",
        "#         title: 'HiFi reads'",
        "#         namespace: 'PacBio'",
        "#         description: 'HiFi (Q20+) reads, from the instrument''s barcode report (unbarcoded cell: its CCS report)'",
        "#         format: '{:,.0f}'",
        "#     hifi_yield:",
        "#         title: 'HiFi yield'",
        "#         namespace: 'PacBio'",
        "#         description: 'HiFi bases'",
        "#         suffix: ' Mb'",
        "#         format: '{:,.1f}'",
        "#     mean_length:",
        "#         title: 'Mean length'",
        "#         namespace: 'PacBio'",
        "#         description: 'Mean HiFi read length'",
        "#         suffix: ' bp'",
        "#         format: '{:,.0f}'",
        "#     median_qv:",
        "#         title: 'Median QV'",
        "#         namespace: 'PacBio'",
        "#         description: 'Median HiFi read quality'",
        "#         format: '{:,.0f}'",
        "#     bc_quality:",
        "#         title: 'Barcode quality'",
        "#         namespace: 'PacBio'",
        "#         description: 'Mean lima barcode score (0-100)'",
        "#         min: 0",
        "#         max: 100",
        "#         format: '{:,.1f}'",
        "#     fail_reads:",
        "#         title: 'Failed reads'",
        "#         namespace: 'PacBio'",
        "#         description: 'Reads below Q20 assigned to this barcode (fail_reads lima counts)'",
        "#         format: '{:,.0f}'",
        "#         hidden: true",
        "Sample\thifi_reads\thifi_yield\tmean_length\tmedian_qv\tbc_quality\tfail_reads",
    ]
    def lines = rows.collect { r -> r.values().join('\t') }
    return ( header + lines ).join('\n') + '\n'
}

// The instrument software that produced a cell's reads, from
// metadata/<movie>.metadata.xml, as [ process, tool, version ] for the
// Software Versions section. Revio ICS is the instrument release; SMRT Link
// is the release of the on-instrument ccs/lima/pbreports that wrote every
// number in this report. Empty when the file is missing.
def revioSoftwareVersions( List samples ) {
    def wanted = [ ics: 'Revio ICS', smrtlink: 'SMRT Link' ]
    return samples.findAll { meta, bam -> meta.movie && revioCell( bam ) }
        .collect { meta, bam -> revioCell( bam ).resolve("metadata/${meta.movie}.metadata.xml") }
        .unique()
        .findAll { xml -> xml.exists() }
        .collectMany { xml ->
            def doc = new groovy.xml.XmlSlurper( false, false ).parse( xml.toFile() )
            doc.depthFirst().findAll { n -> localName( n ) == 'VersionInfo' }.collect { n ->
                def a = n.attributes()
                wanted[ a['Name']?.toString() ] && a['Version'] ? tuple( 'PACBIO_QC', wanted[ a['Name'].toString() ], a['Version'].toString() ) : null
            }.findAll { v -> v != null }
        }
        .unique()
}

// A productivity count from the loading report as a share of all ZMWs, or ''
// (an empty table cell) when the report lacks it.
def zmwShare( Map loading, String key ) {
    def zmws = loading['loading_xml_report.productive_zmws']
    return zmws && loading[key] != null ? 100.0 * loading[key] / zmws : ''
}

// One SMRT Cell's metrics, each n/a-able on its own, from summary.json,
// ccs_report.json and the loading report inside reports.zip. Loading is
// P0 (empty), P1 (productive) and P2 (other) as a share of all ZMWs.
def revioCellMetrics( cell, String movie, String well_sample ) {
    def stats   = cell.resolve('statistics')
    def summary = stats.resolve("${movie}.summary.json")
    def ccs     = stats.resolve("${movie}.ccs_report.json")
    def zip     = stats.resolve("${movie}.reports.zip")

    def s = summary.exists() ? new groovy.json.JsonSlurper().parse( summary.toFile() ) : [:]
    def c = ccs.exists() ? pbAttributes( ccs.text ) : [:]
    def loading = zipEntryText( zip, 'loading.report.json' )
    def l = loading ? pbAttributes( loading ) : [:]
    // N50 from pbreports' ccs.report.json (what SMRT Link shows), never from
    // ccs_report.json: on SMRT Link 25.3 its ccs_readlength_n50_q20 is a copy
    // of the median (2572 and 2572 on r84334_20260817, where pbreports and
    // SMRT Link give 7544).
    def pbccs = zipEntryText( zip, 'ccs.report.json' )
    def p = pbccs ? pbAttributes( pbccs ) : [:]
    def bases = c['ccs_processing.total_number_of_ccs_bases_q20']
    return [
        well_sample: well_sample ?: '',
        movie:       movie,
        p0:          zmwShare( l, 'loading_xml_report.productivity_0_n' ),
        p1:          zmwShare( l, 'loading_xml_report.productivity_1_n' ),
        p2:          zmwShare( l, 'loading_xml_report.productivity_2_n' ),
        hifi_reads:  s.number_hifi_reads ?: c['ccs_processing.number_of_ccs_reads_q20'] ?: '',
        hifi_yield:  bases ? bases / 1e9 : '',
        mean_length: s.hifi_mean_read_length ?: '',
        n50:         p['ccs2.ccs_readlength_n50'] ?: '',
        mean_qv:     s.hifi_mean_qv ?: '',
        q30:         s.q30_base_percentage ?: '',
    ]
}

// One row per SMRT Cell, as MultiQC custom content. The run-level slot that
// the InterOp table fills for Illumina. No loading flag yet: P1 is shown, and
// a threshold waits for the lab to say what under-loaded means on our Revio.
def revioMetricsTable( List samples ) {
    def cells = samples.findAll { meta, bam -> meta.movie && revioCell( bam ) }
        .collect { meta, bam -> [ well: meta.well, cell: revioCell( bam ), movie: meta.movie ] }
        .unique { c -> [ c.well, c.movie ] }
        .sort { c -> c.well }
    if ( !cells ) {
        return null
    }
    def header = [
        "# id: 'pacbio_run_metrics'",
        "# section_name: 'PacBio Run Metrics'",
        "# description: 'One row per SMRT Cell, from the instrument reports (summary.json, ccs_report.json, reports.zip). P0 empty, P1 productive and P2 other ZMWs, as a share of all ZMWs; HiFi reads are Q20+.'",
        "# plot_type: 'table'",
        "# pconfig:",
        "#     id: 'pacbio_run_metrics_table'",
        "#     namespace: 'PacBio'",
        "#     col1_header: 'SMRT Cell'",
        "# headers:",
        "#     well_sample:",
        "#         title: 'Well sample'",
        "#         scale: false",
        "#     movie:",
        "#         title: 'Movie'",
        "#         scale: false",
        "#         hidden: true",
        "#     p0:",
        "#         title: 'P0'",
        "#         description: 'Empty ZMWs'",
        "#         suffix: '%'",
        "#         min: 0",
        "#         max: 100",
        "#         format: '{:,.1f}'",
        "#     p1:",
        "#         title: 'P1'",
        "#         description: 'Productive ZMWs: one polymerase, usable reads'",
        "#         suffix: '%'",
        "#         min: 0",
        "#         max: 100",
        "#         format: '{:,.1f}'",
        "#     p2:",
        "#         title: 'P2'",
        "#         description: 'Other ZMWs, e.g. more than one polymerase'",
        "#         suffix: '%'",
        "#         min: 0",
        "#         max: 100",
        "#         format: '{:,.1f}'",
        "#     hifi_reads:",
        "#         title: 'HiFi reads'",
        "#         format: '{:,.0f}'",
        "#     hifi_yield:",
        "#         title: 'HiFi yield'",
        "#         suffix: ' Gb'",
        "#         format: '{:,.2f}'",
        "#     mean_length:",
        "#         title: 'Mean length'",
        "#         suffix: ' bp'",
        "#         format: '{:,.0f}'",
        "#     n50:",
        "#         title: 'Read N50'",
        "#         description: 'HiFi read length N50'",
        "#         suffix: ' bp'",
        "#         format: '{:,.0f}'",
        "#     mean_qv:",
        "#         title: 'Mean QV'",
        "#         format: '{:,.1f}'",
        "#     q30:",
        "#         title: 'Q30 bases'",
        "#         suffix: '%'",
        "#         min: 0",
        "#         max: 100",
        "#         format: '{:,.1f}'",
        "Cell\twell_sample\tmovie\tp0\tp1\tp2\thifi_reads\thifi_yield\tmean_length\tn50\tmean_qv\tq30",
    ]
    def lines = cells.collect { c ->
        def m = revioCellMetrics( c.cell, c.movie, revioSampleNames( c.cell, c.movie ).well_sample )
        ( [ c.well ] + m.values() ).join('\t')
    }
    return ( header + lines ).join('\n') + '\n'
}

workflow PACBIO_QC {

    take:
    ch_run           // channel: [ val(meta), path(run_dir) ]
    bam_dir          // file or null: --bam_dir, where the HiFi BAM live if not in the run dir
    max_undetermined // number: --max_undetermined, % of reads Unassigned

    main:
    // One BAM per SMRT Cell, or per barcode on a barcoded cell. Collected as
    // one list first, so sample names can be checked for clashes across cells.
    // The BAM are only listed, never read: they name the samples, and each
    // sample's numbers come from the instrument's reports. Unassigned reads
    // are counted in the Demultiplexing table, failed reads per barcode in a
    // (hidden) General Statistics column.
    ch_samples = ch_run.map { meta, run_dir ->
        def root = bam_dir ?: run_dir
        def bams = ( files("${root}/**/*.hifi_reads*.bam") + files("${root}/*.hifi_reads*.bam") )
            .findAll { f -> f.exists() && !f.name.contains('.unassigned.') }
            .unique()
        if ( !bams ) {
            error "PACBIO_QC: no *.hifi_reads*.bam found under ${root}. Point --bam_dir at the SMRT Cell output if it lives elsewhere."
        }
        revioSamples( meta, bams.sort() )
    }

    ch_bam = ch_samples.flatMap { samples -> samples }

    // Reads lima could not assign to a barcode: delivered with the samples as
    // unassigned_<well>.hifi_reads.bam, but not a sample and not QC'd.
    ch_unassigned = ch_run.flatMap { meta, run_dir ->
        def root = bam_dir ?: run_dir
        ( files("${root}/**/*.hifi_reads.unassigned.bam") + files("${root}/*.hifi_reads.unassigned.bam") )
            .findAll { f -> f.exists() }
            .unique()
            .sort()
            .collect { bam -> tuple( meta + [ id: "unassigned_${revioCell( bam )?.name ?: revioBamName( bam )?.movie ?: 'cell'}".toString(), unassigned: true ], bam ) }
    }

    // flatMap, not map + filter: ch_samples is a value channel (see ILLUMINA_QC).
    ch_demux = ch_samples.flatMap { samples ->
        def demux = limaReads( samples )
        demux ? [ tuple( samples[0][0], demux ) ] : []
    }
    ch_read_share = ch_demux
        .map { _meta, demux -> readShareTable( demux, 'Unassigned' ) }
        .collectFile( name: 'demux_read_share_mqc.tsv' )
    ch_demux_summary = ch_demux
        .map { meta, demux -> demuxSummaryTable( meta.run, demux, max_undetermined, 'lima_counts.txt', 'Unassigned', 'PACBIO_QC' ) }
        .collectFile( name: 'demux_summary_mqc.tsv' )

    ch_run_metrics = ch_samples
        .flatMap { samples ->
            def table = revioMetricsTable( samples )
            table ? [ table ] : []
        }
        .collectFile( name: 'pacbio_run_metrics_mqc.tsv' )

    ch_sample_metrics = ch_samples
        .flatMap { samples ->
            def table = revioSampleTable( samples )
            table ? [ table ] : []
        }
        .collectFile( name: 'pacbio_sample_metrics_mqc.tsv' )

    ch_qc = channel.empty()
        .mix( ch_sample_metrics )
        .mix( ch_run_metrics )
        .mix( ch_read_share )
        .mix( ch_demux_summary )

    ch_versions = ch_samples.flatMap { samples -> revioSoftwareVersions( samples ) }

    emit:
    reads    = ch_bam.mix( ch_unassigned ) // HiFi uBAM delivered unchanged, unassigned included
    qc       = ch_qc
    versions = ch_versions
}
