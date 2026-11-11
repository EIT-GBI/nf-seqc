//
// Subworkflow: NANOPORE_QC
//
// A PromethION run folder holds one acquisition directory per flow cell (and
// per restart), each written by MinKNOW:
//
//   <experiment>/<sample>/<date>_<time>_<device>_<flowcell>_<run>/
//       final_summary_*.txt      written last: flow cell, sample, kit (protocol=)
//       sample_sheet_*.csv       barcode aliases, when the run had a sample sheet
//       pod5_pass/barcodeNN/     POD5 MinKNOW classified to one barcode
//       pod5_pass/mixed/         POD5 holding reads of several barcodes
//       pod5_pass/unclassified/  no barcode found
//       pod5_fail/ ...           failed reads (excluded unless --include_failed)
//
// MinKNOW's per-barcode folders are not a complete split: mixed/ holds a large
// share of each barcode's reads (34% of barcode01's on 20260515_Plasmids_A_B).
// So every POD5 folder is basecalled with --kit-name, which tags each read with
// its barcode, and DORADO_DEMUX then splits all of a flow cell's reads by that
// tag. Unbarcoded flow cells skip the demux: each folder's BAM is the sample.
//

include { DORADO_BASECALLER } from '../../modules/dorado/basecaller/main.nf'
include { DORADO_DEMUX      } from '../../modules/dorado/demux/main.nf'
include { SAMTOOLS_FASTQ    } from '../../modules/samtools/fastq/main.nf'
include { NANOPLOT          } from '../../modules/local/nanoplot/main.nf'
include { demuxSummary; readShareTable; demuxSummaryTable } from './demux_report.nf'

// key=value lines of a MinKNOW final_summary_*.txt.
def keyValues( summary ) {
    return summary.readLines()
        .findAll { l -> l.contains('=') }
        .collectEntries { l -> [ ( l.substring( 0, l.indexOf('=') ) ): l.substring( l.indexOf('=') + 1 ) ] }
}

// The barcoding kit from final_summary's protocol=, e.g.
// sequencing/sequencing_MIN114_DNA_e8_2_400K:FLO-MIN114:SQK-RBK114-24:400.
// An expansion kit (EXP-NBD114 with SQK-LSK114) is the one carrying the
// barcodes, so it wins over the sequencing kit.
def protocolKit( String protocol ) {
    def kits = ( protocol ?: '' ).tokenize(':').findAll { t -> t ==~ /^(SQK|EXP|VSK)-.+/ }
    return kits.find { k -> k.startsWith('EXP-') } ?: kits ? kits[-1] : null
}

// Barcode aliases from the acquisition's sample_sheet_*.csv, when it has
// barcode and alias columns (MinKNOW's sample sheet without them only names
// the flow cell and kit).
def sampleSheetAliases( dir ) {
    def sheets = files("${dir}/sample_sheet_*.csv").findAll { f -> f.exists() }
    if ( !sheets ) {
        return [:]
    }
    def rows = sheets[0].splitCsv( header: true )
    return rows.findAll { r -> r.barcode && r.alias }.collectEntries { r -> [ ( r.barcode ): r.alias ] }
}

// One record per POD5 folder to basecall: [ meta, folder ]. meta.flowcell,
// meta.kit and meta.aliases come from the acquisition's final_summary and
// sample sheet. A POD5 tree without final_summary (--pod5_dir pointing at bare
// POD5) is one acquisition named after its directory, and needs --barcode_kit
// if barcoded.
def ontPod5Folders( meta, root, boolean include_failed, barcode_kit ) {
    def pod5s = ( files("${root}/**/*.pod5") + files("${root}/*.pod5") )
        .findAll { f -> f.exists() }
        .unique()
    def folders = pod5s.collect { f -> f.parent }.unique()
    def records = folders.collect { dir ->
        // pod5_pass/barcode01/x.pod5 or pod5_pass/x.pod5 (unbarcoded)
        def status_dir = dir.name.startsWith('pod5') ? dir : dir.parent
        def acq        = status_dir.name.startsWith('pod5') ? status_dir.parent : dir
        def status     = status_dir.name.startsWith('pod5') ? status_dir.name : 'pod5'
        def summary    = files("${acq}/final_summary_*.txt").find { f -> f.exists() }
        def kv         = summary ? keyValues( summary ) : [:]
        def barcoded   = dir.name ==~ /barcode\d+|mixed|unclassified/
        def flowcell   = kv.flow_cell_id ?: acq.name
        def run8       = ( kv.protocol_run_id ?: acq.name ).take(8)
        [
            meta: meta + [
                id:       "${flowcell}_${run8}_${status}_${dir.name}".toString(),
                flowcell: flowcell,
                sample:   kv.sample_id ?: flowcell,
                kit:      barcoded ? ( barcode_kit ?: protocolKit( kv.protocol ) ) : null,
                aliases:  sampleSheetAliases( acq ),
                barcoded: barcoded,
                status:   status,
            ],
            dir: dir,
        ]
    }
    return records
        .findAll { r -> include_failed || r.meta.status != 'pod5_fail' }
        .collect { r -> tuple( r.meta, r.dir ) }
        .sort { t -> t[0].id }
}

// Acquisition directories that have POD5 but no final_summary_*.txt yet:
// MinKNOW writes final_summary last. Checked from the entry workflow, before
// any task runs.
def nanoporeUnfinished( root ) {
    def status_dirs = ( files("${root}/**/pod5*", type: 'dir') + files("${root}/pod5*", type: 'dir') )
        .findAll { d -> d.exists() && files("${d}/**/*.pod5") + files("${d}/*.pod5") }
    return status_dirs
        .collect { d -> d.parent }
        .unique()
        .findAll { acq -> !files("${acq}/final_summary_*.txt").find { f -> f.exists() } }
        .collect { acq -> acq.name }
}

// barcode01 from dorado demux's <kit>_barcode01.bam, or unclassified.
def demuxBarcode( bam ) {
    def m = bam.name =~ /(barcode\d+|unclassified)\.bam$/
    return m.find() ? m.group(1) : bam.name - ~/\.bam$/
}

// Sample ids end up in file and directory names (NanoPlot --prefix/--outdir).
def safeId( String name ) {
    return name.replaceAll(/[^A-Za-z0-9._-]+/, '_')
}

// Per-sample [ meta, bam ] from every flow cell's demux output, named by the
// sample sheet alias, else the barcode. A name used on two flow cells (barcode01
// on each) is suffixed with the flow cell id.
def ontSamples( List demuxed ) {
    def rows = demuxed.collectMany { meta, bams ->
        ( bams instanceof List ? bams : [ bams ] ).collect { bam ->
            def barcode = demuxBarcode( bam )
            def name    = meta.aliases[ barcode ] ?: barcode
            [ meta: meta + [ id: safeId( name ), barcode: barcode ], bam: bam ]
        }
    }
    def counts = rows.findAll { r -> r.meta.barcode != 'unclassified' }.countBy { r -> r.meta.id }
    return rows.collect { r ->
        def id = r.meta.barcode == 'unclassified' || counts[ r.meta.id ] > 1 ? "${r.meta.id}_${r.meta.flowcell}" : r.meta.id
        tuple( r.meta + [ id: id.toString() ], r.bam )
    }
}

// Reads per sample and unclassified, from each flow cell's
// barcoding_summary.txt (one row per read, streamed: a PromethION flow cell
// has millions), as [ meta, reads: [ id: n ], unclassified: n ]. reads is
// empty when there is no summary to count from (no barcoded flow cell).
def barcodeCounts( List samples, List summaries ) {
    def ids        = samples.collectEntries { meta, _bam -> [ ( "${meta.flowcell}/${meta.barcode}".toString() ): meta.id ] }
    def reads      = [:]
    def unassigned = 0L
    summaries.each { meta, summary ->
        def column = -1
        summary.eachLine { line, n ->
            def fields = line.split('\t')
            if ( n == 1 ) {
                column = fields.findIndexOf { f -> f == 'barcode' }
                return
            }
            if ( column < 0 || column >= fields.size() ) {
                return
            }
            def m       = fields[column] =~ /(barcode\d+|unclassified)$/
            def barcode = m.find() ? m.group(1) : fields[column]
            if ( barcode == 'unclassified' ) {
                unassigned += 1
            }
            else {
                def id = ids[ "${meta.flowcell}/${barcode}".toString() ]
                if ( id ) {
                    reads[id] = ( reads[id] ?: 0L ) + 1
                }
            }
        }
    }
    return [ meta: samples ? samples[0][0] : [:], reads: reads, unclassified: unassigned ]
}

// A barcode with fewer than `min_reads` reads is not a sample: on a real run
// a dozen unused barcodes each pick up a read or two from misclassification.
// Their reads count as unclassified in the Demultiplexing table. null when no
// barcode has enough reads.
def demuxAboveCutoff( Map counts, long min_reads ) {
    def low  = counts.reads.findAll { _id, n -> n < min_reads }
    def kept = counts.reads.findAll { _id, n -> n >= min_reads }
    if ( low ) {
        log.info "NANOPORE_QC: ${low.size()} barcode(s) below --min_barcode_reads ${min_reads}, counted as unclassified, not reported as samples: ${low.collect { id, n -> "${id} (${n})" }.join(', ')}"
    }
    return kept ? demuxSummary( kept, counts.unclassified + ( low.values().sum() ?: 0L ), [] ) : null
}

// A JSON number MinKNOW wrote as a string ("328774"), or null.
def asLong( value ) {
    return value != null && value.toString().isNumber() ? value.toString() as long : null
}

// Hours between two MinKNOW timestamps (RFC 3339, nanoseconds), or null.
def hoursBetween( String start, String end ) {
    if ( !start || !end ) {
        return null
    }
    def seconds = java.time.Duration.between( java.time.Instant.parse( start ), java.time.Instant.parse( end ) ).seconds
    return seconds / 3600.0
}

// MinKNOW's STOPPED_PROTOCOL_ENDED as "protocol ended".
def stopReason( String reason ) {
    return reason ? ( reason - ~/^STOPPED_/ ).toLowerCase().replace('_', ' ') : ''
}

// One flow cell acquisition's run metrics, measured live by MinKNOW, from its
// report_*.json. The report holds calibration acquisitions first; the metrics
// come from the last SEQUENCING one. Returns null when there is none (a run
// stopped during the flow cell check). Every field is n/a-able on its own.
def ontReportMetrics( report ) {
    def d   = new groovy.json.JsonSlurper().parse( report.toFile() )
    def acq = ( d.acquisitions ?: [] ).findAll { a -> a.acquisition_run_info?.config_summary?.purpose == 'SEQUENCING' }
    if ( !acq ) {
        return null
    }
    def a        = acq[-1]
    def info     = a.acquisition_run_info ?: [:]
    def yields   = info.yield_summary ?: [:]
    def run      = d.protocol_run_info ?: [:]
    def flowcell = run.flow_cell?.flow_cell_id ?: run.user_info?.user_specified_flow_cell_id ?: report.parent.name
    // N50 of the bases MinKNOW basecalled live, else of its estimate from the
    // raw signal. The pipeline's own (sup) N50 per sample is in NanoStat.
    def lengths  = ( a.read_length_histogram ?: [] ).collectEntries { h -> [ ( h.plot?.read_length_type ?: h.read_length_type ?: 'raw' ): h.plot?.histogram_data?.getAt(0)?.n50 ] }
    // Pore scans run every 1.5 h; the first is at the start of sequencing.
    def scans    = ( d.user_messages ?: [] ).findAll { m -> m.identifier == 'mux_scan_result' && m.extra_data }
    def occ      = ( a.duty_time?.getAt(0)?.pore_occupancy ?: [] ).findAll { v -> v != null }
    def passed   = asLong( yields.basecalled_pass_bases )
    def estimate = asLong( yields.estimated_selected_bases )
    return [
        flowcell:    flowcell,
        protocol:    ( run.run_id ?: '' ).take(8),
        sample:      run.user_info?.sample_id ?: '',
        position:    run.device?.device_id ?: '',
        flowcell_type: run.flow_cell?.product_code ?: '',
        kit:         run.user_info?.kit_info?.sequencing_kit ?: '',
        hours:       hoursBetween( info.data_read_start_time, info.data_read_end_time ) ?: '',
        stop:        stopReason( info.stop_reason ),
        pores_start: scans ? ( asLong( scans[0].extra_data.total_pores ) ?: '' ) : '',
        pores_ready: scans ? ( asLong( scans[0].extra_data.num_pores ) ?: '' ) : '',
        pores_end:   scans ? ( asLong( scans[-1].extra_data.total_pores ) ?: '' ) : '',
        occupancy:   occ ? 100.0 * occ.sum() / occ.size() : '',
        reads:       asLong( yields.read_count ) ?: '',
        pass_bases:  passed != null ? passed / 1e9 : '',
        est_bases:   estimate != null ? estimate / 1e9 : '',
        n50:         lengths.BasecalledBases ?: lengths.EstimatedBases ?: '',
        model:       a.basecaller_config?.model_names?.simplex_model ?: '',
        minknow:     run.software_versions?.minknow?.full ?: '',
    ]
}

// One row per flow cell acquisition, as MultiQC custom content: the run-level
// slot the InterOp table fills for Illumina. A flow cell restarted into a
// second acquisition gets a row per acquisition, suffixed with the protocol
// run id. No flags yet: what counts as low occupancy or pore count waits for
// the lab, as for PacBio P1.
def ontMetricsTable( root ) {
    def reports = ( files("${root}/**/report_*.json") + files("${root}/report_*.json") )
        .findAll { f -> f.exists() }
        .unique()
        .sort()
    def rows = reports.collect { r -> ontReportMetrics( r ) }.findAll { m -> m }
    if ( !rows ) {
        return null
    }
    def counts = rows.countBy { m -> m.flowcell }
    def header = [
        "# id: 'nanopore_run_metrics'",
        "# section_name: 'Nanopore Run Metrics'",
        "# description: 'One row per flow cell, from MinKNOW''s report_*.json: what MinKNOW measured live during sequencing. Pores are from the pore scans MinKNOW runs every 1.5 h (first and last scan); occupancy is the share of available pores sequencing, averaged over the run. Yield and N50 are MinKNOW''s own live basecalls; the pipeline''s re-basecalled reads per sample are under NanoStat.'",
        "# plot_type: 'table'",
        "# pconfig:",
        "#     id: 'nanopore_run_metrics_table'",
        "#     namespace: 'Nanopore'",
        "#     col1_header: 'Flow cell'",
        "# headers:",
        "#     sample:",
        "#         title: 'Sample'",
        "#         description: 'MinKNOW sample id'",
        "#         scale: false",
        "#     position:",
        "#         title: 'Position'",
        "#         scale: false",
        "#     flowcell_type:",
        "#         title: 'Flow cell type'",
        "#         scale: false",
        "#         hidden: true",
        "#     kit:",
        "#         title: 'Kit'",
        "#         scale: false",
        "#     hours:",
        "#         title: 'Run time'",
        "#         suffix: ' h'",
        "#         format: '{:,.1f}'",
        "#         scale: false",
        "#     stop:",
        "#         title: 'Stopped'",
        "#         description: 'Why MinKNOW stopped acquiring'",
        "#         scale: false",
        "#     pores_start:",
        "#         title: 'Pores, start'",
        "#         description: 'Pores found by the first pore scan'",
        "#         format: '{:,.0f}'",
        "#     pores_ready:",
        "#         title: 'Ready, start'",
        "#         description: 'Pores available for immediate sequencing at the first pore scan'",
        "#         format: '{:,.0f}'",
        "#         hidden: true",
        "#     pores_end:",
        "#         title: 'Pores, end'",
        "#         description: 'Pores found by the last pore scan'",
        "#         format: '{:,.0f}'",
        "#     occupancy:",
        "#         title: 'Pore occupancy'",
        "#         description: 'Mean share of available pores that were sequencing'",
        "#         suffix: '%'",
        "#         min: 0",
        "#         max: 100",
        "#         format: '{:,.1f}'",
        "#     reads:",
        "#         title: 'Reads'",
        "#         description: 'Reads MinKNOW generated'",
        "#         format: '{:,.0f}'",
        "#     pass_bases:",
        "#         title: 'Pass yield'",
        "#         description: 'Bases in reads passing MinKNOW''s live basecalling'",
        "#         suffix: ' Gb'",
        "#         format: '{:,.2f}'",
        "#     est_bases:",
        "#         title: 'Estimated yield'",
        "#         description: 'Bases estimated from the raw signal'",
        "#         suffix: ' Gb'",
        "#         format: '{:,.2f}'",
        "#         hidden: true",
        "#     n50:",
        "#         title: 'Read N50'",
        "#         description: 'Read length N50 of MinKNOW''s live basecalls'",
        "#         suffix: ' bp'",
        "#         format: '{:,.0f}'",
        "#     model:",
        "#         title: 'Live model'",
        "#         description: 'Model MinKNOW basecalled with live; the pipeline re-basecalls with --basecall.model'",
        "#         scale: false",
        "#         hidden: true",
        "#     minknow:",
        "#         title: 'MinKNOW'",
        "#         scale: false",
        "#         hidden: true",
        "Flow cell\tsample\tposition\tflowcell_type\tkit\thours\tstop\tpores_start\tpores_ready\tpores_end\toccupancy\treads\tpass_bases\test_bases\tn50\tmodel\tminknow",
    ]
    def lines = rows.collect { m ->
        def id = counts[ m.flowcell ] > 1 ? "${m.flowcell}_${m.protocol}" : m.flowcell
        ( [ id ] + m.subMap( [ 'sample', 'position', 'flowcell_type', 'kit', 'hours', 'stop', 'pores_start', 'pores_ready', 'pores_end', 'occupancy', 'reads', 'pass_bases', 'est_bases', 'n50', 'model', 'minknow' ] ).values() ).join('\t')
    }
    return ( header + lines ).join('\n') + '\n'
}

workflow NANOPORE_QC {

    take:
    ch_run           // channel: [ val(meta), path(run_dir) ]
    pod5_dir         // file or null: --pod5_dir, where the POD5 live if not in the run dir
    include_failed   // boolean: --include_failed (pod5_fail, and the unclassified reads)
    barcode_kit      // string or null: --barcode_kit, overrides the kit in final_summary
    model            // string: --basecall.model
    max_undetermined // number: --max_undetermined, % of reads unclassified
    min_reads        // number: --min_barcode_reads, fewest reads for a barcode to be a sample

    main:
    ch_versions = channel.empty()

    ch_pod5 = ch_run.flatMap { meta, run_dir ->
        def folders = ontPod5Folders( meta, pod5_dir ?: run_dir, include_failed, barcode_kit )
        if ( !folders ) {
            error "NANOPORE_QC: no *.pod5 found under ${pod5_dir ?: run_dir}. Point --pod5_dir at the POD5 directory if they live elsewhere."
        }
        def no_kit = folders.findAll { meta_, _dir -> meta_.barcoded && !meta_.kit }
        if ( no_kit ) {
            error "NANOPORE_QC: ${no_kit[0][0].flowcell} is barcoded but its final_summary names no barcoding kit. Pass --barcode_kit (e.g. SQK-RBK114-24)."
        }
        folders
    }

    DORADO_BASECALLER( ch_pod5, model )

    ch_called = DORADO_BASECALLER.out.ubam.branch { meta, _bam ->
        barcoded: meta.barcoded
        plain:    true
    }

    // Every barcoded folder of a flow cell (barcodeNN, mixed, unclassified) is
    // demultiplexed together, so each barcode ends up in one file.
    ch_to_demux = ch_called.barcoded
        .map { meta, bam -> tuple( meta.flowcell, meta, bam ) }
        .groupTuple()
        .map { flowcell, metas, bams ->
            tuple( metas[0].subMap( [ 'run', 'platform', 'flowcell', 'sample', 'aliases' ] ) + [ id: flowcell ], bams )
        }

    DORADO_DEMUX( ch_to_demux )

    // toList: names are checked for clashes across flow cells.
    ch_demuxed = DORADO_DEMUX.out.ubam
        .toList()
        .flatMap { demuxed -> ontSamples( demuxed ) }

    // An unbarcoded folder's BAM is its sample, named after MinKNOW's sample id.
    ch_plain = ch_called.plain
        .map { meta, bam -> tuple( meta + [ id: safeId( "${meta.sample}_${meta.flowcell}".toString() ), barcode: null ], bam ) }

    ch_all = ch_demuxed.mix( ch_plain )

    // Each list wrapped in a one-element tuple: combine would otherwise spread
    // a list's elements into the tuple.
    ch_counts = ch_all.toList().map { s -> [ s ] }
        .combine( DORADO_DEMUX.out.summary.toList().map { s -> [ s ] } )
        .map { samples, summaries -> barcodeCounts( samples, summaries ) }

    // A barcode becomes a sample only with --min_barcode_reads reads; an
    // unbarcoded flow cell's sample always does, and so does every barcode
    // when there are no counts to judge by (no barcoding_summary.txt).
    ch_bam = ch_all
        .combine( ch_counts )
        .filter { meta, _bam, counts ->
            meta.barcode == 'unclassified' ? include_failed
                : !meta.barcode || !counts.reads || ( counts.reads[ meta.id ] ?: 0L ) >= min_reads
        }
        .map { meta, bam, _counts -> tuple( meta, bam ) }

    SAMTOOLS_FASTQ( ch_bam )
    NANOPLOT( ch_bam )

    ch_demux = ch_counts.flatMap { counts ->
        def demux = demuxAboveCutoff( counts, min_reads as long )
        demux ? [ tuple( counts.meta, demux ) ] : []
    }
    ch_read_share = ch_demux
        .map { _meta, demux -> readShareTable( demux, 'Unclassified' ) }
        .collectFile( name: 'demux_read_share_mqc.tsv' )
    ch_demux_summary = ch_demux
        .map { meta, demux -> demuxSummaryTable( meta.run, demux, max_undetermined, 'dorado demux barcoding_summary.txt', 'Unclassified', 'NANOPORE_QC' ) }
        .collectFile( name: 'demux_summary_mqc.tsv' )

    // From the run folder, not --pod5_dir: MinKNOW writes its reports next to
    // final_summary, whichever disk the POD5 were moved to.
    ch_run_metrics = ch_run
        .flatMap { _meta, run_dir ->
            def table = ontMetricsTable( run_dir )
            table ? [ table ] : []
        }
        .collectFile( name: 'nanopore_run_metrics_mqc.tsv' )

    ch_qc = channel.empty()
        .mix( NANOPLOT.out.stats.map { _meta, stats -> stats } )
        .mix( ch_run_metrics )
        .mix( ch_read_share )
        .mix( ch_demux_summary )

    ch_versions = ch_versions
        .mix( DORADO_BASECALLER.out.versions_dorado )
        .mix( DORADO_DEMUX.out.versions_dorado )
        .mix( NANOPLOT.out.versions_nanoplot )

    emit:
    reads    = SAMTOOLS_FASTQ.out.reads // FASTQ for scientists
    bam      = ch_bam                   // per-sample uBAM, carries modified-base tags
    qc       = ch_qc
    versions = ch_versions
}
