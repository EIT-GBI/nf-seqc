//
// Checks how Nanopore BAMs are named from dorado demux's real file names, which
// the -stub run cannot reach (the module's stub writes barcode01.bam).
// No task runs and no file is read.
//
//     nextflow run tests/ont_samples.nf
//
// A flow cell that ran two MinKNOW protocols (oxford_nanopore_promethion_PCA100516
// / 20260928_Wei_samples, PBO58669) gets one demux BAM per barcode per
// acquisition: each must be delivered under its own name, with one sample row
// per barcode.
//

include { ontSamples; barcodeCounts; ontAcquisitionRuns } from '../subworkflows/local/nanopore_qc.nf'

def check( boolean ok, String what ) {
    if ( !ok ) {
        error "ont_samples: ${what}"
    }
    log.info "ok: ${what}"
}

workflow {
    main:
    def acq_a = '27de4d4f10a94db01d309f53c8d0e6f8ba72544a'
    def acq_b = '1bf26ea348dff2da197b98ce46285119e50d3366'
    // final_summary as MinKNOW writes it (the keys used), one per acquisition.
    def fixture = file( "${workflow.workDir}/ont_samples_run" )
    [ [ 'Wei_samples/20260928_1319_1G_PBO58669_6481f7b6', acq_a, '6481f7b6-cc05-4ebf-9489-67f43681cf8f', '2026-09-28T13:25:30.667351+00:00' ],
      [ 'no_sample_id/20260929_1010_1G_PBO58669_73d69264', acq_b, '73d69264-6c5d-456d-8c5c-9c8829906141', '2026-09-29T10:16:06.331727+00:00' ] ].each { dir, acq, prot, started ->
        def d = fixture.resolve( dir )
        d.mkdirs()
        d.resolve( "final_summary_PBO58669_${prot.take(8)}_${acq.take(8)}.txt" ).text = "instrument=1G\nflow_cell_id=PBO58669\nprotocol_run_id=${prot}\nacquisition_run_id=${acq}\nstarted=${started}\n"
    }
    def runs  = ontAcquisitionRuns( [ fixture, null ] )
    check( runs.collectEntries { a, v -> [ ( a ): v.run ] } == [ ( acq_a ): '6481f7b6', ( acq_b ): '73d69264' ], "protocol runs from final_summary: ${runs}" )
    def kit   = 'SQK-NBD114-96'
    def meta  = [ run: 'r', platform: 'nanopore', flowcell: 'PBO58669', sample: 'Wei_samples', aliases: [:], status: 'pod5_pass' ]
    def bams  = [ "${acq_a}_${kit}_barcode09", "${acq_b}_${kit}_barcode09", "${acq_a}_unclassified", "${acq_b}_unclassified", "${acq_b}_${kit}_barcode15" ]
        .collect { n -> file( "/demux/${n}.bam" ) }
    def fail  = meta + [ status: 'pod5_fail' ]
    def fbams = [ "${acq_a}_${kit}_barcode09", "${acq_b}_${kit}_barcode09" ].collect { n -> file( "/demux_fail/${n}.bam" ) }

    def out = ontSamples( [ [ meta, bams ], [ fail, fbams ] ], runs )
    def ids = out.collect { m, _b -> m.id }
    check( ids.size() == ids.unique( false ).size(), "names are unique: ${ids}" )
    check( ids.containsAll( [ 'barcode09_PBO58669_6481f7b6', 'barcode09_PBO58669_73d69264', 'barcode15_PBO58669_73d69264' ] ), 'pass BAMs named after the protocol run' )
    check( ids.containsAll( [ 'unclassified_PBO58669_6481f7b6', 'unclassified_PBO58669_73d69264' ] ), 'unclassified BAMs named after the protocol run' )
    check( ids.containsAll( [ 'barcode09_PBO58669_6481f7b6_fail', 'barcode09_PBO58669_73d69264_fail' ] ), 'failed BAMs named after the protocol run' )

    def rows = out.findAll { m, _b -> m.report && m.status == 'pod5_pass' && m.barcode != 'unclassified' }.collect { m, _b -> m.report_id }.sort()
    check( rows == [ 'barcode09', 'barcode15' ], "one sample row per barcode, a barcode of the second protocol only included: ${rows}" )
    def row09 = out.find { m, _b -> m.report && m.barcode == 'barcode09' && m.status == 'pod5_pass' }[0].id
    check( row09 == 'barcode09_PBO58669_6481f7b6', "the row is on the protocol that started first: ${row09}" )

    // Two acquisitions of one protocol run: named by acquisition instead.
    def same = ontSamples( [ [ meta, bams.take(2) ] ], [ ( acq_a ): [ run: '6481f7b6', started: '' ], ( acq_b ): [ run: '6481f7b6', started: '' ] ] ).collect { m, _b -> m.id }.sort()
    check( same == [ 'barcode09_PBO58669_1bf26ea3', 'barcode09_PBO58669_27de4d4f' ], "one protocol run, two acquisitions: ${same}" )

    // Read counts go to the shared sample row, whichever BAM the row is on.
    def summary = file( "${workflow.workDir}/ont_samples_barcoding_summary.txt" )
    summary.parent.mkdirs()
    def barcodes = [ 'SQK-NBD114-96_barcode09' ] * 3 + [ 'SQK-NBD114-96_barcode15', 'unclassified', 'unclassified' ]
    summary.text = ( [ 'read_id\tbarcode' ] + barcodes.withIndex().collect { b, i -> "read${i}\t${b}" } ).join('\n') + '\n'
    def counts = barcodeCounts( out, [ [ meta, summary ] ] )
    check( counts.reads == [ barcode09: 3L, barcode15: 1L ] && counts.unclassified == 2L, "reads counted per sample row: ${counts.reads}, ${counts.unclassified} unclassified" )

    // Unchanged when the flow cell ran one protocol: no protocol run in names.
    def single = ontSamples( [ [ meta, bams.findAll { b -> b.name.startsWith( acq_a ) } ] ], runs ).collect { m, _b -> m.id }.sort()
    check( single == [ 'barcode09', 'unclassified_PBO58669' ], "one protocol keeps the old names: ${single}" )

    // And with the module's stub names (no acquisition prefix).
    def stub = ontSamples( [ [ meta, [ file('/demux/barcode01.bam'), file('/demux/unclassified.bam') ] ] ] ).collect { m, _b -> m.id }.sort()
    check( stub == [ 'barcode01', 'unclassified_PBO58669' ], "stub names: ${stub}" )

    // Without a final_summary, the acquisition id's first 8 characters.
    def bare = ontSamples( [ [ meta, bams.take(2) ] ], [:] ).collect { m, _b -> m.id }.sort()
    check( bare == [ 'barcode09_PBO58669_1bf26ea3', 'barcode09_PBO58669_27de4d4f' ], "no final_summary: ${bare}" )
}
