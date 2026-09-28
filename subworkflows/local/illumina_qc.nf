//
// Subworkflow: ILLUMINA_QC
//
// NextSeq demultiplexes on-instrument, so FASTQ already exist inside the run
// folder. This branch therefore only reads: per-sample FASTQ quality (FastQC)
// and run-level metrics from the binary InterOp files (%PF, %occupied, Q30,
// phasing, PhiX error rate) — none of which reach the FASTQ at all.
//
// The run-level half mirrors EIT-GBI/Automate-Seq-Run-Metrics-Collection: same
// metrics, same column names, same %Occupied vs %PF plot. Two behaviours from
// that tool are handled differently here:
//
//   * Its `RTAComplete.txt` + `CopyComplete.txt` gate is kept, but lives in
//     `requireRunComplete()` in main.nf: raising it there, where run_dir is a
//     plain value, gives one clean error instead of an abort inside a channel
//     operator followed by a spurious "empty List" message.
//   * Its `monitoring_complete.txt` marker exists to stop cron reprocessing a
//     run. Nextflow's -resume covers that, so no marker is written; nothing in
//     the run directory is modified.
//

include { FASTQC_FASTQC        } from '../../modules/fastqc/fastqc/main.nf'
include { INTEROP_RUN_METRICS  } from '../../modules/local/interop/main.nf'

workflow ILLUMINA_QC {

    take:
    ch_run // channel: [ val(meta), path(run_dir) ]

    main:
    ch_versions = channel.empty()

    // FASTQ location varies by layout: Analysis/<n>/Data/... on DRAGEN-enabled
    // NextSeq 1000/2000, Alignment_*/.../Fastq on older setups. Glob instead of
    // hardcoding; --fastq_dir overrides when the FASTQ live outside the run dir.
    ch_pairs = ch_run.flatMap { meta, run_dir ->
        def root = params.fastq_dir ? file(params.fastq_dir, checkIfExists: true) : run_dir
        def r1s  = ( [] + file("${root}/**/*_R1_001.fastq.gz") + file("${root}/*_R1_001.fastq.gz") )
            .findAll { f -> f.exists() }
            .unique()

        if ( !params.include_undetermined ) {
            r1s = r1s.findAll { f -> !f.name.startsWith('Undetermined') }
        }
        if ( !r1s ) {
            error "ILLUMINA_QC: no *_R1_001.fastq.gz found under ${root}. Has the run been demultiplexed? Point --fastq_dir at the FASTQ directory if they live elsewhere."
        }

        r1s.sort().collect { r1 ->
            def r2  = file( r1.toString().replace('_R1_001.fastq.gz', '_R2_001.fastq.gz') )
            def sid = r1.name
                .replaceAll(/_S\d+(_L\d{3})?_R1_001\.fastq\.gz$/, '')
                .replaceAll(/_R1_001\.fastq\.gz$/, '')
            tuple( meta + [ id: sid, single_end: !r2.exists() ], r1, r2.exists() ? r2 : [] )
        }
    }

    FASTQC_FASTQC( ch_pairs )
    INTEROP_RUN_METRICS( ch_run )

    // The metrics table and occupancy plot go to MultiQC as custom content;
    // the CSV itself is the machine-readable log row and is published separately.
    ch_qc = channel.empty()
        .mix( FASTQC_FASTQC.out.zip.map { _meta, zip -> zip } )
        .mix( INTEROP_RUN_METRICS.out.mqc_table.map { _meta, tsv -> tsv } )
        .mix( INTEROP_RUN_METRICS.out.plot.map { _meta, png -> png } )

    ch_versions = ch_versions
        .mix( FASTQC_FASTQC.out.versions_fastqc )
        .mix( INTEROP_RUN_METRICS.out.versions_interop )

    emit:
    reads    = ch_pairs.map { meta, r1, r2 -> tuple( meta, [ r1, r2 ].flatten().findAll { f -> f } ) }
    qc       = ch_qc                          // channel: path — flat files for MultiQC
    metrics  = INTEROP_RUN_METRICS.out.metrics // channel: [ meta, run_metrics.csv ]
    versions = ch_versions                     // channel: version records collected from modules
}
