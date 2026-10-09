//
// nf-seQC — batch-level QC for one sequencing run.
//
// Input is the complete output directory of a sequencing run. Layout and
// metadata files differ per instrument, so the platform is detected first and
// the run is dispatched to a platform-specific subworkflow. Each branch
// normalises reads for delivery and emits QC files in a shape MultiQC can read,
// which BATCH_REPORT aggregates into one report for the whole run.
//

include { ILLUMINA_QC  } from './subworkflows/local/illumina_qc.nf'
include { AVITI_QC     } from './subworkflows/local/aviti_qc.nf'
include { NANOPORE_QC  } from './subworkflows/local/nanopore_qc.nf'
include { PACBIO_QC    } from './subworkflows/local/pacbio_qc.nf'
include { BATCH_REPORT } from './subworkflows/local/batch_report.nf'

// Entries under a run directory matching a glob. `**/` requires at least one
// directory level, so root-level files need their own pattern alongside it.
// Kept as a named function rather than a closure: the Nextflow language server
// does not resolve closures assigned to local variables.
def matchingFiles( run_dir, String pattern ) {
    return files("${run_dir}/${pattern}").findAll { f -> f.exists() }
}

// Identify the instrument from the marker files each writes at the top of its
// run folder. AVITI is checked before Illumina because both write a file called
// RunParameters, differing only in extension (.json vs .xml).
def detectPlatform( run_dir ) {
    if ( file("${run_dir}/RunParameters.json").exists() ) {
        return 'aviti'
    }
    if ( file("${run_dir}/RunInfo.xml").exists() ) {
        return 'illumina'
    }
    if ( matchingFiles(run_dir, 'final_summary_*.txt')
      || matchingFiles(run_dir, '*.pod5')
      || matchingFiles(run_dir, '**/*.pod5') ) {
        return 'nanopore'
    }
    // Revio writes its dataset manifests at the top of the run directory.
    if ( matchingFiles(run_dir, '*.consensusreadset.xml')
      || matchingFiles(run_dir, '**/*.consensusreadset.xml')
      || matchingFiles(run_dir, '*.hifi_reads*.bam')
      || matchingFiles(run_dir, '**/*.hifi_reads*.bam') ) {
        return 'pacbio'
    }

    error "Could not identify the sequencing platform of ${run_dir}. Expected one of: RunParameters.json (AVITI), RunInfo.xml (Illumina), final_summary_*.txt or *.pod5 (Nanopore), *.consensusreadset.xml or *.hifi_reads.bam (PacBio). Override with --platform <illumina|aviti|nanopore|pacbio>."
}

// Refuse to report on a run the instrument has not finished writing. Mirrors the
// RTAComplete.txt + CopyComplete.txt gate in
// EIT-GBI/Automate-Seq-Run-Metrics-Collection: InterOp metrics read from a
// partial run are misleading rather than merely incomplete.
def requireRunComplete( run_dir, String run_id, List<String> markers ) {
    def missing = markers.findAll { m -> !file("${run_dir}/${m}").exists() }
    if ( !missing ) {
        return
    }
    if ( !params.allow_incomplete ) {
        error "Run ${run_id} looks unfinished — missing ${missing.join(', ')}. Metrics from a partial run are misleading. Pass --allow_incomplete to override."
    }
    log.warn "Run ${run_id} is missing ${missing.join(', ')}; metrics may be incomplete (--allow_incomplete was set)."
}

workflow {

    main:
    if ( !params.input ) {
        error "No input run directory given. Use --input <run_dir> (one pipeline run = one sequencing run)."
    }

    def run_dir  = file( params.input, checkIfExists: true )
    def platform = params.platform ?: detectPlatform( run_dir )
    def run_id   = params.run_id ?: run_dir.name

    log.info "nf-seQC: run '${run_id}' on platform '${platform}'"

    ch_run = channel.value( tuple( [ id: run_id, run: run_id, platform: platform ], run_dir ) )

    ch_reads    = channel.empty()
    ch_bam      = channel.empty()
    ch_qc       = channel.empty()
    ch_metrics  = channel.empty()
    ch_versions = channel.empty()

    if ( platform == 'illumina' ) {
        requireRunComplete( run_dir, run_id, [ 'RTAComplete.txt', 'CopyComplete.txt' ] )
        ILLUMINA_QC( ch_run )
        ch_reads    = ILLUMINA_QC.out.reads
        ch_qc       = ILLUMINA_QC.out.qc
        ch_metrics  = ILLUMINA_QC.out.metrics
        ch_versions = ILLUMINA_QC.out.versions
    }
    else if ( platform == 'aviti' ) {
        AVITI_QC( ch_run )
        ch_reads    = AVITI_QC.out.reads
        ch_qc       = AVITI_QC.out.qc
        ch_versions = AVITI_QC.out.versions
    }
    else if ( platform == 'nanopore' ) {
        NANOPORE_QC( ch_run )
        ch_reads    = NANOPORE_QC.out.reads
        ch_bam      = NANOPORE_QC.out.bam
        ch_qc       = NANOPORE_QC.out.qc
        ch_versions = NANOPORE_QC.out.versions
    }
    else if ( platform == 'pacbio' ) {
        PACBIO_QC( ch_run )
        ch_reads    = PACBIO_QC.out.reads
        ch_qc       = PACBIO_QC.out.qc
        ch_versions = PACBIO_QC.out.versions
    }
    else {
        error "Unknown --platform '${platform}'. Expected one of: illumina, aviti, nanopore, pacbio."
    }

    BATCH_REPORT( ch_qc, ch_versions, channel.value(run_id) )

    // Publish the [meta, files] channels as-is. Two traps to be aware of, both
    // of which fail silently rather than erroring:
    //
    //  - Do NOT reshape a channel with flatMap before publishing. That loses the
    //    task provenance the output block needs and publishes nothing. `map` is
    //    fine, as are tuple channels straight off a process.
    //  - Only task outputs can be published. `reads` therefore covers the
    //    AVITI (bases2fastq) and Nanopore (samtools fastq) branches only —
    //    Illumina FASTQ and Revio HiFi uBAM are written by the instrument and
    //    stay in the run folder rather than being copied into the results.
    publish:
    reads    = ch_reads
    bam      = ch_bam
    qc       = ch_qc
    metrics  = ch_metrics
    report   = BATCH_REPORT.out.report
    versions = BATCH_REPORT.out.versions
}

output {
    reads {
        path 'reads'
    }
    bam {
        path 'bam'
    }
    qc {
        path 'qc'
    }
    // One-row-per-run CSV, column-compatible with the log kept by
    // EIT-GBI/Automate-Seq-Run-Metrics-Collection so rows can be concatenated.
    metrics {
        path 'metrics'
    }
    report {
        path 'multiqc'
    }
    versions {
        path 'pipeline_info'
    }
}
