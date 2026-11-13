//
// nf-seqc — batch-level QC for one sequencing run.
//
// Input is the complete output directory of a sequencing run. Layout and
// metadata files differ per instrument, so the platform is detected first and
// the run is dispatched to a platform-specific subworkflow. Each branch
// normalises reads for delivery and emits QC files in a shape MultiQC can read,
// which BATCH_REPORT aggregates into one report for the whole run.
//

include { ILLUMINA_QC  } from './subworkflows/local/illumina_qc.nf'
include { AVITI_QC     } from './subworkflows/local/aviti_qc.nf'
include { NANOPORE_QC; nanoporeUnfinished } from './subworkflows/local/nanopore_qc.nf'
include { PACBIO_QC; revioUnfinishedCells } from './subworkflows/local/pacbio_qc.nf'
include { BATCH_REPORT } from './subworkflows/local/batch_report.nf'
include { COPY_READS   } from './modules/local/copy_reads/main.nf'

// Instrument directory names start with the brand: <brand>_<model>_<machine_ID>.
// The platform each brand's run folder should be detected as; null for a brand
// not listed here.
def instrumentPlatform( String instrument ) {
    def brands = [ illumina: 'illumina', element: 'aviti', oxford_nanopore: 'nanopore', pacbio: 'pacbio' ]
    return brands.find { brand, _platform -> instrument.startsWith( "${brand}_" ) }?.value
}

// The run directory named by --instrument/--run. nextflow.config has already
// looked under every <prep_root>/<instrument>/<YYYY>/ for it (it builds the
// output path from the year) and left every match in params.run_matches,
// separated by commas; here that becomes exactly one directory or a clear error.
def resolveRun( String instrument, String run, String prep_root, String found ) {
    def machine_dir = file( "${prep_root}/${instrument}" )
    if ( !machine_dir.isDirectory() ) {
        def known = files( "${prep_root}/*", type: 'dir' ).collect { d -> d.name }.sort()
        error "No instrument directory ${machine_dir}. Instruments under ${prep_root}: ${known ? known.join(', ') : 'none'}."
    }
    def matches = found ? found.tokenize(',') : []
    if ( !matches ) {
        error "Run '${run}' not found under ${machine_dir}/<YYYY>/. Check --run against the run folder name."
    }
    if ( matches.size() > 1 ) {
        error "Run '${run}' found in more than one year under ${machine_dir}: ${matches.join(', ')}. Pass the one wanted with --input and --outdir."
    }
    return file( matches[0] )
}

// A run is published once. The first launch claims the output directory by
// creating qc/nf-seqc_session.txt, holding its session ID, before any
// task runs. A later launch into a directory with files in it is refused unless
// it is a -resume of that same session (Nextflow keeps the session ID on
// -resume; a bare -resume takes the last launch's, which after a refused
// launch is the wrong one, so the error names it) or --overwrite is given.
// The claim file is created atomically, so of two launches that both find the
// directory empty, only one proceeds.
//
// A folder published before the option C layout (reads/, multiqc/, ...) is
// refused outright, even with -resume or --overwrite: either would add the new
// layout beside the old one and leave every read file there twice.
def claimOutput( out, String session, boolean resuming, boolean overwrite ) {
    def old_layout = [ 'reads', 'bam', 'multiqc', 'metrics', 'pipeline_info' ].findAll { d -> out.resolve( d ).isDirectory() }
    if ( old_layout ) {
        error "Output directory ${out} was published in the old layout (${old_layout.collect { d -> "${d}/" }.join(', ')}). Move or remove it before publishing this run again."
    }
    def claim = out.resolve( 'qc/nf-seqc_session.txt' )
    def owner = claim.exists() ? claim.text.trim() : null
    if ( resuming && owner == session ) {
        return
    }
    if ( out.exists() && out.list() ) {
        if ( !overwrite ) {
            def resume = owner ? "-resume ${owner} to continue the session that wrote them" : '-resume <session> to continue the session that wrote them'
            error "Output directory ${out} already has results. Pass ${resume}, or --overwrite to write over them."
        }
        log.warn "Writing over the results already in ${out} (--overwrite was set)."
        claim.parent.mkdirs()
        claim.text = "${session}\n"
        return
    }
    claim.parent.mkdirs()
    try {
        java.nio.file.Files.createFile( claim )
    }
    catch ( java.nio.file.FileAlreadyExistsException _e ) {
        error "Output directory ${out} was claimed by another launch while this one was starting (${claim}). Is the same run being processed twice?"
    }
    claim.text = "${session}\n"
}

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

// Command-line values arrive as Strings: "--include_failed false" gives "false",
// which Groovy truth (and `as boolean`) reads as true. Config defaults arrive
// typed, so only strings need parsing.
def flag( value ) {
    return value instanceof CharSequence ? value.toString().toBoolean() : value as boolean
}

// Refuse to report on a run the instrument has not finished writing. Mirrors the
// RTAComplete.txt + CopyComplete.txt gate in
// EIT-GBI/Automate-Seq-Run-Metrics-Collection: InterOp metrics read from a
// partial run are misleading rather than merely incomplete.
//
// `missing` names what is absent: marker files (Illumina) or unfinished SMRT
// Cells (PacBio).
def requireRunComplete( String run_id, List<String> missing, boolean allow_incomplete ) {
    if ( !missing ) {
        return
    }
    if ( !allow_incomplete ) {
        error "Run ${run_id} looks unfinished — missing ${missing.join(', ')}. Metrics from a partial run are misleading. Pass --allow_incomplete to override."
    }
    log.warn "Run ${run_id} is missing ${missing.join(', ')}; metrics may be incomplete (--allow_incomplete was set)."
}

workflow {

    main:
    def by_name = params.instrument || params.run
    if ( by_name && !( params.instrument && params.run ) ) {
        error "--instrument and --run go together: --instrument <brand>_<model>_<machine_ID> --run <run folder>."
    }
    if ( by_name && params.input ) {
        error "Give either --input or --instrument/--run, not both."
    }
    if ( !by_name && !params.input ) {
        error "No run given. Use --instrument <brand>_<model>_<machine_ID> --run <run folder>, or --input <run_dir> (one pipeline run = one sequencing run)."
    }

    def run_dir = by_name
        ? resolveRun( params.instrument, params.run, params.prep_root, params.run_matches )
        : file( params.input, checkIfExists: true )
    def platform = params.platform ?: detectPlatform( run_dir )
    def run_id   = params.run_id ?: run_dir.name

    // The brand in the instrument name and the files in its run folder must
    // agree; a mismatch means the wrong run, or a run filed under the wrong
    // instrument, and results would be published under that instrument.
    if ( by_name ) {
        def expected = instrumentPlatform( params.instrument )
        if ( !expected ) {
            log.warn "Unknown brand in instrument '${params.instrument}'; cannot check it against the detected platform '${platform}'."
        }
        else if ( expected != platform ) {
            error "Instrument '${params.instrument}' should hold ${expected} runs, but run ${run_dir} looks like ${platform}."
        }
    }

    claimOutput( workflow.outputDir, workflow.sessionId.toString(), workflow.resume, flag( params.overwrite ) )

    log.info "nf-seqc: run '${run_id}' on platform '${platform}'"
    log.info "nf-seqc: reading ${run_dir}"
    log.info "nf-seqc: publishing to ${workflow.outputDir}"

    // params are read here, in the entry workflow, and handed to the
    // subworkflows as inputs: reading params anywhere else is deprecated.
    // Directory overrides are checked now, so a typo fails before any task.
    def fastq_dir = params.fastq_dir ? file( params.fastq_dir, checkIfExists: true ) : null
    def pod5_dir  = params.pod5_dir  ? file( params.pod5_dir,  checkIfExists: true ) : null
    def bam_dir   = params.bam_dir   ? file( params.bam_dir,   checkIfExists: true ) : null
    def multiqc_config = params.multiqc_config ? file( params.multiqc_config, checkIfExists: true ) : null
    // Likewise a number: "2" compares with a BigDecimal without error, and
    // always false
    def max_undetermined = params.max_undetermined as BigDecimal

    ch_run = channel.value( tuple( [ id: run_id, run: run_id, platform: platform ], run_dir ) )

    ch_reads    = channel.empty()
    ch_bam      = channel.empty()
    ch_reports  = channel.empty() // instrument vendors' own HTML reports, made by the pipeline
    ch_qc       = channel.empty()
    ch_metrics  = channel.empty()
    ch_versions = channel.empty()

    if ( platform == 'illumina' ) {
        def missing = [ 'RTAComplete.txt', 'CopyComplete.txt' ].findAll { m -> !file("${run_dir}/${m}").exists() }
        requireRunComplete( run_id, missing, flag( params.allow_incomplete ) )
        ILLUMINA_QC( ch_run, fastq_dir, flag( params.include_undetermined ), max_undetermined )
        // The FASTQ are the instrument's, so copy them to make them publishable.
        COPY_READS( ILLUMINA_QC.out.reads )
        ch_reads    = COPY_READS.out.reads
        ch_qc       = ILLUMINA_QC.out.qc
        ch_metrics  = ILLUMINA_QC.out.metrics
        ch_versions = ILLUMINA_QC.out.versions
    }
    else if ( platform == 'aviti' ) {
        AVITI_QC( ch_run, flag( params.include_undetermined ), max_undetermined )
        ch_reads    = AVITI_QC.out.reads
        ch_reports  = AVITI_QC.out.report
        ch_qc       = AVITI_QC.out.qc
        ch_versions = AVITI_QC.out.versions
    }
    else if ( platform == 'nanopore' ) {
        def unfinished = nanoporeUnfinished( pod5_dir ?: run_dir ).collect { acq -> "${acq}/final_summary_*.txt".toString() }
        requireRunComplete( run_id, unfinished, flag( params.allow_incomplete ) )
        NANOPORE_QC( ch_run, pod5_dir, flag( params.include_failed ), params.barcode_kit, params.basecall.model, max_undetermined, params.min_barcode_reads as long )
        ch_bam      = NANOPORE_QC.out.bam
        ch_qc       = NANOPORE_QC.out.qc
        ch_versions = NANOPORE_QC.out.versions
    }
    else if ( platform == 'pacbio' ) {
        def unfinished = revioUnfinishedCells( bam_dir ?: run_dir ).collect { cell -> "${cell}/metadata/*.transferdone".toString() }
        requireRunComplete( run_id, unfinished, flag( params.allow_incomplete ) )
        PACBIO_QC( ch_run, bam_dir, max_undetermined )
        // The HiFi BAM are the instrument's, so copy them to make them
        // publishable, with their PacBio index (.bam.pbi) where there is one.
        COPY_READS( PACBIO_QC.out.reads.map { meta, bam ->
            def pbi = file( "${bam}.pbi" )
            tuple( meta, pbi.exists() ? [ bam, pbi ] : [ bam ] )
        } )
        ch_reads    = COPY_READS.out.reads
        ch_qc       = PACBIO_QC.out.qc
        ch_versions = PACBIO_QC.out.versions
    }
    else {
        error "Unknown --platform '${platform}'. Expected one of: illumina, aviti, nanopore, pacbio."
    }

    BATCH_REPORT( ch_qc, ch_versions, channel.value(run_id), multiqc_config )

    // Publish the [meta, files] channels as-is. Two traps to be aware of, both
    // of which fail silently rather than erroring:
    //
    //  - Do NOT reshape a channel with flatMap before publishing. That loses the
    //    task provenance the output block needs and publishes nothing. `map` is
    //    fine, as are tuple channels straight off a process.
    //  - Only task outputs can be published. Illumina FASTQ and Revio HiFi
    //    uBAM are written by the instrument, so COPY_READS copies them into a
    //    task first; AVITI (bases2fastq) FASTQ and Nanopore (dorado demux)
    //    uBAM are task outputs already.
    publish:
    reads    = ch_reads
    bam      = ch_bam
    qc       = ch_qc
    metrics  = ch_metrics
    report   = BATCH_REPORT.out.report.mix( ch_reports )
    versions = BATCH_REPORT.out.versions
}

// Run folder layout ("option C", agreed with the team on 2026-11-13): what
// scientists open sits at the top of the run folder (the reports and the
// reads), and everything else (FastQC, MultiQC inputs, run metrics, versions,
// the publish claim) goes in qc/.
output {
    reads {
        path '.'
    }
    // Nanopore per-sample uBAM. Renamed after the sample: dorado demux names
    // them <kit>_barcodeNN.bam under demux/.
    bam {
        path { meta, bam -> bam >> "${meta.id}.bam" }
    }
    // The MultiQC report, and on AVITI Element's own <run>_QC.html.
    report {
        path '.'
    }
    qc {
        path 'qc'
    }
    // One-row-per-run CSV, column-compatible with the log kept by
    // EIT-GBI/Automate-Seq-Run-Metrics-Collection so rows can be concatenated.
    metrics {
        path 'qc'
    }
    versions {
        path 'qc'
    }
}
