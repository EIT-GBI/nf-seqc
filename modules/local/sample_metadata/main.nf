// The run's sample sheet for scientists: one row per sample with its name,
// barcode/index and the run-level details, as <run>_samplesheet.csv at the
// top of the results folder. Written by the platform's extractor from the run
// folder (sample sheet / run manifest / metadata XML / MinKNOW sample sheet);
// no read is opened. The PacBio one is Melissa's original table.
//
// The scripts are also standalone (Python 3.8+, standard library only); see
// README.md here for running them by hand. A script that finds no sample
// (e.g. a run folder without its sample sheet) stops with a "No ... found" /
// "No sample records were extracted" message. Here that exits 3 and the task
// is ignored, so the run goes on without a sample sheet and Nextflow prints a
// note naming the task. Any other failure (a crash on a malformed sheet, a
// time limit, a missing image) stops the run as usual.
//
// Local module — pipeline-shaped; promote to an nf-mod repo only if another
// pipeline needs it. The scripts ship as module binaries (resources/usr/bin/),
// so editing one does not invalidate -resume.

// A value single-quoted for the shell, any ' in it closed and escaped.
def sampleMetadataQuote( value ) {
    return "'" + value.toString().replace("'", "'\\''") + "'"
}

process SAMPLE_METADATA {
    tag "${meta.id}"
    label 'process_low'
    errorStrategy { task.exitStatus == 3 ? 'ignore' : 'terminate' }

    input:
    tuple val(meta), path(run_dir)

    output:
    tuple val(meta), path("${meta.run}_samplesheet.csv"), emit: samplesheet
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/^Python //"'), emit: versions_python, topic: versions

    script:
    def out = sampleMetadataQuote( "${meta.run}_samplesheet.csv" )
    """
    status=0
    ${meta.platform}_sample_metadata.py ${sampleMetadataQuote( run_dir )} -o ${out} 2> extract.err || status=\$?
    cat extract.err >&2
    if [ "\$status" -ne 0 ]; then
        if grep -qE '^No .*(found|extracted)' extract.err; then
            exit 3
        fi
        exit "\$status"
    fi
    """

    stub:
    """
    echo 'Name' > ${sampleMetadataQuote( "${meta.run}_samplesheet.csv" )}
    """
}
