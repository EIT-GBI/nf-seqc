// The run's sample sheet for scientists: one row per sample with its name,
// barcode/index and the run-level details, as <run>_samplesheet.csv at the
// top of the results folder. Written by the platform's extractor from the run
// folder (sample sheet / run manifest / metadata XML / MinKNOW sample sheet);
// no read is opened. The PacBio one is Melissa's original table.
//
// The scripts are also standalone (Python 3.8+, standard library only); see
// README.md here for running them by hand. A script that finds no sample
// (e.g. a run folder without its sample sheet) exits with an error; here that
// leaves the run without a sample sheet and a warning in the Nextflow log
// (errorStrategy 'ignore') rather than stopping the QC.
//
// Local module — pipeline-shaped; promote to an nf-mod repo only if another
// pipeline needs it. The scripts ship as module binaries (resources/usr/bin/),
// so editing one does not invalidate -resume.

process SAMPLE_METADATA {
    tag "${meta.id}"
    label 'process_low'
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(run_dir)

    output:
    tuple val(meta), path("${meta.run}_samplesheet.csv"), emit: samplesheet
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/^Python //"'), emit: versions_python, topic: versions

    script:
    """
    ${meta.platform}_sample_metadata.py ${run_dir} -o '${meta.run}_samplesheet.csv'
    """

    stub:
    """
    echo 'Name' > '${meta.run}_samplesheet.csv'
    """
}
