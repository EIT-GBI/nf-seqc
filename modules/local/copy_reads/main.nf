// Copy reads the instrument already wrote (Illumina FASTQ, Revio HiFi BAM) so
// they can be published. The `output {}` block only publishes task outputs and
// silently skips files that came from the run directory, so without this task
// those reads never reach the results tree.
//
// Inputs are staged under in/ so the copies can take the original names at the
// top of the task directory: published files keep the relative path they had
// in the task, so nothing is emitted from a subdirectory. `type: 'file'` keeps
// the in/ directory itself out of the output glob.
//
// No container: plain cp, run on the host. Local module, not tied to a
// platform; promote to an nf-mod repo only if another pipeline needs to
// publish its inputs.

process COPY_READS {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(reads, stageAs: 'in/*')

    output:
    tuple val(meta), path('*', type: 'file'), emit: reads

    script:
    // -L: the staged inputs are symlinks into the run directory.
    """
    cp -L in/* .
    """

    stub:
    """
    for f in in/*; do touch "\$(basename "\$f")"; done
    """
}
