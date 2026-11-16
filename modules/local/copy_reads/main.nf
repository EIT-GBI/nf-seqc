// Copy reads the instrument already wrote (Illumina FASTQ, Revio HiFi BAM) so
// they can be published. The `output {}` block only publishes task outputs and
// silently skips files that came from the run directory, so without this task
// those reads never reach the results tree.
//
// `names` gives each copy its published name, in the same order as `reads`
// (Revio BAMs are renamed after their sample; Illumina FASTQ keep theirs).
//
// Inputs are staged under in/ so the copies sit at the top of the task
// directory: published files keep the relative path they had in the task, so
// nothing is emitted from a subdirectory. `type: 'file'` keeps the in/
// directory itself out of the output glob.
//
// No container: plain cat, run on the host. Local module, not tied to a
// platform; promote to an nf-mod repo only if another pipeline needs to
// publish its inputs.

process COPY_READS {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(reads, stageAs: 'in/*'), val(names)

    output:
    tuple val(meta), path('*', type: 'file'), emit: reads

    script:
    // cat, not cp: the staged inputs are symlinks into the run directory, and
    // cat reads through them as a plain byte copy, which is what we want (other
    // filesystem in production). It works with every cat, GNU or macOS, and
    // nothing here needs cp's extras. Unlike cp, the copy takes the default
    // mode (umask), not the source's, so a read-only source is not delivered
    // read-only; delivered permissions are still to be agreed with scientific
    // computing. Reads from the sandpit-tokyo Lustre intermittently fail with
    // "No data available" (since 2026-11-13); conf/base.config retries the
    // task after a wait, see the CLAUDE.md gotcha.
    def copies = copyPairs( reads, names ).collect { src, dest -> "cat -- '${src}' > '${dest}'" }.join('\n')
    """
    ${copies}
    """

    stub:
    def touches = copyPairs( reads, names ).collect { _src, dest -> "touch '${dest}'" }.join('\n')
    """
    ${touches}
    """
}

// [ staged path (in/<file>), published name ] for each read file. A single
// staged file arrives as a Path, which iterates as its directory components,
// so wrap it. The staged path already carries the in/ prefix.
def copyPairs( reads, names ) {
    def files = reads instanceof List ? reads : [ reads ]
    if ( files.size() != names.size() ) {
        error "COPY_READS: ${files.size()} files but ${names.size()} names"
    }
    return [ files.collect { f -> f.toString() }, names ].transpose()
}
