//
// Subworkflow: BATCH_REPORT
//
// QC is reported per sequencing run, so every per-sample and run-level metric
// file from the platform branch is collected into one MultiQC report.
//

include { MULTIQC } from '../../modules/local/multiqc/main.nf'

workflow BATCH_REPORT {

    take:
    ch_qc       // channel: path — flat QC files from the platform subworkflow
    ch_versions // channel: [ val(process), val(tool), val(version) ]
    ch_run_id   // channel: val(run_id)

    main:
    ch_software_versions = ch_versions
        .map { _process, tool, version -> "${tool}: ${version}" }
        .unique()
        .collectFile( name: 'software_versions.yml', newLine: true, sort: true )

    ch_config = params.multiqc_config
        ? channel.value( file(params.multiqc_config, checkIfExists: true) )
        : channel.value( file("${projectDir}/assets/NO_FILE") )

    MULTIQC( ch_qc.collect(), ch_run_id, ch_config, ch_software_versions )

    emit:
    report   = MULTIQC.out.report
    data     = MULTIQC.out.data
    versions = ch_software_versions
}
