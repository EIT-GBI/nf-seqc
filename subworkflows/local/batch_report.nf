//
// Subworkflow: BATCH_REPORT
//
// QC is reported per sequencing run, so every per-sample and run-level metric
// file from the platform branch is collected into one MultiQC report.
//

include { MULTIQC } from '../../modules/local/multiqc/main.nf'

workflow BATCH_REPORT {

    take:
    ch_qc          // channel: path — flat QC files from the platform subworkflow
    ch_versions    // channel: [ val(process), val(tool), val(version) ]
    ch_run_id      // channel: val(run_id)
    multiqc_config // file or null: --multiqc_config, applied after the default

    main:
    // The pipeline's own version always heads the list. It also keeps the file
    // from being empty, whatever ran before MultiQC: with no
    // software_versions.yml MultiQC would never start.
    ch_software_versions = ch_versions
        .mix( channel.of( tuple( 'NF_SEQC', 'nf-seqc', workflow.manifest.version ) ) )
        .map { _process, tool, version -> "${tool}: ${version}" }
        .unique()
        .collectFile( name: 'software_versions.yml', newLine: true, sort: true )

    // The pipeline's own config always applies; --multiqc_config is layered on
    // top of it, so it can override anything there.
    def configs = [ file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true) ]
    if ( multiqc_config ) {
        configs << multiqc_config
    }
    ch_config = channel.value( configs )

    MULTIQC( ch_qc.collect(), ch_run_id, ch_config, ch_software_versions )

    emit:
    report   = MULTIQC.out.report
    data     = MULTIQC.out.data
    versions = ch_software_versions
}
