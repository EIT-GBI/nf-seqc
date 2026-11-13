// The PacBio QC plots and tables of the MultiQC report for one SMRT Cell, from
// what SMRT Link's run QC report is drawn from: the cell's reports.zip
// (pbreports) and zmw_metrics.csv.gz. No BAM is read.
//
// meta.id names the cell (its well, as the Run Metrics table names it);
// meta.samples ([ barcode: sample id ]) names the per-barcode series after
// the pipeline's samples. Either input file may be absent ([]); the sections
// it feeds are then left out.
//
// Local module — Revio-specific and pipeline-shaped; promote to an nf-mod repo
// only if another pipeline needs it. The script ships as a module binary
// (resources/usr/bin/), so editing it does not invalidate -resume.

process PACBIO_CELL_QC {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(reports_zip, stageAs: 'reports/*'), path(zmw_metrics, stageAs: 'zmw/*')

    output:
    path "mqc/*_mqc.json", emit: mqc, optional: true
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/^Python //"'), emit: versions_python, topic: versions

    script:
    // Single-quoted for the shell, so a ' in a sample id is closed and escaped.
    def samples = groovy.json.JsonOutput.toJson( meta.samples ?: [:] ).replace("'", "'\\''")
    def zip     = reports_zip ? "--reports-zip ${reports_zip}" : ''
    def zmw     = zmw_metrics ? "--zmw-metrics ${zmw_metrics}" : ''
    """
    pacbio_cell_qc.py \\
        --cell '${meta.id}' \\
        --samples '${samples}' \\
        --mqc-dir mqc \\
        ${zip} ${zmw}
    """

    stub:
    """
    mkdir mqc
    echo '{"id": "pacbio_loading", "plot_type": "bargraph", "data": {"${meta.id}": {"p0": 10, "p1": 80, "p2": 10}}}' > mqc/${meta.id}.pacbio_loading_mqc.json
    """
}
