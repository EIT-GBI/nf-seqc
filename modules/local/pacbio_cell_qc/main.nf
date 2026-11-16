// The PacBio QC plots and tables of the MultiQC report for every SMRT Cell of
// a run, from what SMRT Link's run QC report is drawn from: each cell's
// reports.zip (pbreports) and zmw_metrics.csv.gz. No BAM is read.
//
// meta.cells lists the cells: [ cell: <well>, reports_zip: <file name>,
// zmw_metrics: <file name>, samples: [ barcode: sample id ] ], a missing file
// being null; the sections it feeds are then left out for that cell. The
// files themselves are staged from reports_zips / zmw_metrics (movie-named,
// so unique). All cells go in one task because the distributions are bar
// charts with one bar per cell in each bin row, and MultiQC merges separate
// files by row name, which would drop all but one cell.
//
// Local module — Revio-specific and pipeline-shaped; promote to an nf-mod repo
// only if another pipeline needs it. The script ships as a module binary
// (resources/usr/bin/), so editing it does not invalidate -resume.

process PACBIO_CELL_QC {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(reports_zips, stageAs: 'reports/*'), path(zmw_metrics, stageAs: 'zmw/*')

    output:
    path "*_mqc.json", emit: mqc, optional: true
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/^Python //"'), emit: versions_python, topic: versions

    script:
    // File names as staged, then the whole list single-quoted for the shell,
    // so a ' in a sample id is closed and escaped.
    def cells = ( meta.cells ?: [] ).collect { c ->
        c + [
            reports_zip: c.reports_zip ? "reports/${c.reports_zip}" : null,
            zmw_metrics: c.zmw_metrics ? "zmw/${c.zmw_metrics}" : null,
        ]
    }
    def json = groovy.json.JsonOutput.toJson( cells ).replace("'", "'\\''")
    """
    pacbio_cell_qc.py \\
        --cells '${json}' \\
        --mqc-dir .
    """

    stub:
    def cell = ( meta.cells ?: [ [ cell: 'cell' ] ] )[0].cell
    """
    echo '{"id": "pacbio_loading", "plot_type": "bargraph", "data": {"${cell}": {"p0": 10, "p1": 80, "p2": 10}}}' > pacbio_loading_mqc.json
    """
}
