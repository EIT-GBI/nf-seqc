#!/usr/bin/env python3
"""Collect run-level QC metrics for one completed Illumina NextSeq 2000 run.

Port of utils.collect_run_metrics / plot_occupied_vs_pf from
EIT-GBI/Automate-Seq-Run-Metrics-Collection, which this deliberately mirrors:
the metric set, the column names and the plot are kept identical so rows from
either tool line up in the same log.

Differences from the original, both forced by running under Nextflow:

  * That tool is a cron collector that APPENDS to one central CSV across many
    runs. Here one pipeline run is one sequencing run, so this writes a
    single-row CSV per run (header + one line) and lets the caller concatenate.
    Appending to a shared file from a task would break -resume and race when
    two runs are processed at once.
  * Missing optional inputs degrade to "n/a" instead of raising. A QC report
    that is missing PrimaryAnalysisMetrics.csv is still worth having, and the
    original would traceback on it.

Also emits MultiQC custom-content files so the table and plot appear in the
batch report: *_mqc.tsv and *_mqc.png.
"""

import argparse
import csv
import math
import os
import re
import sys
import xml.etree.ElementTree as ET

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from interop import py_interop_run_metrics, py_interop_summary

# Column order is the contract with the colleague's log CSV — keep in sync.
LOG_COLUMNS = [
    'RunID', 'Run_nr', 'RunStart', 'RunLength', 'Flowcell', 'Cycles',
    'TotalReadsPF (M)', 'Yield (Gb)', 'Q30 (%)', 'PF (%)',
    'PF_median (%)', 'PF_std', 'Occupied (%)', 'Occupied_median (%)', 'Occupied_std',
    'LoadingConc (pM)', 'PhixAligned (%)', 'PhixError (%)',
    'Phasing (%)', 'Prephasing (%)', 'FirstCycleIntensity', 'Concerns',
    'SampleSheet', 'SerialNumber',
]

NA = 'n/a'


def _nan_to_na(value):
    try:
        return NA if math.isnan(value) else value
    except TypeError:
        return value


def _safe(fn, default=NA):
    """Evaluate fn(), returning `default` if it raises or yields NaN.

    Every metric is individually guarded: a patterned-flowcell-only field such
    as percent_occupied is absent on other flowcells, and one missing field
    should not cost the whole row.
    """
    try:
        return _nan_to_na(fn())
    except Exception as exc:  # noqa: BLE001 - deliberately broad, see docstring
        print(f'  warning: could not read metric ({exc})', file=sys.stderr)
        return default


def _parse_run_elapsed_time(text):
    """Parse RunElapsedTime to a DD:HH:MM:SS string.

    The NextSeq 2000 stores this as "DD:HH:MM:SS:subseconds" (five colon-
    separated parts). Earlier/alternative formats use fewer parts or a trailing
    fractional-seconds field separated by '.'. Keep at most four parts and
    strip any fractional tail.
    """
    s = text.strip().split(' ')[0]
    parts = s.split(':')
    if len(parts) > 4:
        parts = parts[:4]
    if parts:
        parts[-1] = parts[-1].split('.')[0]
    return ':'.join(parts)


def _parse_flowcell(cartridge_mode_text):
    """Extract a concise flowcell label from the CartridgeMode string.

      "NextSeq 2000 P2 (100 Cycles)"       -> "P2 100"
      "NextSeq 2000 P2 Reagents (100 ...)" -> "P2 100"

    Cycle count is rounded to the nearest 100 at 100 or above, otherwise 50
    (round(50, -2) is 0, hence the floor).
    Falls back to the raw text if the expected pattern is not found.

    NOTE: the original uses `> 100` here, which sends a 100-cycle cartridge —
    the commonest P2 kit — down the else branch and labels it "P2 50",
    contradicting its own docstring example. Fixed to `>= 100`. Worth fixing
    upstream too, otherwise the two logs disagree on every 100-cycle run.
    """
    match = re.search(r'(P\d)\b.*?\((\d+)\s+Cycles?\)', cartridge_mode_text, re.IGNORECASE)
    if match:
        kit = match.group(1).upper()
        n_cycles = int(match.group(2))
        n_cycles = round(n_cycles, -2) if n_cycles >= 100 else 50
        return f'{kit} {n_cycles}'
    return cartridge_mode_text


def _read_total_reads(run_dir):
    """Total reads PF from PrimaryAnalysisMetrics.csv.

    Positional lookup kept identical to the original's pandas `iloc[2, 2]`:
    third data row, third column, header excluded.
    """
    pam = os.path.join(run_dir, 'PrimaryAnalysisMetrics', 'PrimaryAnalysisMetrics.csv')
    if not os.path.isfile(pam):
        print(f'  warning: {pam} not found, TotalReadsPF unavailable', file=sys.stderr)
        return NA
    with open(pam, newline='') as fh:
        rows = list(csv.reader(fh))
    data = rows[1:]
    try:
        return data[2][2]
    except IndexError:
        print(f'  warning: {pam} has an unexpected shape, TotalReadsPF unavailable', file=sys.stderr)
        return NA


def _read_run_parameters(run_dir):
    """Pull the RunParameters.xml fields, each independently optional."""
    out = {
        'SampleSheet': NA, 'SerialNumber': NA, 'Cycles': NA,
        'Flowcell': NA, 'RunStart': NA, 'RunLength': NA, 'Run_nr': NA,
    }
    xml_path = os.path.join(run_dir, 'RunParameters.xml')
    if not os.path.isfile(xml_path):
        print(f'  warning: {xml_path} not found, run parameters unavailable', file=sys.stderr)
        return out

    root = ET.parse(xml_path).getroot()

    def text_of(tag):
        node = root.find(tag)
        return node.text if node is not None and node.text else None

    ss_raw = text_of('SampleSheetFilePath')
    if ss_raw:
        # Split on both separators: the original uses os.sep, which leaves a
        # Windows-style path from RunParameters.xml completely unsplit when the
        # pipeline runs on Linux, putting the full path in the log.
        parts = [p for p in re.split(r'[\\/]+', ss_raw.strip()) if p]
        # Keep the last two components (project/filename)
        out['SampleSheet'] = '/'.join(parts[-2:]) if len(parts) >= 2 else (parts[-1] if parts else NA)

    out['SerialNumber'] = text_of('FlowCellSerialNumber') or NA

    completed = root.find('CompletedCycles')
    if completed is not None:
        out['Cycles'] = ' '.join(c.text for c in completed if c.text)

    cartridge = text_of('CartridgeMode')
    if cartridge:
        out['Flowcell'] = _parse_flowcell(cartridge)

    start = text_of('RunStartTime')
    if start:
        out['RunStart'] = start.split('.')[0]

    elapsed = text_of('RunElapsedTime')
    if elapsed:
        out['RunLength'] = _parse_run_elapsed_time(elapsed)

    output_folder = text_of('OutputFolder')
    if output_folder:
        # Folder format: YYYYMMDD_<instrument>_<RunNumber>_<FlowcellID>
        name = os.path.basename(os.path.normpath(output_folder))
        bits = name.split('_')
        if len(bits) >= 3:
            out['Run_nr'] = bits[2]

    return out


def collect_run_metrics(run_dir, run_id):
    """Collect QC metrics for a completed NextSeq run. Returns a LOG_COLUMNS dict."""
    run_dir = os.path.abspath(run_dir)

    run_metrics = py_interop_run_metrics.run_metrics()
    run_metrics.read(run_dir)

    summary = py_interop_summary.run_summary()
    py_interop_summary.summarize_run_metrics(run_metrics, summary)

    # total_summary covers all reads+lanes; nonindex_summary excludes index reads
    total = summary.total_summary()
    nonidx = summary.nonindex_summary()
    # Per-read/per-lane metrics: Read 1 (index 0), Lane 1 (index 0)
    read1_lane = summary.at(0).at(0)

    # Quality and yield — nonindex so short index cycles don't dilute Q30
    percent_q30 = _safe(lambda: round(nonidx.percent_gt_q30(), 2))
    yield_g = _safe(lambda: round(nonidx.yield_g(), 2))

    # Loading / occupancy (patterned nanowell flowcell specific)
    percent_occupied = _safe(lambda: round(total.percent_occupied(), 2))
    loading_conc = _safe(lambda: round(total.percent_occupancy_proxy(), 2))

    # PhiX spike-in metrics
    phix_aligned = _safe(lambda: round(total.percent_aligned(), 2))
    phix_error_rate = _safe(lambda: round(total.error_rate(), 2))

    concerns = ''
    if isinstance(percent_occupied, (int, float)) and percent_occupied > 90:
        concerns = 'potentially overloaded'

    params = _read_run_parameters(run_dir)

    return {
        'RunID':               run_id,
        'Run_nr':              params['Run_nr'],
        'RunStart':            params['RunStart'],
        'RunLength':           params['RunLength'],
        'Flowcell':            params['Flowcell'],
        'Cycles':              params['Cycles'],
        'TotalReadsPF (M)':    _read_total_reads(run_dir),
        'Yield (Gb)':          yield_g,
        'Q30 (%)':             percent_q30,
        'PF (%)':              _safe(lambda: round(read1_lane.percent_pf().mean(), 2)),
        'PF_median (%)':       _safe(lambda: round(read1_lane.percent_pf().median(), 2)),
        'PF_std':              _safe(lambda: round(read1_lane.percent_pf().stddev(), 2)),
        'Occupied (%)':        percent_occupied,
        'Occupied_median (%)': _safe(lambda: round(read1_lane.percent_occupied().median(), 2)),
        'Occupied_std':        _safe(lambda: round(read1_lane.percent_occupied().stddev(), 2)),
        'LoadingConc (pM)':    loading_conc,
        'PhixAligned (%)':     phix_aligned,
        'PhixError (%)':       phix_error_rate,
        'Phasing (%)':         _safe(lambda: round(read1_lane.phasing().mean(), 4)),
        'Prephasing (%)':      _safe(lambda: round(read1_lane.prephasing().mean(), 4)),
        'FirstCycleIntensity': _safe(lambda: round(read1_lane.first_cycle_intensity().mean(), 0)),
        'Concerns':            concerns,
        'SampleSheet':         params['SampleSheet'],
        'SerialNumber':        params['SerialNumber'],
    }


def plot_occupied_vs_pf(run_dir, run_id, out_path):
    """Per-tile % Occupied vs % Pass Filter scatter.

    Each point is one tile. The cloud shape reveals loading status:
      - Bottom-left diagonal: underloaded
      - Central cloud: optimal
      - Near-vertical at high %Occ: overloaded

    See https://knowledge.illumina.com/instrumentation/general/instrumentation-general-troubleshooting-list/000002308
    Returns the written path, or None when the run has no tile data (which is
    the case for any non-patterned flowcell).
    """
    run_dir = os.path.abspath(run_dir)

    run_metrics = py_interop_run_metrics.run_metrics()
    run_metrics.read(run_dir)

    # per-tile %PF from TileMetricSet
    pf_by_tile = {}
    tile_vec = run_metrics.tile_metric_set().metrics()
    for i in range(tile_vec.size()):
        m = tile_vec[i]
        pf_by_tile[(m.lane(), m.tile())] = m.percent_pf()

    # per-tile %Occupied from ExtendedTileMetricSet
    occ_by_tile = {}
    ext_vec = run_metrics.extended_tile_metric_set().metrics()
    for i in range(ext_vec.size()):
        m = ext_vec[i]
        occ_by_tile[(m.lane(), m.tile())] = m.percent_occupied()

    # Join on (lane, tile) and drop NaN / non-finite values
    common = set(pf_by_tile) & set(occ_by_tile)
    occ = np.array([occ_by_tile[k] for k in common], dtype=float)
    pf = np.array([pf_by_tile[k] for k in common], dtype=float)

    valid = np.isfinite(occ) & np.isfinite(pf)
    occ, pf = occ[valid], pf[valid]

    if len(occ) == 0:
        print(f'  no valid tile data for {run_id}, skipping occupancy plot', file=sys.stderr)
        return None

    med_occ = float(np.median(occ))
    med_pf = float(np.median(pf))

    fig, ax = plt.subplots(figsize=(7, 6))

    # Background zones (Illumina loading guidance for patterned flow cells)
    ax.axvspan(0, 60, alpha=0.07, color='steelblue', zorder=0)
    ax.axvspan(60, 85, alpha=0.07, color='green', zorder=0)
    ax.axvspan(85, 100, alpha=0.07, color='firebrick', zorder=0)

    for x, label in [(30, 'Underloaded'), (72.5, 'Optimal'), (92.5, 'Overloaded')]:
        ax.text(x, 99.5, label, ha='center', va='top', fontsize=7.5,
                color='dimgray', style='italic')

    ax.scatter(occ, pf, s=5, alpha=0.45, color='steelblue', linewidths=0, zorder=2)

    ax.axvline(med_occ, color='black', lw=1, ls='--', alpha=0.65, zorder=3)
    ax.axhline(med_pf, color='black', lw=1, ls='--', alpha=0.65, zorder=3)

    ax.text(0.02, 0.7,
            f'Median % Occupied : {med_occ:.1f}%\n'
            f'Median % PF       : {med_pf:.1f}%\n'
            f'Tiles             : {len(occ)}',
            transform=ax.transAxes, va='top', ha='left', fontsize=8.5,
            fontfamily='monospace',
            bbox=dict(boxstyle='round,pad=0.45', facecolor='white',
                      alpha=0.85, edgecolor='lightgray'))

    ax.set_xlabel('% Occupied', fontsize=11)
    ax.set_ylabel('% Pass Filter', fontsize=11)
    ax.set_title(f'% Occupied vs % Pass Filter — {run_id}', fontsize=12, fontweight='bold')
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.grid(True, alpha=0.25, zorder=1)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    return out_path


def write_metrics_csv(metrics, path):
    with open(path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=LOG_COLUMNS)
        writer.writeheader()
        writer.writerow(metrics)


def write_multiqc_table(metrics, path, run_id):
    """MultiQC custom content: one-row table keyed on the run id.

    MultiQC has no module for these metrics, so they are injected as custom
    content rather than parsed.
    """
    header = [
        '# id: illumina_run_metrics',
        "# section_name: 'Illumina Run Metrics'",
        "# description: 'Run-level metrics from InterOp, RunParameters.xml and"
        " PrimaryAnalysisMetrics.csv. Mirrors the Automate-Seq-Run-Metrics-Collection log.'",
        '# plot_type: table',
        '# pconfig:',
        '#     id: illumina_run_metrics_table',
        "#     namespace: 'nf-seQC'",
    ]
    cols = [c for c in LOG_COLUMNS if c != 'RunID']
    with open(path, 'w', newline='') as fh:
        fh.write('\n'.join(header) + '\n')
        fh.write('Run\t' + '\t'.join(cols) + '\n')
        fh.write(run_id + '\t' + '\t'.join(str(metrics[c]) for c in cols) + '\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--run-dir', required=True, help='Illumina run output directory')
    ap.add_argument('--run-id', required=True, help='Run identifier used in output filenames')
    ap.add_argument('--out-prefix', required=True, help='Prefix for the emitted files')
    ap.add_argument('--no-plot', action='store_true', help='Skip the %%Occupied vs %%PF plot')
    args = ap.parse_args()

    metrics = collect_run_metrics(args.run_dir, args.run_id)

    write_metrics_csv(metrics, f'{args.out_prefix}.run_metrics.csv')
    write_multiqc_table(metrics, f'{args.out_prefix}_run_metrics_mqc.tsv', args.run_id)

    if not args.no_plot:
        plot_occupied_vs_pf(args.run_dir, args.run_id,
                            f'{args.out_prefix}_occupied_vs_pf_mqc.png')

    for key in LOG_COLUMNS:
        print(f'{key:22}: {metrics[key]}')


if __name__ == '__main__':
    main()
