#!/usr/bin/env python3
"""MultiQC sections for the PacBio Revio SMRT Cells of a run, from the instrument's reports.

Reads what SMRT Link's run QC report is drawn from, in the cell's statistics/
directory; never the reads:

* <movie>.reports.zip: pbreports' JSON reports (raw_data, loading, control,
  ccs, methylation, barcodes) and the binned per-barcode data behind its
  barcode heat maps.
* <movie>.zmw_metrics.csv.gz: per-ZMW read length, polymerase length, read
  quality, HiFi flag and number of passes, for the subset of ZMWs the
  instrument writes it for. Distributions are shares, so the subset is enough.

Writes MultiQC custom-content JSON, one file per section, into --mqc-dir.
Cell-level series and rows are named after the cell (its well, as the Run
Metrics table names it); per-barcode series after the pipeline's samples.
--cells lists every cell of the run, as JSON:

    [{"cell": "1_A01", "reports_zip": "<file>", "zmw_metrics": "<file>",
      "samples": {"<barcode>": "<sample id>"}}, ...]

(either file may be null). All cells go in one call because the
distributions are bar charts with one row per bin and one bar per cell or
sample: MultiQC merges custom-content files by row name, so cells written to
separate files would overwrite each other's rows. MultiQC draws bars
horizontally only, so the bins run down the page.
"""

import argparse
import collections
import concurrent.futures
import csv
import gzip
import json
import math
import os
import re
import zipfile

LEN_BIN = 1000       # HiFi read length, bp
LEN_MAX = 40000      # longer HiFi reads are pooled into the last bin
POL_BIN = 10000      # polymerase read length, bp
POL_MAX = 200000     # longer polymerase reads are pooled into the last bin
HIFI_QV = 20         # HiFi: a CCS read at QV 20 or above
MAX_QV = 60          # pbreports caps QV at 60
QV_BIN = 2           # read quality bins, QV
PASS_BIN = 2         # number-of-passes bins
PASS_MAX = 61        # more passes are pooled into the last bin


def safe_id(name):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name)


def barcode_key(name):
    """bc2053--bc2053 -> bc2053, as barcodeKey() in pacbio_qc.nf."""
    parts = name.split("--")
    return parts[0] if len(parts) == 2 and parts[0] == parts[1] else name


# Two report groups (MultiQC modules), each section in one. Their order, and
# the sections' order within them, is set in assets/multiqc_config.yml.
PARENTS = {
    "pacbio_cell": (
        "PacBio SMRT Cell",
        "How each SMRT Cell ran, from pbreports (reports.zip, what SMRT Link's run QC report shows): loading, "
        "polymerase reads, CCS details, the sequencing control and base modifications.",
    ),
    "pacbio_reads": (
        "PacBio HiFi reads",
        "HiFi read length, quality and passes per SMRT Cell (zmw_metrics.csv.gz), and per barcoded sample "
        "(pbreports' binned barcode data).",
    ),
}
SECTION_PARENT = {
    "pacbio_cell_details": "pacbio_cell", "pacbio_loading": "pacbio_cell",
    "pacbio_polymerase_length": "pacbio_cell", "pacbio_methylation": "pacbio_cell",
}


def section(sid, name, description, plot_type, data, pconfig=None, **extra):
    parent = SECTION_PARENT.get(sid, "pacbio_reads")
    d = {
        "id": sid,
        "parent_id": parent,
        "parent_name": PARENTS[parent][0],
        "parent_description": PARENTS[parent][1],
        "section_name": name,
        "description": description,
        "plot_type": plot_type,
        "pconfig": {"id": f"{sid}_plot", "title": f"PacBio: {name}", **(pconfig or {})},
        "data": data,
    }
    d.update(extra)
    return d



def as_pairs(section_dict):
    """Line-graph series as sorted [x, y] pairs. MultiQC orders a series given
    as a JSON object by its keys as text ("1000" before "126"), which draws a
    zigzag; pairs keep numeric order."""
    if section_dict["plot_type"] == "linegraph":
        section_dict["data"] = {
            name: sorted([float(x), y] for x, y in series.items() if y is not None)
            for name, series in section_dict["data"].items()
        }
    return section_dict

def kb(bp):
    return f"{bp / 1000:g}"


def len_label(b, step, cap):
    return f"{kb(cap)}+ kb" if b >= cap else f"{kb(b)}-{kb(b + step)} kb"


def qv_label(b):
    return f"QV {MAX_QV}" if b >= MAX_QV else f"QV {b}-{b + QV_BIN - 1}"


def pass_label(b):
    return f"{PASS_MAX}+" if b >= PASS_MAX else f"{b}-{b + PASS_BIN - 1}"


def bars(series, label, step=None, start=None, missing=0):
    """Bar-chart data with one row per bin: {row label: {category: value}}.

    series is {category: {bin: value}}, a category being a cell or a sample.
    Rows run from the smallest bin to the largest; with step, every bin in
    between gets a row even where no category has a value, so a gap shows as
    a gap. start is the first bin to show (default the smallest with data).
    missing is a category's value in a bin it has none for: 0 for a share of
    reads; None (no bar) for a value that is not a count, such as a mean."""
    bins = sorted({b for values in series.values() for b in values})
    if not bins:
        return {}
    if step:
        lo = bins[0] if start is None else min(start, bins[0])
        bins = list(range(int(lo), int(bins[-1]) + 1, step))
    return {
        label(b): {c: values.get(b, missing) for c, values in series.items() if b in values or missing is not None}
        for b in bins
    }


def bar_pconfig(xlab, ylab, **extra):
    """A distribution as bars: one row per bin, one bar per cell or sample side
    by side, in bin order (no sample sorting), values as given (no % switch).
    MultiQC labels a bar chart's value axis ylab, which is the horizontal one."""
    return {"ylab": ylab, "xlab": xlab, "cpswitch": False, "stacking": "group", "sort_samples": False, "ymin": 0, **extra}


def pct(counter):
    total = sum(counter.values())
    return {k: 100.0 * counter[k] / total for k in sorted(counter)} if total else {}


def read_zip(path):
    """{entry name: parsed JSON} for the reports.zip entries used here."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    with zipfile.ZipFile(path) as z:
        for name in z.namelist():
            base = os.path.basename(name)
            if "/" in name.strip("/"):
                continue  # per-collection copies in sub-folders repeat the top level
            if base.endswith(".report.json"):
                out[base] = json.loads(z.read(name))
            elif base.endswith(".json.gz"):
                out[base] = json.loads(gzip.decompress(z.read(name)))
    return out


def attributes(report):
    return {a["id"].split(".")[-1]: a.get("value") for a in (report or {}).get("attributes") or []}


def table(report, table_id):
    for t in (report or {}).get("tables") or []:
        if t.get("id") == table_id:
            cols = {c["id"].split(".")[-1]: c.get("values") or [] for c in t.get("columns") or []}
            n = max((len(v) for v in cols.values()), default=0)
            return [{k: (v[i] if i < len(v) else None) for k, v in cols.items()} for i in range(n)]
    return []


def qv(value):
    """A read quality as QV: pbreports writes accuracy (0-1) or a QV."""
    if value is None or math.isnan(value) or value < 0:
        return None
    if value <= 1:
        q = MAX_QV if value >= 1 else min(MAX_QV, -10 * math.log10(1 - value))
    else:
        q = min(MAX_QV, value)
    # Rounded, so accuracy 0.99 is QV 20 (not 19.999...) and counts as HiFi,
    # and boundary values land in their own QV bin.
    return round(q, 6)


def number(value):
    """A float, or None for anything missing, unparsable or NaN."""
    try:
        f = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(f) else f


def zmw_distributions(path):
    """Per-ZMW distributions of one cell from zmw_metrics.csv.gz, or None."""
    if not path or not os.path.exists(path):
        return None
    hifi_len, hifi_bases, pol_len = collections.Counter(), collections.Counter(), collections.Counter()
    quals, passes = collections.Counter(), collections.Counter()
    q_sum, q_n = collections.Counter(), collections.Counter()
    with gzip.open(path, "rt", newline="") as handle:
        for row in csv.DictReader(handle):
            length = int(number(row.get("read_length")) or 0)
            pol = int(number(row.get("pol_length")) or 0)
            q = qv(number(row.get("read_qual")))
            # HiFi by quality, as pbreports counts it; the file's own hifi
            # column is set for ZMWs without a read too.
            hifi = q is not None and q >= HIFI_QV
            n = int(number(row.get("num_passes")) or 0)
            if pol > 0:
                pol_len[min(pol, POL_MAX) // POL_BIN * POL_BIN] += 1
            if q is not None and length > 0:
                quals[min(int(q), MAX_QV) // QV_BIN * QV_BIN] += 1
            if hifi and length > 0:
                b = min(length, LEN_MAX) // LEN_BIN * LEN_BIN
                hifi_len[b] += 1
                hifi_bases[b] += length
                if n > 0:
                    passes[(min(n, PASS_MAX) - 1) // PASS_BIN * PASS_BIN + 1] += 1
                if q is not None:
                    q_sum[b] += q
                    q_n[b] += 1

    total_bases = sum(hifi_bases.values())
    yield_ge, seen = {}, 0
    for b in sorted(hifi_bases, reverse=True):
        seen += hifi_bases[b]
        yield_ge[b] = 100.0 * seen / total_bases
    return {
        "hifi_len": pct(hifi_len),
        "yield_ge": dict(sorted(yield_ge.items())),
        "qv_by_len": {b: q_sum[b] / q_n[b] for b in sorted(q_n) if q_n[b] >= 20},
        "quals": pct(quals),
        "passes": pct(passes),
        "pol_len": pct(pol_len),
    }


def zmw_sections(per_cell):
    """Sections from every cell's zmw_metrics.csv.gz: {cell: zmw_distributions()}."""
    def series(key):
        return {cell: d[key] for cell, d in per_cell.items() if d and d[key]}

    note = ("From zmw_metrics.csv.gz, which the instrument writes for a subset of the cell's ZMWs; "
            "shares, so the subset stands for the cell.")
    out = []
    hifi_len = series("hifi_len")
    if hifi_len:
        out.append(section(
            "pacbio_hifi_length", "HiFi read length distribution",
            f"Share of each cell's HiFi reads (QV {HIFI_QV}+) by read length ({LEN_BIN // 1000} kb bins; "
            f"{LEN_MAX // 1000} kb and longer pooled in the last bin). " + note,
            "bargraph", bars(hifi_len, lambda b: len_label(b, LEN_BIN, LEN_MAX), LEN_BIN, 0),
            bar_pconfig("HiFi read length", "% of HiFi reads"),
        ))
    yield_ge = series("yield_ge")
    if yield_ge:
        # A cumulative curve, as SMRT Link draws it: a line, not bars.
        out.append(section(
            "pacbio_hifi_yield_by_length", "HiFi yield by read length",
            "Share of each cell's HiFi yield in reads at least this long, as SMRT Link's HiFi read length "
            "yield plot. " + note,
            "linegraph", yield_ge,
            {"xlab": "Minimum HiFi read length (bp)", "ylab": "% of HiFi yield", "ymin": 0, "ymax": 100, "xmin": 0},
        ))
    qv_by_len = series("qv_by_len")
    if qv_by_len:
        out.append(section(
            "pacbio_qv_by_length", "HiFi quality by read length",
            f"Mean QV of each cell's HiFi reads in each read-length bin ({LEN_BIN // 1000} kb bins; bins with "
            "fewer than 20 reads left out). " + note,
            "bargraph", bars(qv_by_len, lambda b: len_label(b, LEN_BIN, LEN_MAX), missing=None),
            bar_pconfig("HiFi read length", "Mean QV"),
        ))
    quals = series("quals")
    if quals:
        out.append(section(
            "pacbio_read_quality", "Read quality distribution",
            f"Share of each cell's CCS reads by read quality ({QV_BIN}-QV bins, capped at QV {MAX_QV}), HiFi and "
            "below. HiFi is QV 20 and above. " + note,
            "bargraph", bars(quals, qv_label, QV_BIN),
            bar_pconfig("Read quality", "% of reads"),
        ))
    passes = series("passes")
    if passes:
        out.append(section(
            "pacbio_npasses", "Number of passes",
            f"Share of each cell's HiFi reads by the number of full polymerase passes over the insert "
            f"({PASS_BIN}-pass bins, {PASS_MAX}+ pooled). More passes, higher accuracy; few passes with long "
            "inserts is expected. " + note,
            "bargraph", bars(passes, pass_label, PASS_BIN, 1),
            bar_pconfig("Passes", "% of HiFi reads"),
        ))
    pol_len = series("pol_len")
    if pol_len:
        out.append(section(
            "pacbio_polymerase_length", "Polymerase read length",
            f"Share of each cell's productive ZMWs by polymerase read length ({POL_BIN // 1000} kb bins, "
            f"{POL_MAX // 1000} kb and longer pooled): how far the polymerase got. " + note,
            "bargraph", bars(pol_len, lambda b: len_label(b, POL_BIN, POL_MAX), POL_BIN, 0),
            bar_pconfig("Polymerase read length", "% of ZMWs"),
        ))
    return out


def report_sections(reports, cell):
    out = []
    raw = attributes(reports.get("raw_data.report.json"))
    load = attributes(reports.get("loading.report.json"))
    ctrl = attributes(reports.get("control.report.json"))
    ccs = attributes(reports.get("ccs.report.json"))

    if load:
        out.append(section(
            "pacbio_loading", "ZMW loading",
            "ZMWs per cell by loading, from loading.report.json: P0 empty, P1 productive (one polymerase), "
            "P2 other (e.g. more than one polymerase).",
            "bargraph",
            {cell: {k: load.get(f"productivity_{i}_n") or 0 for i, k in enumerate(["p0", "p1", "p2"])}},
            {"ylab": "ZMWs", "cpswitch_c_active": False},
            categories={
                "p1": {"name": "P1 productive", "color": "#00a651"},
                "p2": {"name": "P2 other", "color": "#f7b538"},
                "p0": {"name": "P0 empty", "color": "#d1d1d1"},
            },
        ))

    row = {}
    if raw:
        row.update({
            "pol_reads": raw.get("nreads"), "pol_bases": (raw.get("nbases") or 0) / 1e9 or None,
            "pol_length": raw.get("read_length"), "pol_n50": raw.get("read_n50"),
            "insert_length": raw.get("insert_length"), "insert_n50": raw.get("insert_n50"),
            "umy": (raw.get("unique_molecular_yield") or 0) / 1e9 or None,
            "base_rate": raw.get("local_base_rate_median"),
        })
    if ccs:
        q30 = number(ccs.get("percent_ccs_bases_q30"))
        adapters = number(ccs.get("percent_missing_adapters"))
        row.update({
            "median_length": ccs.get("median_ccs_readlength"),
            "mean_npasses": ccs.get("mean_npasses"),
            "median_qv": ccs.get("median_qv"),
            # pbreports gives these as fractions despite the names
            "q30": q30 * 100 if q30 is not None else None,
            "missing_adapters": adapters * 100 if adapters is not None else None,
        })
    if ctrl:
        concordance = number(ctrl.get("concordance_mean"))
        row.update({
            "control_reads": ctrl.get("reads_n"), "control_length": ctrl.get("readlength_mean"),
            "control_concordance": concordance * 100 if concordance is not None else None,
        })
    row = {k: v for k, v in row.items() if v is not None}
    if row:
        def col(title, desc, fmt="{:,.0f}", suffix="", hidden=False, **kw):
            h = {"title": title, "description": desc, "format": fmt, "scale": False}
            if suffix:
                h["suffix"] = suffix
            if hidden:
                h["hidden"] = True
            h.update(kw)
            return h
        out.append(section(
            "pacbio_cell_details", "SMRT Cell details",
            "Per-cell figures from pbreports (reports.zip) that the Run Metrics table leaves out: polymerase "
            "(raw) reads, insert sizes, CCS details and the sequencing control. Hidden columns via Configure "
            "columns.",
            "table", {cell: row},
            {"col1_header": "SMRT Cell"},
            headers={
                "pol_reads": col("Polymerase reads", "Productive ZMWs' polymerase reads"),
                "pol_bases": col("Polymerase yield", "Bases in polymerase reads", "{:,.1f}", " Gb"),
                "pol_length": col("Polymerase mean", "Mean polymerase read length", suffix=" bp"),
                "pol_n50": col("Polymerase N50", "Polymerase read length N50", suffix=" bp"),
                "insert_length": col("Insert mean", "Mean insert (subread) length", suffix=" bp"),
                "insert_n50": col("Insert N50", "Insert length N50", suffix=" bp"),
                "umy": col("Unique molecular yield", "Sum of the longest subread of each ZMW", "{:,.1f}", " Gb", hidden=True),
                "base_rate": col("Base rate", "Median local base rate", "{:,.2f}", " b/s", hidden=True),
                "median_length": col("CCS median length", "Median CCS read length", suffix=" bp", hidden=True),
                "mean_npasses": col("Mean passes", "Mean number of passes per CCS read", "{:,.0f}"),
                "median_qv": col("Median QV", "Median CCS read quality", "{:,.0f}"),
                "q30": col("Q30 bases", "Share of CCS bases at Q30 or above", "{:,.1f}", "%"),
                "missing_adapters": col("Missing adapters", "Share of CCS reads with a missing adapter", "{:,.2f}", "%"),
                "control_reads": col("Control reads", "Sequencing control (spike-in) reads", hidden=True),
                "control_length": col("Control length", "Mean control polymerase read length", suffix=" bp", hidden=True),
                "control_concordance": col("Control concordance", "Mean control read concordance (accuracy)", "{:,.1f}", "%"),
            },
        ))

    meth = table(reports.get("methylation.report.json"), "methylation.summary")
    if meth:
        out.append(section(
            "pacbio_methylation", "Base modifications",
            "Share of motif sites the instrument scored for each modification, and the share of scored sites "
            "called modified, from methylation.report.json (Jasmine).",
            "table",
            {f"{cell} {m.get('modification')} {m.get('motif')}": {
                "scored": 100.0 * (number(m.get("frac_sites_scored")) or 0),
                "modified": 100.0 * (number(m.get("frac_sites_modified")) or 0),
            } for m in meth},
            {"col1_header": "Cell, modification, motif", "sort_rows": False},
            headers={
                "scored": {"title": "Sites scored", "suffix": "%", "format": "{:,.1f}", "min": 0, "max": 100},
                "modified": {"title": "Modified", "suffix": "%", "format": "{:,.2f}", "min": 0, "max": 100},
            },
        ))

    return out


# Per barcode: pbreports' binned length and barcode-quality data (its heat
# maps), as one bar per sample in each bin. pbreports' own bins: 1 kb read
# length, 2-point barcode quality.
BARCODE_PLOTS = [
    ("binned_readlength.json.gz", "pacbio_barcode_length", "Read length by sample", "HiFi read length",
     lambda b: f"{kb(b)} kb",
     "Share of each barcoded sample's HiFi reads by read length, from the binned data behind SMRT Link's "
     "barcode read-length heat map (its 1 kb bins)."),
    ("binned_bcqual.json.gz", "pacbio_barcode_quality", "Barcode quality by sample", "Barcode quality score",
     lambda b: f"{b:g}",
     "Share of each barcoded sample's reads by lima barcode quality score (0-100), from the binned data "
     "behind SMRT Link's barcode quality heat map. Good barcodes pile up near 100."),
]


def barcode_series(reports, samples):
    """{zip entry: {sample: {bin: % of reads}}} for one cell."""
    out = {}
    for key, *_ in BARCODE_PLOTS:
        binned = (reports.get(key) or {}).get("data") or []
        if not binned:
            continue
        b = binned[0]
        for i, bc in enumerate(b.get("x") or []):
            sample = samples.get(barcode_key(str(bc)))
            if not sample:
                continue
            column = collections.Counter({y: (row[i] if i < len(row) else 0) for y, row in zip(b.get("y") or [], b.get("z") or [])})
            if sum(column.values()):
                out.setdefault(key, {})[sample] = pct(column)
    return out


def barcode_sections(per_sample):
    """per_sample: {zip entry: {sample: {bin: %}}}, every cell's samples."""
    out = []
    for key, sid, name, xlab, label, desc in BARCODE_PLOTS:
        if per_sample.get(key):
            out.append(section(sid, name, desc, "bargraph", bars(per_sample[key], label),
                               bar_pconfig(xlab, "% of reads")))
    return out


def merge_sections(sections):
    """One section per id: the cells' table and loading rows are keyed by cell,
    so their data dicts merge without clashes."""
    merged = {}
    for sec in sections:
        if sec["id"] in merged:
            merged[sec["id"]]["data"].update(sec["data"])
        else:
            merged[sec["id"]] = sec
    return list(merged.values())


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--cells", required=True, help="JSON list of the run's cells (see the module docstring)")
    parser.add_argument("--mqc-dir", required=True)
    parser.add_argument("--threads", type=int, default=1, help="cells read in parallel")
    args = parser.parse_args()

    cells = json.loads(args.cells)
    # zmw_metrics.csv.gz is the slow part (a CSV row per ZMW in the subset):
    # read the cells' files in parallel processes.
    with concurrent.futures.ProcessPoolExecutor(max_workers=max(1, args.threads)) as pool:
        distributions = list(pool.map(zmw_distributions, [c.get("zmw_metrics") for c in cells]))

    per_cell_sections, zmw, per_sample = [], {}, {}
    for c, dist in zip(cells, distributions):
        cell = c["cell"]
        samples = {barcode_key(k): v for k, v in (c.get("samples") or {}).items()}
        reports = read_zip(c.get("reports_zip"))
        per_cell_sections += report_sections(reports, cell)
        zmw[cell] = dist
        for key, series in barcode_series(reports, samples).items():
            per_sample.setdefault(key, {}).update(series)

    sections = merge_sections(per_cell_sections) + barcode_sections(per_sample) + zmw_sections(zmw)
    os.makedirs(args.mqc_dir, exist_ok=True)
    for s in sections:
        if s["data"]:
            with open(os.path.join(args.mqc_dir, f"{s['id']}_mqc.json"), "w") as handle:
                json.dump(as_pairs(s), handle)


if __name__ == "__main__":
    main()
