#!/usr/bin/env python3
"""
Extract Oxford Nanopore (MinKNOW / PromethION) sample metadata from one or
more run folders into one CSV.

Companion to pacbio_sample_metadata.py. Standard library only.

MinKNOW writes one folder per acquisition:
    <experiment>/<sample>/<date>_<time>_<position>_<flowcell>_<id>/
holding sample_sheet_*.csv, final_summary_*.txt and pod5_pass/. Any folder
holding one of those counts as an acquisition, so pass the experiment folder,
an acquisition folder, or a glob of either.

One row per barcode:
  - from the sample sheet, when it has a `barcode` column (BarcodeSource =
    sample_sheet; Name = alias, else barcode);
  - otherwise from the pod5_pass/barcodeNN folders MinKNOW split into
    (BarcodeSource = pod5_folder). These are only the barcodes MinKNOW saw
    live; many reads sit in pod5_pass/mixed, so this is not a read count;
  - otherwise one row for the acquisition (BarcodeSource = none).

Acquisitions that were restarted on the same flow cell each get their own rows;
the Acquisition column tells them apart.

Output columns:
Name, Barcode, Alias, SampleId, ExperimentId, FlowCellId, FlowCellProductCode,
Kit, PositionId, Instrument, ProtocolRunId, Protocol, Started, Acquisition,
BarcodeSource, SampleSheet
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import sys
from pathlib import Path

OUTPUT_COLUMNS = [
    "Name",
    "Barcode",
    "Alias",
    "SampleId",
    "ExperimentId",
    "FlowCellId",
    "FlowCellProductCode",
    "Kit",
    "PositionId",
    "Instrument",
    "ProtocolRunId",
    "Protocol",
    "Started",
    "Acquisition",
    "BarcodeSource",
    "SampleSheet",
]

MARKERS = ("sample_sheet_*.csv", "final_summary_*.txt")


def is_acquisition(path: Path) -> bool:
    return any(next(path.glob(m), None) for m in MARKERS) or (path / "pod5_pass").is_dir()


def find_acquisitions(root: Path) -> list[Path]:
    """Walk down from root, stopping at each acquisition: never lists the POD5 folders."""
    found = []
    for dirpath, dirnames, _files in os.walk(root):
        path = Path(dirpath)
        if is_acquisition(path):
            found.append(path)
            dirnames.clear()
        else:
            dirnames.sort()
    return found


def read_final_summary(acq: Path) -> dict[str, str]:
    summary = next(iter(sorted(acq.glob("final_summary_*.txt"))), None)
    if summary is None:
        return {}
    values = {}
    for line in summary.read_text(errors="replace").splitlines():
        key, sep, value = line.partition("=")
        if sep:
            values[key.strip()] = value.strip()
    return values


def kit_from_protocol(protocol: str) -> str:
    """'sequencing/...:FLO-PRO114M:SQK-NBD114-24:400' -> 'SQK-NBD114-24'."""
    return next((p for p in protocol.split(":") if p.startswith(("SQK-", "EXP-", "VSK-"))), "")


def read_sample_sheet(acq: Path) -> tuple[Path | None, list[dict[str, str]]]:
    sheet = next(iter(sorted(acq.glob("sample_sheet_*.csv"))), None)
    if sheet is None:
        return None, []
    with open(sheet, newline="", encoding="utf-8-sig") as fh:
        return sheet, [{k.strip(): (v or "").strip() for k, v in r.items() if k} for r in csv.DictReader(fh)]


def extract_samples(acq: Path) -> list[dict[str, str]]:
    summary = read_final_summary(acq)
    sheet, sheet_rows = read_sample_sheet(acq)
    first = sheet_rows[0] if sheet_rows else {}
    protocol = summary.get("protocol", "")

    base = {
        "SampleId": first.get("sample_id") or summary.get("sample_id", ""),
        "ExperimentId": first.get("experiment_id") or summary.get("protocol_group_id", ""),
        "FlowCellId": first.get("flow_cell_id") or summary.get("flow_cell_id", ""),
        "FlowCellProductCode": first.get("flow_cell_product_code", ""),
        "Kit": first.get("kit") or kit_from_protocol(protocol),
        "PositionId": first.get("position_id") or summary.get("position", ""),
        "Instrument": summary.get("instrument", ""),
        "ProtocolRunId": first.get("protocol_run_id") or summary.get("protocol_run_id", ""),
        "Protocol": protocol,
        "Started": summary.get("started", ""),
        "Acquisition": acq.name,
        "SampleSheet": str(sheet) if sheet else "",
    }

    rows = []
    if any(r.get("barcode") for r in sheet_rows):
        for r in sheet_rows:
            barcode, alias = r.get("barcode", ""), r.get("alias", "")
            rows.append({**base, "Name": alias or barcode, "Barcode": barcode, "Alias": alias,
                         "SampleId": r.get("sample_id") or base["SampleId"],
                         "BarcodeSource": "sample_sheet"})
        return rows

    folders = sorted(p.name for p in (acq / "pod5_pass").glob("barcode*") if p.is_dir())
    for barcode in folders:
        rows.append({**base, "Name": barcode, "Barcode": barcode, "Alias": "", "BarcodeSource": "pod5_folder"})
    if not rows:
        rows.append({**base, "Name": base["SampleId"], "Barcode": "", "Alias": "", "BarcodeSource": "none"})
    return rows


def resolve_inputs(inputs: list[str]) -> list[Path]:
    paths: list[Path] = []
    for item in inputs:
        p = Path(item)
        paths.extend([p] if p.exists() else sorted(Path(x) for x in glob.glob(item, recursive=True)))

    acquisitions, seen = [], set()
    for p in paths:
        if p.is_file():
            p = p.parent
        if not p.is_dir():
            continue
        found = find_acquisitions(p)
        if not found:
            print(f"WARNING: no MinKNOW acquisition found in {p}", file=sys.stderr)
        for acq in found:
            key = str(acq.resolve())
            if key not in seen:
                seen.add(key)
                acquisitions.append(acq)
    return acquisitions


def main():
    parser = argparse.ArgumentParser(
        description="Extract Oxford Nanopore sample metadata from MinKNOW run folders into one CSV."
    )
    parser.add_argument("inputs", nargs="+", help="Experiment or acquisition folders, and/or glob patterns")
    parser.add_argument("-o", "--output", default="nanopore_samples_combined.csv", help="Combined output CSV")
    args = parser.parse_args()

    acquisitions = resolve_inputs(args.inputs)
    if not acquisitions:
        raise SystemExit("No MinKNOW acquisition folders found.")

    rows, failures = [], 0
    for acq in acquisitions:
        try:
            samples = extract_samples(acq)
            rows.extend(samples)
            print(f"{acq}: {len(samples)} sample(s)")
        except (OSError, csv.Error, UnicodeDecodeError) as exc:
            failures += 1
            print(f"WARNING: skipped {acq}: {exc}")

    if not rows:
        raise SystemExit("No sample records were extracted from the input folders.")

    rows.sort(key=lambda r: (r["FlowCellId"], r["Acquisition"]))  # stable within an acquisition
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=OUTPUT_COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nAcquisitions read: {len(acquisitions)}")
    print(f"Samples written: {len(rows)}")
    print(f"Output: {args.output}")
    if failures:
        print(f"Skipped folders: {failures}")


if __name__ == "__main__":
    main()
