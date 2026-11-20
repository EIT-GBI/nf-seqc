#!/usr/bin/env python3
"""
Extract Element AVITI sample metadata from one or more run folders into one
CSV, one row per [SAMPLES] row of RunManifest.csv.

Companion to pacbio_sample_metadata.py. Standard library only.

Per run folder: samples from RunManifest.csv ([SAMPLES] section; lines
starting with # are comments), run-level columns from RunParameters.json.
A sample listed with several index pairs (e.g. PhiX, four pairs) gives one
row per pair, as in the manifest.

FilterMask is included because some runs set their own (e.g. R1-only), which
changes what bases2fastq filters on and makes runs not directly comparable.

Output columns:
Name, Index1, Index2, Lane, Project, RunName, RunFolderName, RunDescription,
Side, FlowcellID, InstrumentName, Date, KitConfiguration, ChemistryVersion,
PreparationWorkflow, LibraryType, ThroughputSelection, ReadCycles, FilterMask,
RunManifest
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import re
import sys
from pathlib import Path

OUTPUT_COLUMNS = [
    "Name",
    "Index1",
    "Index2",
    "Lane",
    "Project",
    "RunName",
    "RunFolderName",
    "RunDescription",
    "Side",
    "FlowcellID",
    "InstrumentName",
    "Date",
    "KitConfiguration",
    "ChemistryVersion",
    "PreparationWorkflow",
    "LibraryType",
    "ThroughputSelection",
    "ReadCycles",
    "FilterMask",
    "RunManifest",
]

RUN_KEYS = [
    "RunName",
    "RunFolderName",
    "RunDescription",
    "Side",
    "FlowcellID",
    "InstrumentName",
    "Date",
    "KitConfiguration",
    "ChemistryVersion",
    "PreparationWorkflow",
    "LibraryType",
    "ThroughputSelection",
    "FilterMask",
]


def read_sections(path: Path) -> dict[str, list[list[str]]]:
    """Split RunManifest.csv into {SECTION: rows}; comments and blank cells dropped."""
    sections: dict[str, list[list[str]]] = {}
    current = None
    with open(path, newline="", encoding="utf-8-sig") as fh:
        for row in csv.reader(fh):
            while row and not row[-1].strip():
                row.pop()
            if not row or row[0].lstrip().startswith("#"):
                continue
            match = re.fullmatch(r"\[(.+)\]", row[0].strip())
            if match:
                current = match.group(1).strip().upper()
                sections.setdefault(current, [])
            elif current is not None:
                sections[current].append([cell.strip() for cell in row])
    return sections


def run_parameters(run_dir: Path) -> dict[str, str]:
    path = run_dir / "RunParameters.json"
    if not path.is_file():
        print(f"WARNING: no RunParameters.json in {run_dir}", file=sys.stderr)
        return {}
    try:
        params = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        print(f"WARNING: could not read {path}: {exc}", file=sys.stderr)
        return {}

    run = {key: str(params.get(key) or "") for key in RUN_KEYS}
    cycles = params.get("Cycles") or {}
    if isinstance(cycles, dict):
        order = [r.strip() for r in str(params.get("ReadOrder", "")).split(",") if r.strip()] or list(cycles)
        run["ReadCycles"] = "-".join(f"{r}:{cycles[r]}" for r in order if r in cycles)
    return run


def extract_samples(manifest: Path) -> list[dict[str, str]]:
    sections = read_sections(manifest)
    table = sections.get("SAMPLES")
    if not table:
        raise ValueError("no [SAMPLES] section")

    header = table[0]
    run = run_parameters(manifest.parent)
    rows = []
    for values in table[1:]:
        sample = dict(zip(header, values + [""] * (len(header) - len(values))))
        name = sample.get("SampleName", "")
        if not name:
            continue
        rows.append(
            {
                "Name": name,
                "Index1": sample.get("Index1", ""),
                "Index2": sample.get("Index2", ""),
                "Lane": sample.get("Lane", ""),
                "Project": sample.get("Project", ""),
                **run,
                "RunManifest": str(manifest),
            }
        )
    return rows


def resolve_inputs(inputs: list[str]) -> list[Path]:
    """Run folders, RunManifest.csv files, or globs of either -> manifests."""
    paths: list[Path] = []
    for item in inputs:
        p = Path(item)
        paths.extend([p] if p.exists() else sorted(Path(x) for x in glob.glob(item, recursive=True)))

    manifests, seen = [], set()
    for p in paths:
        manifest = p / "RunManifest.csv" if p.is_dir() else p
        if not manifest.is_file():
            print(f"WARNING: no RunManifest.csv in {p}", file=sys.stderr)
            continue
        key = str(manifest.resolve())
        if key not in seen:
            seen.add(key)
            manifests.append(manifest)
    return manifests


def main():
    parser = argparse.ArgumentParser(
        description="Extract Element AVITI sample metadata from run folders into one CSV."
    )
    parser.add_argument("inputs", nargs="+", help="Run folders, RunManifest.csv files, and/or glob patterns")
    parser.add_argument("-o", "--output", default="aviti_samples_combined.csv", help="Combined output CSV")
    args = parser.parse_args()

    manifests = resolve_inputs(args.inputs)
    if not manifests:
        raise SystemExit("No RunManifest.csv files found.")

    rows, failures = [], 0
    for manifest in manifests:
        try:
            samples = extract_samples(manifest)
            rows.extend(samples)
            print(f"{manifest}: {len(samples)} sample row(s)")
        except (OSError, ValueError, csv.Error, UnicodeDecodeError) as exc:
            failures += 1
            print(f"WARNING: skipped {manifest}: {exc}")

    if not rows:
        raise SystemExit("No sample records were extracted from the input files.")

    rows.sort(key=lambda r: r.get("RunName", ""))  # stable: keeps manifest order within a run
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=OUTPUT_COLUMNS, restval="", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nManifests read: {len(manifests)}")
    print(f"Sample rows written: {len(rows)}")
    print(f"Output: {args.output}")
    if failures:
        print(f"Skipped files: {failures}")


if __name__ == "__main__":
    main()
