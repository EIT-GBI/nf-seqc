#!/usr/bin/env python3
"""
Extract Illumina (NextSeq 1000/2000) sample metadata from one or more run
folders into one CSV, one row per sample sheet row.

Companion to pacbio_sample_metadata.py. Standard library only.

Where the sample sheet is looked for, per run folder, first match wins:
  1. Analysis/<n>/Data/Reports/SampleSheet.csv   (the sheet BCL Convert used;
                                                   highest <n> = latest re-analysis)
  2. Analysis/<n>/Data/SampleSheet.csv
  3. any *.csv in the run root with a [BCLConvert_Data] or [Data] section
     (the instrument copies the user's sheet there under its own name,
     e.g. SEQN2K000112_samplesheet_2x50.csv)

Run-level columns come from RunInfo.xml and RunParameters.xml in the run folder.
Both v2 ([BCLConvert_Data], with optional [Cloud_Data]) and v1 ([Data])
sample sheets are read.

Output columns:
Name, SampleName, Project, Lane, Index, Index2, I7_Index_ID, I5_Index_ID,
LibraryName, LibraryPrepKit, IndexAdapterKit, ExperimentName, RunId, Flowcell,
FlowcellMode, CartridgeMode, ReadCycles, InstrumentId, InstrumentName,
InstrumentType, RunStartTime, SampleSheet
"""

from __future__ import annotations

import argparse
import csv
import glob
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

OUTPUT_COLUMNS = [
    "Name",
    "SampleName",
    "Project",
    "Lane",
    "Index",
    "Index2",
    "I7_Index_ID",
    "I5_Index_ID",
    "LibraryName",
    "LibraryPrepKit",
    "IndexAdapterKit",
    "ExperimentName",
    "RunId",
    "Flowcell",
    "FlowcellMode",
    "CartridgeMode",
    "ReadCycles",
    "InstrumentId",
    "InstrumentName",
    "InstrumentType",
    "RunStartTime",
    "SampleSheet",
]

DATA_SECTIONS = ("BCLConvert_Data", "Data")


def read_sections(path: Path) -> dict[str, list[list[str]]]:
    """Split an Illumina sample sheet into {section: rows}, trailing empty cells removed."""
    sections: dict[str, list[list[str]]] = {}
    current = None
    with open(path, newline="", encoding="utf-8-sig") as fh:
        for row in csv.reader(fh):
            while row and not row[-1].strip():
                row.pop()
            if not row:
                continue
            match = re.fullmatch(r"\[(.+)\]", row[0].strip())
            if match:
                current = match.group(1).strip()
                sections.setdefault(current, [])
            elif current is not None:
                sections[current].append([cell.strip() for cell in row])
    return sections


def as_records(rows: list[list[str]]) -> list[dict[str, str]]:
    """Table section (first row = header) -> list of dicts."""
    if not rows:
        return []
    header = rows[0]
    return [dict(zip(header, row + [""] * (len(header) - len(row)))) for row in rows[1:]]


def as_keyvalues(rows: list[list[str]]) -> dict[str, str]:
    return {row[0]: (row[1] if len(row) > 1 else "") for row in rows}


def has_data_section(path: Path) -> bool:
    try:
        return any(s in read_sections(path) for s in DATA_SECTIONS)
    except (OSError, UnicodeDecodeError, csv.Error):
        return False


def analysis_number(path: Path) -> int:
    for part in reversed(path.parts):
        if part.isdigit():
            return int(part)
    return -1


def find_sample_sheet(run_dir: Path) -> Path | None:
    for pattern in ("Analysis/*/Data/Reports/SampleSheet.csv", "Analysis/*/Data/SampleSheet.csv"):
        sheets = sorted(run_dir.glob(pattern), key=analysis_number)
        if sheets:
            return sheets[-1]
    root_sheets = [p for p in sorted(run_dir.glob("*.csv")) if has_data_section(p)]
    return root_sheets[0] if root_sheets else None


def find_run_dir(sheet: Path) -> Path | None:
    """Walk up from a sample sheet to the folder that holds RunInfo.xml."""
    for parent in [sheet.parent, *sheet.parents]:
        if (parent / "RunInfo.xml").is_file():
            return parent
    return None


def xml_text(root, tag: str) -> str:
    if root is None:
        return ""
    node = root.find(f".//{tag}")
    return (node.text or "").strip() if node is not None else ""


def parse_xml(path: Path):
    try:
        return ET.parse(path).getroot() if path.is_file() else None
    except ET.ParseError as exc:
        print(f"WARNING: could not parse {path}: {exc}", file=sys.stderr)
        return None


def run_info(run_dir: Path | None) -> dict[str, str]:
    if run_dir is None:
        return {}
    info = parse_xml(run_dir / "RunInfo.xml")
    params = parse_xml(run_dir / "RunParameters.xml")

    run = info.find(".//Run") if info is not None else None
    cycles = []
    if info is not None:
        for read in info.iter("Read"):
            prefix = "I" if read.get("IsIndexedRead") == "Y" else "R"
            cycles.append(f"{prefix}{read.get('NumCycles', '')}")

    return {
        "ExperimentName": xml_text(params, "ExperimentName"),
        "RunId": run.get("Id", "") if run is not None else "",
        "Flowcell": xml_text(info, "Flowcell") or xml_text(params, "FlowCellSerialNumber"),
        "FlowcellMode": xml_text(params, "FlowCellMode"),
        "CartridgeMode": xml_text(params, "CartridgeMode"),
        "ReadCycles": "-".join(cycles),
        "InstrumentId": xml_text(info, "Instrument") or xml_text(params, "InstrumentSerialNumber"),
        "InstrumentName": xml_text(params, "InstrumentName"),
        "InstrumentType": xml_text(params, "InstrumentType"),
        "RunStartTime": xml_text(params, "RunStartTime"),
    }


def extract_samples(sheet: Path, run_dir: Path | None) -> list[dict[str, str]]:
    sections = read_sections(sheet)
    data_section = next((s for s in DATA_SECTIONS if s in sections), None)
    if data_section is None:
        raise ValueError("no [BCLConvert_Data] or [Data] section")

    header = as_keyvalues(sections.get("Header", []))
    cloud = {r.get("Sample_ID", ""): r for r in as_records(sections.get("Cloud_Data", []))}
    run = run_info(run_dir)
    if not run.get("ExperimentName"):
        run["ExperimentName"] = header.get("RunName", "") or header.get("Experiment Name", "")
    if not run.get("InstrumentType"):
        run["InstrumentType"] = header.get("InstrumentPlatform", "") or header.get("InstrumentType", "")

    rows = []
    for sample in as_records(sections[data_section]):
        sample_id = sample.get("Sample_ID", "")
        if not sample_id:
            continue
        extra = cloud.get(sample_id, {})
        rows.append(
            {
                "Name": sample_id,
                "SampleName": sample.get("Sample_Name", ""),
                "Project": sample.get("Sample_Project", "") or extra.get("ProjectName", ""),
                "Lane": sample.get("Lane", ""),
                "Index": sample.get("Index", "") or sample.get("index", ""),
                "Index2": sample.get("Index2", "") or sample.get("index2", ""),
                "I7_Index_ID": sample.get("I7_Index_ID", ""),
                "I5_Index_ID": sample.get("I5_Index_ID", ""),
                "LibraryName": extra.get("LibraryName", ""),
                "LibraryPrepKit": extra.get("LibraryPrepKitName", ""),
                "IndexAdapterKit": extra.get("IndexAdapterKitName", ""),
                **run,
                "SampleSheet": str(sheet),
            }
        )
    return rows


def resolve_inputs(inputs: list[str]) -> list[tuple[Path, Path | None]]:
    """Each input (run folder, sample sheet, or glob of either) -> (sheet, run_dir)."""
    paths: list[Path] = []
    for item in inputs:
        p = Path(item)
        paths.extend([p] if p.exists() else sorted(Path(x) for x in glob.glob(item, recursive=True)))

    pairs, seen = [], set()
    for p in paths:
        if p.is_dir():
            sheet = find_sample_sheet(p)
            if sheet is None:
                print(f"WARNING: no sample sheet found in {p}", file=sys.stderr)
                continue
            run_dir = p if (p / "RunInfo.xml").is_file() else find_run_dir(sheet)
        elif p.is_file():
            sheet, run_dir = p, find_run_dir(p)
        else:
            continue
        key = str(sheet.resolve())
        if key not in seen:
            seen.add(key)
            pairs.append((sheet, run_dir))
    return pairs


def main():
    parser = argparse.ArgumentParser(
        description="Extract Illumina sample metadata from run folders / sample sheets into one CSV."
    )
    parser.add_argument("inputs", nargs="+", help="Run folders, sample sheet CSVs, and/or glob patterns")
    parser.add_argument("-o", "--output", default="illumina_samples_combined.csv", help="Combined output CSV")
    args = parser.parse_args()

    pairs = resolve_inputs(args.inputs)
    if not pairs:
        raise SystemExit("No sample sheets found.")

    rows, failures = [], 0
    for sheet, run_dir in pairs:
        try:
            samples = extract_samples(sheet, run_dir)
            rows.extend(samples)
            print(f"{sheet}: {len(samples)} sample(s)")
        except (OSError, ValueError, csv.Error, UnicodeDecodeError) as exc:
            failures += 1
            print(f"WARNING: skipped {sheet}: {exc}")

    if not rows:
        raise SystemExit("No sample records were extracted from the input files.")

    rows.sort(key=lambda r: (r["RunId"], r["Lane"]))  # stable: keeps sheet order within a run
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=OUTPUT_COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nSample sheets read: {len(pairs)}")
    print(f"Samples written: {len(rows)}")
    print(f"Output: {args.output}")
    if failures:
        print(f"Skipped files: {failures}")


if __name__ == "__main__":
    main()
