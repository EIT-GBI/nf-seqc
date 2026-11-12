#!/usr/bin/env python3
"""
Extract PacBio Revio sample metadata from one or more run folders into one
CSV, one row per BioSample of each SMRT Cell.

Original script by Melissa; reworked to match the other extractors here.
Standard library only.

Revio writes, per SMRT Cell:
    <run>/<cell>/metadata/<movie>.metadata.xml
plus two copies taken while the cell was still acquiring,
<movie>_basic_preview_4hr.metadata.xml and <movie>_full_preview_23hr.metadata.xml,
which list the same samples. Those are skipped when found through a folder or
glob (a file passed by its own path is always read), and rows repeated across
files (e.g. a downloaded `<movie>.metadata(1).xml`) are written once.

Only the first DNABarcode of each BioSample is kept. BarcodeName is its left
barcode (bc2053--bc2053 -> bc2053), as in the original script; BarcodePair is
the full name, which tells asymmetric pairs (bc2001--bc2002) apart.

Output columns:
Name, ReadSet, RunName, MovieName, WellName, ReactionWell, BarcodeName,
application, library_type, loading_concentration, insert_size, is_ccs,
is_isoseq, instrument_type, instrument_id, instrument_name, BarcodePair
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
    "ReadSet",
    "RunName",
    "MovieName",
    "WellName",
    "ReactionWell",
    "BarcodeName",
    "application",
    "library_type",
    "loading_concentration",
    "insert_size",
    "is_ccs",
    "is_isoseq",
    "instrument_type",
    "instrument_id",
    "instrument_name",
    "BarcodePair",
]

# WellSample child element -> output column
WELL_FIELDS = {
    "Application": "application",
    "LibraryType": "library_type",
    "OnPlateLoadingConcentration": "loading_concentration",
    "InsertSize": "insert_size",
    "IsCCS": "is_ccs",
    "IsoSeq": "is_isoseq",
}


# --- XML helpers: Revio files are namespaced, so match on the local tag name ---

def local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def descendants(parent: ET.Element, name: str):
    return (elem for elem in parent.iter() if local_name(elem.tag) == name)


def first_descendant(parent: ET.Element, name: str) -> ET.Element | None:
    return next(descendants(parent, name), None)


def child_text(parent: ET.Element, name: str) -> str:
    for child in parent:
        if local_name(child.tag) == name:
            return (child.text or "").strip()
    return ""


# --- Names ---

def movie_from_filename(filename: str) -> str:
    """m84334_..._s1.metadata.xml / .metadata(1).xml -> m84334_..._s1"""
    return re.sub(r"\.metadata(?:\(\d+\))?\.xml$", "", filename, flags=re.IGNORECASE)


def clean_barcode(barcode_name: str) -> str:
    """bc2033--bc2033 -> bc2033; bc2033 -> bc2033"""
    return barcode_name.split(";", 1)[0].split("--", 1)[0]


def readset_and_reactionwell(xml_file: Path) -> tuple[str, str]:
    """
    <run>/<cell>/metadata/<movie>.metadata.xml -> (<run>, <cell>); else blank.
    absolute(), not resolve(): a run reached through a symlink keeps the name
    it was given, not the link target's.
    """
    path = xml_file.absolute()
    if path.parent.name == "metadata":
        return path.parent.parent.parent.name, path.parent.parent.name
    return "", ""


def is_preview(path: Path) -> bool:
    return "_preview_" in path.name


# --- Extraction ---

def extract_samples(xml_file: Path) -> list[dict[str, str]]:
    root = ET.parse(xml_file).getroot()
    read_set, reaction_well = readset_and_reactionwell(xml_file)

    rows = []
    for run in descendants(root, "Run"):
        for collection in descendants(run, "CollectionMetadata"):
            well = first_descendant(collection, "WellSample")
            biosamples = first_descendant(well, "BioSamples") if well is not None else None
            if biosamples is None:
                continue

            shared = {
                "ReadSet": read_set,
                "RunName": run.get("Name", ""),
                # Context is the movie name; the file name is the fallback.
                "MovieName": collection.get("Context", "") or movie_from_filename(xml_file.name),
                "WellName": well.get("Name", ""),
                "ReactionWell": reaction_well,
                **{column: child_text(well, tag) for tag, column in WELL_FIELDS.items()},
                "instrument_type": run.get("InstrumentType", ""),
                "instrument_id": collection.get("InstrumentId", ""),
                "instrument_name": collection.get("InstrumentName", ""),
            }

            for biosample in biosamples:
                if local_name(biosample.tag) != "BioSample":
                    continue
                barcode = first_descendant(biosample, "DNABarcode")
                pair = barcode.get("Name", "") if barcode is not None else ""
                rows.append(
                    {
                        "Name": biosample.get("Name", ""),
                        "BarcodeName": clean_barcode(pair),
                        "BarcodePair": pair,
                        **shared,
                    }
                )
    return rows


def resolve_inputs(inputs: list[str]) -> list[Path]:
    """Run/cell folders, metadata XML files, or globs of either -> XML files."""
    files: list[Path] = []
    for item in inputs:
        p = Path(item)
        if p.is_file():
            files.append(p)
        elif p.is_dir():
            files.extend(x for x in sorted(p.rglob("*.metadata*.xml")) if not is_preview(x))
        else:
            matches = sorted(Path(x) for x in glob.glob(item, recursive=True))
            files.extend(x for x in matches if x.is_file() and not is_preview(x))

    unique, seen = [], set()
    for f in files:
        key = f.resolve()
        if key not in seen:
            seen.add(key)
            unique.append(f)
    return unique


def main():
    parser = argparse.ArgumentParser(
        description="Extract PacBio Revio sample metadata from metadata.xml files into one CSV."
    )
    parser.add_argument("inputs", nargs="+", help="Run or SMRT Cell folders, metadata XML files, and/or glob patterns")
    parser.add_argument("-o", "--output", default="pacbio_samples_combined.csv", help="Combined output CSV")
    args = parser.parse_args()

    xml_files = resolve_inputs(args.inputs)
    if not xml_files:
        raise SystemExit("No metadata XML files found.")

    rows, seen, failures = [], set(), 0
    for xml_file in xml_files:
        try:
            samples = extract_samples(xml_file)
        except (ET.ParseError, OSError) as exc:
            failures += 1
            print(f"WARNING: skipped {xml_file}: {exc}", file=sys.stderr)
            continue
        new = 0
        for row in samples:
            key = tuple(row[c] for c in OUTPUT_COLUMNS)
            if key not in seen:
                seen.add(key)
                rows.append(row)
                new += 1
        duplicates = f" ({len(samples) - new} already listed)" if new < len(samples) else ""
        print(f"{xml_file}: {len(samples)} sample(s){duplicates}")

    if not rows:
        raise SystemExit("No sample records were extracted from the input files.")

    rows.sort(key=lambda r: (r["ReadSet"], r["ReactionWell"]))  # stable: keeps XML order within a cell
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=OUTPUT_COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nFiles read: {len(xml_files)}")
    print(f"Samples written: {len(rows)}")
    print(f"Output: {args.output}")
    if failures:
        print(f"Skipped files: {failures}")


if __name__ == "__main__":
    main()
