#!/usr/bin/env python3
"""Per-barcode read statistics from MinKNOW sequencing_summary_*.txt files.

MinKNOW writes one row per read: its length, mean Q-score, pass/fail and the
barcode it classified the read to. That is everything NanoPlot computed by
reading every BAM, so a flow cell's QC needs one streaming pass over a text
file instead.

Every summary given is one acquisition of the same flow cell (a restart writes
a second one); their reads are pooled. Memory stays flat whatever the read
count: lengths and Q-scores go into histograms, from which N50 and the median
are exact (Q to 0.1).

Output, tab-separated, one row per barcode plus `all` (every barcode pooled,
for an unbarcoded flow cell):

    barcode  pass_reads  fail_reads  pass_bases  mean_length  n50  median_q

Pass reads only in the length and quality columns, matching pod5_pass, which
is what the pipeline basecalls by default.
"""

import argparse
import collections
import csv
import sys

csv.field_size_limit(sys.maxsize)

COLUMNS = ["barcode", "pass_reads", "fail_reads", "pass_bases", "mean_length", "n50", "median_q"]


class Stats:
    def __init__(self):
        self.pass_reads = 0
        self.fail_reads = 0
        self.bases = 0
        self.lengths = collections.Counter()
        self.quals = collections.Counter()  # Q x 10, as int

    def add(self, passed, length, qual):
        if not passed:
            self.fail_reads += 1
            return
        self.pass_reads += 1
        self.bases += length
        self.lengths[length] += 1
        if qual is not None:
            self.quals[round(qual * 10)] += 1

    def n50(self):
        half, seen = self.bases / 2, 0
        for length in sorted(self.lengths, reverse=True):
            seen += length * self.lengths[length]
            if seen >= half:
                return length
        return ""

    def median_q(self):
        total = sum(self.quals.values())
        if not total:
            return ""
        seen = 0
        for q in sorted(self.quals):
            seen += self.quals[q]
            if seen >= total / 2:
                return f"{q / 10:.1f}"
        return ""

    def row(self, barcode):
        mean = f"{self.bases / self.pass_reads:.0f}" if self.pass_reads else ""
        return [barcode, self.pass_reads, self.fail_reads, self.bases, mean,
                self.n50() if self.pass_reads else "", self.median_q()]


def as_bool(value):
    return value is None or value.strip().upper() in ("TRUE", "1", "")


def as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("summaries", nargs="+", help="sequencing_summary_*.txt of one flow cell")
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    by_barcode = collections.defaultdict(Stats)
    everything = Stats()
    for path in args.summaries:
        with open(path, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for read in reader:
                length = int(read.get("sequence_length_template") or 0)
                qual = as_float(read.get("mean_qscore_template"))
                passed = as_bool(read.get("passes_filtering"))
                # MinKNOW names barcodes barcode01; older versions and some
                # kits write the arrangement as <kit>_barcode01.
                barcode = (read.get("barcode_arrangement") or "unclassified").split("_")[-1]
                by_barcode[barcode].add(passed, length, qual)
                everything.add(passed, length, qual)

    with open(args.output, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(COLUMNS)
        for barcode in sorted(by_barcode):
            writer.writerow(by_barcode[barcode].row(barcode))
        writer.writerow(everything.row("all"))


if __name__ == "__main__":
    main()
