#!/usr/bin/env python3
"""Per-barcode read statistics and MultiQC plots for one ONT flow cell.

Reads MinKNOW's own records of the flow cell, never the reads themselves:

* sequencing_summary_*.txt: one row per read (length, mean Q-score, pass/fail,
  time, channel, end reason, barcode). Every summary given is one acquisition
  of the same flow cell (a restart writes a second one); their reads are
  pooled. Streamed once; memory stays flat whatever the read count, because
  everything goes into histograms (N50 and the median Q are still exact).
* report_*.json: what MinKNOW's own HTML report is drawn from (pore activity,
  pore occupancy, temperature, bias voltage, run messages). Optional.

Writes:

* -o: tab-separated, one row per barcode plus `all` (every barcode pooled, for
  an unbarcoded flow cell). Pass reads only in the length and quality columns,
  matching pod5_pass, which the pipeline basecalls by default.

      barcode  pass_reads  fail_reads  pass_bases  mean_length  n50  median_q

* --mqc-dir: MultiQC custom-content JSON, one file per section. Series are
  named after the pipeline's samples (sample-sheet alias, else the barcode;
  see ontSamples() in nanopore_qc.nf), flow-cell-level series after the flow
  cell. Each section has a fixed id, so the files of several flow cells merge
  into one plot; that is why no plot here has several datasets (MultiQC keeps
  only the last file's datasets).
"""

import argparse
import collections
import csv
import datetime
import json
import math
import os
import re
import sys

csv.field_size_limit(sys.maxsize)

COLUMNS = ["barcode", "pass_reads", "fail_reads", "pass_bases", "mean_length", "n50", "median_q"]

# Read lengths in log-spaced bins, 20 per decade from 10 bp to 10 Mb: fine
# enough for a smooth curve, a few hundred points per sample at most.
LEN_BINS_PER_DECADE = 20
YIELD_BIN_S = 900          # yield over time: 15-minute bins
TIME_BIN_S = 3600          # Q-score and active channels over time: 1 h bins
BIAS_BIN_S = 600           # bias voltage: 10-minute steps

# MinKNOW's channel-state groups and colours, as in its own pore activity plot.
STATE_GROUPS = {
    "sequencing": ("Sequencing", "#00a651"),
    "pore": ("Pore available", "#90c6e7"),
    "recovering": ("Recovering", "#f7b538"),
    "inactive": ("Inactive", "#8e8e8e"),
    "unclassified": ("Unclassified", "#d1d1d1"),
}
DEFAULT_STATE_GROUP = {
    "strand": "sequencing", "adapter": "sequencing", "pore": "pore",
    "unavailable": "recovering", "unblocking": "recovering",
}


def safe_id(name):
    """The pipeline's sample ids, as safeId() in nanopore_qc.nf."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name)


def len_bin(length):
    """Lower edge of a read length's log bin."""
    if length < 10:
        return 0
    return int(math.floor(math.log10(length) * LEN_BINS_PER_DECADE))


def len_bin_centre(b):
    return round(10 ** ((b + 0.5) / LEN_BINS_PER_DECADE))


class Stats:
    """Everything kept per sample: counts and histograms, never reads."""

    def __init__(self):
        self.pass_reads = 0
        self.fail_reads = 0
        self.bases = 0
        self.fail_bases = 0
        self.lengths = collections.Counter()      # exact pass lengths, for N50
        self.quals = collections.Counter()        # pass Q x 10, for the median
        self.q_all = collections.Counter()        # every read's Q x 2 (0.5 bins)
        self.len_reads = collections.Counter()    # log length bin -> pass reads
        self.len_bases = collections.Counter()    # log length bin -> pass bases
        self.len_qsum = collections.Counter()     # log length bin -> sum of Q
        self.len_qn = collections.Counter()
        self.speed = collections.Counter()        # bases/s in 10 b/s bins
        self.end_reasons = collections.Counter()
        self.bc_score = collections.Counter()     # barcode score, integer bins
        self.yield_bins = collections.Counter()   # time bin -> pass bases

    def add(self, row):
        length = row["length"]
        qual = row["qual"]
        if qual is not None:
            self.q_all[int(qual * 2)] += 1
        if row["end_reason"]:
            self.end_reasons[row["end_reason"]] += 1
        if row["bc_score"] is not None:
            self.bc_score[int(row["bc_score"])] += 1
        if not row["passed"]:
            self.fail_reads += 1
            self.fail_bases += length
            return
        self.pass_reads += 1
        self.bases += length
        self.lengths[length] += 1
        b = len_bin(length)
        self.len_reads[b] += 1
        self.len_bases[b] += length
        if qual is not None:
            self.quals[round(qual * 10)] += 1
            self.len_qsum[b] += qual
            self.len_qn[b] += 1
        if row["duration"] and row["duration"] > 0 and length:
            self.speed[int(length / row["duration"] // 10) * 10] += 1
        if row["end"] is not None:
            self.yield_bins[int(row["end"] // YIELD_BIN_S)] += length

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


def barcode_of(read):
    # MinKNOW names barcodes barcode01; older versions and some kits write the
    # arrangement as <kit>_barcode01.
    return (read.get("barcode_arrangement") or "unclassified").split("_")[-1]


# ---------------------------------------------------------------------------
# MultiQC custom content
# ---------------------------------------------------------------------------

# Two report groups (MultiQC modules), each section in one. Their order, and
# the sections' order within them, is set in assets/multiqc_config.yml.
PARENTS = {
    "nanopore_flowcell": (
        "Nanopore flow cell",
        "How each flow cell behaved, from MinKNOW's report_*.json (what its own HTML report shows) and "
        "sequencing_summary: pore activity and occupancy, channels, temperature, bias voltage, warnings.",
    ),
    "nanopore_reads": (
        "Nanopore reads",
        "Read-level QC per sample from MinKNOW's sequencing_summary_*.txt: MinKNOW's live basecalls and "
        "barcoding, one row per read. Samples are named as the pipeline names them; % Reads in General "
        "Statistics comes from the pipeline's own demultiplexing.",
    ),
}
SECTION_PARENT = {
    "ont_pore_activity_total": "nanopore_flowcell", "ont_pore_activity": "nanopore_flowcell",
    "ont_pore_occupancy": "nanopore_flowcell", "ont_active_channels": "nanopore_flowcell",
    "ont_q_over_time": "nanopore_flowcell", "ont_temperature": "nanopore_flowcell",
    "ont_bias_voltage": "nanopore_flowcell", "ont_minknow_messages": "nanopore_flowcell",
}


def section(sid, name, description, plot_type, data, pconfig=None, **extra):
    parent = SECTION_PARENT.get(sid, "nanopore_reads")
    d = {
        "id": sid,
        "parent_id": parent,
        "parent_name": PARENTS[parent][0],
        "parent_description": PARENTS[parent][1],
        "section_name": name,
        "description": description,
        "plot_type": plot_type,
        "pconfig": {"id": f"{sid}_plot", "title": f"Nanopore: {name}", **(pconfig or {})},
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

def pct(counter, total):
    return {k: 100.0 * v / total for k, v in sorted(counter.items())} if total else {}


def summary_sections(samples, per_flowcell, flowcell, start_s):
    """Sections drawn from sequencing_summary. samples: name -> Stats."""
    out = []
    named = {n: s for n, s in samples.items() if s.pass_reads + s.fail_reads}

    out.append(section(
        "ont_pass_fail", "Reads passing MinKNOW's filter",
        "Reads per sample that passed or failed MinKNOW's live Q-score filter, from sequencing_summary. "
        "Not every run applies the filter (then every read passes).",
        "bargraph",
        {n: {"pass": s.pass_reads, "fail": s.fail_reads} for n, s in named.items()},
        {"ylab": "Reads", "cpswitch_counts_label": "Reads"},
        categories={"pass": {"name": "Pass", "color": "#00a651"}, "fail": {"name": "Fail", "color": "#d9534f"}},
    ))

    out.append(section(
        "ont_yield_over_time", "Yield over time",
        "Cumulative pass bases per sample against time since the flow cell started sequencing (read end time), "
        "from sequencing_summary. A curve that flattens early points to pores running out.",
        "linegraph",
        {n: cumulative_gb(s.yield_bins, start_s) for n, s in named.items() if s.yield_bins},
        {"xlab": "Hours", "ylab": "Pass yield (Gb)", "ymin": 0, "xmin": 0},
    ))

    out.append(section(
        "ont_read_length", "Read length distribution",
        "Share of each sample's pass reads by read length (log scale, 20 bins per decade), from sequencing_summary.",
        "linegraph",
        {n: {len_bin_centre(b): v for b, v in sorted(pct(s.len_reads, s.pass_reads).items())}
         for n, s in named.items() if s.pass_reads},
        {"xlab": "Read length (bp)", "ylab": "% of pass reads", "xlog": True, "ymin": 0},
    ))

    out.append(section(
        "ont_read_length_yield", "Yield by read length",
        "Share of each sample's pass bases by read length: where the yield is, rather than where the reads are.",
        "linegraph",
        {n: {len_bin_centre(b): v for b, v in sorted(pct(s.len_bases, s.bases).items())}
         for n, s in named.items() if s.bases},
        {"xlab": "Read length (bp)", "ylab": "% of pass bases", "xlog": True, "ymin": 0},
    ))

    out.append(section(
        "ont_qscore", "Q-score distribution",
        "Share of each sample's reads by mean read Q-score (0.5 bins), pass and fail together, so MinKNOW's "
        "filter threshold shows as the step where failed reads stop.",
        "linegraph",
        {n: {k / 2: v for k, v in pct(s.q_all, sum(s.q_all.values())).items()}
         for n, s in named.items() if s.q_all},
        {"xlab": "Mean read Q-score", "ylab": "% of reads", "ymin": 0},
    ))

    out.append(section(
        "ont_q_by_length", "Q-score by read length",
        "Mean Q-score of each sample's pass reads in each read-length bin. Short reads with low Q are often "
        "adapter or chimeric fragments.",
        "linegraph",
        {n: {len_bin_centre(b): s.len_qsum[b] / s.len_qn[b] for b in sorted(s.len_qn) if s.len_qn[b] >= 5}
         for n, s in named.items() if s.len_qn},
        {"xlab": "Read length (bp)", "ylab": "Mean Q-score", "xlog": True},
    ))

    out.append(section(
        "ont_speed", "Translocation speed",
        "Share of each sample's pass reads by translocation speed (bases per second of template signal, "
        "10 b/s bins). R10.4.1 runs at a nominal 400 b/s; a broad or shifted peak can point to a pore or "
        "temperature problem.",
        "linegraph",
        {n: pct(s.speed, sum(s.speed.values())) for n, s in named.items() if s.speed},
        {"xlab": "Bases per second", "ylab": "% of pass reads", "ymin": 0, "xmin": 0, "xmax": 1000},
    ))

    reasons = sorted({r for s in named.values() for r in s.end_reasons})
    out.append(section(
        "ont_end_reasons", "Read end reasons",
        "Why each read ended, from sequencing_summary. signal_positive is a strand leaving the pore normally; "
        "a large unblock_mux_change share means many reads were cut short by pore unblocking, and "
        "data_service_unblock_mux_change is adaptive sampling rejecting reads.",
        "bargraph",
        {n: dict(s.end_reasons) for n, s in named.items() if s.end_reasons},
        {"ylab": "Reads", "cpswitch_c_active": False},
        categories=reasons,
    ))

    scored = {n: s for n, s in named.items() if sum(s.bc_score.values())}
    if scored:
        out.append(section(
            "ont_barcode_score", "Barcode score",
            "Share of each sample's reads by MinKNOW's barcode classification score (0-100). Barcoded reads "
            "should pile up near 100; unclassified reads sit low.",
            "linegraph",
            {n: pct(s.bc_score, sum(s.bc_score.values())) for n, s in scored.items()},
            {"xlab": "Barcode score", "ylab": "% of reads", "ymin": 0, "xmin": 0, "xmax": 100},
        ))

    fc = per_flowcell
    if fc["q_time"]:
        out.append(section(
            "ont_q_over_time", "Q-score over time",
            "Median mean-read Q-score of each flow cell's reads, per hour of sequencing. A steady decline is "
            "normal pore wear; a sudden drop is worth a look.",
            "linegraph",
            {flowcell: {bin_hours(t, start_s): median_from(c) for t, c in sorted(fc["q_time"].items())}},
            {"xlab": "Hours", "ylab": "Median Q-score", "xmin": 0},
        ))
    if fc["channels"]:
        out.append(section(
            "ont_active_channels", "Active channels over time",
            "Channels that produced at least one read in each hour of sequencing, per flow cell.",
            "linegraph",
            {flowcell: {bin_hours(t, start_s): len(c) for t, c in sorted(fc["channels"].items())}},
            {"xlab": "Hours", "ylab": "Channels", "ymin": 0, "xmin": 0},
        ))
    return out


def cumulative_gb(bins, start_s):
    total, out = 0, {}
    for t in sorted(bins):
        total += bins[t]
        out[round(((t + 1) * YIELD_BIN_S - start_s) / 3600, 3)] = total / 1e9
    return out


def bin_hours(t, start_s):
    """Middle of an hourly bin, in hours since sequencing started."""
    return round(max(0.0, ((t + 0.5) * TIME_BIN_S - start_s) / 3600), 2)


def median_from(counter):
    total, seen = sum(counter.values()), 0
    for q in sorted(counter):
        seen += counter[q]
        if seen >= total / 2:
            return q / 10
    return None


def parse_time(stamp):
    if not stamp:
        return None
    # MinKNOW writes nanoseconds; Python takes microseconds.
    stamp = re.sub(r"(\.\d{6})\d+", r"\1", stamp.replace("Z", "+00:00"))
    try:
        return datetime.datetime.fromisoformat(stamp)
    except ValueError:
        return None


def report_sections(reports, flowcell):
    """Sections drawn from MinKNOW's report_*.json: one series per sequencing
    acquisition (a restarted flow cell has two, suffixed with the protocol run)."""
    acquisitions = []
    for path in reports:
        with open(path) as handle:
            d = json.load(handle)
        run = d.get("protocol_run_info") or {}
        seq = [a for a in d.get("acquisitions") or []
               if ((a.get("acquisition_run_info") or {}).get("config_summary") or {}).get("purpose") == "SEQUENCING"]
        if seq:
            acquisitions.append((run, seq[-1], d.get("user_messages") or []))
    if not acquisitions:
        return []
    multi = len(acquisitions) > 1

    activity, totals, occupancy, temperature, bias, messages = {}, {}, {}, {}, {}, {}
    for run, acq, user_messages in acquisitions:
        name = f"{flowcell}_{(run.get('run_id') or '')[:8]}" if multi else flowcell
        info = acq.get("acquisition_run_info") or {}
        groups = {}
        for g in ((info.get("config_summary") or {}).get("channel_state_info") or {}).get("groups") or []:
            for s in g.get("states") or []:
                groups[s.get("name")] = g.get("name")
        duty = (acq.get("duty_time") or [{}])[0]
        buckets = duty.get("bucket_ranges") or []
        hours = [((b.get("start") or 0) + (b.get("end") or 0)) / 2 / 3600 for b in buckets]

        by_group = collections.defaultdict(lambda: [0.0] * len(buckets))
        for state, v in (duty.get("channel_states") or {}).items():
            group = groups.get(state) or DEFAULT_STATE_GROUP.get(state, "inactive")
            for i, t in enumerate(v.get("state_times") or []):
                if i < len(buckets):
                    by_group[group][i] += float(t or 0)
        bucket_total = [sum(by_group[g][i] for g in by_group) for i in range(len(buckets))]
        run_total = sum(bucket_total)
        totals[name] = {g: 100.0 * sum(by_group[g]) / run_total for g in by_group} if run_total else {}
        for g, values in by_group.items():
            label = STATE_GROUPS.get(g, (g, None))[0]
            activity[f"{name}: {label}"] = {
                round(hours[i], 2): 100.0 * values[i] / bucket_total[i] for i in range(len(buckets)) if bucket_total[i]
            }

        occ = duty.get("pore_occupancy") or []
        occupancy[name] = {round(hours[i], 2): 100.0 * float(v) for i, v in enumerate(occ) if v is not None and i < len(hours)}

        temps = ((acq.get("temperature") or [{}])[0]).get("temperatures") or []
        temperature[name] = {}
        for i, t in enumerate(temps):
            device = next((v for k, v in t.items() if k != "target_temperature" and isinstance(v, dict)), {})
            value = device.get("heatsink_temperature", device.get("flowcell_temperature"))
            if value is not None and i < len(hours):
                temperature[name][round(hours[i], 2)] = value

        # The voltage in force at the end of each 10-minute bin. MinKNOW logs a
        # reading only when it changes (tens of thousands over a three-day run,
        # flicking during every pore scan), so this is a step line.
        last_in_bin = {}
        for v in ((acq.get("bias_voltage") or [{}])[0]).get("bias_voltages") or []:
            if "bias_voltage" in v and "time_seconds" in v:
                last_in_bin[int(float(v["time_seconds"]) // BIAS_BIN_S)] = float(v["bias_voltage"])
        if last_in_bin:
            voltage, steps = None, {}
            for b in range(min(last_in_bin), max(last_in_bin) + 1):
                voltage = last_in_bin.get(b, voltage)
                steps[b] = voltage
            bias[name] = {round((b + 1) * BIAS_BIN_S / 3600, 3): v for b, v in steps.items()}

        start = parse_time(info.get("data_read_start_time") or info.get("start_time"))
        for i, m in enumerate(user_messages):
            severity = (m.get("severity") or "").replace("MESSAGE_SEVERITY_", "").lower()
            if severity in ("warning", "error"):
                when = parse_time(m.get("time"))
                hours_in = (when - start).total_seconds() / 3600 if when and start else None
                messages[f"{name} #{i + 1}"] = {
                    "severity": severity,
                    "hours": hours_in if hours_in is not None else "",
                    "identifier": m.get("identifier") or "",
                    "message": (m.get("user_message") or "").strip(),
                }

    states = sorted({g for t in totals.values() for g in t}, key=lambda g: list(STATE_GROUPS).index(g) if g in STATE_GROUPS else 99)
    out = [
        section(
            "ont_pore_activity_total", "Pore activity",
            "Share of channel time in each MinKNOW channel state over the whole run, per flow cell, from "
            "report_*.json (what MinKNOW's Pore activity plot shows). Sequencing is strand + adapter; "
            "Pore available is an open pore waiting for a strand.",
            "bargraph",
            {n: {g: v for g, v in t.items()} for n, t in totals.items()},
            {"ylab": "% of channel time", "cpswitch": False, "stacking": "normal", "ymax": 100},
            categories={g: {"name": STATE_GROUPS.get(g, (g,))[0], "color": STATE_GROUPS.get(g, (g, None))[1]} for g in states},
        ),
        section(
            "ont_pore_activity", "Pore activity over time",
            "Share of channel time in each MinKNOW channel-state group, per time bucket of the run (MinKNOW's own "
            "buckets, usually 2 h). Sequencing falling while Inactive rises is pores dying; Pore available rising "
            "with Sequencing falling is the library running out.",
            "linegraph",
            activity,
            {"xlab": "Hours", "ylab": "% of channel time", "ymin": 0, "ymax": 100, "xmin": 0,
             # MinKNOW's colours, so the plot reads like its own pore activity plot
             "colors": {n: STATE_GROUPS[g][1] for n in activity for g in STATE_GROUPS
                        if n.endswith(": " + STATE_GROUPS[g][0])}},
        ),
        section(
            "ont_pore_occupancy", "Pore occupancy",
            "Share of available pores that were sequencing (strand + adapter), per time bucket, from report_*.json.",
            "linegraph",
            occupancy,
            {"xlab": "Hours", "ylab": "Pore occupancy (%)", "ymin": 0, "ymax": 100, "xmin": 0},
        ),
    ]
    if any(temperature.values()):
        out.append(section(
            "ont_temperature", "Temperature",
            "Device heatsink temperature per time bucket, from report_*.json. MinKNOW holds it near the target "
            "(35-37 °C); drift means the device struggled to regulate it.",
            "linegraph",
            {n: t for n, t in temperature.items() if t},
            {"xlab": "Hours", "ylab": "°C", "xmin": 0},
        ))
    if bias:
        out.append(section(
            "ont_bias_voltage", "Bias voltage",
            "Applied bias voltage over the run, from report_*.json. MinKNOW steps it up as pores age, to keep "
            "translocation speed constant.",
            "linegraph",
            bias,
            {"xlab": "Hours", "ylab": "mV", "xmin": 0},
        ))
    if messages:
        out.append(section(
            "ont_minknow_messages", "MinKNOW warnings",
            "Warnings and errors MinKNOW raised during the run, from report_*.json. Hours are from the start of "
            "sequencing.",
            "table",
            messages,
            {"col1_header": "Flow cell", "sort_rows": False},
            headers={
                "severity": {"title": "Severity", "scale": False},
                "hours": {"title": "Hours", "format": "{:,.1f}", "scale": False},
                "identifier": {"title": "Type", "scale": False},
                "message": {"title": "Message", "scale": False},
            },
        ))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("summaries", nargs="+", help="sequencing_summary_*.txt of one flow cell")
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--report", nargs="*", default=[], help="report_*.json of the same flow cell")
    parser.add_argument("--flowcell", default="flowcell", help="flow cell id, names flow-cell-level series")
    parser.add_argument("--plain-name", help="sample name when the flow cell is not barcoded")
    parser.add_argument("--suffix-flowcell", action="store_true",
                        help="append _<flowcell> to sample names (the run has several flow cells)")
    parser.add_argument("--aliases", help="JSON {barcode: alias} from the sample sheet")
    parser.add_argument("--mqc-dir", help="write MultiQC custom-content JSON here")
    parser.add_argument("--min-reads", type=int, default=0,
                        help="leave barcodes with fewer reads out of the plots, as the pipeline leaves them out of its samples")
    args = parser.parse_args()

    aliases = json.loads(args.aliases) if args.aliases else {}
    by_barcode = collections.defaultdict(Stats)
    everything = Stats()
    q_time = collections.defaultdict(collections.Counter)
    channels = collections.defaultdict(set)
    start_s = None
    for path in args.summaries:
        with open(path, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for read in reader:
                start = as_float(read.get("start_time"))
                duration = as_float(read.get("duration"))
                row = {
                    "length": int(read.get("sequence_length_template") or 0),
                    "qual": as_float(read.get("mean_qscore_template")),
                    "passed": as_bool(read.get("passes_filtering")),
                    "end_reason": read.get("end_reason") or "",
                    "bc_score": as_float(read.get("barcode_score")) if read.get("barcode_arrangement") else None,
                    "duration": as_float(read.get("template_duration")) or duration,
                    "end": start + (duration or 0) if start is not None else None,
                }
                by_barcode[barcode_of(read)].add(row)
                everything.add(row)
                if start is not None:
                    start_s = start if start_s is None else min(start_s, start)
                    t = int(start // TIME_BIN_S)
                    if row["qual"] is not None:
                        q_time[t][round(row["qual"] * 10)] += 1
                    if read.get("channel"):
                        channels[t].add(read["channel"])

    with open(args.output, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(COLUMNS)
        for barcode in sorted(by_barcode):
            writer.writerow(by_barcode[barcode].row(barcode))
        writer.writerow(everything.row("all"))

    if not args.mqc_dir:
        return
    os.makedirs(args.mqc_dir, exist_ok=True)

    # Name series as the pipeline names its samples.
    barcoded = any(b != "unclassified" for b in by_barcode)
    if barcoded:
        samples = {}
        for barcode, s in by_barcode.items():
            if barcode != "unclassified" and s.pass_reads + s.fail_reads < args.min_reads:
                continue
            if barcode == "unclassified":
                name = f"unclassified_{args.flowcell}"
            else:
                name = safe_id(aliases.get(barcode) or barcode)
                if args.suffix_flowcell:
                    name = f"{name}_{args.flowcell}"
            samples[name] = s
    else:
        samples = {args.plain_name or args.flowcell: everything}

    sections = summary_sections(samples, {"q_time": q_time, "channels": channels}, args.flowcell, start_s or 0)
    sections += report_sections(args.report, args.flowcell)
    for s in sections:
        if not s["data"]:
            continue
        with open(os.path.join(args.mqc_dir, f"{safe_id(args.flowcell)}.{s['id']}_mqc.json"), "w") as handle:
            json.dump(as_pairs(s), handle)


if __name__ == "__main__":
    main()
