# nf-seqc

Batch-level QC for a sequencing run. Point it at the **complete output directory
of one sequencing run** and it detects which instrument produced it, normalises
the reads, and reports QC for the whole batch in a single MultiQC report.

One pipeline run = one sequencing run.

## Quickstart

```bash
# Clone with submodules — the nf-mod-* modules are git submodules
git clone --recurse-submodules https://github.com/EIT-GBI/nf-seqc.git

# Smoke-test the wiring on the bundled marker-file fixtures (no real data)
nextflow run . -profile test -stub
nextflow run . -profile test_aviti -stub
nextflow run . -profile test_nanopore -stub
nextflow run . -profile test_pacbio -stub
nextflow run . -profile test_publish -stub   # --instrument/--run route

# Real run (GBI cluster: see below)
nextflow run . -profile cluster --instrument illumina_nextseq_VH02427 --run <run folder>
```

### Where runs are read from and published to

Instruments write to
`/mnt/instrument-prep-zone/sequencing/<brand>_<model>_<machine_ID>/<YYYY>/<run>`.
Name a run by its instrument and folder:

```bash
nextflow run . -profile cluster --instrument pacbio_revio_84334 --run r84334_20260817_141703
```

The pipeline finds the year directory itself and publishes the MultiQC report,
QC files and reads to the same path under `/mnt/instrument-data/sequencing/`,
where scientists read them:
`/mnt/instrument-data/sequencing/pacbio_revio_84334/2026/r84334_20260817_141703/`.
The roots are `--prep_root` and `--publish_root`; `--outdir` overrides the
destination. `--input <run_dir>` still takes any run folder, publishing to
`--outdir` (default `results`).

Inside the run folder, what scientists open sits at the top and everything
else is in `qc/`:

```
<run>/
├── <run>_multiqc_report.html     the QC report for the whole run
├── <run>_QC.html                 AVITI only: Element's own run report
├── <sample>…fastq.gz             Illumina, AVITI
├── <sample>.bam                  Nanopore (uBAM, keeps modified-base tags)
├── <movie>.hifi_reads.<bc>.bam   Revio, with its .pbi index
└── qc/                           FastQC, MultiQC inputs, run metrics CSV,
                                  software versions, publish record
```

- The brand in `--instrument` (`illumina`, `element`, `oxford_nanopore`,
  `pacbio`) must match the platform detected from the run's files, or the run
  stops before any task.
- A run is published once. The first launch records its session ID in
  `qc/nf-seqc_session.txt`. If the destination already has files,
  the run stops, unless it is `-resume <that session ID>` (the error names it)
  or `--overwrite` is given. Two launches of the same run at once: only one
  gets past the start. `--overwrite` writes over files of the same name; it
  does not remove files from an earlier run.
- Results are **copied**, since the destination is another filesystem from
  `work/`. Illumina FASTQ and Revio HiFi BAM (with `.pbi`) are copied into a
  task first (`COPY_READS`), because only task outputs can be published, so
  they briefly exist three times: run folder, `work/`, destination. Clean
  `work/` once a run is published.

### Running on the GBI cluster

`-profile cluster` targets the sandpit SLURM clusters and switches Apptainer on
by itself; don't add `singularity` to it. Run folders are under
`/mnt/instrument-prep-zone/sequencing/<instrument>/<YYYY>/`.

On **sandpit-tokyo**, Apptainer exists only on the compute nodes, so the
Nextflow head job must itself run as a SLURM job. Nextflow (with its Java)
comes from the shared module, `module load nextflow`; a batch script needs a
login shell (`#!/bin/bash -l`) for `module` to exist:

```bash
#!/bin/bash -l
#SBATCH -p cpu -c 2 --mem 8G -t 24:00:00 -J nf-seqc -o nf-seqc-%j.out
set -euo pipefail
module load nextflow
cd ~/nf-seqc
nextflow run . -profile cluster \
    --instrument element_aviti_AV252906 --run <run> \
    -w /mnt/lustre/users/$USER/nextflow-work
```

Keep `-w` on Lustre, not the NFS home: bases2fastq writes the full FASTQ set
into `work/`. Images are cached in the shared
`/mnt/lustre/containers/.staging/nf/cache/`.

The head job pulls any image missing from that cache, and converting the
11 GB dorado image to SIF needs more memory than an 8 GB head job has
(`mksquashfs` is OOM-killed). Before the first Nanopore run on a new dorado
image, pull it once in a job of its own, e.g. `-c 2 --mem 24G` (it peaked at
~20 GB), with
`APPTAINER_TMPDIR` on Lustre, then move the file into the cache under the name
Nextflow expects (`ghcr.io-eit-gbi-nf-mod-dorado-<tag>.img`).

A bare `-resume` picks up the **most recent** run launched from that directory,
whatever run it was. If you have launched another run from the same clone
since, it silently re-runs everything, and for an AVITI run whose data has
dropped out of the Alluxio cache that is about 2 h of bases2fastq. Resume a
specific run by its session ID from `nextflow log`
(`-resume <session-id>`), or launch each run from its own directory.

### Real-data test runs

`-stub` exercises only the wiring. The `test_real_*` profiles run the tools on
real, miniaturised runs kept on sandpit-tokyo, so they work only there, next to
`cluster`:

```bash
nextflow run . -profile cluster,test_real_illumina   # NextSeq 2000 subset
nextflow run . -profile cluster,test_real_aviti      # AVITI, bases2fastq on 2 tiles
nextflow run . -profile cluster,test_real_nanopore   # P24 subset, needs a GPU node
nextflow run . -profile cluster,test_real_pacbio     # Revio, one SMRT Cell
```

Submit them with `sbatch` as above. `conf/test_real.config` lists what each
fixture holds.

## Supported platforms

The platform is detected from the marker files each instrument writes; override
with `--platform illumina|aviti|nanopore|pacbio`.

| Platform | Detected by | Reads delivered | QC |
| --- | --- | --- | --- |
| Illumina NextSeq | `RunInfo.xml` | already demultiplexed on-instrument, left in place | FastQC + InterOp run metrics (see below) + BCL Convert demultiplexing reports |
| Element AVITI | `RunParameters.json` | `bases2fastq` writes FASTQ in the pipeline | FastQC + RunStats/Metrics |
| ONT PromethION P24 | `final_summary_*.txt`, `*.pod5` | `dorado` basecalls on **GPU** to uBAM and demultiplexes; uBAM delivered, no FASTQ | per-sample reads, yield, N50 and median Q from MinKNOW's `sequencing_summary_*.txt` (its live basecalls, no BAM is read) + per-flow-cell run metrics from MinKNOW's `report_*.json` (pore scans, occupancy, live yield and N50) + barcode counts from the pipeline's demux |
| PacBio Revio | `<well>/pb_formats/*.consensusreadset.xml`, `<well>/hifi_reads/*.hifi_reads*.bam` | HiFi uBAM already final, left in place | read from the instrument's own reports, no BAM is read: per-sample HiFi reads, yield, length, QV and barcode quality (`barcodes.report.json`) + per-SMRT-Cell run metrics (P0/P1/P2, HiFi yield, N50, QV) + lima barcode counts |

MinKNOW does not split POD5 cleanly by barcode: `pod5_pass/mixed/` holds part
of each barcode's reads. So every POD5 folder is basecalled with the run's
barcoding kit (`--kit-name`, read from `final_summary_*.txt`), and
`dorado demux` then splits each flow cell's reads by barcode. Samples are named
after the sample sheet's `alias`, or the barcode. Unclassified reads are counted
in the Demultiplexing table but not delivered unless `--include_failed`.

### Illumina run metrics

The Illumina branch mirrors
[Automate-Seq-Run-Metrics-Collection](https://github.com/EIT-GBI/Automate-Seq-Run-Metrics-Collection):
same metric set, same 24 column names, same %Occupied vs %PF tile plot, so rows
from either tool concatenate into one log. It reads InterOp via Illumina's Python
bindings plus `RunParameters.xml` and `PrimaryAnalysisMetrics.csv` (the
`INTEROP_RUN_METRICS` process from
[nf-mod-interop](https://github.com/EIT-GBI/nf-mod-interop)), and emits:

- `metrics/<run>.run_metrics.csv` — header + one row, the log row for this run
- `qc/<run>_run_metrics_mqc.tsv` — the same table as a MultiQC section
- `qc/<run>_occupied_vs_pf_mqc.png` — the loading plot, embedded in the report

Differences, both forced by running under Nextflow rather than cron:

- That tool appends to one central CSV across runs; here one pipeline run is one
  sequencing run, so a single-row CSV is emitted per run. Appending to a shared
  file from a task would race and break `-resume`.
- Its `monitoring_complete.txt` marker exists to stop cron reprocessing a run.
  `-resume` covers that, so **nothing in the run directory is modified**.

The `RTAComplete.txt` + `CopyComplete.txt` completeness gate is kept — metrics
from a partial run are misleading. Override with `--allow_incomplete`.

## Parameters

| Param | Default | Purpose |
| --- | --- | --- |
| `--instrument` | — | instrument directory, `<brand>_<model>_<machine_ID>`; with `--run` |
| `--run` | — | run folder name, looked up under `<prep_root>/<instrument>/<YYYY>/` |
| `--prep_root` | `/mnt/instrument-prep-zone/sequencing` | where instruments write runs |
| `--publish_root` | `/mnt/instrument-data/sequencing` | where `--instrument`/`--run` results are published |
| `--input` | — | run directory, instead of `--instrument`/`--run` |
| `--outdir` | `<publish_root>/<instrument>/<YYYY>/<run>`, or `results` with `--input` | output directory |
| `--overwrite` | `false` | publish into a destination that already has results |
| `--platform` | auto-detect | force a platform |
| `--run_id` | run dir name | overrides the run label |
| `--fastq_dir` | — | Illumina: FASTQ live outside the run folder |
| `--pod5_dir` | — | Nanopore: POD5 live outside the run folder |
| `--barcode_kit` | — | Nanopore: barcoding kit (e.g. `SQK-RBK114-24`) when `final_summary` does not name one |
| `--min_barcode_reads` | `100` | Nanopore: fewest reads for a barcode to be reported as a sample; reads of smaller barcodes count as unclassified |
| `--bam_dir` | — | PacBio: HiFi BAM live outside the run folder |
| `--include_undetermined` | `false` | QC the Undetermined/Unassigned reads too |
| `--include_failed` | `false` | Nanopore: include `pod5_fail` and the unclassified reads. PacBio failed reads are always counted per barcode (hidden General Statistics column) |
| `--allow_incomplete` | `false` | Report even on an unfinished run: Illumina without `RTAComplete.txt`/`CopyComplete.txt`, PacBio with a SMRT Cell lacking `metadata/*.transferdone` |
| `--max_undetermined` | `15` | Illumina, AVITI and barcoded PacBio: flag the run as high in the Demultiplexing table, plus a log warning, above this % of Undetermined (Illumina) or Unassigned (AVITI, PacBio Not Barcoded) reads |
| `--basecall.model` | baked-in `sup@v5.0.0` | dorado model path inside the container |
| `--basecall.modified_bases` | `null` | e.g. `5mCG_5hmCG` to emit MM/ML tags |
| `--multiqc_config` | — | extra MultiQC config, applied on top of `assets/multiqc_config.yml` (report order, and hides bases2fastq plots FastQC already shows) |

## Layout

- `main.nf` — platform detection, dispatch, and the `output {}` publish targets
- `subworkflows/local/` — one subworkflow per platform, plus `batch_report.nf`
- `modules/<tool>/` — git submodules from `nf-mod-*` repos
- `modules/local/` — tools with no `nf-mod-*` repo yet (bases2fastq, multiqc, ont_summary_stats)
- `conf/` — `base.config` (resources), `test.config`. Container images are set in
  `nextflow.config`: each released module's `conf/module.config` is included, and
  local-module images are pinned in its `process {}` block. Nothing uses `:latest`
  except `-profile dev`
- `assets/test/` — marker-file fixtures for `-stub` smoke runs

## Adding a module

```bash
git submodule add https://github.com/EIT-GBI/nf-mod-<tool>.git modules/<tool>
git -C modules/<tool> checkout v<version>
```

Record the version in `modules.versions` (the source of truth — the sync
workflow moves the submodule pointer to match), add
`includeConfig 'modules/<tool>/conf/module.config'` to the module-containers
section of `nextflow.config` plus a `:latest` line for it in the `dev` profile,
then include it in the relevant subworkflow:

```groovy
include { TOOL_SUBCOMMAND } from '../../modules/<tool>/<subcommand>/main.nf'
```
