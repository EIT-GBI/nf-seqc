# nf-seQC

Batch-level QC for a sequencing run. Point it at the **complete output directory
of one sequencing run** and it detects which instrument produced it, normalises
the reads, and reports QC for the whole batch in a single MultiQC report.

One pipeline run = one sequencing run.

## Quickstart

```bash
# Clone with submodules — the nf-mod-* modules are git submodules
git clone --recurse-submodules https://github.com/EIT-GBI/nf-seQC.git

# Smoke-test the wiring on the bundled marker-file fixtures (no real data)
nextflow run . -profile test -stub
nextflow run . -profile test_aviti -stub
nextflow run . -profile test_nanopore -stub
nextflow run . -profile test_pacbio -stub

# Real run (GBI cluster: see below)
nextflow run . -profile cluster --input /path/to/run_dir --outdir results
```

### Running on the GBI cluster

`-profile cluster` targets the sandpit SLURM clusters and switches Apptainer on
by itself; don't add `singularity` to it. Run folders are under
`/mnt/instrument-data/<instrument>/` (read-only).

On **sandpit-tokyo**, Apptainer exists only on the compute nodes, so the
Nextflow head job must itself run as a SLURM job. A batch job does not load
your interactive shell, so put Java and Nextflow on `PATH` explicitly (here
from SDKMAN and `~/.local/bin`):

```bash
sbatch -p cpu -c 2 --mem 8G -t 24:00:00 -J nf-seqc -o nf-seqc-%j.out --wrap="
  export JAVA_HOME=\$HOME/.sdkman/candidates/java/current
  export PATH=\$JAVA_HOME/bin:\$HOME/.local/bin:\$PATH
  cd ~/nf-seQC && nextflow run . -profile cluster \
    --input /mnt/instrument-data/element_aviti_1/<run> \
    --outdir /mnt/lustre/users/\$USER/nf-seqc/<run> \
    -w /mnt/lustre/users/\$USER/nextflow-work"
```

Keep `-w` on Lustre, not the NFS home: bases2fastq writes the full FASTQ set
into `work/`, and `--outdir` must be on the same filesystem, because results
are published as hard links. Images are cached in the shared
`/mnt/lustre/containers/.staging/nf/cache/`.

A bare `-resume` picks up the **most recent** run launched from that directory,
whatever run it was. If you have launched another run from the same clone
since, it silently re-runs everything, and for an AVITI run whose data has
dropped out of the Alluxio cache that is about 2 h of bases2fastq. Resume a
specific run by its session ID from `nextflow log`
(`-resume <session-id>`), or launch each run from its own directory.

## Supported platforms

The platform is detected from the marker files each instrument writes; override
with `--platform illumina|aviti|nanopore|pacbio`.

| Platform | Detected by | Reads delivered | QC |
| --- | --- | --- | --- |
| Illumina NextSeq | `RunInfo.xml` | already demultiplexed on-instrument, left in place | FastQC + InterOp run metrics (see below) |
| Element AVITI | `RunParameters.json` | `bases2fastq` writes FASTQ to `<outdir>/reads` | FastQC + RunStats/Metrics |
| ONT PromethION P24 | `final_summary_*.txt`, `*.pod5` | `dorado` basecalls on **GPU** to uBAM, then FASTQ | NanoPlot |
| PacBio Revio | `*.consensusreadset.xml`, `*.hifi_reads*.bam` | HiFi uBAM already final, left in place | NanoPlot + samtools flagstat |

POD5 are expected to be already split per barcode: each directory containing
`*.pod5` is treated as one sample, named after that directory. There is no
`dorado demux` step.

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
| `--input` | — | run directory (required) |
| `--outdir` | `results` | output directory |
| `--platform` | auto-detect | force a platform |
| `--run_id` | run dir name | overrides the run label |
| `--fastq_dir` | — | Illumina: FASTQ live outside the run folder |
| `--pod5_dir` | — | Nanopore: POD5 live outside the run folder |
| `--bam_dir` | — | PacBio: HiFi BAM live outside the run folder |
| `--include_undetermined` | `false` | QC the Undetermined/Unassigned reads too |
| `--include_failed` | `false` | include `pod5_fail` / unassigned PacBio BAM |
| `--allow_incomplete` | `false` | Illumina: report even without `RTAComplete.txt`/`CopyComplete.txt` |
| `--basecall.model` | baked-in `sup@v5.0.0` | dorado model path inside the container |
| `--basecall.modified_bases` | `null` | e.g. `5mCG_5hmCG` to emit MM/ML tags |
| `--multiqc_config` | — | extra MultiQC config, applied on top of `assets/multiqc_config.yml` (which puts run-level metrics first) |

## Layout

- `main.nf` — platform detection, dispatch, and the `output {}` publish targets
- `subworkflows/local/` — one subworkflow per platform, plus `batch_report.nf`
- `modules/<tool>/` — git submodules from `nf-mod-*` repos
- `modules/local/` — tools with no `nf-mod-*` repo yet (bases2fastq, nanoplot, multiqc)
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
