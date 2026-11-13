# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this pipeline is for

nf-seqc takes the **complete output directory of one sequencing run** and reports QC for the **batch/run**, across four instruments. One pipeline run = one sequencing run.

| Platform | Native run output | Reads | QC source |
| --- | --- | --- | --- |
| Illumina NextSeq | `RunInfo.xml`, `RunParameters.xml`, `InterOp/`, `PrimaryAnalysisMetrics/` | demultiplexed on-instrument, left in place | FastQC + InterOp (%PF, %occupied, Q30, phasing, PhiX error — none of which reach the FASTQ) + on-instrument BCL Convert `Reports/` |
| Element AVITI | `RunParameters.json` | none in the run folder — `bases2fastq` is part of this pipeline's job | FastQC + `RunStats.json`/`Metrics.csv` |
| ONT PromethION P24 | `<sample>/<acquisition>/pod5_pass/{barcodeNN,mixed,unclassified}/`, `final_summary_*.txt` | `dorado basecaller --kit-name` (GPU) per POD5 folder → `dorado demux` per flow cell → per-sample uBAM (delivered as is; no FASTQ, by the team's choice) | MinKNOW `sequencing_summary_*.txt` + `report_*.json` → per-sample stats and the Nanopore flow cell / reads plots (`ONT_SUMMARY_STATS`) + run metrics + barcode counts |
| PacBio Revio | one dir per SMRT Cell (`1_A01/`): HiFi uBAM, `pb_formats/*.consensusreadset.xml` (sample names), `statistics/` (CCS, loading, lima) | already final, delivered as uBAM | the instrument's own reports only: per-barcode `barcodes.report.json`, per-cell metrics, lima counts, and the SMRT Cell / HiFi reads plots from `reports.zip` + `zmw_metrics.csv.gz` (`PACBIO_CELL_QC`); no task reads the BAM |

Only two of the four platforms hand you FASTQ, so read normalisation is per-platform and happens *inside* the pipeline. **POD5 are *not* cleanly split per barcode.** MinKNOW's `pod5_pass/mixed/` holds a large share of each barcode's reads (34% of barcode01's on `20260515_Plasmids_A_B`), so a folder is not a sample. Every POD5 folder is basecalled with `--kit-name` (kit from `final_summary`'s `protocol=`, or `--barcode_kit`), which tags reads with their barcode, and `DORADO_DEMUX --no-classify` splits each flow cell's reads on that tag. Samples are named by the sample sheet `alias`, else the barcode; `unclassified` is the Unassigned share.

**Dorado is the only genuinely GPU-accelerated tool here.** bases2fastq, FastQC, MultiQC and `ONT_SUMMARY_STATS` are all CPU-only; putting them in CUDA images would cost image size for nothing. GPU resources are set by name in `conf/base.config` (`accelerator = 1`) and `containerOptions` per profile in `nextflow.config` (`--gpus all` for Docker, `--nv` for Singularity/Apptainer); the Singularity and `cluster` profiles also run every container with `--contain -B /tmp -B /var/tmp` (see Gotchas). If reference-based QC is ever added, `nf-mod-parabricks` is the GPU path.

## Commands

```bash
# Smoke-test the wiring on the bundled marker-file fixtures (see conf/test.config)
nextflow run . -profile test -stub            # illumina
nextflow run . -profile test_aviti -stub
nextflow run . -profile test_nanopore -stub
nextflow run . -profile test_pacbio -stub

# Real miniaturised runs, sandpit-tokyo only (sbatch the head job)
nextflow run . -profile cluster,test_real_nanopore   # or _illumina, _aviti, _pacbio

nextflow run . -profile test_publish -stub    # --instrument/--run → mirrored outdir
nextflow run tests/copy_reads.nf --src <dir> # COPY_READS for real (see the copy-reads CI job)

# Real run
nextflow run . -profile cluster --instrument <brand>_<model>_<machine_ID> --run <run>   # head job via sbatch; see README

# Config / syntax checks — no linter or test framework is configured
nextflow config . -profile test
nextflow run . -resume

# Lint. Reports ERRORS ONLY — the deprecation/style warnings come from the
# language server (IDE), so a clean `nextflow lint` does not mean warning-free.
nextflow lint main.nf subworkflows/local/ modules/local/ nextflow.config conf/
```

Requires Nextflow `>=25.04.0` (manifest-enforced); developed against 25.10.2; the cluster's shared module is 26.04.6 (Nanopore and PacBio ran on it on 2026-11-05; Illumina and AVITI not yet).

Style rules the language server enforces, all of which this repo already follows — keep new code consistent:

- **`channel.*`, not `Channel.*`** — the capitalised factory namespace is deprecated.
- **Declare explicit closure parameters**; implicit `it` is deprecated. Prefix deliberately unused ones with `_` (`.map { _meta, zip -> zip }`).
- **Don't assign closures to local variables inside a function** — the language server does not resolve them and reports `X is not defined` as a hard error. Use a named top-level function instead (see `matchingFiles` in `main.nf`).

`assets/test/` holds marker files and zero-byte reads: enough to exercise platform detection, sample discovery, channel wiring and publishing under `-stub`, but the tools themselves are never run. Real miniaturised runs live on sandpit-tokyo only: `-profile cluster,test_real_<platform>` (fixtures listed in `conf/test_real.config`). They are not in CI.

## Architecture

`main.nf` detects the platform from marker files, dispatches to one of four subworkflows in `subworkflows/local/`, and declares the publish targets. Every platform subworkflow emits `qc` (flat files MultiQC can read) and `versions`, which is what makes `BATCH_REPORT` uniform. Reads differ by platform: Illumina, AVITI and PacBio emit `reads`, Nanopore emits `bam` instead (no FASTQ), and AVITI also emits `report` (Element's `<run>_QC.html`). Files for scientists that MultiQC does not read go out as `files` (Illumina/AVITI: FastQC HTML; AVITI: bases2fastq's `IndexAssignment.csv`, `UnassignedSequences.csv`, log, `RunManifestErrors.json`) or `barcoding` (Nanopore: dorado's `barcoding_summary.txt` per flow cell, renamed `<flowcell>.barcoding_summary.txt`), all published in `qc/` with the MultiQC data folder. AVITI's PhiX spike-in is not a sample: no FASTQ delivered, no FastQC. Revio BAMs are published renamed after the sample (`<sample>.hifi_reads.bam` + `.pbi`; `COPY_READS` takes the published names), the instrument names staying in the prep zone. Undetermined (Illumina), Unassigned (AVITI), unclassified (Nanopore) and unassigned (Revio, as `unassigned_<well>.hifi_reads.bam`, never QC'd) reads are always delivered (team decision 2026-11-13); `--include_undetermined` / `--include_failed` only add them to the QC (FastQC, sample rows). `meta` carries `[id, run, platform]`, where `id` is the sample (barcode, SMRT Cell, or Illumina/AVITI sample name) and `run` is constant for the whole pipeline run.

**Where runs come from and go.** Instruments write to `/mnt/instrument-prep-zone/sequencing/<brand>_<model>_<machine_ID>/<YYYY>/<run>` (`--prep_root`); results are published to the same relative path under `/mnt/instrument-data/sequencing` (`--publish_root`) for scientists. With `--instrument`/`--run`, `nextflow.config` globs the year directory itself (it has to: `outputDir` is fixed before the workflow starts) and leaves every match in `params.input`, comma-separated; `resolveRun()` in `main.nf` turns that into one directory or an error. The brand prefix must match the detected platform, and a destination that already has files is refused unless it is a `-resume` of the session recorded in `qc/nf-seqc_session.txt` (created atomically, so concurrent launches can't both pass) or `--overwrite`. The lookup goes into the internal `params.run_matches`, never `params.input`, so a user `--input` alongside `--instrument` is detectable. Publishing is `copy` (another filesystem from `work/`), and Illumina FASTQ and Revio BAM go through `COPY_READS` so they are task outputs. Inside the run folder the layout is "option C" (team decision 2026-11-13): reports (MultiQC, and Element's `<run>_QC.html` on AVITI) and reads at the top, everything else under `qc/`; see the `output {}` block. Scientific computing (Christiaan) owns these storage conventions; permissions/group ownership of the destination are still to be agreed with them.

Platform detection order matters: AVITI is checked before Illumina because both write a file called `RunParameters`, differing only in extension.

### Illumina run metrics mirror an existing internal tool

`INTEROP_RUN_METRICS` ([`nf-mod-interop`](https://github.com/EIT-GBI/nf-mod-interop), submodule `modules/interop/`) is a deliberate port of [EIT-GBI/Automate-Seq-Run-Metrics-Collection](https://github.com/EIT-GBI/Automate-Seq-Run-Metrics-Collection) (a cron-based collector for NextSeq 2000). **Keep `LOG_COLUMNS` in `collect_illumina_metrics.py` in sync with that repo's `utils.py`** — the whole point is that rows from either tool concatenate into one log. The metric set, column names and %Occupied vs %PF plot are intentionally identical, and the image pins InterOp/numpy/matplotlib to that repo's `uv.lock`. The script ships as a **module binary** (`run_metrics/resources/usr/bin/`), which is why `nextflow.config` sets `nextflow.enable.moduleBinaries = true` — remove that and the task fails with `command not found`.

Deliberate divergences, all documented in the script's module docstring:

- **One-row CSV per run, not an append to a central log.** One pipeline run is one sequencing run; appending to a shared file from a task would race and break `-resume`.
- **Nothing in the run directory is written.** The upstream `monitoring_complete.txt` marker exists to stop cron reprocessing; `-resume` covers that here.
- **Missing optional inputs degrade to `n/a`** rather than raising. Every metric is individually wrapped in `_safe()`, because patterned-flowcell-only fields (`percent_occupied`) are absent on other flowcells and one missing field should not cost the whole row.
- **Three upstream bugs are fixed here** — if they are fixed upstream, drop the notes: `_parse_flowcell` used `> 100`, labelling a 100-cycle P2 cartridge "P2 50" and contradicting its own docstring; and the sample-sheet path and `OutputFolder` were split on `os.sep`, which leaves a Windows path from `RunParameters.xml` unsplit when running on Linux (whole path in `SampleSheet`, `<dir>\<run>` in `Run_nr`).
- The MultiQC custom-content namespace is `Illumina InterOp`, not the pipeline name, since the module is reusable.

**Composition via git submodules.** Per-tool processes live in their own `EIT-GBI/nf-mod-*` repos under `modules/<tool>/`, one subdirectory per subcommand (`modules/fastqc/fastqc/main.nf` → `FASTQC_FASTQC`, `modules/dorado/basecaller/main.nf` → `DORADO_BASECALLER`). Always check whether an `nf-mod-*` repo already exists in the org before writing a process. `modules/local/` holds tools not yet usable from a repo — `bases2fastq`, `multiqc` (no repo), `ont_summary_stats` and `pacbio_cell_qc` (platform-specific; ship their scripts as module binaries), `bam_check` (samtools quickcheck + modified-base tags of every delivered long-read BAM, on the `nf-mod-samtools` image; a candidate subcommand for that repo) and `copy_reads` (plain `cp`, no container) — and each is marked for promotion. NanoPlot is no longer used: both long-read platforms take per-sample stats from the instrument's own reports. Whether a tool gets its own repo is decided case by case: anything reusable by another pipeline does.

**`modules.versions` is the source of truth for module versions, not the submodule pointers.** The flow is:

1. Renovate (config inherited from `EIT-GBI/nf-renovate-config`) opens a PR editing a `nf-mod-<tool>=<version>` line. The shared preset suppresses PATCH propagation and routes MAJOR bumps to the dashboard, so only MINOR bumps land automatically.
2. `sync-submodules.yml` fires on push to `main` touching `modules.versions` **only when the head commit message starts with `feat(deps)`**. It maps repo name → submodule path via `.gitmodules`, checks out tag `v<version>`, and commits `chore(submodules) sync pointers to modules.versions`.
3. `auto-release.yml` fires on that `chore(submodules)` commit and cuts a MINOR release.
4. `flag-major-updates.yml` runs weekday mornings and opens a deduplicated tracking issue per available major.

Never hand-edit a submodule pointer to change a version — edit `modules.versions` and let the sync workflow move it, or the two drift.

**Containers: never `:latest` outside `-profile dev`** (same scheme as `nf-dnaseq`). The `nf-mod-*` processes carry no `container` directive. Released repos ship `conf/module.config` assigning `ghcr.io/eit-gbi/nf-mod-<tool>:v<version>`, kept in step with the release by that repo's `release.yml`; `nextflow.config` includes those before `profiles` (currently fastqc, interop, dorado; nf-mod-samtools was dropped when Nanopore stopped making FASTQ), so moving the submodule moves the image. An nf-mod with no release yet would be pinned in the same section to the `sha-<short>` image of the exact commit its submodule points at. The `dev` profile points every nf-mod image back at `:latest` for testing unreleased module changes — add a line there for each new module. Local-module images (bases2fastq, multiqc, ont_summary_stats, pacbio_cell_qc, bam_check) are pinned in the `process {}` block after `profiles`.

**Config layering.** `nextflow.config` (manifest, params) → `conf/base.config` (resources) → `modules/*/conf/module.config` (images) → profiles (`dev` swaps nf-mod images to `:latest`, `docker`/`singularity` add GPU container options, `cluster` sets SLURM queues, `test*` set fixtures) → local-module images → `outputDir` last.

## Gotchas that cost time

These were all found empirically and all fail **silently**:

- **Publishing only works for task outputs.** Files that came from the input run directory cannot be published by the `output {}` block. That is why Illumina FASTQ and Revio HiFi BAM go through `COPY_READS` first, and why the PacBio CCS reports never appear under `--outdir`. Adding a publish target for them does nothing rather than erroring.
- **Be wary of `flatMap` before publishing.** A `flatMap`-derived channel once lost the task provenance the output block needs and published nothing at all. It is not universal: AVITI's `qc` channel flattens `BASES2FASTQ.out.stats` with `flatMap` and its files do publish (stub, 2026-11-13). Prefer `map`; if a published channel has to come from `flatMap`, check in a stub run that its files appear.
- **`outputDir = params.outdir` must come after the `profiles` block.** Profiles are merged as the config is parsed, so assigning it earlier captures the default and silently ignores any profile that overrides `params.outdir`. The template shipped it in the wrong place.
- **Process output globs keep their relative subpaths when published.** A tool that writes nested output (bases2fastq → `Samples/<project>/<sample>/`) will leak that structure into the results tree, so flatten inside the process before emitting.
- **`**/` requires at least one directory level**, so root-level files need their own glob pattern alongside it.
- **A module that hardcodes `publishDir` publishes twice** — once from the module, once from the `output {}` tree. `nf-mod-fastqc` v3 and `nf-mod-samtools` v1 removed it; `nf-mod-dorado` v1 removed it; `nf-mod-nanoplot` still has it (not used here). CI's `No publishDir` step fails on any uncommented `publishDir` in a `.nf` file, submodules included. Check for it before adopting any module.
- **`eval()` version commands run even under `-stub`, and a non-zero exit fails the whole task.** The existing modules get away with `tool --version | sed ...` because the pipeline's exit status is `sed`'s. Anything that can exit non-zero (e.g. `python3 -c "import x"`) needs `2>/dev/null || echo unknown`.
- **`error()` inside a channel operator prints your message and then a spurious `Cannot access first() element from an empty List`.** Raise from the workflow body where the value is plain, not from inside a `map`/`flatMap` — see `requireRunComplete()` in `main.nf`.
- **bases2fastq 2.4.0 ignores `--include-tile`.** It accepts the flag and still logs `240 tiles submitted for processing`. Subset an AVITI test run with `--exclude-tile` instead: `--exclude-tile '^(?!L[12]R06C01S1$).*$'` keeps one tile per lane.
- **Two `Demultiplex_Stats.csv` on a NextSeq 2000.** The one under `Analysis/<n>/Data/RunInstrumentAnalyticsMetrics/` has samples renamed `0001`… and indices `ELIDED`; only `Reports/` has real names. MultiQC's bclconvert module also needs a `RunInfo.xml` in the *same directory* as `Demultiplex_Stats.csv` (it logs an error and skips both otherwise) — `Reports/` has none, so the run root's is staged alongside. It does not read `Adapter_Metrics.csv` (1.35), and its "Clusters/run" is a share of *assigned* reads; the `% Reads` General Statistics column counts Undetermined too.
- **`-resume` does not see changes to a module binary.** Editing `modules/*/resources/usr/bin/*` (e.g. `collect_illumina_metrics.py`) leaves the task cached with its old output; the cache key covers the process script, not the files it calls. After changing one, run without `-resume` or delete that task's work directory.
- **Apptainer/Singularity containers share the node's `/dev/shm` unless `--contain` is set.** Nextflow always adds `--pid`, so a tool has the same PID in every container, and Intel OpenMP (dorado's bundled `libiomp5`) registers itself as `/dev/shm/__KMP_REGISTERED_LIB_<pid>_<uid>`. Dorado tasks starting together on one node raced on that file: `OMP: Error #15 ... unknown library already initialized` (exit 134) or a bus error (135), at random, and every abort left stale files that made the next run likelier to fail (EIT-GBI/nf-seqc#10). `runOptions = '--contain -B /tmp -B /var/tmp'` in the singularity and `cluster` profiles fixes it for every process; don't drop it, and keep the `/tmp` binds, because `--contain` alone gives a 64 MB tmpfs `/tmp`. The `--pid` default can only be turned off with `NXF_APPTAINER_NEW_PID_NAMESPACE=false` in the head job's environment, not from config.
- **`withLabel` beats a bare `cpus =`.** `conf/test.config` therefore uses `resourceLimits` to cap resources, because setting `process.cpus = 1` would not override `withLabel: process_high`.

## Known gaps

- `INTEROP_RUN_METRICS` has been run on real InterOp data only for NovaSeq and MiSeq folders (nf-mod-interop's tests). The **NextSeq 2000-only fields** (`Flowcell`, `Cycles`, `TotalReadsPF`, `RunStart`, `LoadingConc`) have never been read from a real NextSeq run.
- MultiQC parses FastQC, bases2fastq (`RunStats.json`), and BCL Convert (`Analysis/<n>/Data/Reports/`) output; the InterOp table and plot, the PacBio and Nanopore Run Metrics tables, both long-read platforms' per-sample General Statistics columns, and their flow cell / SMRT Cell and reads plots arrive as custom content (`*_mqc.*`). The long-read plots are one `*_mqc.json` per section per flow cell or cell, written by `ont_summary_stats.py` and `pacbio_cell_qc.py`; files sharing a section `id` merge by sample name, but MultiQC keeps only the last file's data for multi-dataset (`data_labels`) plots and heatmaps, so those scripts use neither. Sections are grouped with `parent_id` and ordered in `assets/multiqc_config.yml`. Nanopore's come from MinKNOW's `sequencing_summary_*.txt` via `ONT_SUMMARY_STATS`: MinKNOW's **live** basecalls and barcoding (hac on some runs), not the pipeline's re-basecalled reads, whose demux still gives the sample rows and % Reads. Not every run applies MinKNOW's Q filter (`20260515_Plasmids_A_B` has every read `passes_filtering=TRUE`). **The bases2fastq module only exists from MultiQC 1.29** — on 1.25.1 the AVITI report silently had no bases2fastq section — so never pin MultiQC below that. It has **no module for PacBio CCS reports**, so `pacbio_qc.nf` reads `summary.json`, `ccs_report.json`, and `loading.report.json` and `barcodes.report.json` (inside `reports.zip`) itself, and writes the PacBio Run Metrics table and the per-sample General Statistics columns. NanoPlot is not run on PacBio: the instrument's per-barcode report already has reads, yield, mean length and median QV. Software versions come from `metadata/<movie>.metadata.xml` (Revio ICS, SMRT Link), and `BATCH_REPORT` always adds nf-seqc's own, because MultiQC never starts if `software_versions.yml` is empty.
- The `cluster` profile targets the GBI sandpit SLURM clusters (checked on sandpit-tokyo, 2026-10-09: `cpu` and `gpu` partitions, GPU nodes `gpu:h100:8`). It switches Apptainer on itself and uses the shared cache `/mnt/lustre/containers/.staging/nf/cache/`. On tokyo, Apptainer exists only on the compute nodes, so the head job must be `sbatch`ed, from a `#!/bin/bash -l` script that runs `module load nextflow` (Java comes with the module).

## Release conventions

`manifest.version` in `nextflow.config` is bumped by CI, never by hand.

- Every PR to `main` must carry **exactly one** of `bump:patch`, `bump:minor`, `bump:major`; `pr-label-check.yml` fails it otherwise.
- On merge, `release.yml` bumps the version, commits `chore: bump version to X.Y.Z [skip ci]`, tags, and releases.
- Commit prefixes `feat(deps)` and `chore(submodules)` are load-bearing for the automation above — preserve them when squashing or rewriting.

`nf-illumina-qc` in the org was an earlier testing pipeline and may be deleted; nf-seqc supersedes it.
