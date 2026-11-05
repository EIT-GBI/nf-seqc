# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this pipeline is for

nf-seqc takes the **complete output directory of one sequencing run** and reports QC for the **batch/run**, across four instruments. One pipeline run = one sequencing run.

| Platform | Native run output | Reads | QC source |
| --- | --- | --- | --- |
| Illumina NextSeq | `RunInfo.xml`, `RunParameters.xml`, `InterOp/`, `PrimaryAnalysisMetrics/` | demultiplexed on-instrument, left in place | FastQC + InterOp (%PF, %occupied, Q30, phasing, PhiX error — none of which reach the FASTQ) + on-instrument BCL Convert `Reports/` |
| Element AVITI | `RunParameters.json` | none in the run folder — `bases2fastq` is part of this pipeline's job | FastQC + `RunStats.json`/`Metrics.csv` |
| ONT PromethION P24 | `<sample>/<acquisition>/pod5_pass/{barcodeNN,mixed,unclassified}/`, `final_summary_*.txt` | `dorado basecaller --kit-name` (GPU) per POD5 folder → `dorado demux` per flow cell → uBAM → `samtools fastq` | NanoPlot + MinKNOW `report_*.json` run metrics + barcode counts |
| PacBio Revio | one dir per SMRT Cell (`1_A01/`): HiFi uBAM, `pb_formats/*.consensusreadset.xml` (sample names), `statistics/` (CCS, loading, lima) | already final, delivered as uBAM | NanoPlot + per-cell run metrics + lima barcode counts |

Only two of the four platforms hand you FASTQ, so read normalisation is per-platform and happens *inside* the pipeline. **POD5 are *not* cleanly split per barcode.** MinKNOW's `pod5_pass/mixed/` holds a large share of each barcode's reads (34% of barcode01's on `20260515_Plasmids_A_B`), so a folder is not a sample. Every POD5 folder is basecalled with `--kit-name` (kit from `final_summary`'s `protocol=`, or `--barcode_kit`), which tags reads with their barcode, and `DORADO_DEMUX --no-classify` splits each flow cell's reads on that tag. Samples are named by the sample sheet `alias`, else the barcode; `unclassified` is the Unassigned share.

**Dorado is the only genuinely GPU-accelerated tool here.** bases2fastq, FastQC, NanoPlot, samtools and MultiQC are all CPU-only; putting them in CUDA images would cost image size for nothing. GPU resources are set by name in `conf/base.config` (`accelerator = 1`) and `containerOptions` per profile in `nextflow.config` (`--gpus all` for Docker, `--nv` for Singularity). If reference-based QC is ever added, `nf-mod-parabricks` is the GPU path.

## Commands

```bash
# Smoke-test the wiring on the bundled marker-file fixtures (see conf/test.config)
nextflow run . -profile test -stub            # illumina
nextflow run . -profile test_aviti -stub
nextflow run . -profile test_nanopore -stub
nextflow run . -profile test_pacbio -stub

# Real miniaturised runs, sandpit-tokyo only (sbatch the head job)
nextflow run . -profile cluster,test_real_nanopore   # or _illumina, _aviti, _pacbio

# Real run
nextflow run . -profile cluster --input <run_dir> --outdir <dir>   # head job via sbatch; see README

# Config / syntax checks — no linter or test framework is configured
nextflow config . -profile test
nextflow run . -resume

# Lint. Reports ERRORS ONLY — the deprecation/style warnings come from the
# language server (IDE), so a clean `nextflow lint` does not mean warning-free.
nextflow lint main.nf subworkflows/local/ modules/local/ nextflow.config conf/
```

Requires Nextflow `>=25.04.0` (manifest-enforced); developed against 25.10.2.

Style rules the language server enforces, all of which this repo already follows — keep new code consistent:

- **`channel.*`, not `Channel.*`** — the capitalised factory namespace is deprecated.
- **Declare explicit closure parameters**; implicit `it` is deprecated. Prefix deliberately unused ones with `_` (`.map { _meta, zip -> zip }`).
- **Don't assign closures to local variables inside a function** — the language server does not resolve them and reports `X is not defined` as a hard error. Use a named top-level function instead (see `matchingFiles` in `main.nf`).

`assets/test/` holds marker files and zero-byte reads: enough to exercise platform detection, sample discovery, channel wiring and publishing under `-stub`, but the tools themselves are never run. Real miniaturised runs live on sandpit-tokyo only: `-profile cluster,test_real_<platform>` (fixtures listed in `conf/test_real.config`). They are not in CI.

## Architecture

`main.nf` detects the platform from marker files, dispatches to one of four subworkflows in `subworkflows/local/`, and declares the publish targets. Each platform subworkflow emits the same three channels — `reads`, `qc` (flat files MultiQC can read), `versions` — which is what makes `BATCH_REPORT` uniform. `meta` carries `[id, run, platform]`, where `id` is the sample (barcode, SMRT Cell, or Illumina/AVITI sample name) and `run` is constant for the whole pipeline run.

Platform detection order matters: AVITI is checked before Illumina because both write a file called `RunParameters`, differing only in extension.

### Illumina run metrics mirror an existing internal tool

`INTEROP_RUN_METRICS` ([`nf-mod-interop`](https://github.com/EIT-GBI/nf-mod-interop), submodule `modules/interop/`) is a deliberate port of [EIT-GBI/Automate-Seq-Run-Metrics-Collection](https://github.com/EIT-GBI/Automate-Seq-Run-Metrics-Collection) (a cron-based collector for NextSeq 2000). **Keep `LOG_COLUMNS` in `collect_illumina_metrics.py` in sync with that repo's `utils.py`** — the whole point is that rows from either tool concatenate into one log. The metric set, column names and %Occupied vs %PF plot are intentionally identical, and the image pins InterOp/numpy/matplotlib to that repo's `uv.lock`. The script ships as a **module binary** (`run_metrics/resources/usr/bin/`), which is why `nextflow.config` sets `nextflow.enable.moduleBinaries = true` — remove that and the task fails with `command not found`.

Deliberate divergences, all documented in the script's module docstring:

- **One-row CSV per run, not an append to a central log.** One pipeline run is one sequencing run; appending to a shared file from a task would race and break `-resume`.
- **Nothing in the run directory is written.** The upstream `monitoring_complete.txt` marker exists to stop cron reprocessing; `-resume` covers that here.
- **Missing optional inputs degrade to `n/a`** rather than raising. Every metric is individually wrapped in `_safe()`, because patterned-flowcell-only fields (`percent_occupied`) are absent on other flowcells and one missing field should not cost the whole row.
- **Three upstream bugs are fixed here** — if they are fixed upstream, drop the notes: `_parse_flowcell` used `> 100`, labelling a 100-cycle P2 cartridge "P2 50" and contradicting its own docstring; and the sample-sheet path and `OutputFolder` were split on `os.sep`, which leaves a Windows path from `RunParameters.xml` unsplit when running on Linux (whole path in `SampleSheet`, `<dir>\<run>` in `Run_nr`).
- The MultiQC custom-content namespace is `Illumina InterOp`, not the pipeline name, since the module is reusable.

**Composition via git submodules.** Per-tool processes live in their own `EIT-GBI/nf-mod-*` repos under `modules/<tool>/`, one subdirectory per subcommand (`modules/fastqc/fastqc/main.nf` → `FASTQC_FASTQC`, `modules/dorado/basecaller/main.nf` → `DORADO_BASECALLER`). Always check whether an `nf-mod-*` repo already exists in the org before writing a process. `modules/local/` holds tools not yet usable from a repo — `bases2fastq`, `multiqc` (no repo) and `nanoplot` (`nf-mod-nanoplot` exists but only takes aligned BAM and has no release) — and each is marked for promotion. Whether a tool gets its own repo is decided case by case: anything reusable by another pipeline does.

**`modules.versions` is the source of truth for module versions, not the submodule pointers.** The flow is:

1. Renovate (config inherited from `EIT-GBI/nf-renovate-config`) opens a PR editing a `nf-mod-<tool>=<version>` line. The shared preset suppresses PATCH propagation and routes MAJOR bumps to the dashboard, so only MINOR bumps land automatically.
2. `sync-submodules.yml` fires on push to `main` touching `modules.versions` **only when the head commit message starts with `feat(deps)`**. It maps repo name → submodule path via `.gitmodules`, checks out tag `v<version>`, and commits `chore(submodules) sync pointers to modules.versions`.
3. `auto-release.yml` fires on that `chore(submodules)` commit and cuts a MINOR release.
4. `flag-major-updates.yml` runs weekday mornings and opens a deduplicated tracking issue per available major.

Never hand-edit a submodule pointer to change a version — edit `modules.versions` and let the sync workflow move it, or the two drift.

**Containers: never `:latest` outside `-profile dev`** (same scheme as `nf-dnaseq`). The `nf-mod-*` processes carry no `container` directive. Released repos ship `conf/module.config` assigning `ghcr.io/eit-gbi/nf-mod-<tool>:v<version>`, kept in step with the release by that repo's `release.yml`; `nextflow.config` includes those before `profiles` (currently fastqc, samtools, interop), so moving the submodule moves the image. An unreleased nf-mod (dorado) is pinned in the same section to the `sha-<short>` image of the exact commit the submodule points at. The `dev` profile points every nf-mod image back at `:latest` for testing unreleased module changes — add a line there for each new module. Local-module images (bases2fastq, nanoplot, multiqc) are pinned in the `process {}` block after `profiles`.

**Config layering.** `nextflow.config` (manifest, params) → `conf/base.config` (resources) → `modules/*/conf/module.config` + the unreleased-module pins (images) → profiles (`dev` swaps nf-mod images to `:latest`, `docker`/`singularity` add GPU container options, `cluster` sets SLURM queues, `test*` set fixtures) → local-module images → `outputDir` last.

## Gotchas that cost time

These were all found empirically and all fail **silently**:

- **Publishing only works for task outputs.** Files that came from the input run directory cannot be published by the `output {}` block. That is why Illumina FASTQ, Revio HiFi BAM and the PacBio CCS reports never appear under `--outdir` — they are instrument-produced and stay in the run folder. Adding a publish target for them does nothing rather than erroring.
- **Never `flatMap` a channel before publishing it.** A `flatMap`-derived channel loses the task provenance the output block needs and publishes nothing at all. `map` and raw tuple channels are fine.
- **`outputDir = params.outdir` must come after the `profiles` block.** Profiles are merged as the config is parsed, so assigning it earlier captures the default and silently ignores any profile that overrides `params.outdir`. The template shipped it in the wrong place.
- **Process output globs keep their relative subpaths when published.** A tool that writes nested output (bases2fastq → `Samples/<project>/<sample>/`) will leak that structure into the results tree, so flatten inside the process before emitting.
- **`**/` requires at least one directory level**, so root-level files need their own glob pattern alongside it.
- **A module that hardcodes `publishDir` publishes twice** — once from the module, once from the `output {}` tree. `nf-mod-fastqc` v3 and `nf-mod-samtools` v1 removed it; `nf-mod-nanoplot` still has it (`nf-mod-dorado` drops it in its v1 PR). Check for it before adopting any module.
- **`eval()` version commands run even under `-stub`, and a non-zero exit fails the whole task.** The existing modules get away with `tool --version | sed ...` because the pipeline's exit status is `sed`'s. Anything that can exit non-zero (e.g. `python3 -c "import x"`) needs `2>/dev/null || echo unknown`.
- **`error()` inside a channel operator prints your message and then a spurious `Cannot access first() element from an empty List`.** Raise from the workflow body where the value is plain, not from inside a `map`/`flatMap` — see `requireRunComplete()` in `main.nf`.
- **bases2fastq 2.4.0 ignores `--include-tile`.** It accepts the flag and still logs `240 tiles submitted for processing`. Subset an AVITI test run with `--exclude-tile` instead: `--exclude-tile '^(?!L[12]R06C01S1$).*$'` keeps one tile per lane.
- **Two `Demultiplex_Stats.csv` on a NextSeq 2000.** The one under `Analysis/<n>/Data/RunInstrumentAnalyticsMetrics/` has samples renamed `0001`… and indices `ELIDED`; only `Reports/` has real names. MultiQC's bclconvert module also needs a `RunInfo.xml` in the *same directory* as `Demultiplex_Stats.csv` (it logs an error and skips both otherwise) — `Reports/` has none, so the run root's is staged alongside. It does not read `Adapter_Metrics.csv` (1.35), and its "Clusters/run" is a share of *assigned* reads; the `% Reads` General Statistics column counts Undetermined too.
- **`-resume` does not see changes to a module binary.** Editing `modules/*/resources/usr/bin/*` (e.g. `collect_illumina_metrics.py`) leaves the task cached with its old output; the cache key covers the process script, not the files it calls. After changing one, run without `-resume` or delete that task's work directory.
- **`withLabel` beats a bare `cpus =`.** `conf/test.config` therefore uses `resourceLimits` to cap resources, because setting `process.cpus = 1` would not override `withLabel: process_high`.

## Known gaps

- `nf-mod-dorado` has **no git tag or release** yet. Its v1 is [EIT-GBI/nf-mod-dorado#1](https://github.com/EIT-GBI/nf-mod-dorado/pull/1), opened from the `paulinefx` fork because the user's account is read-only there; Cristian maintains it. The submodule points at that PR's commit, so it cannot be pinned in `modules.versions`, and the image stays on `sha-9f93e96` (the last build of `main`: same dorado 0.8.3 and sup v5.0.0) because the PR's `conf/module.config` pins a `v0.0.0` that does not exist until the release. Once released: move the submodule to the tag, `includeConfig 'modules/dorado/conf/module.config'`, drop the sha pin, and add `nf-mod-dorado=<version>` to `modules.versions`.
- `INTEROP_RUN_METRICS` has been run on real InterOp data only for NovaSeq and MiSeq folders (nf-mod-interop's tests). The **NextSeq 2000-only fields** (`Flowcell`, `Cycles`, `TotalReadsPF`, `RunStart`, `LoadingConc`) have never been read from a real NextSeq run.
- MultiQC parses FastQC, bases2fastq (`RunStats.json`), BCL Convert (`Analysis/<n>/Data/Reports/`) and NanoStat output; the InterOp table and plot, and the PacBio and Nanopore Run Metrics tables, arrive as custom content (`*_mqc.*`). **The bases2fastq module only exists from MultiQC 1.29** — on 1.25.1 the AVITI report silently had no bases2fastq section — so never pin MultiQC below that. It has **no module for PacBio CCS reports**, so `pacbio_qc.nf` reads `summary.json`, `ccs_report.json` and `loading.report.json` (inside `reports.zip`) itself and writes the PacBio Run Metrics table.
- The `cluster` profile targets the GBI sandpit SLURM clusters (checked on sandpit-tokyo, 2026-10-09: `cpu` and `gpu` partitions, GPU nodes `gpu:h100:8`). It switches Apptainer on itself and uses the shared cache `/mnt/lustre/containers/.staging/nf/cache/`. On tokyo, Apptainer exists only on the compute nodes, so the head job must be `sbatch`ed, with Java put on `PATH` by hand (SDKMAN is only loaded by interactive shells).

## Release conventions

`manifest.version` in `nextflow.config` is bumped by CI, never by hand.

- Every PR to `main` must carry **exactly one** of `bump:patch`, `bump:minor`, `bump:major`; `pr-label-check.yml` fails it otherwise.
- On merge, `release.yml` bumps the version, commits `chore: bump version to X.Y.Z [skip ci]`, tags, and releases.
- Commit prefixes `feat(deps)` and `chore(submodules)` are load-bearing for the automation above — preserve them when squashing or rewriting.

`nf-illumina-qc` in the org was an earlier testing pipeline and may be deleted; nf-seqc supersedes it.
