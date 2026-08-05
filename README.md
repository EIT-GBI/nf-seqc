# nf-seqc


Nextflow pipeline for seqc.

## Quickstart

```bash
nextflow run eit-gbi/nf-seqc -profile test,docker
```

## Layout

- `main.nf` — entry workflow
- `nextflow.config` — manifest, params, profiles
- `conf/` — split configuration (`base.config`, `test.config`)
- `modules/` — per-tool processes pulled in as git submodules from `nf-mod-*` repos
- `bin/` — helper scripts added to `$PATH` inside processes
- `assets/` — small static inputs the pipeline ships with

## Adding a module

```bash
git submodule add https://github.com/eit-gbi/nf-mod-<tool>.git modules/<tool>
```

Then in `main.nf`:

```groovy
include { TOOL } from './modules/<tool>/main.nf'
```
