# nf-__NAME__

<!-- TEMPLATE:START -->
> **Template repo.** Create a new pipeline via GitHub's "Use this template", naming it `nf-<name>`. The `init-from-template` workflow auto-runs on the first push: it derives `<name>` from the repo name, replaces `__NAME__` / `__NAME_UPPER__` tokens, strips this block, and removes itself.
>
> Manual fallback (e.g. running locally before pushing, or if the repo name doesn't start with `nf-`): `./scripts/init.sh <name>`.
>
> Afterwards, edit the remaining `TODO`s in `main.nf`, `nextflow.config`, and `conf/test.config` for the pipeline's actual inputs, processes, and parameters.
<!-- TEMPLATE:END -->

Nextflow pipeline for __NAME__.

## Quickstart

```bash
nextflow run eit-gbi/nf-__NAME__ -profile test,docker
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
