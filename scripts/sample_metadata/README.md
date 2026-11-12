# Sample metadata extractors

Standalone scripts that read the run folders the instruments write and produce
one CSV per batch with **one row per sample**: name, barcode/index and the
run-level details. Not part of the Nextflow pipeline; run them by hand.
Python 3.8+, standard library only.

| Script | Platform | Reads |
| --- | --- | --- |
| `pacbio_sample_metadata.py` (after Melissa's original) | PacBio Revio | `<cell>/metadata/<movie>.metadata.xml` |
| `illumina_sample_metadata.py` | Illumina NextSeq 1000/2000 | sample sheet (`Analysis/<n>/Data/Reports/SampleSheet.csv`, else the sheet in the run root), `RunInfo.xml`, `RunParameters.xml` |
| `aviti_sample_metadata.py` | Element AVITI | `RunManifest.csv` `[SAMPLES]`, `RunParameters.json` |
| `nanopore_sample_metadata.py` | ONT MinKNOW (P24) | per acquisition: `sample_sheet_*.csv`, `final_summary_*.txt`, `pod5_pass/barcode*/` |

All four take any mix of run folders, files and glob patterns, and `-o` for the
output CSV:

```bash
python3 illumina_sample_metadata.py /mnt/instrument-data/illumina_nextseq_1/261112_* -o nextseq.csv
python3 aviti_sample_metadata.py   /mnt/instrument-data/element_aviti_1/20260924_* -o aviti.csv
python3 nanopore_sample_metadata.py <experiment_dir> -o nanopore.csv
python3 pacbio_sample_metadata.py  <revio_run_dir> -o revio.csv
```

Platform notes:

- **PacBio.** A Revio `metadata/` folder also holds
  `*_basic_preview_4hr.metadata.xml` and `*_full_preview_23hr.metadata.xml`
  with the same samples. These are skipped when found through a folder or glob
  pattern, so each sample is listed once. A preview file passed by its full
  path is still read. Identical rows from several files (e.g. a downloaded
  `<movie>.metadata(1).xml`) are written once. The columns are Melissa's
  original ones plus `BarcodePair` at the end: `BarcodeName` keeps only the left
  barcode, so `BarcodePair` is what tells asymmetric pairs (`bc2001--bc2002`)
  apart.
- **Illumina.** v2 (`[BCLConvert_Data]`, plus `[Cloud_Data]` when present) and v1
  (`[Data]`) sheets are both read. With several `Analysis/<n>` re-runs, the
  highest `<n>` wins. Before the on-instrument analysis has run, the sheet is the
  one in the run folder that `RunParameters.xml`'s `SampleSheetFilePath` names.
- **AVITI.** A sample with several index pairs (PhiX has four) gives one row per
  pair, as in the manifest. `FilterMask` is shown because some runs set their
  own (e.g. R1-only), which makes runs hard to compare.
- **Nanopore.** Barcodes come from the sample sheet's `barcode`/`alias` columns.
  When the sheet has none, they come from the `barcodeNN` folders MinKNOW created
  under `pod5_pass/` (or `pod5/` when pass/fail splitting is off). Those show the
  barcodes MinKNOW saw during the run, not the final demux (`BarcodeSource` says
  which). Each restarted acquisition on a flow cell gets its own rows.
  Acquisitions with no reads and no `final_summary` (a start that was aborted
  straight away, like 1359 and 1405 on `20260515_Plasmids_A_B`) are skipped. An
  acquisition with reads but no barcode folders gets a single `none` row.

Tested 2026-11-12 on sandpit-tokyo: NextSeq `261006_VH02427_108_AAJJW55M5` and
`261112_VH02427_122_222H77WNX`, AVITI `SEQEBA000015`, P24
`20260515_Plasmids_A_B` (mini), Revio `r84334_20260817_141703` (1_A01).
