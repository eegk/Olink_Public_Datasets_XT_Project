# Olink master object (not normalized)

These files combine the 18 per-sample datasets from `processed/`: 2,144 samples, 838 patients and 1,494 proteins. Coffey 2024 (row 62) was added from the supplementary file you downloaded. Every value is copied unchanged from the per-study matrices, which in turn match the raw downloads.

**Nothing is normalized.** Scales differ between datasets (see `value_scale`), so check the following before combining:
- Most datasets are NPX (log2).
- Gao (row 32) is centred NPX as deposited.
- Yu (row 103) is on a linear scale.
- Batch/plate effects between studies are untouched.

The summary-only datasets (rows 44, 78 and 96) are not included because they have no per-sample values. They remain in `processed/`.

## Files

| File | Contents |
|---|---|
| `master_object.h5ad` | AnnData object: `X` = samples x proteins (raw values, NaN = not measured), `obs` = master metadata, `var` = protein features, `layers['measured']` = 0/1 mask. Open in Python with `anndata.read_h5ad`, or in R with `zellkonverter::readH5AD` or `anndata::read_h5ad`. |
| `master_matrix.csv.gz` | Proteins (rows, UniProt IDs) x samples (columns, `sample_uid`). |
| `master_metadata.csv` | One row per sample, in the same order as the matrix columns, with harmonized core columns (below). |
| `master_metadata_full.csv` | The core columns plus every original metadata column from each study (prefixed `orig__`). Same-named columns from different studies share one column. |
| `master_features.csv` | One row per protein: name, number of datasets, number of samples measured, `in_IO_core_91`, panels, and how many alternative assays were left out of the matrix. |
| `master_long.csv.gz` | Lossless long table of every assay value (including alternative assays). `used_in_master_matrix` marks the assay used in the wide matrix. |

## One assay per protein in the matrix

A few studies measured the same protein with more than one assay (several panels, or PRINCE's unlabelled duplicates). The wide matrix uses one assay per protein per dataset, chosen in this order: Immuno-Oncology panel > Inflammation > Oncology > Immune Response > other, then file order. For PRINCE, the matrix uses measurement 1. The other assays (6,977 values, mainly CITN-10 and Johnson multi-panel duplicates) are kept in `master_long.csv.gz`.

## Harmonized metadata columns

- **`sample_uid` / `patient_uid`**: `<dataset prefix>__<original id>`, so they are unique across studies.
- **`group`**: `patient` or `healthy control` (Yu row 103 has 10 healthy donors).
- **`disease`, `disease_detail`, `treatment`, `treatment_detail`**: study-level values, plus per-patient detail where the study gives it.
- **`sample_type`**: plasma or serum, as stated in the paper or deposit. Coffey (row 62) is bone marrow plasma.
- **`olink_panels`, `value_scale`**: which panels were run and what scale the values are on.
- **`timepoint_original`, `timepoint_std`**: the study's own timepoint label, and a harmonized version (`baseline`, `on-treatment`, `on/post-treatment`, `end of treatment`, `unknown`, `healthy control`).
- **`is_baseline`, `is_primary_baseline`**: every pre-treatment sample, and one per patient (BL/Pre > Lead-in Day 1 > Cycle 1 Day 1).
- **`response_std`, `response_definition`**: R/NR using each study's own response definition. Definitions differ between studies, so always check `response_definition`. Coffey uses IMWG myeloma criteria (best response of PR or better = R); `bor_std` there is written as "IMWG <category>" and `orr_recist` is left blank.
- **`bor_std`, `orr_recist`**: best overall response (CR/PR/SD/PD/NE; uPR = unconfirmed PR) where reported, and R/NR derived from it (CR or PR = R; NE left blank).
- **`pfs_months`, `pfs_event`, `pfs_type`, `os_months`, `os_event`**: survival in months with event = 1. PORTER reports radiographic PFS (rPFS). PRINCE days were divided by 30.4375 (raw days remain in `orig__clinical_observation_*`).
- **`irae`**: immune-related adverse event (Nunez only).
- **`age`, `sex`**: as given. Gao's age and sex are coded categories.
- **`qc_warning`, `n_assays_measured`, `n_proteins_in_master_matrix`**: QC flag from the source and how many values each sample has.

## Assumptions to be aware of

- **Johnson (row 12):** B1W1 (block 1, week 1) is treated as baseline.
- **Coffey (row 62):** Screening is treated as baseline; Day 90 and Day 365 are after transplant. Per-patient PFS was not released.
- **CITN-13 (row 86):** C01 is treated as on-treatment.
- **PORTER (row 105):**
  - Cohorts A and B: both BL and Cycle1.Day1 count as baseline.
  - Cohort C: Lead-in.Day1 is the baseline.
- **Gao (row 32):** only the validation cohort (the "24R" samples) is known to be pre-treatment.
- **Nunez (row 70):** three patients appear twice in the source, with identical T1/T2 values. They have one `patient_uid`; the series number is kept in `orig__series`.

Rebuild with `scripts/build_master.py`. Running it needs pandas, plus `anndata` for the .h5ad file.
