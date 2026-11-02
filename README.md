# Processed Olink datasets

The master folder contains all aggregated datasets.

- `expression_matrix.csv` - rows = markers (harmonized Olink assay names), columns = samples. An assay measured on more than one panel keeps one row per panel (`name|panel`). PRINCE (row072) releases IL-10, IL-6 and LAMP3 twice per sample; these are kept as `|measurement 1 of 2` / `|measurement 2 of 2`, not averaged.
- `sample_metadata.csv` - one row per matrix column, same order (`sample_id` = column name): patient, timepoint, treatment and the outcome variables that were released (response, PFS/OS, irAE, etc.).
- `marker_annotation.csv` - original assay name, harmonized UniProt accession, panel, OlinkID.
- `README.txt` - scale, caveats and how the IDs/outcomes were derived.

`dataset_inventory.csv` lists every dataset with marker/sample/patient counts and which outcome columns exist.

`_overlap/` holds the cross-study comparison. Markers are matched by UniProt accession (names differ between sources, e.g. AGRG1 = ADGRG1, IL8 = CXCL8):
- `overlap_summary.txt` gives the headline numbers.
- `marker_presence_by_dataset.csv` is a 0/1 matrix of 1,494 proteins x 21 datasets.
- `pairwise_shared_markers.csv` and `pairwise_shared_markers_heatmap.png` show how many proteins each pair of datasets shares.
- `IO_core_91_markers.csv` lists the 91 IO-panel proteins shared by all 10 per-sample IO datasets.

Summary-only datasets (no per-sample values were released): row044 (median NPX by week), row078 (PFS hazard ratios), row096 (figure-level delta-NPX).

`scripts/` contains the scripts that regenerate everything: `python3 build_processed.py` then `python3 overlap_summary.py`. `verify_independent.py` re-reads every raw file without the processing code and compares every value with the matrices.

Accession merges: older Olink annotations use superseded UniProt IDs for CD6 (Q8WWJ7), FGF5 (Q8NF90) and IGLC2 (P0CG05); these are mapped to P30203, P12034 and P0DOY2. NT-proBNP (no accession in Olink files) is mapped to P16860.

`_master/` holds all 18 per-sample datasets combined, with values unchanged (not normalized): master_object.h5ad, master_matrix.csv.gz, master_metadata.csv, master_metadata_full.csv, master_features.csv and master_long.csv.gz. See _master/README.md.
