# Soil–air thermal contrasts and water-use efficiency in Northern Hemisphere vegetation

This package supports the analyses of vegetation north of 20°N during 2001–2024. It contains spatial-block sufficient statistics for the primary model and temporal-allocation comparisons, complete result tables, figure-source tables, and dataset citations.

## Reproduce the estimates

Run from the repository directory:

```bash
python -m pip install -r requirements.txt
python code/reproduce_models.py
```

The script refits the joint GPP, T and ET models from sufficient statistics, propagates covariance to WUE and its components, aligns spatial blocks for phase and allocation contrasts, and applies Benjamini–Hochberg adjustment within the original families. It reproduces all 60 primary-family estimates and 252 temporal-allocation/support estimates and compares coefficients, standard errors, confidence intervals and P values with the reference tables. Both initially examined exposures remain in the primary family. Computed tables and numerical checks are written to `output/`.

The tested environment used Python 3.14.6, NumPy 2.5.2, pandas 3.0.5 and SciPy 1.18.1. Run the additional mathematical checks with `python -m unittest discover -s tests`.

## Plotting

```bash
python code/plot_figures.py
```

The default command reproduces Figure 5 from the included source tables as PNG, PDF, SVG and TIFF files in `output/figures/`. The tested plotting environment used matplotlib 3.11.1 and Pillow 12.3.0. Large spatial arrays and generated figures are kept outside this lightweight repository. To reproduce Figure S1 with the separate spatial supplement:

```bash
python code/plot_figures.py --gis-classes /path/to/GIS_retained_classes.npz
```

The two optional spatial arrays and their checksums are listed in `metadata/external_spatial_inputs.csv`. They are not needed for the 312 model estimates. Figure S2's retention map uses the separate `spatial_retention.npz`; its summary tables are included here. Figure S5 can be reproduced by supplying the paired monthly observations:

```bash
python code/plot_figures.py --monthly-pairs /path/to/monthly_pairs.csv
```

The optional table requires `endpoint`, `site`, `tower` and `product` columns. It contains one row per matched site-month, with GPP in g C m⁻² d⁻¹ and ET in mm d⁻¹. The plotting function checks its site-balanced correlations and RMSE against the distributed aggregate metrics before drawing. Paired tower observations are not bundled; aggregate monthly GPP and ET comparison metrics are in `data/source_tables/tower_metrics.csv`.

Figure 1's descriptive means, area fractions and latitude profiles are provided in `data/source_tables/Figure_1_thermal_summary.csv` and `Figure_1_latitude_profile.csv`. Their weighting and source identities are specified in `metadata/Figure_1_summary_definition.json`. Regenerating those summaries from cell-level means requires all three options: `--thermal-context`, `--cell-metadata` and `--gis-classes`.

## Package contents

- `data/model_inputs/`: block-level sufficient-statistic matrices for the joint three-response fits.
- `data/reference/`: the complete primary and allocation reference results.
- `data/source_tables/`: source tables for the manuscript and supporting figures; figure identifiers are mapped in `figure_data_map.csv`.
- `data/complete_results/`: full historical model families and direct-contrast results, including non-significant estimates.
- `data/exploratory_tower/`: Supplementary Data E1, the previously completed phenological-phase and WUE diagnostics, with methods in its README. The manuscript's primary tower comparison concerns GPP and ET; the exploratory diagnostics are retained here as a separate record.
- `metadata/`: dataset access routes, station citations, matrix definitions and column descriptions.
- `file_manifest.csv`: SHA-256 checksums for the distributed files.

Three larger CSV tables use lossless gzip compression. Read them directly with `pandas.read_csv('filename.csv.gz')` or decompress them with an archive utility. No rows or numerical precision have been removed. Their uncompressed sizes and SHA-256 checksums are recorded in `metadata/compressed_tables.csv`.

GPP/ET is ecosystem water-use efficiency, GPP/T is transpiration-based efficiency, and T/ET is the transpiration fraction. Coefficients are on the log-response scale. Figures display `100 × beta`; exact percentage changes use `100 × [exp(beta) − 1]`. Ordinary 95% spatial-cluster intervals and multiplicity-adjusted P values answer different inferential questions. Empty numeric CSV fields indicate missing values.

## Reproduction scope

The model inputs already contain the cell-within-transformation, area weights and year-indicator cross-products. They permit re-estimation of the listed models without distributing the large gridded or tower archives. Rebuilding these matrices from original products is a separate preprocessing workflow described in the Supporting Information. The package does not contain raw satellite, reanalysis or tower observations. Access routes and product citations are given in `metadata/dataset_access.csv` and `metadata/DATA_CREDITS.md`.

The analysis grid is 0.05°; harmonisation does not increase the native resolution of coarser products. Spatial clusters, not repeated gridded observations, govern the uncertainty estimates. Moisture-background labels describe groups of precipitation divided by absolute potential evaporation and do not denote independently validated climatic zones.

## Rights and citation

The existing repository software licence remains applicable to its code. This package does not apply a new licence to third-party datasets. Dataset citations and source licences must accompany reuse of the corresponding data. Cite the article and the exact archived repository version used when a persistent release identifier is available.
