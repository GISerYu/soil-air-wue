# Supplementary Data E1: exploratory phase-specific tower comparisons

These tables preserve the previously completed phase-specific GPP, ET and WUE comparisons. They are separate from the monthly GPP and ET evaluation reported in Figure S5. All phase-by-variable combinations are retained; their inclusion does not enlarge the scope of the primary tower comparison.

## Methods

The comparisons used the same GOSIF-defined seasonal thirds as the gridded models. GPP and ET shared valid days, with at least 80% day-weighted phase coverage and three paired years per site. Identical duplicate daily records were collapsed and conflicting records excluded. The matched sample comprised 109 sites, 770 site-years and 2,310 phase records.

Original-scale anomalies were departures from the paired-year site-by-phase mean. Log-WUE anomalies were calculated by taking the natural logarithm of GPP/ET for positive phase fluxes before centring. Site-balanced Pearson correlations and 95% percentile intervals used 1,000 whole-site bootstrap resamples, with seed 20261010. Spearman correlations used conventional pooled ranks and whole-site bootstrapping. These diagnostics describe tower–product agreement, not the conditional soil–air thermal regression coefficient.

## Files

- `phase_anomalies.csv`: compact correlation summary on original and logarithmic scales.
- `tower_phase_metrics.csv`: original phase-level and phase-anomaly metrics.
- `tower_diagnostics_metrics.csv`: full raw/log, Pearson/rank, spread and slope comparisons.
- `tower_temporal_scale_metrics.csv`: corresponding comparisons by time scale.
- `tower_phase_distributions.csv`: phase-distribution summaries.
- `tower_leave_one_site_summary.csv`: leave-one-site influence summaries.
- `tower_site_metrics.csv`: site-level aggregate comparison metrics.

GPP is expressed as g C m⁻² d⁻¹, ET as mm d⁻¹, and WUE as g C kg⁻¹ H₂O. Log-scale values are natural logarithms of positive ratios or their centred anomalies. Empty numeric fields denote missing values. The summary statistics do not redistribute the underlying daily or monthly tower observations.
