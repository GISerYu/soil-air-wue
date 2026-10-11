# Dataset credits and access

The article's data-source section and `dataset_access.csv` identify the products used. The supplied numerical tables and sufficient statistics are derived outputs, not copies of the underlying raster or tower archives.

## Gridded products

GOSIF: Li, X. and Xiao, J. (2019), Remote Sensing 11, 517, https://doi.org/10.3390/rs11050517.

PML-V2.2c: https://doi.org/10.11888/Terre.tpdc.303314. Xu et al. (2026), Earth System Science Data 18, 5663–5695, https://doi.org/10.5194/essd-18-5663-2026; Zhang et al. (2019), Remote Sensing of Environment 222, 165–182, https://doi.org/10.1016/j.rse.2018.12.031.

GLEAM4: Miralles et al. (2025), Scientific Data 12, 416, https://doi.org/10.1038/s41597-025-04610-y; Miralles et al. (2011), Hydrology and Earth System Sciences 15, 453–469, https://doi.org/10.5194/hess-15-453-2011. Module references: Zhong et al. (2022), https://doi.org/10.5194/hess-26-5647-2022; Hulsman et al. (2023), https://doi.org/10.1029/2022WR033731; Koppa et al. (2022), https://doi.org/10.1038/s41467-022-29543-7. GLEAM data are distributed at https://www.gleam.eu/ under the provider's registration and use conditions.

ERA5-Land monthly averaged data: https://doi.org/10.24381/cds.68d2bb30. Muñoz-Sabater et al. (2021), Earth System Science Data 13, 4349–4383, https://doi.org/10.5194/essd-13-4349-2021.

MODIS Collection 6.1: MCD12C1, https://doi.org/10.5067/MODIS/MCD12C1.061; MCD12Q2, https://doi.org/10.5067/MODIS/MCD12Q2.061.

## Eddy-covariance products

FLUXNET data products were produced and harmonized by eddy covariance regional networks and data processing centers, including AmeriFlux, ChinaFlux, European Fluxes Database, ICOS, JapanFlux, KoFlux, OzFlux, SAEON, and TERN. These products also include a modified version of ERA5 hourly data provided by the Copernicus Climate Change Service.

ChinaFlux data are from the National ecosystem science data center https://www.nesdc.org.cn. Relevant network references are Yu et al. (2006), https://doi.org/10.1016/j.agrformet.2006.02.011, and Ueyama et al. (2025), https://doi.org/10.5194/essd-17-3807-2025. The ONEFlux pipeline is described by Pastorello et al. (2020), https://doi.org/10.1038/s41597-020-0534-3. The product archives used here are ONEFlux v1.3 r1/r2, not the FLUXNET2015 release.

Station-specific archive names, identifiers, provider citations and acknowledgments are provided in `tower_dataset_citations.csv` and `tower_dataset_citations.bib`. These products use CC BY 4.0: https://creativecommons.org/licenses/by/4.0/. Provider-required attribution applies to reuse. The current public access portal is https://data.icos-cp.eu/portal/; this access route is not a record of the original download event.

One archive identifier needs special care: the public record associated with `ICOS_DE-Kli_FLUXNET_2004-2025_v1.3_r1.zip` identifies DE-RuW rather than the local archive's DE-Kli station. The mismatch is retained explicitly in `tower_dataset_citations.csv` rather than attributing those observations to the wrong site. Empty author fields reproduce the available provider metadata and are not filled with guessed names.

## Scope of tower results

Monthly GPP and ET provide the primary tower consistency comparisons. Previously completed phenological-phase and WUE diagnostics are retained in `data/exploratory_tower/` as exploratory results. This separation changes their presentation, not the underlying observations or estimates. Site measurements and native raster inputs are accessed through their providers rather than redistributed in this package.
