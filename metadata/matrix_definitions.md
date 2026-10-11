# Model matrices and inference

Each NumPy archive contains `s`, `counts` and `ids`. Optional `b20` entries in the original primary files describe coarser spatial aggregation and are not used by this reproduction script.

`s[g]` is the weighted cross-product matrix for spatial block `ids[g]`. Its variable order is eleven scientific predictors, 23 year-indicator columns (2002–2024), and three responses. The three responses are natural logarithms of GPP, T and ET. Row and column dimensions are 37 × 37. The predictor order is:

1. VPD / 0.1 kPa;
2. negative soil water content / 0.01 m³ m⁻³;
3. air temperature (°C);
4. log shortwave radiation;
5. log(1 + precipitation);
6. wind index / 0.5;
7. net radiation / 10;
8. soil temperature minus air temperature (°C);
9. start of season / 10 days;
10. growing-season length / 10 days;
11. interaction of within-cell-centred VPD and soil-dryness predictors.

All design and response columns, including year indicators, have undergone within-cell centring. The interaction is constructed from the individually centred predictors before this design transformation. `counts[g,0]` is the number of retained cell-years; `counts[g,1]` is the number of distinct retained cells. Empty blocks are removed before fitting.

The joint coefficient matrix is `(X'WX)⁻¹X'WY`. A diagonal scaling is applied before matrix inversion. Cluster scores are `X_g'W_gY_g − X_g'W_gX_g beta`. The finite-sample multiplier for cluster covariance is `G/(G−1) × (N−1)/(N−Ncell−K)`, where K counts active design columns. Degrees of freedom are G−1; contrasts involving more than one fit use the smaller contributing degrees of freedom. Cross-response covariance is retained by combining the three response-specific cluster scores before calculating each standard error. Shared block IDs are aligned for cross-phase and cross-model contrasts.

`Primary` and `Overlap_full` represent the full day-overlap sample; `Midpoint_full` uses each method's own valid support. The two `common` configurations use identical valid cell-years. `Allocation_common` is midpoint minus overlap on common support. `Support_overlap` and `Support_midpoint` isolate the corresponding support changes. Each phase is one equal-duration third of the GOSIF-derived season.

The primary family retains 60 tests: two exposures × six responses × three phases, plus two phase contrasts for each exposure and response. Allocation estimates use the 72-, 36- and 72-member families recorded in the reference table. Reproduced method-specific full-support estimates do not receive a new, smaller correction family.

Numerical comparisons allow absolute error of 2×10⁻¹¹ for beta, standard error and interval endpoints, or 2×10⁻¹⁰ for P and adjusted P, plus relative error 2×10⁻⁹. These tolerances accommodate floating-point and CSV precision; they are fixed before the package is run. The scaled design condition number must be below 10¹⁰.
