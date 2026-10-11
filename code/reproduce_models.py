"""Refit frozen spatial-block sufficient statistics and verify published estimates.

Run from any directory: python code/reproduce_models.py --output output
"""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.stats import t

ROOT = Path(__file__).resolve().parents[1]
PHASES = ("Early", "Middle", "Late")
RESPONSES = {"GPP": (1, 0, 0), "T": (0, 1, 0), "ET": (0, 0, 1),
             "WUE": (1, 0, -1), "GPP_T": (1, -1, 0), "T_ET": (0, 1, -1)}


def bh(p):
    """Benjamini-Hochberg adjusted P values; missing tests remain in family size."""
    p = np.asarray(p, dtype=float)
    finite = np.isfinite(p)
    order = np.argsort(np.where(finite, p, 1.0), kind="stable")
    q = np.minimum.accumulate((np.where(finite, p, 1.0)[order] * len(p) /
                               np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.minimum(q, 1.0)
    out[~finite] = np.nan
    return out


def fit(path):
    """Joint three-response within-cell/year FE fit with cluster covariance."""
    with np.load(path, allow_pickle=False) as z:
        counts = z["counts"].copy()
        keep = counts[:, 0] > 0
        blocks, counts, ids = z["s"][keep], counts[keep], z["ids"][keep]
    total = blocks.sum(axis=0)
    k = total.shape[0] - 3
    xx, xy = total[:k, :k], total[:k, k:]
    active = np.diag(xx) > 1e-9
    year_columns = np.flatnonzero(active & (np.arange(k) >= 11))
    year_block = xx[np.ix_(year_columns, year_columns)]
    if len(year_columns) > 1 and abs(year_block.sum()) < 1e-8 * np.trace(year_block):
        active[year_columns[0]] = False
    scale = np.sqrt(np.diag(xx)[active])
    scaled = xx[np.ix_(active, active)] / scale[:, None] / scale[None, :]
    condition = np.linalg.cond(scaled)
    if condition >= 1e10:
        raise ValueError(f"Ill-conditioned active design: {condition}")
    inv = np.zeros_like(xx)
    inv[np.ix_(active, active)] = np.linalg.inv(scaled) / scale[:, None] / scale[None, :]
    beta = inv @ xy
    score = blocks[:, :k, k:] - blocks[:, :k, :k] @ beta
    n, cells = (int(counts[:, j].sum()) for j in (0, 1))
    g = len(ids)
    correction = g / (g - 1) * (n - 1) / (n - cells - int(active.sum()))
    influence = np.einsum("ab,gbc->gac", inv, score) * np.sqrt(correction)
    return dict(beta=beta, influence=influence, ids=ids, df=g-1,
                n_clusters=g, n_cells=cells, n_cell_years=n, condition=condition)


def estimate(fits, terms, response, column=7):
    """Terms are (fit key, weight); block IDs are aligned before summation."""
    ids = np.unique(np.concatenate([fits[key]["ids"] for key, _ in terms]))
    score = np.zeros(len(ids), dtype=float)
    beta = 0.0
    co = np.asarray(RESPONSES[response])
    for key, weight in terms:
        f = fits[key]
        beta += weight * float(f["beta"][column] @ co)
        score[np.searchsorted(ids, f["ids"])] += weight * (f["influence"][:, column, :] @ co)
    se = float(np.linalg.norm(score))
    df = min(fits[key]["df"] for key, _ in terms)
    critical = t.ppf(.975, df)
    return dict(beta=beta, se=se, ci_low=beta-critical*se, ci_high=beta+critical*se,
                p=float(2*t.sf(abs(beta/se), df)) if se else np.nan,
                df=df, n_clusters=len(ids))


def primary():
    fits = {p: fit(ROOT / "data/model_inputs" / f"primary_{p.lower()}.npz") for p in PHASES}
    rows = []
    comparisons = {p: [(p, 1)] for p in PHASES}
    comparisons.update({"Early minus Middle": [("Early", 1), ("Middle", -1)],
                        "Late minus Middle": [("Late", 1), ("Middle", -1)]})
    for phase, terms in comparisons.items():
        for topic, column in (("D1_Compound_dryness", 10), ("D2_Soil_air_thermal", 7)):
            for endpoint in RESPONSES:
                rows.append(dict(setting="Primary", phase=phase, topic=topic, endpoint=endpoint,
                                 family_id="Primary_60", **estimate(fits, terms, endpoint, column)))
    result = pd.DataFrame(rows)
    result["q"] = bh(result.p)
    return result


def allocation():
    models = ("Overlap_full", "Midpoint_full", "Overlap_common", "Midpoint_common")
    fits = {(m, p): fit(ROOT / "data/model_inputs" / f"allocation_{m}_{p}.npz")
            for m in models for p in PHASES}
    phases = {p: {p: 1} for p in PHASES}
    phases.update({"Early minus Middle": {"Early": 1, "Middle": -1},
                   "Late minus Middle": {"Late": 1, "Middle": -1},
                   "Early minus Late": {"Early": 1, "Late": -1}})
    definitions = {m: {m: 1} for m in models}
    definitions.update({"Allocation_common": {"Midpoint_common": 1, "Overlap_common": -1},
                        "Support_overlap": {"Overlap_common": 1, "Overlap_full": -1},
                        "Support_midpoint": {"Midpoint_common": 1, "Midpoint_full": -1}})
    families = {"Overlap_full": "FULL_REPRODUCTION_NO_NEW_Q", "Midpoint_full": "FULL_REPRODUCTION_NO_NEW_Q",
                "Overlap_common": "G2_common_support_72", "Midpoint_common": "G2_common_support_72",
                "Allocation_common": "G2_allocation_differences_36",
                "Support_overlap": "G2_support_differences_72", "Support_midpoint": "G2_support_differences_72"}
    rows = []
    for model, model_weights in definitions.items():
        for contrast, phase_weights in phases.items():
            terms = [((m, p), u*v) for m, u in model_weights.items() for p, v in phase_weights.items()]
            for endpoint in RESPONSES:
                first = fits[terms[0][0]]
                rows.append(dict(model_id=model, contrast=contrast, endpoint=endpoint,
                                 family_id=families[model], n_cells=first["n_cells"],
                                 n_cell_years=first["n_cell_years"], **estimate(fits, terms, endpoint)))
    result = pd.DataFrame(rows)
    result["q"] = np.nan
    for family, indices in result.groupby("family_id").groups.items():
        if family != "FULL_REPRODUCTION_NO_NEW_Q":
            result.loc[indices, "q"] = bh(result.loc[indices, "p"])
    return result


def verify(actual, reference, keys, name):
    joined = actual.merge(reference, on=keys, suffixes=("_new", "_ref"), validate="one_to_one")
    if len(joined) != len(reference) or len(actual) != len(reference):
        raise AssertionError(f"Incomplete {name} result membership")
    rows = []
    # Absolute tolerances account for summation and CSV precision, not inference changes.
    # Inverse fitting uses a scaled matrix; all designs enforce condition < 1e10.
    for field in ("beta", "se", "ci_low", "ci_high", "p", "q"):
        a, b = joined[field+"_new"].to_numpy(), joined[field+"_ref"].to_numpy()
        finite = np.isfinite(a) & np.isfinite(b)
        missing_matches = np.array_equal(np.isnan(a), np.isnan(b))
        error = float(np.max(abs(a[finite]-b[finite]))) if finite.any() else 0.0
        tol = 2e-10 if field in ("p", "q") else 2e-11
        passed = missing_matches and np.allclose(a, b, rtol=2e-9, atol=tol, equal_nan=True)
        rows.append(dict(test=name, quantity=field, rows=len(joined), max_abs_error=error,
                         absolute_tolerance=tol, relative_tolerance=2e-9, passed=bool(passed)))
    for item in rows:
        if not item["passed"]:
            raise AssertionError(item)
    return rows


def checksum_test():
    manifest = pd.read_csv(ROOT / "file_manifest.csv")
    for row in manifest.itertuples():
        content = (ROOT / row.path).read_bytes()
        if len(content) != row.bytes or hashlib.sha256(content).hexdigest() != row.sha256:
            raise AssertionError(f"Checksum mismatch: {row.path}")
    return dict(test="package_checksums", files=len(manifest), passed=True)


def run(output):
    output.mkdir(parents=True, exist_ok=True)
    checks = [checksum_test()]
    p, a = primary(), allocation()
    checks += verify(p, pd.read_csv(ROOT / "data/reference/primary_family.csv"),
                     ["phase", "topic", "endpoint"], "primary_refit")
    checks += verify(a, pd.read_csv(ROOT / "data/reference/allocation_results.csv"),
                     ["model_id", "contrast", "endpoint"], "allocation_refit")
    p.to_csv(output / "primary_family.csv", index=False)
    a.to_csv(output / "allocation_results.csv", index=False)
    early = p.query("topic == 'D2_Soil_air_thermal' and phase == 'Early'").set_index("endpoint")
    identity_error = abs(early.loc["WUE", "beta"] - early.loc["GPP_T", "beta"] - early.loc["T_ET", "beta"])
    checks.append(dict(test="log_component_identity", max_abs_error=float(identity_error),
                       passed=bool(identity_error < 1e-12)))
    effect = early.loc["WUE"]
    summary = dict(primary_estimates=len(p), allocation_estimates=len(a),
                   early_wue_percent=float(100*np.expm1(effect.beta)),
                   early_wue_ci_percent=[float(100*np.expm1(effect.ci_low)), float(100*np.expm1(effect.ci_high))],
                   early_wue_adjusted_P=float(effect.q), checks=checks)
    (output / "verification.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    if not all(c["passed"] for c in checks):
        raise AssertionError("One or more verification checks failed")
    print(json.dumps({k: v for k, v in summary.items() if k != "checks"}, indent=2))
    print(f"All {len(checks)} verification checks passed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "output")
    run(parser.parse_args().output.resolve())
