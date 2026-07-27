#!/usr/bin/env python3
"""
Bounded asymmetry analysis of curated DisProt disorder fractions and
AlphaFold very-low-confidence fractions.

All quantities are computed from the four DisProt JSON exports via the
protein-level table written by 01_build_tables.py. No network retrieval is
performed (all required databases were blocked at the egress proxy; see
results/diagnostics/alphafold_retrieval_manifest.csv).

Portability note (2026-07-27 refactor): absolute paths replaced with
project-relative pathlib paths plus DISPROT_OUT override, and outputs routed
into results/tables and results/diagnostics. The random seed, replicate count,
grid, estimators, statistics and loop order are unchanged, so every reported
value is bit-identical to the archived deposit.
"""
import itertools
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# ------------------------------------------------------------------- paths
ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = Path(os.environ.get("DISPROT_OUT", ROOT / "results"))
TABLES_DIR = RESULTS_DIR / "tables"
DIAG_DIR = RESULTS_DIR / "diagnostics"
for _d in (TABLES_DIR, DIAG_DIR):
    _d.mkdir(parents=True, exist_ok=True)

rng = np.random.default_rng(20260726)

GRID = np.linspace(0.0, 1.0, 512)
DX = GRID[1] - GRID[0]
EPS_DEFAULT = 1e-3
NBOOT = 2000
GROUPS = ["Human", "Fungal", "Bacterial", "Nematode"]

_prot_path = TABLES_DIR / "disprot_protein_level.csv"
if not _prot_path.exists():
    raise FileNotFoundError(
        f"{_prot_path} not found. Run 01_build_tables.py first, or use "
        f"run_all.py to execute the pipeline in order."
    )
prot = pd.read_csv(_prot_path)


# ---------------------------------------------------------------- densities
def kde_density(x, bw_scale=1.0, grid=GRID):
    """Gaussian KDE evaluated on [0,1], renormalised to unit integral.

    Renormalisation on the closed domain is used instead of reflection so
    that the same estimator applies to every group; it is a boundary
    approximation, and the histogram estimator is reported alongside.
    """
    x = np.asarray(x, float)
    if len(x) < 3 or np.std(x) == 0:
        return None
    k = stats.gaussian_kde(x, bw_method="scott")
    k.set_bandwidth(k.factor * bw_scale)
    d = k(grid)
    return d / np.trapezoid(d, grid)


def hist_density(x, nbins=20, grid=GRID):
    x = np.asarray(x, float)
    h, edges = np.histogram(x, bins=nbins, range=(0, 1), density=True)
    idx = np.clip(np.digitize(grid, edges) - 1, 0, nbins - 1)
    d = h[idx]
    s = np.trapezoid(d, grid)
    return d / s if s > 0 else d


# ------------------------------------------------------------------ indices
def asym_profile(pA, pB, eps=EPS_DEFAULT):
    return (pA - pB) / (pA + pB + 2 * eps)


def indices(S, grid=GRID):
    """Domain-normalised A-bar, B-bar and directional consistency C."""
    span = grid[-1] - grid[0]
    A = np.trapezoid(np.abs(S), grid) / span
    B = np.trapezoid(S, grid) / span
    C = abs(B) / A if A > 0 else 0.0
    return A, B, C


def analyse_pair(xA, xB, eps=EPS_DEFAULT, bw=1.0, estimator="kde",
                 nbins=20):
    if estimator == "kde":
        pA, pB = kde_density(xA, bw), kde_density(xB, bw)
    else:
        pA, pB = hist_density(xA, nbins), hist_density(xB, nbins)
    if pA is None or pB is None:
        return None
    S = asym_profile(pA, pB, eps)
    A, B, C = indices(S)
    return {"S": S, "pA": pA, "pB": pB, "A": A, "B": B, "C": C}


# ------------------------------------------------------- feature definitions
FEATURES = {
    "D_exp": "Curated DisProt disorder fraction",
    "af_very_low_content": "AlphaFold very-low-confidence fraction",
}

rows_summary, rows_pairs, profile_store, boot_store = [], [], {}, []

for feat in FEATURES:
    sub = prot.dropna(subset=[feat])
    for g in GROUPS:
        v = sub.loc[sub.group == g, feat].to_numpy(float)
        rows_summary.append({
            "feature": feat, "group": g, "n": len(v),
            "mean": v.mean(), "median": np.median(v),
            "q25": np.percentile(v, 25), "q75": np.percentile(v, 75),
            "iqr": np.percentile(v, 75) - np.percentile(v, 25),
            "min": v.min(), "max": v.max(),
            "frac_eq_0": float((v == 0).mean()),
            "frac_eq_1": float((v == 1).mean()),
        })

    for gA, gB in itertools.combinations(GROUPS, 2):
        xA = sub.loc[sub.group == gA, feat].to_numpy(float)
        xB = sub.loc[sub.group == gB, feat].to_numpy(float)
        res = analyse_pair(xA, xB)
        if res is None:
            continue
        profile_store[(feat, gA, gB)] = res

        ks = stats.ks_2samp(xA, xB)
        mw = stats.mannwhitneyu(xA, xB, alternative="two-sided")
        # rank-biserial effect size from U
        rb = 2 * mw.statistic / (len(xA) * len(xB)) - 1

        # protein-level bootstrap of the three indices
        bA, bB, bC = [], [], []
        for _ in range(NBOOT):
            ra = rng.choice(xA, len(xA), replace=True)
            rb_ = rng.choice(xB, len(xB), replace=True)
            r = analyse_pair(ra, rb_)
            if r is None:
                continue
            bA.append(r["A"]); bB.append(r["B"]); bC.append(r["C"])
        bA, bB, bC = map(np.asarray, (bA, bB, bC))

        rows_pairs.append({
            "feature": feat, "group_A": gA, "group_B": gB,
            "n_A": len(xA), "n_B": len(xB),
            "A_bar": res["A"],
            "A_bar_lo": np.percentile(bA, 2.5),
            "A_bar_hi": np.percentile(bA, 97.5),
            "B_bar": res["B"],
            "B_bar_lo": np.percentile(bB, 2.5),
            "B_bar_hi": np.percentile(bB, 97.5),
            "B_bar_excludes_zero": bool(
                np.percentile(bB, 2.5) > 0 or np.percentile(bB, 97.5) < 0),
            "C": res["C"],
            "C_lo": np.percentile(bC, 2.5),
            "C_hi": np.percentile(bC, 97.5),
            "triangle_ok": bool(abs(res["B"]) <= res["A"] + 1e-12),
            "S_max_abs": float(np.max(np.abs(res["S"]))),
            "S_in_open_unit": bool(np.all(np.abs(res["S"]) < 1)),
            "n_zero_crossings": int(
                np.sum(np.diff(np.sign(res["S"])) != 0)),
            "ks_stat": ks.statistic, "ks_p": ks.pvalue,
            "mw_U": mw.statistic, "mw_p": mw.pvalue,
            "rank_biserial": rb,
            "median_A": float(np.median(xA)),
            "median_B": float(np.median(xB)),
        })

pd.DataFrame(rows_summary).to_csv(
    TABLES_DIR / "group_feature_summary.csv", index=False)
pairs = pd.DataFrame(rows_pairs)
pairs.to_csv(TABLES_DIR / "thales_curated_sample_results.csv", index=False)

# ------------------------------------------------------------- robustness
rob = []
for feat in FEATURES:
    sub = prot.dropna(subset=[feat])
    for gA, gB in itertools.combinations(GROUPS, 2):
        xA = sub.loc[sub.group == gA, feat].to_numpy(float)
        xB = sub.loc[sub.group == gB, feat].to_numpy(float)
        specs = ([("kde", bw, 20, EPS_DEFAULT) for bw in (0.75, 1.0, 1.25)]
                 + [("hist", 1.0, nb, EPS_DEFAULT) for nb in (10, 20, 40)]
                 + [("kde", 1.0, 20, e) for e in (1e-6, 1e-2, 1e-1)])
        for est, bw, nb, eps in specs:
            r = analyse_pair(xA, xB, eps=eps, bw=bw, estimator=est, nbins=nb)
            if r is None:
                continue
            rob.append({"feature": feat, "group_A": gA, "group_B": gB,
                        "estimator": est, "bw_scale": bw, "nbins": nb,
                        "eps": eps, "A_bar": r["A"], "B_bar": r["B"],
                        "C": r["C"],
                        "sign_B": int(np.sign(r["B"]))})
rob = pd.DataFrame(rob)
rob.to_csv(DIAG_DIR / "robustness_sensitivity.csv", index=False)

# sign stability of B_bar per comparison
sign_tab = (rob.groupby(["feature", "group_A", "group_B"])["sign_B"]
            .agg(n_specs="size", n_pos=lambda s: int((s > 0).sum()),
                 n_neg=lambda s: int((s < 0).sum()))
            .reset_index())
sign_tab["sign_stable"] = (
    (sign_tab.n_pos == sign_tab.n_specs) | (sign_tab.n_neg == sign_tab.n_specs))
sign_tab.to_csv(DIAG_DIR / "sign_stability.csv", index=False)

# ------------------------- experimental vs AlphaFold concordance per group
conc = []
paired = prot.dropna(subset=["D_exp", "af_very_low_content"])
for g in GROUPS:
    s = paired[paired.group == g]
    x = s["D_exp"].to_numpy(float)
    y = s["af_very_low_content"].to_numpy(float)
    rho, p = stats.spearmanr(x, y)
    bs = []
    for _ in range(NBOOT):
        i = rng.integers(0, len(x), len(x))
        r_, _ = stats.spearmanr(x[i], y[i])
        if np.isfinite(r_):
            bs.append(r_)
    bs = np.asarray(bs)
    conc.append({
        "group": g, "n_paired": len(s), "spearman_rho": rho, "p": p,
        "rho_lo": np.percentile(bs, 2.5), "rho_hi": np.percentile(bs, 97.5),
        "median_D_exp": float(np.median(x)),
        "median_af": float(np.median(y)),
        "median_signed_gap": float(np.median(y - x)),
        "mean_signed_gap": float(np.mean(y - x)),
        "wilcoxon_p": stats.wilcoxon(y, x).pvalue,
    })
conc = pd.DataFrame(conc)
conc.to_csv(TABLES_DIR / "cross_taxon_calibration_summary.csv", index=False)

# pairwise test of whether the AF-minus-experimental gap differs by taxon
gap_rows = []
for gA, gB in itertools.combinations(GROUPS, 2):
    a = (paired[paired.group == gA]["af_very_low_content"]
         - paired[paired.group == gA]["D_exp"]).to_numpy(float)
    b = (paired[paired.group == gB]["af_very_low_content"]
         - paired[paired.group == gB]["D_exp"]).to_numpy(float)
    mw = stats.mannwhitneyu(a, b, alternative="two-sided")
    gap_rows.append({"group_A": gA, "group_B": gB,
                     "median_gap_A": float(np.median(a)),
                     "median_gap_B": float(np.median(b)),
                     "mw_p": mw.pvalue,
                     "rank_biserial": 2 * mw.statistic / (len(a) * len(b)) - 1})
pd.DataFrame(gap_rows).to_csv(DIAG_DIR / "gap_by_taxon_tests.csv", index=False)

np.savez_compressed(DIAG_DIR / "asymmetry_profiles.npz",
                    grid=GRID,
                    **{f"{f}|{a}|{b}": v["S"]
                       for (f, a, b), v in profile_store.items()})

pd.set_option("display.width", 200)
print("== group summaries ==")
print(pd.DataFrame(rows_summary).round(4).to_string(index=False))
print("\n== pairwise indices ==")
print(pairs[["feature", "group_A", "group_B", "A_bar", "B_bar", "B_bar_lo",
             "B_bar_hi", "C", "n_zero_crossings", "ks_p", "mw_p",
             "rank_biserial"]].round(4).to_string(index=False))
print("\n== sign stability ==")
print(sign_tab.to_string(index=False))
print("\n== AF vs experimental concordance ==")
print(conc.round(4).to_string(index=False))
print("\n== gap differences by taxon ==")
print(pd.DataFrame(gap_rows).round(4).to_string(index=False))
