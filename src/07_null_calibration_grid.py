#!/usr/bin/env python3
"""
Synthetic null-calibration audit of A_bar vs. A_w.

STATISTICAL-METHOD DIAGNOSTIC ONLY. No biological data, no DisProt input, no
network access. This script exists to test one narrow question: does the
density-weighted asymmetry magnitude A_w (defined in
revision_addendum/ADDENDUM_permutation_calibration.pdf, section 4) genuinely
reduce the finite-sample null pathology of the unweighted magnitude A_bar
across heterogeneous synthetic sampling regimes, or was the apparent
improvement specific to the twelve real DisProt comparisons previously
examined?

REPRODUCIBILITY GAP THIS SCRIPT ADDRESSES: the addendum describes a script
named `profile_inference.py` implementing permutation tests, A_w, T_max,
Holm correction, simultaneous bands and mass trimming. No such file, and no
A_w implementation of any kind, exists anywhere in this repository's
committed history (verified by exhaustive grep across every .py/.md/.csv
file prior to writing this module). This script implements the smallest
transparent reference sufficient for the null-calibration question only: it
does NOT implement T_max, simultaneous bands, or mass trimming, which are not
needed to answer this question and are explicitly out of scope.

DEFINITIONS (reconstructed from src/02_asymmetry.py and cross-checked against
the addendum; NOT silently changed from the production pipeline):

    S(x)   = (p_A(x) - p_B(x)) / (p_A(x) + p_B(x) + 2*eps)
    A_bar  = (1/|Omega|) * integral |S(x)| dx
    B_bar  = (1/|Omega|) * integral S(x) dx
    A_w    = integral |S(x)| * p_bar(x) dx / integral p_bar(x) dx,
             p_bar = (p_A + p_B) / 2

Domain: [0, 1] (matches production). Epsilon: EPS_DEFAULT = 1e-3 (matches
src/02_asymmetry.py exactly; the addendum's eps = delta/|Omega|
reparameterization is a no-op on this unit domain per the addendum's own
text, so no numerical difference arises). Bandwidth: 1-D Scott's rule,
h = sigma_hat(ddof=1) * n**(-1/5), which is exactly what
scipy.stats.gaussian_kde(x, bw_method="scott") computes in 1-D (verified by
a cross-check test against scipy's own implementation in
tests/test_null_calibration_grid.py) -- reimplemented directly here (rather
than imported) for two reasons: src/02_asymmetry.py's filename is not a
valid Python module identifier, and this reimplementation is what makes the
matched-bandwidth diagnostic control (Condition 2) possible at all, since
scipy.stats.gaussian_kde does not expose a "shared absolute bandwidth across
two independently-sized samples" option directly.

Grid: PRIMARY experiments use a 257-point grid on [0,1] (production uses
512); this is a deliberate, documented reduction for computational
tractability across the ~56,000+ replicate/permutation evaluations this
audit requires, not a silent change. Insensitivity to this choice is
verified directly in the epsilon/grid-resolution robustness check
(run_epsilon_grid_sensitivity), which includes grid=512 (the production
value) as one of its tested resolutions.

Renormalization: closed-domain (divide the raw kernel-density evaluation by
its own trapezoidal integral over [0,1]), matching
src/02_asymmetry.py's kde_density docstring verbatim ("Renormalisation on
the closed domain is used instead of reflection").

Labels used throughout this module's outputs follow the project convention:
EXACT ALGEBRA, EXACT COMPUTATION, MONTE CARLO RESULT, EMPIRICAL SIMULATION
RESULT, HEURISTIC, OPEN, NEGATIVE RESULT. Simulation findings are never
called theorems; no claim of general statistical optimality, and no claim
that A_w or the total-variation identity is novel.
"""
import itertools
import json
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# ------------------------------------------------------------------- trapz
# np.trapezoid does not exist before NumPy 2.0; this environment's installed
# NumPy is 1.26.4 (see PRE-COMMIT notes), so src/02_asymmetry.py's use of
# np.trapezoid would itself raise AttributeError here. Using np.trapz
# (available in both 1.x and 2.x) makes this diagnostic run regardless.
_trapz = getattr(np, "trapezoid", None) or np.trapz

# ------------------------------------------------------------------- paths
ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = Path(os.environ.get("DISPROT_NULLCAL_OUT", ROOT / "results" / "null_calibration"))

MASTER_SEED = 20260906  # documented once; every condition's RNG is spawned
                          # deterministically from this single seed sequence.

DOMAIN = (0.0, 1.0)
EPS_DEFAULT = 1e-3  # matches src/02_asymmetry.py EPS_DEFAULT exactly

# ================================================================= Section 1
# Core statistics: pure functions, matching src/02_asymmetry.py exactly.
# ================================================================= Section 1

def scott_bandwidth(n, sigma):
    """1-D Scott's rule: h = sigma * n**(-1/5).
    EXACT ALGEBRA -- this is exactly what scipy.stats.gaussian_kde computes
    internally for bw_method='scott' in one dimension (covariance_factor =
    n**(-1/(d+4)) with d=1); cross-checked numerically in the test suite."""
    return sigma * n ** (-1.0 / 5.0)


def gaussian_kde_eval(x, grid, h):
    """Vectorised Gaussian KDE evaluation at bandwidth h (a kernel std, not
    a relative factor). EXACT COMPUTATION of the standard KDE formula."""
    x = np.asarray(x, dtype=np.float64)
    diffs = (grid[None, :] - x[:, None]) / h
    dens = np.exp(-0.5 * diffs * diffs).sum(axis=0)
    dens /= x.size * h * np.sqrt(2.0 * np.pi)
    return dens


def kde_density(x, grid, bw_scale=1.0, h_override=None):
    """Closed-domain-renormalised Gaussian KDE, matching
    src/02_asymmetry.py's kde_density. If h_override is given, that absolute
    bandwidth is used (Condition 2: matched/shared bandwidth); otherwise the
    sample's own Scott bandwidth is used (Condition 1: real pipeline)."""
    x = np.asarray(x, dtype=np.float64)
    n = x.size
    sigma = np.std(x, ddof=1) if n > 1 else 0.0
    if sigma <= 0:
        sigma = 1e-6
    h = h_override if h_override is not None else scott_bandwidth(n, sigma) * bw_scale
    d = gaussian_kde_eval(x, grid, h)
    area = _trapz(d, grid)
    return d / area if area > 0 else d


def asym_profile(pA, pB, eps=EPS_DEFAULT):
    """EXACT ALGEBRA -- identical to src/02_asymmetry.py."""
    return (pA - pB) / (pA + pB + 2.0 * eps)


def indices(S, grid):
    """A_bar, B_bar -- identical construction to src/02_asymmetry.py's
    `indices` (domain span here is always 1.0 since grid spans [0,1])."""
    span = grid[-1] - grid[0]
    A = _trapz(np.abs(S), grid) / span
    B = _trapz(S, grid) / span
    return A, B


def a_w(pA, pB, S, grid):
    """Density-weighted magnitude, exactly as defined in
    revision_addendum/ADDENDUM_permutation_calibration.pdf section 4:
        A_w = integral |S| * p_bar dx / integral p_bar dx,  p_bar=(pA+pB)/2
    EXACT ALGEBRA given pA, pB, S."""
    pbar = 0.5 * (pA + pB)
    num = _trapz(np.abs(S) * pbar, grid)
    den = _trapz(pbar, grid)
    return num / den if den > 0 else 0.0


def compute_all_statistics(xA, xB, grid, eps=EPS_DEFAULT, shared_bandwidth=False):
    """One pair -> (A_bar, A_w, B_bar, ks_stat, ks_p). `shared_bandwidth`
    selects Condition 2 (matched pooled-Scott bandwidth) vs. Condition 1
    (each sample's own independent Scott bandwidth, the real-pipeline
    default)."""
    if shared_bandwidth:
        pooled = np.concatenate([xA, xB])
        sigma_pooled = np.std(pooled, ddof=1)
        h_shared = scott_bandwidth(pooled.size, sigma_pooled if sigma_pooled > 0 else 1e-6)
        pA = kde_density(xA, grid, h_override=h_shared)
        pB = kde_density(xB, grid, h_override=h_shared)
    else:
        pA = kde_density(xA, grid)
        pB = kde_density(xB, grid)
    S = asym_profile(pA, pB, eps)
    A, B = indices(S, grid)
    Aw = a_w(pA, pB, S, grid)
    ks = stats.ks_2samp(xA, xB)
    return A, Aw, B, ks.statistic, ks.pvalue


# ================================================================= Section 4
# Synthetic null distribution families. FIXED, PRE-DECLARED, on [0,1].
# ================================================================= Section 4

FAMILIES = {
    "A_interior_unimodal": {"kind": "beta", "params": (5.0, 5.0)},
    "B_broad_interior":    {"kind": "beta", "params": (2.0, 2.0)},
    "C_boundary_heavy":    {"kind": "beta", "params": (0.7, 3.0)},
    "D_bimodal":           {"kind": "mixture",
                             "components": [(0.5, 8.0, 2.0), (0.5, 2.0, 8.0)]},
}


def sample_family(rng, family_key, n):
    spec = FAMILIES[family_key]
    if spec["kind"] == "beta":
        a, b = spec["params"]
        return rng.beta(a, b, size=n)
    elif spec["kind"] == "mixture":
        comps = spec["components"]
        weights = np.array([c[0] for c in comps])
        counts = rng.multinomial(n, weights)
        parts = [rng.beta(a, b, size=c) for (w, a, b), c in zip(comps, counts)]
        out = np.concatenate(parts)
        rng.shuffle(out)
        return out
    raise ValueError(family_key)


# ================================================================= Section 5
# Sample-size grid. FIXED, PRE-DECLARED.
# ================================================================= Section 5

SIZE_CONFIGS = [
    (50, 50), (200, 200), (800, 800),          # equal-size controls
    (50, 100), (50, 250), (50, 500),           # unequal-size controls
    (150, 1350),                                # severe imbalance (repo-scale)
]

BANDWIDTH_CONDITIONS = ["independent", "shared_pooled"]


# ============================================================ RNG machinery
def make_condition_rngs():
    """Deterministically spawn one independent RNG per
    (family, size_config, bandwidth_condition) from a single master seed, in
    a fixed enumeration order, so re-running this module reproduces every
    condition's synthetic draws exactly regardless of run order or which
    subset of conditions is executed."""
    ss = np.random.SeedSequence(MASTER_SEED)
    keys = [
        (fam, cfg, bwcond)
        for fam in FAMILIES
        for cfg in SIZE_CONFIGS
        for bwcond in BANDWIDTH_CONDITIONS
    ]
    children = ss.spawn(len(keys))
    return {key: np.random.default_rng(child) for key, child in zip(keys, children)}


# ================================================================= Section 7
# Run one condition: N_outer null replicates -> distribution of statistics.
# ================================================================= Section 7

def run_condition(rng, family_key, n_A, n_B, grid, n_outer, shared_bandwidth, eps=EPS_DEFAULT):
    A_arr = np.empty(n_outer)
    Aw_arr = np.empty(n_outer)
    B_arr = np.empty(n_outer)
    ks_arr = np.empty(n_outer)
    ksp_arr = np.empty(n_outer)
    for i in range(n_outer):
        xA = sample_family(rng, family_key, n_A)
        xB = sample_family(rng, family_key, n_B)  # SAME family -> F_A = F_B
        A, Aw, B, ksd, ksp = compute_all_statistics(
            xA, xB, grid, eps=eps, shared_bandwidth=shared_bandwidth)
        A_arr[i], Aw_arr[i], B_arr[i] = A, Aw, B
        ks_arr[i], ksp_arr[i] = ksd, ksp
    return {"A_bar": A_arr, "A_w": Aw_arr, "B_bar": B_arr,
            "ks_stat": ks_arr, "ks_p": ksp_arr}


def summarize(arr):
    return {
        "median": float(np.median(arr)), "mean": float(np.mean(arr)),
        "std": float(np.std(arr, ddof=1)),
        "q90": float(np.percentile(arr, 90)),
        "q95": float(np.percentile(arr, 95)),
        "q99": float(np.percentile(arr, 99)),
    }


# ================================================================ Section 10
# Full-pipeline permutation calibration (bandwidth reselected every replicate)
# ================================================================ Section 10

def permutation_pvalue(xA, xB, grid, n_perm, rng, statistic="A_bar",
                        shared_bandwidth=False, eps=EPS_DEFAULT):
    """One permutation test. Group sizes preserved; the ENTIRE estimation
    pipeline (bandwidth selection included) is recomputed inside every
    replicate -- never frozen from the observed data -- matching the
    addendum's stated design requirement."""
    nA, nB = len(xA), len(xB)
    pooled = np.concatenate([xA, xB])
    A_obs, Aw_obs, B_obs, ks_obs, _ = compute_all_statistics(
        xA, xB, grid, eps=eps, shared_bandwidth=shared_bandwidth)
    obs = {"A_bar": abs(A_obs), "A_w": abs(Aw_obs), "ks_stat": ks_obs}[statistic]
    count_ge = 0
    for _ in range(n_perm):
        perm = rng.permutation(pooled)
        pA, pB = perm[:nA], perm[nA:]
        A_p, Aw_p, B_p, ks_p, _ = compute_all_statistics(
            pA, pB, grid, eps=eps, shared_bandwidth=shared_bandwidth)
        val = {"A_bar": abs(A_p), "A_w": abs(Aw_p), "ks_stat": ks_p}[statistic]
        if val >= obs:
            count_ge += 1
    return (count_ge + 1) / (n_perm + 1)


def run_type1_error_check(rng, family_key, n_A, n_B, grid, n_outer_calib, n_perm,
                           statistic, shared_bandwidth=False, alpha=0.05,
                           eps=EPS_DEFAULT):
    """Empirical type-I error of the permutation-calibrated statistic under
    a TRUE null (F_A = F_B by construction). MONTE CARLO RESULT."""
    rejections = 0
    for _ in range(n_outer_calib):
        xA = sample_family(rng, family_key, n_A)
        xB = sample_family(rng, family_key, n_B)
        p = permutation_pvalue(xA, xB, grid, n_perm, rng, statistic=statistic,
                                shared_bandwidth=shared_bandwidth, eps=eps)
        if p < alpha:
            rejections += 1
    rate = rejections / n_outer_calib
    se = np.sqrt(alpha * (1 - alpha) / n_outer_calib)
    return rate, se


# ================================================================= Section 9
# Planted alternatives (predeclared, NOT tuned for A_w).
# ================================================================= Section 9

ALTERNATIVES = {
    "alt1_location_shift": {"A": ("beta", (5.0, 5.0)), "B": ("beta", (6.0, 4.0))},
    "alt2_variance_change": {"A": ("beta", (5.0, 5.0)), "B": ("beta", (2.0, 2.0))},
    "alt3_localized_mixture": {
        "A": ("beta", (5.0, 5.0)),
        "B": ("mixture", [(0.85, 5.0, 5.0), (0.15, 20.0, 2.0)]),
    },
}


def _sample_alt_spec(rng, spec, n):
    kind = spec[0]
    if kind == "beta":
        a, b = spec[1]
        return rng.beta(a, b, size=n)
    elif kind == "mixture":
        comps = spec[1]
        weights = np.array([c[0] for c in comps])
        counts = rng.multinomial(n, weights)
        parts = [rng.beta(a, b, size=c) for (w, a, b), c in zip(comps, counts)]
        out = np.concatenate(parts)
        rng.shuffle(out)
        return out
    raise ValueError(kind)


def run_power_check(rng, alt_key, n_A, n_B, grid, n_outer, n_perm, alpha=0.05,
                     eps=EPS_DEFAULT):
    spec = ALTERNATIVES[alt_key]
    rej_A = rej_Aw = rej_ks = 0
    for _ in range(n_outer):
        xA = _sample_alt_spec(rng, spec["A"], n_A)
        xB = _sample_alt_spec(rng, spec["B"], n_B)
        pA_val = permutation_pvalue(xA, xB, grid, n_perm, rng, statistic="A_bar", eps=eps)
        pAw_val = permutation_pvalue(xA, xB, grid, n_perm, rng, statistic="A_w", eps=eps)
        ks = stats.ks_2samp(xA, xB)
        rej_A += pA_val < alpha
        rej_Aw += pAw_val < alpha
        rej_ks += ks.pvalue < alpha
    return {"power_A_bar_perm": rej_A / n_outer,
            "power_A_w_perm": rej_Aw / n_outer,
            "power_KS": rej_ks / n_outer}


# =============================================================== Section 13
# Epsilon / grid-resolution sensitivity (secondary, bounded check).
# =============================================================== Section 13

def run_epsilon_grid_sensitivity(rng, family_key, n_A, n_B, n_outer,
                                  eps_values=(1e-4, 1e-3, 1e-2),
                                  grid_sizes=(129, 257, 512)):
    rows = []
    for eps in eps_values:
        grid = np.linspace(*DOMAIN, 257)
        arrs = run_condition(rng, family_key, n_A, n_B, grid, n_outer,
                              shared_bandwidth=False, eps=eps)
        rows.append({"axis": "epsilon", "value": eps,
                     "median_A_bar": np.median(arrs["A_bar"]),
                     "median_A_w": np.median(arrs["A_w"])})
    for gs in grid_sizes:
        grid = np.linspace(*DOMAIN, gs)
        arrs = run_condition(rng, family_key, n_A, n_B, grid, n_outer,
                              shared_bandwidth=False, eps=EPS_DEFAULT)
        rows.append({"axis": "grid_size", "value": gs,
                     "median_A_bar": np.median(arrs["A_bar"]),
                     "median_A_w": np.median(arrs["A_w"])})
    return pd.DataFrame(rows)


# ===================================================================== main
def main():
    t0 = time.time()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    grid = np.linspace(*DOMAIN, 257)  # primary grid; see module docstring
    N_OUTER = int(os.environ.get("NULLCAL_N_OUTER", "1000"))

    rngs = make_condition_rngs()
    summary_rows = []
    for family_key in FAMILIES:
        for (n_A, n_B) in SIZE_CONFIGS:
            for bwcond in BANDWIDTH_CONDITIONS:
                key = (family_key, (n_A, n_B), bwcond)
                rng = rngs[key]
                shared = bwcond == "shared_pooled"
                res = run_condition(rng, family_key, n_A, n_B, grid, N_OUTER, shared)
                row = {"family": family_key, "n_A": n_A, "n_B": n_B,
                       "n_eff": 2.0 / (1.0 / n_A + 1.0 / n_B),
                       "imbalance_ratio": max(n_A, n_B) / min(n_A, n_B),
                       "bandwidth_condition": bwcond, "n_outer": N_OUTER}
                for stat_name in ("A_bar", "A_w", "B_bar", "ks_stat"):
                    for k, v in summarize(res[stat_name]).items():
                        row[f"{stat_name}_{k}"] = v
                row["abs_B_bar_median"] = float(np.median(np.abs(res["B_bar"])))
                summary_rows.append(row)
                print(f"[{time.time()-t0:7.1f}s] {family_key:24s} "
                      f"n=({n_A:5d},{n_B:5d}) {bwcond:15s} "
                      f"median A_bar={row['A_bar_median']:.4f} "
                      f"A_w={row['A_w_median']:.4f}")

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT_DIR / "null_calibration_summary.csv", index=False)
    print(f"[{time.time()-t0:7.1f}s] wrote null_calibration_summary.csv "
          f"({len(summary)} rows)")

    # ---- Section 10: permutation type-I error, representative subset -----
    calib_rows = []
    N_OUTER_CALIB = int(os.environ.get("NULLCAL_N_OUTER_CALIB", "200"))
    N_PERM = int(os.environ.get("NULLCAL_N_PERM", "200"))
    calib_conditions = [
        ("A_interior_unimodal", 200, 200, "equal_interior"),
        ("A_interior_unimodal", 150, 1350, "severe_imbalance_interior"),
        ("C_boundary_heavy", 150, 1350, "severe_imbalance_boundary"),
    ]
    calib_ss = np.random.SeedSequence(MASTER_SEED + 1)
    calib_children = calib_ss.spawn(len(calib_conditions) * 3)
    ci = 0
    for family_key, n_A, n_B, label in calib_conditions:
        for statistic in ("A_bar", "A_w", "ks_stat"):
            rng = np.random.default_rng(calib_children[ci]); ci += 1
            rate, se = run_type1_error_check(
                rng, family_key, n_A, n_B, grid, N_OUTER_CALIB, N_PERM, statistic)
            calib_rows.append({"condition": label, "family": family_key,
                                "n_A": n_A, "n_B": n_B, "statistic": statistic,
                                "n_outer_calib": N_OUTER_CALIB, "n_perm": N_PERM,
                                "empirical_type1_rate": rate, "monte_carlo_se": se})
            print(f"[{time.time()-t0:7.1f}s] type-I check {label:28s} "
                  f"{statistic:10s} rate={rate:.3f} (se={se:.3f})")
    calib_df = pd.DataFrame(calib_rows)
    calib_df.to_csv(OUT_DIR / "permutation_type1_error.csv", index=False)

    # ---- Section 9: planted-alternative power sanity check ---------------
    power_rows = []
    N_OUTER_POWER = int(os.environ.get("NULLCAL_N_OUTER_POWER", "150"))
    N_PERM_POWER = int(os.environ.get("NULLCAL_N_PERM_POWER", "150"))
    power_conditions = [(200, 200, "equal_interior"), (150, 1350, "severe_imbalance")]
    power_ss = np.random.SeedSequence(MASTER_SEED + 2)
    power_children = power_ss.spawn(len(ALTERNATIVES) * len(power_conditions))
    pi = 0
    for alt_key in ALTERNATIVES:
        for (n_A, n_B, label) in power_conditions:
            rng = np.random.default_rng(power_children[pi]); pi += 1
            res = run_power_check(rng, alt_key, n_A, n_B, grid,
                                   N_OUTER_POWER, N_PERM_POWER)
            res.update({"alternative": alt_key, "n_A": n_A, "n_B": n_B,
                        "condition": label, "n_outer": N_OUTER_POWER,
                        "n_perm": N_PERM_POWER})
            power_rows.append(res)
            print(f"[{time.time()-t0:7.1f}s] power {alt_key:24s} {label:18s} "
                  f"A_bar={res['power_A_bar_perm']:.2f} "
                  f"A_w={res['power_A_w_perm']:.2f} KS={res['power_KS']:.2f}")
    power_df = pd.DataFrame(power_rows)
    power_df.to_csv(OUT_DIR / "alternative_power_summary.csv", index=False)

    # ---- Section 13: epsilon / grid-resolution sensitivity ---------------
    sens_rng = np.random.default_rng(np.random.SeedSequence(MASTER_SEED + 3))
    N_OUTER_SENS = int(os.environ.get("NULLCAL_N_OUTER_SENS", "200"))
    sens_df = run_epsilon_grid_sensitivity(
        sens_rng, "A_interior_unimodal", 150, 1350, N_OUTER_SENS)
    sens_df.to_csv(OUT_DIR / "epsilon_grid_sensitivity.csv", index=False)
    print(f"[{time.time()-t0:7.1f}s] wrote epsilon_grid_sensitivity.csv")

    manifest = {
        "master_seed": MASTER_SEED, "grid_size_primary": 257,
        "n_outer": N_OUTER, "n_outer_calib": N_OUTER_CALIB, "n_perm": N_PERM,
        "n_outer_power": N_OUTER_POWER, "n_perm_power": N_PERM_POWER,
        "n_outer_sens": N_OUTER_SENS,
        "families": {k: v for k, v in FAMILIES.items()},
        "size_configs": SIZE_CONFIGS,
        "bandwidth_conditions": BANDWIDTH_CONDITIONS,
        "alternatives": {k: str(v) for k, v in ALTERNATIVES.items()},
        "total_runtime_seconds": time.time() - t0,
    }
    with open(OUT_DIR / "run_manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"[{time.time()-t0:7.1f}s] DONE. wrote run_manifest.json")


if __name__ == "__main__":
    main()
