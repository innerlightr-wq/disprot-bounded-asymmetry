# Synthetic null-calibration audit of A_bar vs. A_w

STATISTICAL-METHOD DIAGNOSTIC. Documents the experiment already executed by
`src/07_null_calibration_grid.py` (run manifest: `run_manifest.json`, master
seed `20260906`, total wall-clock 915.2s, exit code 0). No DisProt data, no
biological claim, no network access. This file documents results already
produced; it does not add, rerun, or alter any simulation.

## Research question

Does the density-weighted asymmetry magnitude `A_w` (defined in
`revision_addendum/ADDENDUM_permutation_calibration.pdf`, section 4)
genuinely reduce the finite-sample null pathology of the unweighted
magnitude `A_bar` across heterogeneous synthetic sampling regimes, or was
the apparent improvement specific to the twelve real DisProt comparisons
previously examined?

## Preregistered simulation design

**Four density families**, fixed on `[0,1]`, not changed after seeing
results:

| Key | Distribution |
|---|---|
| `A_interior_unimodal` | Beta(5, 5) |
| `B_broad_interior` | Beta(2, 2) |
| `C_boundary_heavy` | Beta(0.7, 3) |
| `D_bimodal` | 0.5·Beta(8, 2) + 0.5·Beta(2, 8) |

**Seven sample-size configurations**: (50,50), (200,200), (800,800) —
equal-size controls; (50,100), (50,250), (50,500) — unequal-size controls;
(150,1350) — severe imbalance, approximating the repository's real
Human-vs-*E. coli* scale.

**Two bandwidth conditions**: `independent` (each sample gets its own
1-D Scott bandwidth, `h = sigma_hat(ddof=1) * n**(-1/5)` — the real-pipeline
default, matching `src/02_asymmetry.py` exactly) and `shared_pooled` (one
bandwidth derived from the pooled sample's Scott rule, applied to both
groups — a mechanistic control only, never proposed as a replacement
estimator).

**`N_outer = 1000`** independent null datasets per (family × size × bandwidth
condition) cell = 56 cells = 56,000 null datasets for the main grid.
Additional, smaller replicate counts for the permutation-calibration subset
(`N_outer_calib=200`, `N_perm=200`, 3 representative conditions), the
planted-alternative power check (`N_outer=150`, `N_perm=150`, 3 alternatives
× 2 conditions), and the epsilon/grid-resolution sensitivity check
(`N_outer=200`, 6 settings). **Master RNG seed `20260906`**, deterministically
spawned per condition via `numpy.random.SeedSequence`, so every condition is
independently reproducible regardless of run order.

Primary grid resolution: 257 points on `[0,1]` (production uses 512 — a
documented, tractability-driven reduction; grid-resolution sensitivity was
checked directly, including grid=512, and found small: ~2–3% relative
variation for both statistics — see `epsilon_grid_sensitivity.csv`).
`EPS_DEFAULT = 1e-3`, matching `src/02_asymmetry.py` exactly.

## Output files

| File | Contents |
|---|---|
| `null_calibration_summary.csv` | 56 rows: median/mean/SD/q90/q95/q99 of A_bar, A_w, B_bar, KS statistic per condition |
| `permutation_type1_error.csv` | 9 rows: empirical type-I error at alpha=0.05 for permutation-calibrated A_bar, A_w, KS at 3 representative conditions |
| `alternative_power_summary.csv` | 6 rows: rejection rate under 3 predeclared alternatives x 2 conditions |
| `epsilon_grid_sensitivity.csv` | 6 rows: median A_bar/A_w across 3 epsilon values and 3 grid resolutions |
| `run_manifest.json` | Exact configuration and total runtime for this run |

## Reimplementation status

**FORMULA-LEVEL REPRODUCTION: YES.** `src/07_null_calibration_grid.py`
implements `A_bar`, `B_bar`, and `A_w` exactly as defined in
`ADDENDUM_permutation_calibration.pdf` section 4, and its Scott-bandwidth
KDE was cross-checked numerically against `scipy.stats.gaussian_kde(x,
bw_method="scott")` (`atol=1e-6`, see `tests/test_null_calibration_grid.py`).

**ORIGINAL-CODE REPRODUCTION: NO / UNAVAILABLE.** The addendum describes a
script, `profile_inference.py`, that reportedly implemented permutation
tests, `A_w`, `T_max`, Holm correction, simultaneous bands, and mass
trimming. No such file, and no prior `A_w` implementation of any kind,
exists anywhere in this repository's committed history (confirmed by
exhaustive `grep`/`find` across every tracked `.py`/`.md`/`.csv` file before
writing this module). `src/07_null_calibration_grid.py` is therefore a
**new reimplementation from the published formula**, not a recovery of the
original missing script, and must not be described as such. `T_max`,
simultaneous bands, and mass trimming are out of scope for this experiment
and were not implemented.

## Environment / dependency note

The active environment's NumPy is **1.26.4**, while this repository's
`pyproject.toml`/`requirements.txt` declare NumPy>=2.0. This is an
**ENVIRONMENT/DEPENDENCY MISMATCH**, not a defect in `src/02_asymmetry.py`:
that script's use of `np.trapezoid` (added in NumPy 2.0) would raise
`AttributeError` in this specific environment, but that is a property of
the environment relative to the repository's declared requirement, not
evidence the shipped code is wrong for its intended target. This diagnostic
uses `np.trapz` instead (available in both NumPy 1.x and 2.x); its
trapezoidal integration was independently validated numerically — it
matches the closed-form integral of `x^2` on `[0,1]` to `1.7e-11` and an
independently-written manual trapezoidal-rule sum to machine precision
(`5.6e-17`).

## Test status

19/19 tests pass in `tests/test_null_calibration_grid.py`, covering:
degenerate-case zeros, swap symmetry, bounds, the TV-distance limit as
epsilon shrinks, the scipy Scott-bandwidth cross-check, shared-bandwidth
swap-symmetry, exact seeded reproducibility, and permutation-machinery
sanity (group-size preservation, p-values in `[0,1]`).

---

## Results summary

### A_w verdict: C — CONTEXT-DEPENDENT REWEIGHTING

Not A (not a universally strong replacement), not D (the reduction is real
and substantial in most tested shapes), not B or E. Specifically:

- **Strong raw-null reduction** for narrow-unimodal (`A_interior_unimodal`,
  median ratio A_w/A_bar ≈0.42) and boundary-heavy (`C_boundary_heavy`,
  ≈0.38) densities.
- **Modest reduction** for broad-interior (`B_broad_interior`, ≈0.82).
- **Negligible reduction** for bimodal densities (`D_bimodal`, ≈0.93; at the
  smallest equal sample size the q95 ratio is slightly *above* 1 — A_w's
  worst-case tail can be no better than A_bar's for this shape).
- **No general protection against the strongest imbalance pathology**: the
  one family with a genuine imbalance-specific null-magnitude penalty beyond
  what effective sample size alone predicts (`D_bimodal`, excess ratio up to
  ~1.7 at severe imbalance) hits `A_bar` and `A_w` almost equally (mean
  excess 1.42 vs. 1.35) — `A_w` does not shield against it.
- **Approximately calibrated permutation inference**: at all 3 representative
  conditions tested, permutation-calibrated `A_w`'s empirical type-I error
  was within Monte Carlo uncertainty of the nominal 5% (0.015–0.065 across
  conditions, SE≈0.015).
- **No systematic power loss** relative to `A_bar` across the 3 predeclared
  alternatives tested (location shift, variance/shape change, localized
  mixture change) — both statistics reached ≥92% power in every cell, and
  both clearly outperformed KS on the variance/shape-change alternative at
  equal sample sizes (100% vs. 79%).
- **Improved epsilon robustness**: at the representative severe-imbalance
  condition, `A_w`'s null median varied ~8.9% across the tested epsilon
  range vs. `A_bar`'s ~27.7%.
- **No claim of universal superiority is supported by this experiment.** The
  benefit is real but shape-conditional, and the one clearest failure mode
  (bimodal densities, worst-case imbalance) is exactly where a user would
  most want protection.

### A_bar classification: DESCRIPTIVE + CALIBRATED INFERENTIAL USE

Raw `A_bar` carries a large, shape- and size-dependent positive null
magnitude (medians 0.031–0.279 across tested conditions) and must not be
interpreted inferentially on its own. However, **raw positive null magnitude
is not the same thing as invalid permutation-calibrated inference**: at
every one of the 3 representative conditions tested, full-pipeline
permutation calibration (bandwidth reselected inside every replicate, never
frozen from the observed data) produced an empirical type-I error within
Monte Carlo uncertainty of the nominal 5% (0.020–0.030) for `A_bar`
specifically. Raw `A_bar` is therefore descriptive only; properly calibrated
`A_bar` is usable for inference.

### Bandwidth-mechanism result — H4 FAILED UNDER THE POOLED-SCOTT CONTROL USED HERE

This experiment does **not** establish that differential bandwidth selection
is irrelevant to the null-inflation mechanism, and that claim must not be
made from these results. What was tested and found:

Forcing one Scott bandwidth, derived from the **pooled** sample, onto both
groups did not "nearly remove the pathology" as H4 predicted. In the
majority of conditions across families A, B, and C, the shared-pooled
bandwidth **increased** rather than decreased both `A_bar`'s and `A_w`'s
null medians relative to independent per-sample bandwidths, and the
`A_bar - A_w` gap did not collapse under it.

The reason this is not a clean test of "differential smoothing" as a
mechanism: the pooled-Scott control changes two things at once — (1) it
equalizes the bandwidth between the two groups, and (2) it changes the
*absolute* smoothing strength, since a bandwidth computed from `n_A + n_B`
(the pooled size) is systematically narrower than what the minority sample
would select on its own under severe imbalance. Under this specific control,
the minority sample is effectively under-smoothed relative to its own
appropriate scale — a confound distinct from "differential smoothing between
groups" per se. **The mechanism question therefore remains OPEN.** A cleaner
common-bandwidth intervention — for example, using the larger of the two
individual per-group Scott bandwidths, or their average, rather than one
derived from the pooled sample — is a candidate future experiment, not
attempted here.

### General methodological result (restrained)

Supported by this experiment: unsigned KDE contrast, aggregated uniformly
over the domain, can carry a substantial, shape-dependent positive
finite-sample null magnitude. Pooled-density weighting can materially reduce
that magnitude for some density shapes, but not universally, and not for
bimodal densities in particular. `A_w` is not shown to be a debiased or
unbiased estimator, and is not shown to be a universal replacement for
`A_bar` — neither claim is made here.

**STRUCTURAL INTERPRETATION (not a theorem):** a useful way to think about
the difference is that `A_w` is not merely a corrected `A_bar` — it is a
mass-weighted view of the density discrepancy, whereas `A_bar` weights
domain location uniformly regardless of how much probability mass is there.
This framing is offered as interpretive scaffolding, not as an established
mathematical result beyond the exact TV-distance limit already verified
(module docstring, `tests/test_null_calibration_grid.py`).

### Scaling analysis (equal-n conditions, independent bandwidth)

Convention: fitting `median(A) ~ C * n^{-alpha}` by log-log linear
regression, `alpha` is reported as a **positive** decay exponent (so a
*larger* `alpha` means the null magnitude shrinks *faster* as `n` grows).
The raw regression coefficient (the log-log slope) is the negative of
`alpha`; both are shown below to keep the two conventions unambiguous.

| Family | alpha (A_bar) | log-log slope (A_bar) | alpha (A_w) | log-log slope (A_w) |
|---|---|---|---|---|
| A_interior_unimodal | 0.283 | -0.283 | 0.370 | -0.370 |
| B_broad_interior | 0.364 | -0.364 | 0.362 | -0.362 |
| C_boundary_heavy | 0.335 | -0.335 | 0.386 | -0.386 |
| D_bimodal | 0.344 | -0.344 | 0.373 | -0.373 |

Fit uses exactly 3 points per family (n=50,200,800); this is a descriptive
summary only, not an asymptotic theorem. `A_w`'s exponent is modestly larger
than `A_bar`'s in 3 of 4 families (a genuine, if small, faster-decay
property, not just a smaller multiplicative constant); for
`B_broad_interior` the exponents are effectively identical and the benefit
is purely a smaller constant.

## What this milestone does not claim

No claim that `A_w` is a general solution to the null-magnitude problem, no
claim of uniform power superiority over `A_bar` or KS, no claim that `A_w`
or the total-variation identity is novel, no claim of original-code
equivalence with the missing `profile_inference.py`, no claim that
differential bandwidth selection is or is not the mechanism in general, and
no generalization beyond the four density families and seven sample-size
configurations actually tested.
