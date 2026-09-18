# Novelty and provenance audit

**Date:** 2026-09-18 · **Repository state audited:** `88ca56f` (= `origin/main`), tag `v1.0.0` untouched
· **Bibliography:** [`references.bib`](../references.bib), 38 entries, every DOI resolved
· **Deposit audited alongside it:** `DeJesus2026disprot` (Zenodo v8; concept DOI `10.5281/zenodo.21613280`)

This is a literature, provenance, statistical-methodology and biological-interpretation audit. It was
asked to answer four questions, and the short answers are:

1. **Is the bounded asymmetry statistic mathematically distinct?** No. `S(x)` is a bounded, strictly
   monotone transform of the **density ratio**, its unweighted summary `Ā` is a continuous
   **Canberra**-type distance between the two density estimates, and its density-weighted summary
   `A_w` is **exactly the total-variation distance** at `ε = 0`. See §3–§5.
2. **Does the density-weighted form reduce to a standard distance?** Yes, exactly, not merely
   asymptotically. `A_w(ε=0) = TV`. See §5.
3. **Are the negative/calibration findings scientifically interesting and publishable?** Yes, and they
   are the strongest thing here — but the *biological* negative result is stronger than the
   statistical one, and it is stronger than the repository currently states. See §9 and §14.
4. **How should the biological conclusions be framed?** As a **protein-level calibration failure with
   taxon-specific offsets**, explicitly not as a verdict on residue-level pLDDT-based disorder
   prediction, where the published benchmarks are good. See §8–§11.

Three factual errors were found in the repository's own documentation. They are in §15. None of them
is an error in a computed result: **every number in the repository reproduces** (§13).

---

## 1. Claim inventory

Classification of every substantive claim, before any literature comparison.

| # | Claim | Where | Initial class |
|---|---|---|---|
| C1 | `S(x) = (p_A−p_B)/(p_A+p_B+2ε)` is bounded in (−1,1), antisymmetric, zero where estimates agree | README, manuscript §2 | mathematical identity (trivial) |
| C2 | `Ā`, `B̄`, `C = |B̄|/Ā` summarise the profile | README, `src/02_asymmetry.py` | proposed statistic |
| C3 | `A_w = ∫|S|p̄ / ∫p̄ → TV` as `ε→0` | README, addendum §4 | known statistical object |
| C4 | `Ā` inflates with histogram bin count (0.181 → 0.333) | manuscript §3.4 | estimator-dependent empirical quantity |
| C5 | `B̄` changes sign in 5/12 comparisons under histogram estimation | `sign_stability.csv` | estimator-dependent empirical quantity |
| C6 | `Ā`'s null median is 0.089–0.398 across the 12 comparisons | addendum | permutation-calibrated finding |
| C7 | Human–*E. coli* localisation is real | **original manuscript headline** | **withdrawn** |
| C8 | Two zero crossings mark real sign-change regions | original manuscript | **withdrawn as inferential** |
| C9 | The profile adds sensitivity beyond scalar two-sample tests | original manuscript | **reversed** |
| C10 | Fungal enrichment in curated disorder vs human and vs *E. coli* | manuscript §3.1 | permutation-calibrated finding (retained) |
| C11 | AlphaFold very-low-confidence fraction is a poor protein-level proxy (ρ = 0.07–0.24) | README finding 1 | negative result / proxy-validation result (retained, promoted to lead) |
| C12 | The proxy's offset vs curated disorder is taxon-dependent | README finding 1 | apparently distinct empirical result |
| C13 | `A_w` reduces `Ā`'s null magnitude for some shapes but not bimodal | `results/null_calibration/` | computational finding |
| C14 | `A_w` is not established as unbiased or superior | `results/null_calibration/` | negative result |
| C15 | Full-pipeline permutation calibration gives ≈nominal type-I error for both `Ā` and `A_w` | `permutation_type1_error.csv` | computational finding |
| C16 | 13 CSVs and 13 profile arrays value-identical to the pre-refactor run; 7 PNGs byte-identical | `results/README.md` | computational finding |

No claim in the repository is classified `claimed novelty`: the README makes **no novelty claim at
all** for the statistic. That is already the correct posture, and this audit confirms it is the
required one.

**Superseded claims are not mixed with current ones.** C7, C8 and C9 appear in the repository only
inside the revision notice, the addendum's disposition table, and this document — with one exception,
which is §15.2.

---

## 2. The exact statistic, as implemented

Pinned from `src/02_asymmetry.py` (not from the prose):

| Property | Value in the committed pipeline |
|---|---|
| `p_A`, `p_B` | `scipy.stats.gaussian_kde(x, bw_method="scott")`, evaluated on the grid, then **divided by their own trapezoidal integral** — so they are normalised densities on the domain |
| Bandwidth | **Scott's rule, selected independently within each group**; in 1-D `h = σ̂(ddof=1)·n^(−1/5)` |
| `ε` | **fixed at `1e-3`**, not data-dependent; sensitivity checked at `1e-6, 1e-2, 1e-1` |
| Domain | `[0,1]` — the feature is a fraction, so this is the natural support, and it is closed |
| Discretisation | `np.linspace(0, 1, 512)`, uniform; trapezoidal integration |
| Boundary handling | closed-domain renormalisation (not reflection) — documented in the code as an approximation |
| `Ā` | `∫|S| dx / (domain length)` → with length 1, the **uniform-measure mean of `|S|`** |
| `B̄` | `∫S dx /` (domain length) |
| `C` | `|B̄|/Ā` |
| `A_w` | `∫|S| p̄ dx / ∫p̄ dx`, `p̄ = (p_A+p_B)/2` — **defined only in the addendum and in `src/07`, not in `src/02`** |
| Uncertainty | 2,000-replicate protein-level bootstrap, percentile CIs, bandwidth re-selected inside each replicate |
| Comparators | `scipy.stats.ks_2samp`, `mannwhitneyu`, rank-biserial effect size |

Two facts matter for everything below: the bandwidths are **per-group**, and `Ā` weights the
**domain**, not the probability mass.

---

## 3. `S(x)` is a bounded density-ratio transform

Verified symbolically (SymPy; the script is quoted in §16).

With `R(x) = p_A(x)/p_B(x)` and `ε = 0`:

```
S = (p_A − p_B)/(p_A + p_B) = (R − 1)/(R + 1) = tanh( ½ log R )
```

Both equalities are exact. `d/dR[(R−1)/(R+1)] = 2/(R+1)² > 0`, so `S` is a **strictly increasing,
bounded reparameterisation of the density ratio** — equivalently, of the log density ratio squashed
through `tanh`. Where both densities are positive, `S(x)` and `R(x)` carry *identical* information:
neither can localise anything the other cannot.

**Consequence for novelty.** The pointwise profile is not a new object. It is the density ratio in
bounded coordinates. The literature on estimating exactly this object is large and explicit
(`SugiyamaEtAl2012` for density-ratio estimation; `SugiyamaEtAl2013` for direct density-*difference*
estimation, which also argues against forming a ratio of two separately-estimated densities — the
very thing done here).

**One precision point the repository should add.** With `ε > 0` the identity breaks:

```
S = q(R − 1) / (qR + q + 2ε)          (q = p_B)
```

which depends on `q` as well as `R`. So the regularised profile is **not** a function of the density
ratio alone; `ε` makes it sensitive to the absolute density level. This is not a defect — it is what
stops `|S| → 1` in empty regions — but it means "bounded density-ratio transform" is exact only in
the `ε → 0` limit, and the repository should say so if it adopts that framing.

---

## 4. `Ā` is a continuous Canberra-type distance

At `ε = 0`, `|S| = |p_A − p_B|/(p_A + p_B)`. The **Canberra metric** of `LanceWilliams1966` is

```
d_Canberra(a,b) = Σ_i |a_i − b_i| / (|a_i| + |b_i|)
```

so `Ā`, evaluated as it is on a uniform 512-point grid, is exactly the **grid mean of the Canberra
distance between the two evaluated density vectors**, with `ε` playing the role of the standard
zero-protection constant that every practical Canberra implementation needs. The `(p−q)/(p+q)`
normalisation itself is the Lance–Williams "relative difference" family (`LanceWilliams1967`), of
which `BrayCurtis1957` dissimilarity is the summed-numerator, summed-denominator member.

**Classification: `KNOWN / REPARAMETERIZED`.** Not a new statistic; a continuous, regularised,
uniformly-weighted Canberra distance between two kernel density estimates.

### 4.1 `Ā` has no probabilistic interpretation, and this is why it misbehaves

`Ā = E_U[|S|]` where `U` is the **uniform** measure on the domain — not either data distribution.
Three consequences follow directly, and all three were observed empirically by the repository:

1. **It is domain-dependent.** Widen the window, and `Ā` changes, because the uniform weight
   redistributes. It is not invariant under reparameterisation of `x`.
2. **It is dominated by low-mass regions.** In the tails, both densities approach `ε`, so `|S|`
   becomes the ratio of two small noisy numbers and drifts toward 1 — exactly the repository's own
   observation that empty histogram bins "drive `|S| → 1`". This is the *mechanism* of the
   finite-sample null inflation reported in the addendum, and it is structural rather than a
   small-sample accident.
3. **It is not an f-divergence between `p_A` and `p_B`.** It is a fixed-measure integral of a bounded
   contrast. So no standard divergence theory — no data-processing inequality, no known null
   distribution, no calibration results — transfers to it.

No prior work was located that proposes this specific uniformly-weighted integral as a two-sample
statistic. That is not a novelty finding in the repository's favour: the closest relatives
(`AndersonHallTitterington1994`, `Duong2013`) weight by density or test locally *because* the
unweighted version has the defects above.

---

## 5. `A_w` is exactly total variation (steps 6 and 45)

Symbolically, for normalised `p_A`, `p_B` on a domain of length 1:

```
|S_{ε=0}| · p̄ = [ |p_A − p_B| / (p_A + p_B) ] · (p_A + p_B)/2 = |p_A − p_B| / 2      (exact)
∫ p̄ dx = 1
  ⇒   A_w(ε = 0) = ½ ∫ |p_A − p_B| dx = TV(p_A, p_B)                                  (exact)
```

For `ε > 0` the integrand carries a pointwise shrinkage factor

```
(p_A + p_B) / (p_A + p_B + 2ε)  ∈ (0, 1]
```

whose derivative in `ε` is negative, so `A_w(ε)` **increases monotonically to `TV` as `ε ↓ 0`** and is
a strict under-estimate of `TV` for every `ε > 0`. This is stronger and cleaner than the
repository's "converges exactly to the total-variation distance as ε→0": it is an identity at
`ε = 0` plus an explicit, signed, one-sided regularisation error.

Note the scope of the identity: it is `TV` between the two **kernel density estimates**, not between
the underlying distributions. `A_w` is therefore a plug-in `TV` estimator, and plug-in `L¹`
functionals of KDEs are known to be biased upward under the null — which is exactly the pattern the
repository's simulation found. `GibbsSu2002` is the standard map of how `TV` relates to the other
metrics used here.

**Classification: `KNOWN / REPARAMETERIZED`.** `A_w` is a regularised plug-in total-variation
distance. Equivalently, it is the continuous **Bray–Curtis** dissimilarity (`BrayCurtis1957`) of the
two densities, which for normalised inputs coincides with `TV`.

### 5.1 An exact identity linking the two summaries

This audit's one mathematical contribution, and it is elementary. On a unit-length domain with both
densities normalised, `E_U[p̄] = 1`, so

```
A_w = E_U[|S| p̄] = Cov_U(|S|, p̄) + E_U[|S|]·E_U[p̄]

      ⇒   A_w = Ā + Cov_U(|S|, p̄)                                   (exact)
```

Verified numerically on all twelve real comparisons: residual ≤ `1.1e-16`.

This turns the repository's "STRUCTURAL INTERPRETATION (not a theorem)" in
`results/null_calibration/README.md` into an exact statement, and it explains the simulation results
that were left OPEN there:

* `A_w < Ā` **exactly when** `|S|` is negatively correlated with the pooled density — i.e. when the
  contrast is largest where there is least probability mass. `Cov_U < 0` in **all twelve** real
  comparisons, which is why `A_w < Ā` in all twelve.
* `A_w = Ā` **exactly when `p̄` is uniform** on the domain. This predicts the simulation's shape
  ordering without any further experiment: the benefit is large for peaked densities
  (Beta(5,5): ratio ≈ 0.42; Beta(0.7,3): ≈ 0.38) and small for flat-ish ones (Beta(2,2): ≈ 0.82;
  the bimodal mixture, whose `p̄` is nearly flat across `[0,1]`: ≈ 0.93). The reported "negligible
  reduction for bimodal densities" is therefore not a puzzle and not a failure of `A_w` — it is the
  identity evaluated at `Cov_U ≈ 0`.

**Exact conditions** (re-derived and re-checked 2026-09-18). Writing `U` for the uniform measure on
the domain `Ω` and `E` for expectation under it, the general form is

```
A_w = Ā + |Ω| · Cov_U(|S|, p̄)        whenever  ∫_Ω p_A = ∫_Ω p_B = 1
```

The identity requires exactly two things: both density estimates integrate to 1 over `Ω` (which the
pipeline enforces by dividing each KDE by its own trapezoidal integral), and `Ā` is defined with the
domain-normalised measure `∫|S|dx / |Ω|` (which `src/02_asymmetry.py` does). On this repository's unit
domain `|Ω| = 1` and `E_U[p̄] = 1` exactly, so it reduces to `A_w = Ā + Cov_U(|S|, p̄)`.

It does **not** require `ε = 0`: `S` enters only through `|S|`, whatever regularisation produced it.
Verified numerically for the Human–Bacterial `D_exp` pair at `ε = 0, 1e-6, 1e-3, 1e-2, 1e-1`, with
residuals of at most `1.4e-17` at every value, and across all twelve real comparisons at the
production `ε = 1e-3` (residual ≤ `1.1e-16`). The total-variation identity of §5, by contrast, *does*
require `ε = 0`, and was re-confirmed for the same pair: `A_w(0) = 0.0857335688 = ½∫|p_A − p_B|`.

`Classification: ELEMENTARY CONSEQUENCE.` It should be stated as a one-line algebraic remark, never
as a theorem.

---

## 6. Triangular discrimination (step 46)

At `ε = 0`, `S²·(p_A+p_B) = (p_A−p_B)²/(p_A+p_B)` exactly, so

```
Δ(p_A,p_B)  :=  ∫ (p_A−p_B)²/(p_A+p_B) dx  =  ∫ S² (p_A+p_B) dx  =  2 ∫ S² p̄ dx
```

`Δ` is **triangular discrimination**, a.k.a. the Vincze–Le Cam distance (`LeCam1986`, `Topsoe2000`),
an f-divergence with generator `f(t) = (t−1)²/(t+1)` — verified convex (`f'' = 8/(1+t)³ > 0`) with
`f(1) = 0` (`LieseVajda2006`).

So the family is fully accounted for by two standard divergences:

| Weighting of the pointwise contrast `S` | Resulting quantity |
|---|---|
| `∫ |S| · 1 dx` (uniform) | `Ā` — regularised continuous Canberra distance; **no divergence interpretation** |
| `∫ |S| · p̄ dx` (mass) | `A_w` — **exactly total variation** |
| `∫ S² · 2p̄ dx` (mass, squared) | **exactly triangular discrimination** |

The repository does not compute the squared form. It is worth recording as the reason `S` "feels"
like a divergence integrand: it is one — just not in the unweighted `L¹` combination actually used.

---

## 7. Ancestry of the bounded asymmetry statistic

| Step | Object | Standard? | Estimator-dependent? | Repo-specific? | Prior source |
|---|---|---|---|---|---|
| 1 | Two samples of a protein-level fraction on `[0,1]` | yes | — | no | — |
| 2 | Per-group Gaussian KDE, Scott's rule, closed-domain renormalisation | yes | **yes** — bandwidth, boundary rule | choice of renormalisation over reflection | `Scott1992`, `SheatherJones1991` |
| 3 | Pointwise normalised contrast `(p_A−p_B)/(p_A+p_B)` | yes | inherits 2 | no | `LanceWilliams1966`, `LanceWilliams1967` |
| 4 | `ε`-regularisation `+2ε` in the denominator | yes (standard Canberra practice) | **yes** — `ε` is a tuning constant | the specific value `1e-3` | `LanceWilliams1966` |
| 5 | Bounded-ratio reading `S = tanh(½ log R)` | yes | inherits 2 | no | `SugiyamaEtAl2012` |
| 6 | Unsigned uniform integral `Ā` | **not located as a two-sample statistic** | **yes, strongly** | **yes** | — (nearest: `AndersonHallTitterington1994`) |
| 7 | Signed uniform integral `B̄`, ratio `C` | no prior located | yes | **yes** | — |
| 8 | Mass-weighted integral `A_w` | yes — **is** `TV` | yes (plug-in bias) | no | `BrayCurtis1957`, `GibbsSu2002` |
| 9 | Squared mass-weighted integral | yes — **is** `Δ` | yes | not computed | `LeCam1986`, `Topsoe2000` |
| 10 | Full-pipeline permutation calibration with bandwidth re-selection | yes | resolves 2 and 6 | correct application | `HemerikGoeman2018`, `Holm1979` |
| 11 | Localisation via the profile shape | yes | yes | no | `Duong2013`, `Gretton2012` |

Steps 6 and 7 are the only repo-specific constructions, and step 6 is the one the calibration showed
to be problematic. **The novelty is concentrated exactly where the statistical behaviour is worst.**

---

## 8. Localisation versus power (steps 11–13)

### 8.1 What each comparator detects

| Method | Detects | Localises | Source |
|---|---|---|---|
| Kolmogorov–Smirnov | sup CDF gap | one point only | used in repo |
| Cramér–von Mises / Anderson–Darling | integrated CDF gap, tail-weighted | no | `ScholzStephens1987` |
| Energy distance | all differences, distribution-free | no | `SzekelyRizzo2013` |
| MMD (kernel two-sample) | all differences in the RKHS | **yes — via the witness function** | `Gretton2012` |
| Wasserstein | quantile displacement | partially | `GibbsSu2002` |
| Total variation | max event-probability gap | no (but `= A_w`) | `GibbsSu2002` |
| Local KDE two-sample test | density difference, pointwise, with multiplicity control | **yes, with confidence statements** | `Duong2013` |
| `S(x)` + `Ā`/`B̄` | bounded density contrast, uniform weight | **yes, descriptively** | this repo |

**`S(x)` is a normalised MMD-style witness function.** The MMD witness is the smoothed density
difference `p_A − p_B`; `S` is that difference divided by the pooled density. The standard
terminology for what the repository is doing is therefore **"witness function"** (`Gretton2012`) or
**"local two-sample testing"** (`Duong2013`), and the repository should adopt it — not because the
name matters, but because `Duong2013` already supplies what the repository lacks: pointwise
significance with multiplicity control for exactly this KDE-difference setting.

### 8.2 Does it add power? Reproduced answer: no

Recomputed from `revision_addendum/permutation_results.csv` and
`results/tables/thales_curated_sample_results.csv`. The addendum's Holm families are of size 6 (six
comparisons within each of the two features — confirmed: smallest raw `p_A` = 4.998e-4, its Holm
value 2.9985e-3, ratio exactly 6.00). Applying the **same** per-feature Holm correction to KS, for a
like-for-like comparison:

| Statistic | Holm-significant comparisons at α = 0.05 |
|---|---|
| `Ā` | **5** |
| `A_w` | **5** — *the identical set* |
| `T_max` | 3 |
| KS | **6** |
| KS, uncorrected | 7 |

`Ā`'s rejection set is a **strict subset** of KS's, and the subset relation survives whether KS is
Holm-corrected (5 ⊂ 6) or not (5 ⊂ 7). The comparison KS finds and `Ā` misses under matched
correction is **`D_exp`, Fungal–Nematode**. So the reversal the repository reports is not an artefact
of an unfair correction — it is robust. **`Ā` adds no power over KS in this dataset.**

### 8.3 Does it add localisation value? Qualified yes

Localisation is a real and separate service, and the repository is right to keep the profile as a
descriptive instrument. But the surviving value is weaker than "localisation" unqualified:

* The profile's localisation is **not calibrated** in `src/`. The simultaneous band and `T_max` live
  only in the addendum, and `A_w`, `T_max`, the bands and mass trimming are **not implemented in the
  committed pipeline** at all (the repository states this; `src/07` is a from-formula
  reimplementation of a subset, not a recovery of the missing `profile_inference.py`).
* For the one comparison the original paper used to *demonstrate* localisation (human–*E. coli*),
  the simultaneous band excludes zero **nowhere**.
* `Duong2013` already provides calibrated localisation for KDE differences.

Recommended framing: **a descriptive visualisation of a witness function, to be read only after a
calibrated global test has established that a difference exists** — which is what the README now
says, and it is the defensible position.

---

## 9. Bandwidth and sample-size sensitivity (steps 7–9)

The repository's finding — that `Ā`'s null baseline is materially nonzero and varies with shape and
sample size — is **expected**, not anomalous, and the literature says why:

* Plug-in `L¹`/`L²` functionals of kernel density estimates carry a positive bias under the null
  that decays slowly with `n`; `AndersonHallTitterington1994` analyses precisely this for two-sample
  KDE discrepancies and is the closest prior art for the whole construction.
* Bandwidth choice, not sample size alone, drives the bias. `Scott1992`'s rule is optimal for a
  single Gaussian-ish density under MISE — it is not designed to make *two* estimates comparable.
  `SheatherJones1991` is the standard better selector, and neither is a two-sample selector.
* Under `n_A ≠ n_B`, `h ∝ σ̂ n^(−1/5)` gives the two groups **different** smoothing. The repository's
  real comparisons run to 1339 vs 58 — a 23× imbalance, so a ≈ 1.9× bandwidth ratio from the `n`
  factor alone, before any difference in `σ̂`. Differential smoothing perturbs the two estimates
  differently, and `Ā`'s unsigned, uniformly-weighted integral **accumulates** those perturbations
  instead of cancelling them. The repository's own diagnosis is correct.
* `SugiyamaEtAl2013` is the direct argument against this whole pipeline shape: estimate the
  difference (or ratio) *directly* rather than subtracting two separately-tuned density estimates.

**On the repository's own H4 test.** `results/null_calibration/README.md` reports that forcing a
single pooled-Scott bandwidth did *not* remove the pathology, and correctly declines to conclude that
differential smoothing is therefore not the mechanism, because the pooled rule also changes the
absolute smoothing level. That reasoning is sound, and the proposed cleaner interventions (max or
mean of the two per-group bandwidths) are the right next experiment. This audit adds one more
candidate from the literature: a **two-sample-specific** bandwidth chosen to minimise the variance of
the *difference*, rather than the MISE of either density.

---

## 10. Permutation calibration (step 10)

| Item | Value |
|---|---|
| Original claim | human–*E. coli* curated-disorder localisation (paper's methodological headline) |
| Null | labels permuted, group sizes preserved, pooled sample |
| Pipeline inside each replicate | **entire** estimation recomputed, bandwidth re-selected — never frozen |
| Replicates | 2,001 permutations (smallest attainable raw p = 4.998e-4) |
| Statistics | `Ā`, `A_w`, `T_max` |
| Correction | Holm, family size 6 (per feature) |
| Type-I error under a true null | 0.020–0.030 for `Ā`, 0.015–0.065 for `A_w`, at 3 representative conditions (SE ≈ 0.015) |
| Comparator | KS, Mann–Whitney |
| Outcome for the headline claim | `Ā = 0.1431` vs null median 0.1107, 95th pct 0.1742, raw p = 0.178 → **inside its own null** |

The procedure is correct, and the distinction the repository draws is the right one and worth
preserving explicitly: **a large positive null magnitude does not invalidate permutation inference.**
Raw `Ā` is uninterpretable on its own; permutation-calibrated `Ā` has approximately nominal size.
Those are different statements and the repository keeps them apart.

`Classification: the calibration itself is a REPLICATION / EXTENSION of standard practice
(HemerikGoeman2018, Holm1979) — correctly applied, not novel.`

---

## 11. Biological input → repository interpretation (step 37)

| Input / proxy | Source | What it actually measures | Repository use | Limitations | Ground truth? |
|---|---|---|---|---|---|
| DisProt consensus "Structural state" = `D` | `Quaglia2022`, `Aspromonte2024`, `Piovesan2017` | residues with **published experimental evidence** of disorder, manually curated | `D_exp = n_D / length` per protein | annotation is **partial by construction** — only studied regions are annotated, so `D_exp` is a lower bound on true disorder for most proteins | **No.** "curated disorder annotation", never biological truth |
| `alphafold_very_low_content` | field as deposited in the DisProt export; the underlying models are `Jumper2021` as served by `Varadi2022` | **presumed** fraction of residues with pLDDT < 50 (the AlphaFold DB "very low" confidence band, `Varadi2022`) — the repository states this reading is **unverified** | the proxy whose concordance is tested | never retrieved from AlphaFold DB (HTTP 403 at the egress proxy, logged); threshold, model version, and residue mapping all unconfirmed | **No.** Provenance is second-hand |
| Taxon label | `ncbi_taxon_id` inside each record | strain-level NCBI taxon (9606, 559292, 83333, 6239) | group assignment; the script refuses to run if a file's records disagree with its name | two of the four are strain-level, so "Fungal"/"Bacterial" overstate scope | **Yes for the label**, no for the clade |
| Protein-level disorder fraction | derived | ratio of annotated disordered residues to chain length | the feature compared across taxa | aggregation destroys residue-level structure; composition and domain architecture confound it | **No** |

### 11.1 Do not equate low pLDDT with disorder (step 17)

The repository's language is already correct — it says "AlphaFold-derived very-low-confidence
fraction" and "surrogate", never "low pLDDT = disorder". The literature supports keeping it that way:
`RuffPappu2021` and `WilsonEtAl2022` set out why low confidence is not the same as disorder (missing
homologous information, conditional folding, flexible-but-structured regions, model limitations), and
`PiovesanEtAl2022` is the quantitative treatment of disorder and **conditional folding** in
AlphaFold DB. `Akdel2022` documents the community-scale caveats. One wording recommendation only:
prefer "**very-low-confidence fraction**" over "proxy" in the title and abstract, and reserve
"proxy" for places where the surrogate relationship is the actual subject.

---

## 12. The biological negative result is stronger than stated

Reproduced exactly (§13). Per-taxon concordance between curated `D_exp` and the very-low-confidence
fraction:

| Taxon | n | Spearman ρ | 95% bootstrap CI | OLS slope (af on `D_exp`) | median signed gap |
|---|---|---|---|---|---|
| Human | 1246 | 0.0708 | 0.0105 – 0.1328 | 0.101 (SE 0.021) | +0.0011 |
| Fungal | 219 | 0.2356 | 0.0893 – 0.3691 | 0.128 (SE 0.041) | −0.0405 |
| Bacterial | 143 | 0.1767 | 0.0054 – 0.3449 | 0.079 (SE 0.024) | −0.0833 |
| Nematode | 54 | 0.1763 | −0.1146 – 0.4663 | 0.302 (SE 0.117) | +0.0195 |

Pooled ρ = 0.102 (n = 1662) — so the weak association is *within* taxa, not a Simpson artefact.

### 12.1 Range compression, not just weak correlation (steps 21, 23)

A decile calibration of the human data, which the repository does not report, shows the proxy is
**severely compressed** rather than merely noisy:

| `D_exp` decile (median) | 0.017 | 0.035 | 0.052 | 0.074 | 0.106 | 0.147 | 0.203 | 0.290 | 0.461 | 0.910 |
|---|---|---|---|---|---|---|---|---|---|---|
| median proxy | 0.252 | 0.178 | 0.132 | 0.116 | 0.124 | 0.155 | 0.162 | 0.205 | 0.284 | 0.269 |
| signed gap | **+0.235** | +0.143 | +0.080 | +0.041 | +0.018 | +0.008 | −0.041 | −0.084 | −0.177 | **−0.641** |

Across a 53-fold range in curated disorder, the proxy stays inside `0.12–0.28`. The sharpest form of
the same fact:

* 77 proteins whose **entire** annotated length is curated disorder (`D_exp = 1.0`) have median proxy
  **0.148** (human, n = 54; min 0.0).
* 375 proteins with `D_exp < 0.05` have median proxy **0.160**.

**The proxy does not separate fully-disordered proteins from nearly-ordered ones.** That is a
stronger, more interpretable and more useful statement than "ρ = 0.07–0.24", it is a calibration
statement rather than an association statement, and it is exactly what §19 of the brief asks for.
Every OLS slope is far below 1 (0.08–0.30), which is the same finding in regression form.

**Caution on direction.** Part of the positive gap at low `D_exp` is expected and benign: curated
annotation is partial, so `D_exp` under-counts disorder while the proxy covers the whole chain. But
that mechanism cannot explain the negative gap at high `D_exp` — a protein annotated disordered
end-to-end should show a *high* very-low-confidence fraction, not 0.15. Two readings remain open and
the repository cannot currently distinguish them, because no residue-level pLDDT was ever retrieved:
either AlphaFold is confident about many experimentally-disordered chains (conditional folding,
`PiovesanEtAl2022`), or `alphafold_very_low_content` is not the quantity its name suggests
(different threshold, different denominator, or partial residue mapping). **Resolving this requires
retrieving residue-level pLDDT, and it should be stated as the single highest-value next step.**

### 12.2 Protein-level versus residue-level (steps 19, 20, 24, 25)

The distinction is essential and the repository already draws it. The literature is unambiguous that
residue-level performance is good: pLDDT-based disorder prediction is competitive in `Necci2021`
(CAID round 1), `DelConte2023` (CAID 2) and `Mehdiabadi2025` (CAID 3), and `PiovesanEtAl2022`
quantifies the relationship directly. A weak **protein-level fraction** correlation is compatible
with strong **residue-level** discrimination, because:

* aggregation to a fraction discards all positional information and compresses variance;
* the denominators differ — curated fraction is over annotated regions, the proxy over the modelled
  chain;
* conditional folding puts genuinely-disordered residues in the confident band;
* composition and domain architecture differ systematically by taxon.

So the correct scope is: **the very-low-confidence fraction is not a calibrated taxon-neutral
protein-level disorder feature.** Nothing here licenses a claim about residue-level disorder
prediction, and the repository must not be read as contradicting CAID.

### 12.3 Is the weak correlation already known? (steps 24, 28)

Partly. That pLDDT is not interchangeable with disorder is well established
(`RuffPappu2021`, `WilsonEtAl2022`, `PiovesanEtAl2022`, `Akdel2022`), and benchmark context is in
`Necci2021`, `DelConte2023`, `Mehdiabadi2025`, `Piovesan2023`, `Erdos2021`, `Hu2021`. What was **not**
located in the literature is a quantification of the **protein-level fraction** relationship
**stratified by taxon**, with signed offsets. So:

* "low pLDDT ≠ disorder" → `STANDARD BIOINFORMATICS` / `REPLICATION`;
* "protein-level fraction agreement is weak" → `REPLICATION / EXTENSION`;
* "**the offset is taxon-dependent, and the proxy is range-compressed**" →
  `APPARENTLY DISTINCT EMPIRICAL RESULT`, with the caveat that it may be a property of the DisProt
  export's field rather than of AlphaFold.

### 12.4 Taxonomic disorder variation (step 15)

Higher disorder in eukaryotes than bacteria is a long-established pattern
(`Ward2004`, `XueEtAl2012`, `Peng2015`), so the cross-taxon *distribution* differences the profile
describes are **biologically expected** and are not a discovery of this work. Two constraints:

* Sample sizes are 1339 / 223 / 146 / 58 — no rationale for the choice of taxa is documented beyond
  availability, and 58 proteins cannot characterise a nematode proteome.
* DisProt is a **curated, literature-driven** sample, not a random draw from any proteome
  (`Quaglia2024` gives the curation practice; the repository's own `data/README.md` correctly calls
  its ascertainment bias "undocumented rather than small"). So these distributions must never be
  described as proteome-level estimates, and no causal evolutionary explanation follows from them.

---

## 13. Numerical reproduction (steps 29, 44)

Full pipeline re-run from `data/raw` in a clean virtual environment.

| Check | Result |
|---|---|
| `python run_all.py` | exit 0, all six scripts |
| CSV outputs byte-identical to committed | **11 of 13** |
| CSV outputs differing | 2 — `thales_curated_sample_results.csv`, `robustness_sensitivity.csv`, max abs diff **9.0e-16** (max rel 1.4e-14) |
| Stored profile arrays (`asymmetry_profiles.npz`) | 13/13 keys match, max abs diff **1.7e-15** |
| Figure PNGs byte-identical | 0 of 7 — matplotlib 3.11.2 here vs 3.10.8 in the lock file |
| `Ā` values vs the addendum's independent `A` column | **max difference exactly 0.0** across all 12 comparisons |
| Protein / region / paired counts | 1766 / 7535 / 1662 — all match; 103 proteins lack a proxy value |
| Per-taxon counts | 1339 / 223 / 146 / 58 — match |
| Spearman ρ, gaps, KS, Mann–Whitney, bootstrap CIs | match to the digits published |

Deviations are attributable to the environment, not the code: the locked versions
(Python 3.12.3, NumPy 2.4.4, pandas 3.0.2, SciPy 1.17.1, matplotlib 3.10.8) could not be installed on
this machine's Python 3.14.7, so the run used NumPy 2.5.3 / pandas 3.0.6 / SciPy 1.18.1 /
matplotlib 3.11.2. **A 2,000-replicate bootstrap reproducing to 1e-15 across a NumPy, SciPy and
pandas minor-version change is unusually strong evidence of determinism.**

### 13.1 Software import-path audit (step 30) — clean, with two documentation gaps

| Check | Result |
|---|---|
| `conftest.py` anywhere | **none** — no pytest shadowing |
| How tests import the module under test | `importlib.util.spec_from_file_location` on the real `src/07_null_calibration_grid.py` — the shipped file, not a copy |
| Duplicate stale implementations | none stale; `src/02` and `src/07` deliberately hold **two** implementations of `S`, `Ā`, `B̄`, cross-checked against `scipy.stats.gaussian_kde` to `atol=1e-6` in the test suite, and the reason is documented (`src/02`'s numeric filename is not importable) |
| Scripts and notebooks share code paths | no notebooks; `run_all.py` invokes `src/01`–`src/06` with `sys.executable` |
| Generated data paths | explicit, `DISPROT_DATA` / `DISPROT_OUT`, all `pathlib`, no absolute paths |
| Test coverage of the manuscript pipeline | **zero.** 19 tests cover `src/07` (a diagnostic that produces no manuscript number); `src/01`–`src/06`, which produce every table and figure, have **no tests** |

The duplicate-implementation risk is real but currently managed: the two implementations can drift,
and only one of them is tested. The cheapest durable fix is to rename `src/02_asymmetry.py` to an
importable module name and have `src/07` and the tests import it.

---

## 14. What the DisProt project actually adds

| Feature | Known statistically? | Known biologically? | Existing benchmark result? | Repo contribution | Publishable value | Evidence |
|---|---|---|---|---|---|---|
| Pointwise bounded contrast `S(x)` | **yes** — density-ratio transform, Canberra family | — | — | bounded coordinates + `ε` regularisation | low on its own | §3, §4 |
| `Ā` (uniform integral) | no exact prior located, but its defects are why | — | — | the construction, and the finding that it is estimator-dependent | **as a caution, yes** | §4.1, §9 |
| `B̄`, `C` | no prior located | — | — | sign and consistency diagnostics | low–moderate | §7 |
| `A_w` | **yes — exactly TV** | — | — | none mathematically | none as a statistic | §5 |
| `A_w = Ā + Cov_U(|S|,p̄)` | elementary | — | — | explains the simulation's shape-dependence | moderate, as a remark | §5.1 |
| Permutation calibration with bandwidth re-selection | standard | — | — | correct application to this statistic | **yes — the self-correction is the story** | §10 |
| "Profile does not beat KS here" | — | — | — | a clean, documented negative result | **yes** | §8.2 |
| Weak protein-level pLDDT/disorder agreement | — | partly known | CAID, `PiovesanEtAl2022` | replication in a new stratification | **yes** | §12 |
| **Taxon-dependent offset + range compression** | — | **not located** | not located | the strongest empirical claim | **yes — likely the lead** | §12.1 |
| Reproducibility record (13 CSVs, 13 arrays, 7 PNGs) | — | — | — | exemplary | **yes, as an artifact** | §13 |

### 14.1 Which outcome is it? (step 42)

**Outcome 3 and Outcome 4 jointly, with Outcome 5 as the honest packaging.** Not Outcome 1: the
statistic is not distinct (§3–§6). Not Outcome 2 alone: the DisProt application is distinct, but the
distinct part is the *proxy calibration result*, not the statistic. The defensible headline:

> After full-pipeline permutation calibration, a bounded normalised density contrast — which is a
> transform of the density ratio, whose mass-weighted form is exactly total variation — does not
> outperform Kolmogorov–Smirnov on these data, its rejection set being a strict subset of KS's. The
> more consequential result is biological: the AlphaFold-derived very-low-confidence fraction is not
> a calibrated protein-level measure of curated disorder in any of four taxa, showing taxon-specific
> offsets and, more sharply, a compressed dynamic range that fails to separate fully-disordered from
> nearly-ordered proteins. This is a calibration failure at the protein-fraction level and is not a
> statement about residue-level disorder prediction, where published benchmarks are good.

This is **publishable as a methodological-caution / negative-result report**, which is what the
repository's own title already says it is. It is not publishable as a new-statistic paper, and the
repository does not claim to be one.

---

## 15. Errors found in the repository's documentation

None of these changes a computed number. The first two matter.

### 15.1 The `A_w` Holm-rejection claim is wrong

`README.md` (finding 3, note on `A_w`) states that `A_w`

> "in the same table … produces *fewer* Holm-significant rejections than `Ā`, not more."

Recomputed from `revision_addendum/permutation_results.csv`: `p_A_holm < 0.05` in **5** rows,
`p_Aw_holm < 0.05` in **5** rows, and **the two sets are identical**. At uncorrected α = 0.05, `A_w`
has *more* (6 vs 5 — it additionally reaches `p_Aw = 0.034` on `D_exp` Fungal–Nematode). So the
sentence is wrong under both corrections. The surrounding argument — that `A_w` is not established as
superior — is unaffected and still correct; only this comparison is wrong. **Corrected in this PR.**

### 15.2 The "revised manuscript" is not revised

`README.md` and `revision_addendum/README.md` describe
`revision_addendum/bounded-asymmetry-profiles-cross-taxon-protein-disorder-2026.pdf` as

> "the revised manuscript incorporating the additional statistical calibration and related
> clarifications"

and direct readers to it for "the current, authoritative account". It is not that document:

* the string "permut" does **not appear anywhere** in it;
* a normalised text diff against `paper/manuscript.pdf` differs **only** in reference line-wrapping
  and one Zenodo DOI (`…21613281` → `…21627908`);
* its abstract still reads "the profile behaves as intended as a descriptive instrument: for the
  human–*E. coli* comparison … it resolves localised structure (Ā = 0.143, C = 0.57) in a comparison
  that global tests read as null" — which is verbatim the claim the addendum **withdraws**.

So the repository's withdrawn headline claim **does** currently survive in a PDF that the README
presents as authoritative. The only document carrying the calibration is
`ADDENDUM_permutation_calibration.pdf`, whose §5 disposition table is explicit and correct (two
withdrawn, one withdrawn-as-inferential, one rewritten, one reversed, the rest retained).
**This PR corrects the descriptions and points readers at the addendum; it does not modify any PDF.**

### 15.3 Smaller items

| Item | Status |
|---|---|
| `README.md` "Repository layout" and "Results map" omit `src/07_null_calibration_grid.py`, `tests/`, and `results/null_calibration/` — all added in `88ca56f` | **corrected in this PR** |
| `results/null_calibration/README.md` cites a `pyproject.toml` that does not exist in the repository | **corrected in this PR** |
| `requirements.txt` caps `pandas>=2.2,<3` but `requirements-lock.txt` pins `pandas==3.0.2`, which the cap excludes | **corrected in this PR** (cap raised to `<4`; the lock is the version actually verified) |
| Zenodo **concept DOI** TODO in `README.md` and `CITATION.cff` | **resolved**: DataCite records `10.5281/zenodo.21628406` as `IsVersionOf` **`10.5281/zenodo.21613280`**, which is the concept DOI. Filled in. |
| `data/README.md` TODO for exact DisProt query strings and retrieval dates | **left open** — only the author knows these; still flagged |
| Manuscript title calls the field a "proxy" | wording recommendation only (§11.1); no PDF changed |

---

## 16. Claim-strength audit (steps 39, 40)

Repository-wide sweep for `new`, `novel`, `first`, `sensitive`, `more sensitive`, `detects`,
`localizes`, `robust`, `universal`, `disorder proxy`, `predicts disorder`, `validates`, `superior`,
`significant`, `biologically meaningful`.

| Wording | Where | Verdict |
|---|---|---|
| "novel" | `results/null_calibration/README.md` ×2 | **SAFE** — both are denials ("no claim that `A_w` or the total-variation identity is novel") |
| "universal", "superiority" | `results/null_calibration/README.md` ×6 | **SAFE** — all are explicit disclaimers |
| "outperformed KS" | `results/null_calibration/README.md` line 139 | **NEEDS QUALIFICATION** — true, but only on a predeclared *synthetic* variance/shape alternative at equal `n`; it must not be read across to the DisProt comparisons, where the opposite holds (§8.2) |
| "poor protein-level proxy" | `README.md` finding 1 | **SAFE** — correctly scoped to protein level |
| "descriptive instrument, not currently an inferential one" | `README.md` finding 2 | **SAFE** — and the right framing |
| "estimator-dependent" | `README.md` finding 3 | **SAFE** — supported by `sign_stability.csv` and the null grid |
| "resolves localised structure" | **both manuscript PDFs' abstracts** | **REMOVE** — this is the withdrawn claim (§15.2); the repository's own README already contradicts it |
| "proxy" (as a bare noun) | title, abstract, `README.md` ×2 | **NEEDS QUALIFICATION** — prefer "very-low-confidence fraction"; reserve "proxy" for the surrogate relationship itself |
| "validated" | `results/null_calibration/README.md` line 98 | **SAFE** — refers to numerically validating trapezoidal integration, not to validating the proxy |
| "significant" | 4 occurrences | **SAFE** — each accompanied by a test and a correction |
| "ground truth", "predicts disorder", "biologically meaningful", "more sensitive" | — | **absent** |

On the word **proxy** specifically (step 40): the data support a *weak, taxon-dependent, poorly
calibrated surrogate*, not a calibrated estimator, not a validated screening heuristic, and nothing
at all at residue level. The repository's usage is already close to correct; the recommendation is
only to prefer the concrete noun over the evaluative one.

---

## 17. Calibration history (step 38)

This table is the record of the project correcting itself. It should be kept, not hidden — it is a
substantial part of what makes the work publishable.

| Stage | Original claim | Calibration performed | Result | Claim retained? | Current wording |
|---|---|---|---|---|---|
| 1 | The profile resolves structure global tests miss (methodological headline) | permutation, bandwidth re-selected per replicate, 2001 reps | sole support was human–*E. coli*, which does not reject | **No — withdrawn** | "descriptive only"; README revision notice |
| 2 | Human–*E. coli* localisation in curated disorder | same | `Ā = 0.1431` inside its own null (median 0.1107, p95 0.1742, raw p = 0.178) | **No — withdrawn** | quoted only as a descriptive example |
| 3 | Two zero crossings mark sign-change regions | simultaneous null band | band excludes zero nowhere in this comparison | **Withdrawn as inferential**, survives as description | "two zero crossings" reported without inference |
| 4 | Profile adds sensitivity beyond scalar two-sample tests | rejection-set comparison vs KS | `Ā` set (5) ⊊ KS set (6 Holm / 7 raw) | **Reversed** | "the original claim … is reversed, not supported" |
| 5 | Fungal enrichment vs human and vs *E. coli* | permutation + Holm | `Ā` Holm p = 1.5e-2 and 3.0e-3 | **Retained, now supported** | reported as calibrated |
| 6 | AlphaFold very-low fraction is a poor protein-level surrogate | unchanged numerics; three statistics agree | ρ = 0.07–0.24, taxon-dependent offsets | **Retained, promoted to lead result** | README finding 1, "gating result" |
| 7 | `Ā` is estimator-dependent | histogram/bandwidth/`ε` sweep + null grid | inflates 0.181 → 0.333 with bins; null median 0.089–0.398 | **Retained, strengthened** | README finding 3 |
| 8 | `A_w` improves on `Ā` | 56,000-dataset synthetic grid | shape-conditional; negligible for bimodal; no imbalance protection | **Retained only as "context-dependent"** | `results/null_calibration/README.md` |
| 9 | *(this audit)* `A_w` yields fewer Holm rejections than `Ā` | recomputation | identical sets, 5 = 5 | **No — wrong, corrected** | §15.1 |

---

## 18. Novelty and provenance after Zotero review

| Result / claim | Repository location | Closest prior work | Equivalent standard concept | Data/estimator dependence | Status | Classification | Confidence | Recommended wording | Keys |
|---|---|---|---|---|---|---|---|---|---|
| `S(x)` bounded, antisymmetric, zero at agreement | README; `02_asymmetry.py:77` | `LanceWilliams1966` | normalised/relative difference | none (algebraic) | formal | `CLASSICAL` | high | "a bounded normalised density contrast" | `LanceWilliams1966`, `LanceWilliams1967` |
| `S = tanh(½ log R)` at `ε=0` | not stated | `SugiyamaEtAl2012` | density-ratio transform | none | formal, verified | `KNOWN / REPARAMETERIZED` | high | "a bounded monotone transform of the density ratio" | `SugiyamaEtAl2012`, `SugiyamaEtAl2013` |
| `Ā` = uniform mean of `|S|` | `02_asymmetry.py:81` | `LanceWilliams1966`; nearest two-sample: `AndersonHallTitterington1994` | continuous Canberra distance | **strong** — bandwidth, `ε`, grid, domain | empirical | `KNOWN / REPARAMETERIZED` (statistic) + `COMPUTATIONAL FINDING` (its instability) | high | "a regularised continuous Canberra distance between KDEs" | `LanceWilliams1966`, `AndersonHallTitterington1994` |
| `A_w → TV` | README; addendum §4 | `BrayCurtis1957`, `GibbsSu2002` | **total variation** (exactly, at `ε=0`) | plug-in KDE bias | formal, verified | `KNOWN / REPARAMETERIZED` | high | "exactly the total-variation distance between the density estimates at `ε=0`" | `BrayCurtis1957`, `GibbsSu2002` |
| `∫S²(p+q) = Δ` | not stated | `LeCam1986`, `Topsoe2000` | triangular discrimination / Vincze–Le Cam | — | formal, verified | `CLASSICAL` | high | "the mass-weighted square is triangular discrimination" | `LeCam1986`, `Topsoe2000`, `LieseVajda2006` |
| `A_w = Ā + Cov_U(|S|,p̄)` | **this audit** | — | covariance decomposition | — | formal, verified | `ELEMENTARY CONSEQUENCE` | high | "an algebraic identity, not a theorem" | — |
| `Ā` estimator-dependence | manuscript §3.4 | `Scott1992`, `SheatherJones1991`, `AndersonHallTitterington1994` | plug-in `L¹` bias; per-group bandwidth | **strong** | empirical | `COMPUTATIONAL FINDING` | high | keep as-is | `Scott1992`, `SheatherJones1991` |
| Nonzero null baseline, imbalance-driven | addendum; null grid | `AndersonHallTitterington1994`, `SugiyamaEtAl2013` | plug-in bias under unequal smoothing | **strong** | empirical | `COMPUTATIONAL FINDING` / `NEGATIVE RESULT` | high | keep; cite the bias literature | `AndersonHallTitterington1994` |
| Profile ⊄ adds power over KS | addendum; §8.2 here | `Duong2013`, `Gretton2012` | rejection-set comparison | — | empirical, reproduced | `NEGATIVE RESULT` | high | "rejection set is a strict subset of KS's" | `Duong2013`, `Gretton2012` |
| Localisation value | README finding 2 | `Duong2013`, `Gretton2012` | witness function / local two-sample test | — | descriptive | `KNOWN / REPARAMETERIZED` | med-high | "a witness-function visualisation, read after a calibrated global test" | `Duong2013`, `Gretton2012` |
| Permutation calibration procedure | addendum | `HemerikGoeman2018`, `Holm1979` | permutation test + Holm | resolves the above | empirical | `REPLICATION / EXTENSION` | high | keep | `HemerikGoeman2018`, `Holm1979` |
| Weak protein-level proxy agreement | README finding 1 | `PiovesanEtAl2022`, `RuffPappu2021`, `WilsonEtAl2022` | proxy calibration | dataset-dependent | empirical, reproduced | `REPLICATION / EXTENSION` | high | "not a calibrated protein-level measure" | `PiovesanEtAl2022`, `Necci2021` |
| Taxon-dependent offset | README finding 1 | none located | — | dataset-dependent | empirical, reproduced | `APPARENTLY DISTINCT EMPIRICAL RESULT` | med | "taxon-specific offsets; direction differs by taxon" | `Ward2004`, `XueEtAl2012`, `Peng2015` |
| **Range compression of the proxy** | **this audit, §12.1** | none located | calibration / dynamic range | dataset-dependent | empirical | `APPARENTLY DISTINCT EMPIRICAL RESULT` | med | "does not separate fully-disordered from nearly-ordered proteins" | `PiovesanEtAl2022` |
| Cross-taxon disorder differences | manuscript §3.1 | `Ward2004`, `XueEtAl2012`, `Peng2015` | taxonomic disorder variation | curation-biased sample | empirical | `STANDARD BIOINFORMATICS` | high | "consistent with the established eukaryote/bacteria pattern" | `Ward2004`, `Peng2015` |
| Residue-level pLDDT performance | not claimed | `Necci2021`, `DelConte2023`, `Mehdiabadi2025` | CAID benchmarks | — | — | `STANDARD BIOINFORMATICS` | high | "published benchmarks are good; this work says nothing about them" | `Necci2021`, `DelConte2023`, `Mehdiabadi2025` |
| Reproducibility record | `results/README.md` | — | — | — | computational | `COMPUTATIONAL FINDING` | high | keep | — |

---

## 19. Statistical status of the bounded asymmetry method

| Quantity | Known equivalent? | Estimator-dependent? | Permutation-calibrated? | Adds power? | Adds localisation? | Recommended use |
|---|---|---|---|---|---|---|
| `S(x)` | **yes** — `tanh(½ log R)`; normalised witness function | yes (bandwidth, `ε`) | n/a (a profile) | n/a | **yes, descriptively** | plot it, after a calibrated global test; call it a witness function |
| `Ā` | **yes** — continuous Canberra distance | **yes, strongly** | yes, and then ≈ nominal size | **no** — strict subset of KS | indirectly | descriptive; inferential **only** with full-pipeline permutation |
| `B̄` | no exact prior located | moderate (sign flips under histograms) | yes; bootstrap agrees on 12/12 | not tested | direction only | the more trustworthy of the two summaries |
| `C = |B̄|/Ā` | no prior located | inherits both | not separately | no | descriptive consistency check |
| `A_w` | **yes — exactly TV** | yes (plug-in bias) | yes; ≈ nominal size | no (identical rejection set to `Ā`) | no | report as "plug-in TV", not as a new statistic |
| `T_max` | sup-type statistic | yes | yes | fewer rejections (3) | yes, with a band | the right tool for localisation claims — **but not implemented in `src/`** |
| `∫S²·2p̄` | **yes** — triangular discrimination | yes | not computed | unknown | no | optional remark only |

---

## 20. Manuscript revision classification (step 41)

**D — major novelty reduction, statistical method mostly known; with an E component.**

* **D**, because `S(x)` is a density-ratio transform, `Ā` is a Canberra-type distance, `A_w` is
  exactly total variation, and the squared mass-weighted form is triangular discrimination. None of
  the statistical machinery is new, and the manuscript should cite the prior art rather than present
  the construction as a definition without ancestry.
* **E component**, because both PDFs' abstracts still assert the withdrawn localisation claim
  (§15.2). That is a correction, not a reframing, and it is the one thing that must change before
  any further deposit.

Not C, because the required change is more than framing. Not A or B, because the abstract carries a
retracted claim.

The good news is that the *repository* has already done most of the work: the README's three current
findings, its revision notice, and the addendum's disposition table are all correct and
well-calibrated. What remains is to propagate that into the manuscript and to add the literature.

**No Zenodo record was updated, no new manuscript version was prepared, and no PDF was modified by
this audit.** Recommended sequence, for the author to decide on: (1) fix the two documentation errors
(done in this PR); (2) issue a manuscript version whose abstract matches the addendum; (3) add the
prior-art citations from `references.bib`; (4) retrieve residue-level pLDDT to settle §12.1.

---

## 21. Unresolved

1. **Whether `alphafold_very_low_content` means what its name implies.** The pLDDT threshold, model
   version, denominator and residue mapping are all unverified, and §12.1 shows a pattern that is
   hard to explain if the field is simply "fraction pLDDT < 50". AlphaFold DB returned HTTP 403 from
   this machine as well, so this audit could not settle it either. **Highest-value next step.**
2. The exact DisProt query strings and retrieval dates (`data/README.md` TODO) — author-only
   knowledge.
3. Whether a two-sample-specific bandwidth (or the max/mean of the two per-group Scott bandwidths)
   removes `Ā`'s null inflation. The repository's H4 experiment is inconclusive for a reason it
   correctly identifies.
4. Whether `B̄`, which no prior work was located for, has a usable null theory. It is the summary
   that behaved best here, and it is the least examined.
5. Whether the taxon-specific offsets survive in a non-curated sample. Untestable within DisProt.

---

## 22. How this audit was carried out

* Every DOI in `references.bib` was resolved and the entry generated from the registrar's own
  metadata (Crossref `/transform/application/x-bibtex`, DataCite for the Zenodo deposit) rather than
  retyped. Two entries have no DOI (`Gretton2012`, JMLR; `Holm1979`, pre-DOI) and are marked.
* The four identities in §3, §5, §5.1 and §6 were verified **symbolically** with SymPy, not
  numerically — as the brief required. The `A_w = Ā + Cov` identity was additionally checked against
  the twelve real profiles (residual ≤ 1.1e-16).
* All reported repository numbers were recomputed from `data/raw` (§13).
* No alternative statistic was searched for after seeing that `Ā` underperforms. The comparators in
  §8.1 are the standard ones named in the brief, and the only new quantities computed were the decile
  calibration of §12.1 and the Holm/KS rejection sets of §8.2 — both diagnostics of existing claims,
  neither a replacement statistic.

---

# ALPHAFOLD PROXY PROVENANCE

Added 2026-09-18, after a dedicated validation audit. **This section corrects §12.1 above.** The
range-compression claim made there does not survive, and the reason is a protein-length confound.
The claims withdrawn are named explicitly in §31.

## 23. What `alphafold_very_low_content` actually is

**Verdict A — exact direct fraction.** It is the fraction of residues whose AlphaFold pLDDT is
below 50, with the number of modelled residues as denominator.

### 23.1 Provenance chain

| Step | Where | What happens |
|---|---|---|
| final CSV column `af_very_low_content` | `results/tables/disprot_protein_level.csv` | written by `01_build_tables.py` |
| transformation | `src/01_build_tables.py:162` | `pd.to_numeric(..., errors="coerce")` — **type coercion only** |
| assignment | `src/01_build_tables.py:108,127` | `af = r.get("alphafold_very_low_content")` — copied verbatim |
| input field | `data/raw/disprot_taxon_*.json` | DisProt record key `alphafold_very_low_content` |
| provider | **DisProt** release 2026_06, "with ambiguous evidences" per-organism JSON export | a DisProt-computed field, not an AlphaFold DB field |
| upstream quantity | AlphaFold DB per-residue pLDDT | not present anywhere in this repository |

**The repository performs no calculation on this field at all** — no threshold, no denominator, no
rescaling. Its semantics are entirely DisProt's, and the repository was right to call its own reading
"unverified". Two further provenance facts:

* the field is **absent from DisProt's single-record API** (`/api/DP00086` has no such key) and
  present only in the search/export endpoint, i.e. it is an export-layer annotation;
* **DisProt does not document it.** Its `/about`, `/help` and `/statistics` pages and its API return
  no definition, so the name is the only provider-side documentation available. This is why the
  question had to be settled empirically.

### 23.2 The export is reproducible

Re-downloading the *E. coli* export today from
`https://disprot.org/api/search?release=current&show_ambiguous=true&show_obsolete=false&format=json&ncbi_taxon_id=83333`
yields a file whose **MD5 is `df0b07151610e12043d687362cb935c2` — identical to the committed
`data/raw/disprot_taxon_83333_escherichia_coli_k12.json`** and to the hash recorded in
`data/metadata/source_queries.csv`. The export is stable and the committed data are faithful.

### 23.3 The denominator (step 10 of the brief)

Settled two ways, and they agree:

1. **Arithmetically.** For **1652 of 1662** records the stored value times DisProt's `length` is an
   integer to 1e-6 — so the denominator is the full sequence length, not a sub-count.
2. **Empirically.** For all **1661** proteins whose AlphaFold model was retrieved, the number of
   modelled residues equals DisProt's `length` exactly: **coverage = 1.000 for every protein**.

So options A, B and C of the brief coincide here: full sequence length = modelled residues = mapped
residues. There is no denominator mismatch, and denominator mismatch is **not** the explanation for
anything.

### 23.4 The threshold is 50, read off the provider's own data

AlphaFold DB's confidence files carry a `confidenceCategory` array alongside the scores. Reading the
band edges directly out of that array rather than from memory:

| Category | pLDDT range observed | Meaning |
|---|---|---|
| `D` | 32.78 – 49.97 | very low |
| `L` | 50.19 – 68.81 | low |
| `M` | 70.00 – 89.19 | confident |
| `H` | 90.00 – 98.69 | very high |

And across all 1661 reconstructions, the fraction of category-`D` residues equals the fraction with
pLDDT < 50 **exactly** (max difference 0.00e+00). So "very low" is `pLDDT < 50`, on the provider's own
definition.

### 23.5 Stored versus directly reconstructed (step 11)

Reconstruction is by `analysis/validate_alphafold_proxy.py`, which fetches
`AF-<acc>-F1-confidence_v6.json` from AlphaFold DB and counts residues below the threshold.

| Quantity | Value |
|---|---|
| proteins with a stored value | 1662 |
| AlphaFold models retrieved | **1661** (1 failure: `O43236-1`, an isoform accession with no AlphaFold entry) |
| coverage = 1.0 | 1661 / 1661 |
| Spearman(stored, reconstructed f50) | **0.9937** |
| Pearson | **0.9960** |
| median \|stored − f50\| | **0.0047** |
| mean (stored − f50) | **−0.0011** |
| median stored / f50 | **1.0000** |
| exact matches (< 1e-9) | 335 |
| agreement within 0.01 | 1184 / 1661 |
| agreement within 0.05 | 1611 / 1661 |
| max \|stored − f50\| | 0.1418 |

The residual scatter is **AlphaFold model-version drift**, not a semantic difference: DisProt computed
its value against an earlier AlphaFold release, and AlphaFold DB now serves only `v6` (v1–v5 return
404), so the historical values cannot be re-fetched. The 335 exact matches are the entries whose
models have not changed. Nothing about the scatter is systematic — the mean difference is −0.001 and
the median ratio is exactly 1.

**Conclusion.** `alphafold_very_low_content` = `count(pLDDT < 50) / n_residues`. The repository's
reading of the field was correct. Verdict **A**, with the model-version caveat recorded.

## 24. The biological finding does not survive as previously stated

### 24.1 Protein length is a strong confounder, acting in opposite directions

| Association | Spearman |
|---|---|
| length vs curated disorder fraction `D` | **−0.396** |
| length vs low-pLDDT fraction `f50` | **+0.414** |
| `D` vs `f50` (pooled — the published figure) | **+0.104** |
| `D` vs `f50`, **partialling out length** | **+0.321** |

Short DisProt entries have a high annotated fraction; long proteins have more low-confidence
residues. The two associations pull in opposite directions and **suppress** the correlation. Within
length bands the association is far stronger than the pooled figure: ρ = 0.45 (150–300 aa),
0.40 (300–500), 0.27 (500–800).

The published ρ = 0.07–0.24 values are arithmetically correct. Their interpretation as "a poor
surrogate" was not: a large part of the weakness is a length-suppression artefact.

### 24.2 The extreme-group comparison reverses under length matching

The two groups were never comparable. Median length is **136 aa** in the fully-annotated-disordered
group and **676 aa** in the near-zero group; 43 of 77 of the first group are under 150 aa and **none**
of the 375 in the second is under 207 aa.

| Comparison (reconstructed `f50`) | n | medians | AUROC | Cliff's δ | p |
|---|---|---|---|---|---|
| pooled — as previously reported | 77 / 375 | 0.153 vs 0.164 | **0.479** | −0.041 | 0.57 |
| **length-matched, 207–929 aa** | 24 / 265 | **0.426 vs 0.127** | **0.810** | **+0.620** | 5.0e-07 |
| length-matched, human only | 17 / 199 | 0.483 vs 0.164 | 0.781 | +0.561 | 1.2e-04 |

Within the fully-disordered group the low-pLDDT fraction tracks length strongly
(Spearman **+0.747**): median `f50` is 0.031 below 100 aa, 0.144 at 100–200 aa, 0.285 at 200–400 aa
and **0.640** at 400–1000 aa. The previous "the proxy cannot separate the extremes" conclusion was an
artefact of comparing a short-peptide-dominated group against a long-protein group.

The biology of the short-protein failure is known: AlphaFold2 assigns confident predictions to many
genuinely disordered regions that fold **conditionally** — on binding or in context
(`Alderson2023`). A confident model over a short disordered peptide is an expected outcome, not a
prediction error.

### 24.3 Range compression is not supported

| Set | statistic | IQR(P)/IQR(D) | interdecile ratio |
|---|---|---|---|
| all 1661 | `f50` | **1.046** | 0.809 |
| length ≥ 300 aa | `f50` | 1.640 | 1.165 |
| all 1661 | `f70` | 1.401 | 1.054 |
| length ≥ 300 aa | `f70` | 2.016 | 1.397 |

The proxy's spread is **as large as or larger than** the curated fraction's. It is not compressed into
a narrow band. The earlier appearance of compression came from a decile table of `D` whose bins are
effectively a length gradient — median length falls from **882 aa** in the lowest-disorder decile to
**228 aa** in the highest.

### 24.4 The failure is specific to the `< 50` threshold

A pre-specified robustness check at AlphaFold DB's other published band edge:

| Statistic | Spearman with `D` | pooled extreme-group AUROC | length-matched AUROC |
|---|---|---|---|
| `f50` (pLDDT < 50) | +0.104 | 0.479 | 0.810 |
| **`f70` (pLDDT < 70)** | **+0.315** | **0.757** | **0.924** |
| mean pLDDT (exploratory) | −0.257 | — | — |
| median pLDDT (exploratory) | −0.246 | — | — |

So the problem is **not** whole-protein aggregation as such. Aggregating at the "low confidence"
boundary (< 70, i.e. categories `D` + `L`) discriminates the extremes well even without length
matching. The very-low band alone is simply too strict to capture disorder that AlphaFold models with
moderate confidence. The `f70` figures are a robustness check at a documented band edge, not a tuned
threshold, and the mean/median pLDDT rows are marked exploratory.

### 24.5 What does not change: it is not a calibrated estimator

| Metric (all 1661) | `f50` | `f70` |
|---|---|---|
| descriptive slope of P on D | **+0.110** | +0.368 |
| slope, length ≥ 300 aa | +0.326 | — |
| MAE \|P − D\| | 0.201 | 0.201 |
| median absolute deviation | 0.119 | 0.123 |
| median signed (P − D) | −0.011 | +0.061 |

A slope far below 1 with MAE ≈ 0.20 means the low-pLDDT fraction is **not** a one-to-one estimator of
curated disorder content, even where it correlates. That part of the repository's conclusion stands.

### 24.6 Taxon offsets survive, attenuated

Restricted to 300–800 aa, which removes the taxon/length mix (median length 266 aa bacterial vs
508 aa nematode):

| Taxon | n | median `D` | median `f50` | median gap | ρ |
|---|---|---|---|---|---|
| Human | 624 | 0.110 | 0.186 | **+0.017** | +0.280 |
| Fungal | 97 | 0.146 | 0.180 | −0.018 | +0.460 |
| Bacterial | 59 | 0.078 | **0.025** | **−0.038** | +0.223 |
| Nematode | 30 | 0.097 | 0.132 | −0.007 | +0.418 |

The human-versus-bacterial offset difference narrows from 0.084 (pooled) to 0.055 but does not vanish.
Bacterial proteins receive markedly fewer very-low-confidence residues at comparable curated disorder.
**Taxon-dependent offsets are retained as a finding**, with length controlled.

## 25. Residue-level control, and why DisProt limits it

Across 626 proteins with cached models, mapping DisProt's `Structural state` type-`D` segments onto
the pLDDT array:

| Residue class | n | median pLDDT | fraction < 50 | fraction < 70 |
|---|---|---|---|---|
| annotated disordered | 57,866 | **48.4** | 0.524 | 0.720 |
| not annotated | 265,117 | **90.9** | 0.169 | 0.235 |

AUROC for low pLDDT marking an annotated-disordered residue: **0.784**.

**This is not a clean benchmark and must not be reported as one.** DisProt's annotation covers a
median of only **12.7%** of each protein; across the whole dataset there are **2,694 disorder
segments and only 4 structured segments**, and just 84 of 1766 proteins are fully covered. The
negative class above is therefore *unassessed residues*, not experimentally ordered ones, so the
figure is a descriptive contrast, not a validation. Its direction is nonetheless clear and consistent
with the CAID literature (`Necci2021`, `DelConte2023`, `Mehdiabadi2025`): **residue-level signal is
strong while protein-level fraction agreement is moderate at best.**

## 26. The two extreme groups, as they actually are

**Both group definitions were introduced by the September 2026 audit, not by the repository's code or
manuscript**, and both need restating.

* **The "77 fully disordered" group** is `D_exp >= 1.0`, i.e. every residue of the DisProt length is
  covered by a type-`D` consensus segment. That is a genuine full-length disorder annotation, but the
  group is dominated by short entries (median 136 aa, minimum 24 aa) and includes classic short
  disordered peptides such as P62328 (thymosin β-4, 44 aa) which AlphaFold models confidently. It
  should be described as *"proteins whose entire annotated length is curated as disordered, mostly
  short"*.
* **The "375 almost no disorder" group** is `D_exp < 0.05`. These are **not ordered proteins.** They
  are long proteins (median 676 aa) carrying a small amount of disorder annotation — median 2
  annotated regions and 17 annotated residues. Given 12.7% median annotation coverage, most of their
  sequence is simply unassessed. They must be described as *"proteins with little curated disorder
  annotation"*, never as ordered.

This is the deeper problem with the whole protein-level calibration framing: `D_exp` divides annotated
disordered residues by the **full** length, so it is a lower bound on disorder whose tightness varies
per protein — and it varies **with length**, which is exactly the confounder of §24.1.

## 27. Missing values, and one counting correction

| Fact | Value |
|---|---|
| proteins with no usable proxy value | **104**, not 103 |
| of which JSON `null` | 103 |
| of which the literal string `"NaN"` | 1 (`P09651-2`) |
| isoform-suffixed accessions in the dataset | 56 |
| isoform-suffixed accessions lacking a value | **52 of 56** |
| of the 104 missing, isoform-suffixed | 52 |

The README and `failed_accessions.csv` report 103. The 104th is `P09651-2`, whose stored value is the
string `"NaN"`; `pd.to_numeric(..., errors="coerce")` correctly turns it into a missing value, so no
computed result is affected — but the count should read 104. Missingness is concentrated on isoform
accessions, which is coherent: AlphaFold DB has no model for most isoforms.

## 28. Literature: the narrow question

The brief's question was not whether pLDDT relates to disorder, but whether the **fraction of
residues below a pLDDT threshold has been validated as an estimator of the whole-protein fraction of
curated disorder**. Searching OpenAlex for that specific claim
(`protein-level disorder content AlphaFold pLDDT fraction correlation`,
`fraction of residues pLDDT below 50 disorder content DisProt`,
`AlphaFold confidence whole-protein disorder fraction estimator`, and
`pLDDT threshold disorder content per protein agreement curated`) returns **no such validation
study**. The adjacent literature is about residue- and region-level behaviour:

* `Alderson2023` — conditionally folded IDRs predicted confidently by AlphaFold2. The mechanism
  behind the short-protein failure mode; **added to the bibliography by this audit**.
* `PiovesanEtAl2022` — disorder and conditional folding across AlphaFold DB.
* `Necci2021`, `DelConte2023`, `Mehdiabadi2025` — CAID rounds 1–3, residue-level performance.
* `RuffPappu2021`, `WilsonEtAl2022`, `Akdel2022` — why confidence is not disorder.

Keeping the three claim levels of the brief apart:

| Level | Claim | Status here |
|---|---|---|
| **A — residue** | low pLDDT carries information about disorder at individual residues | **supported**, and by the CAID literature. Nothing here contradicts it |
| **B — segment** | low-confidence stretches overlap disordered regions | **supported**; not tested in detail here |
| **C — whole-protein fraction** | the low-pLDDT fraction estimates the curated disordered fraction across proteins | **this is the repository's subject.** Correlated but not calibrated; strongly length-dependent; threshold-sensitive |

## 29. Biological verdict: BIO-B (qualified negative result)

Field semantics are confirmed, the reconstruction agrees, and the mappings are sound — but annotation
coverage, protein length and threshold choice **substantially limit the generality** of the negative
result, and two specific sub-claims from the previous audit are withdrawn outright (§31). BIO-A is
unavailable because the extreme-group result does not survive; BIO-C does not apply because the field
*is* a direct pLDDT fraction; BIO-D would overstate matters because the calibration failure itself is
real; BIO-E does not apply because the question was answered.

## 30. Biological wording now justified

> In this DisProt sample the AlphaFold-derived very-low-confidence fraction (`pLDDT < 50`) is a
> correlated but poorly calibrated protein-level index of curated disorder content: the descriptive
> slope against the curated fraction is 0.11 (0.33 among proteins of at least 300 residues) with a
> mean absolute deviation of 0.20, and the residual offset differs by taxon. Its apparent weakness is
> substantially a protein-length artefact — length is negatively associated with the curated fraction
> and positively with the low-confidence fraction, so the pooled Spearman correlation of 0.10 rises
> to 0.32 when length is held constant, and a comparison of fully-annotated-disordered against
> sparsely-annotated proteins reverses from an AUROC of 0.48 to 0.81 under length matching. The
> failure is also specific to the very-low band: using the low-confidence boundary (`pLDDT < 70`)
> instead raises the same AUROC to 0.76 without length matching. Because DisProt annotates a median
> of 12.7% of each sequence and almost never annotates ordered residues, the curated fraction is a
> per-protein lower bound of varying tightness, and these figures should be read as agreement between
> two partial indices rather than as validation against ground truth. Nothing here bears on
> residue-level disorder prediction, where low pLDDT separates annotated-disordered from unassessed
> residues with an AUROC of 0.78 in the same data and where published benchmarks are good.

## 31. Claims withdrawn by this audit

Withdrawn from §12.1 above, and from the README summary added in the first audit commit:

1. **"The proxy does not separate fully-disordered from nearly-ordered proteins."** Withdrawn. It was
   a length confound; under length matching the AUROC is 0.81. The underlying numbers (medians 0.148
   and 0.160) are correct but not comparable.
2. **"The very-low-confidence fraction is range-compressed."** Withdrawn. IQR ratio is 1.05 pooled and
   1.64 for proteins of at least 300 residues — not compressed.
3. **"Across a 53-fold span of curated disorder the proxy's median moves only within 0.12–0.28."**
   Withdrawn as an interpretation: the deciles of curated disorder are a length gradient (median
   length 882 aa down to 228 aa), so that table measures length, not calibration.

Retained and strengthened: the field's semantics (now verified), the absence of one-to-one
calibration, the taxon-dependent offsets (attenuated but present), and the fact that residue-level
signal coexists with moderate protein-level agreement.
