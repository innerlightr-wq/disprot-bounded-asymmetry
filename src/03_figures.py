#!/usr/bin/env python3
"""Figures for the bounded asymmetry manuscript. Inputs are the CSV/NPZ
outputs of 01_build_tables.py and 02_asymmetry.py.

Portability note (2026-07-27 refactor): absolute paths replaced with
project-relative pathlib paths plus DISPROT_OUT override; figures are written
to results/figures. No plotted quantity, style setting or figure dimension was
changed.
"""
import itertools
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

# ------------------------------------------------------------------- paths
ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = Path(os.environ.get("DISPROT_OUT", ROOT / "results"))
TABLES_DIR = RESULTS_DIR / "tables"
DIAG_DIR = RESULTS_DIR / "diagnostics"
FIG = RESULTS_DIR / "figures"
FIG.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({"font.size": 9, "axes.grid": True,
                     "grid.alpha": 0.25, "figure.dpi": 160,
                     "savefig.bbox": "tight"})

GROUPS = ["Human", "Fungal", "Bacterial", "Nematode"]
COL = {"Human": "#1f4e79", "Fungal": "#c0504d",
       "Bacterial": "#4f8a3d", "Nematode": "#8064a2"}

prot = pd.read_csv(TABLES_DIR / "disprot_protein_level.csv")
pairs = pd.read_csv(TABLES_DIR / "thales_curated_sample_results.csv")
conc = pd.read_csv(TABLES_DIR / "cross_taxon_calibration_summary.csv")
rob = pd.read_csv(DIAG_DIR / "robustness_sensitivity.csv")
npz = np.load(DIAG_DIR / "asymmetry_profiles.npz")
grid = npz["grid"]

# ---- Fig 1: sample composition and AlphaFold availability -------------
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
c = prot.groupby("group").size().reindex(GROUPS)
a = prot.groupby("group")["af_available"].sum().reindex(GROUPS)
ax[0].bar(GROUPS, c, color=[COL[g] for g in GROUPS], alpha=.45,
          label="DisProt proteins")
ax[0].bar(GROUPS, a, color=[COL[g] for g in GROUPS],
          label="with AlphaFold value")
for i, (t, v) in enumerate(zip(c, a)):
    ax[0].text(i, t, f"{t}\n({t - v} missing)", ha="center", va="bottom",
               fontsize=7)
ax[0].set_yscale("log"); ax[0].set_ylabel("proteins (log)")
ax[0].legend(fontsize=7); ax[0].set_title("a  Curated sample size")
ax[0].set_ylim(top=c.max() * 3)

for g in GROUPS:
    v = np.sort(prot.loc[prot.group == g, "length"])
    ax[1].plot(v, np.linspace(0, 1, len(v)), color=COL[g], label=g)
ax[1].set_xscale("log"); ax[1].set_xlabel("protein length (residues)")
ax[1].set_ylabel("ECDF"); ax[1].legend(fontsize=7)
ax[1].set_title("b  Length distribution")
fig.savefig(FIG / "fig1_composition.pdf")
fig.savefig(FIG / "fig1_composition.png"); plt.close(fig)

# ---- Fig 2: feature distributions -------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9), sharey=True)
for j, (feat, lab) in enumerate([
        ("D_exp", "curated DisProt disorder fraction"),
        ("af_very_low_content", "AlphaFold very-low-confidence fraction")]):
    for g in GROUPS:
        v = np.sort(prot.loc[prot.group == g, feat].dropna())
        ax[j].plot(v, np.linspace(0, 1, len(v)), color=COL[g],
                   label=f"{g} (n={len(v)})")
    ax[j].set_xlabel(lab); ax[j].set_xlim(0, 1)
    ax[j].legend(fontsize=7, loc="lower right")
ax[0].set_ylabel("ECDF")
ax[0].set_title("a  Experimental"); ax[1].set_title("b  AlphaFold proxy")
fig.savefig(FIG / "fig2_distributions.pdf")
fig.savefig(FIG / "fig2_distributions.png"); plt.close(fig)

# ---- Fig 3: asymmetry profiles ---------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True)
for j, feat in enumerate(["D_exp", "af_very_low_content"]):
    for gA, gB in itertools.combinations(GROUPS, 2):
        k = f"{feat}|{gA}|{gB}"
        if k not in npz:
            continue
        ax[j].plot(grid, npz[k], lw=1.2, label=f"{gA}\u2212{gB}")
    ax[j].axhline(0, color="k", lw=.7)
    ax[j].set_xlabel("feature value $x$"); ax[j].set_ylim(-1, 1)
    ax[j].legend(fontsize=6, ncol=2)
ax[0].set_ylabel("$S(x)$")
ax[0].set_title("a  Curated disorder fraction")
ax[1].set_title("b  AlphaFold very-low fraction")
fig.savefig(FIG / "fig3_asymmetry_profiles.pdf")
fig.savefig(FIG / "fig3_asymmetry_profiles.png"); plt.close(fig)

# ---- Fig 4: bootstrap forest of B_bar --------------------------------
fig, ax = plt.subplots(figsize=(6.2, 3.6))
p = pairs.copy()
p["lab"] = (p.feature.map({"D_exp": "exp", "af_very_low_content": "AF"})
            + ": " + p.group_A + "\u2212" + p.group_B)
p = p.sort_values(["feature", "B_bar"])
y = np.arange(len(p))
ax.errorbar(p.B_bar, y,
            xerr=[p.B_bar - p.B_bar_lo, p.B_bar_hi - p.B_bar],
            fmt="o", ms=4, lw=1, capsize=2,
            color="#1f4e79")
ax.axvline(0, color="k", lw=.8)
ax.set_yticks(y); ax.set_yticklabels(p.lab, fontsize=7)
ax.set_xlabel(r"$\overline{\mathcal{B}}$  (95% protein-level bootstrap CI)")
ax.set_xlim(-1, 1)
fig.savefig(FIG / "fig4_bootstrap_forest.pdf")
fig.savefig(FIG / "fig4_bootstrap_forest.png"); plt.close(fig)

# ---- Fig 5: AF vs experimental concordance ---------------------------
fig, ax = plt.subplots(1, 4, figsize=(7.6, 2.3), sharex=True, sharey=True)
pr = prot.dropna(subset=["D_exp", "af_very_low_content"])
for i, g in enumerate(GROUPS):
    s = pr[pr.group == g]
    ax[i].scatter(s.D_exp, s.af_very_low_content, s=5, alpha=.4,
                  color=COL[g], lw=0)
    ax[i].plot([0, 1], [0, 1], "k--", lw=.7)
    row = conc[conc.group == g].iloc[0]
    ax[i].set_title(f"{g}\n"
                    r"$\rho$=" f"{row.spearman_rho:.2f} "
                    f"[{row.rho_lo:.2f},{row.rho_hi:.2f}]", fontsize=7.5)
    ax[i].set_xlim(0, 1); ax[i].set_ylim(0, 1)
    ax[i].set_xlabel("curated $D$", fontsize=8)
ax[0].set_ylabel("AlphaFold\nvery-low fraction", fontsize=8)
fig.savefig(FIG / "fig5_concordance.pdf")
fig.savefig(FIG / "fig5_concordance.png"); plt.close(fig)

# ---- Fig 6: robustness -----------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0))
r = rob[rob.feature == "D_exp"].copy()
r["spec"] = np.where(r.estimator == "kde",
                     "KDE bw=" + r.bw_scale.astype(str)
                     + " eps=" + r.eps.astype(str),
                     "hist " + r.nbins.astype(str) + " bins")
order = r.spec.unique()
for gA, gB in itertools.combinations(GROUPS, 2):
    s = r[(r.group_A == gA) & (r.group_B == gB)].set_index("spec").loc[order]
    ax[0].plot(range(len(order)), s.B_bar, "o-", ms=3, lw=1,
               label=f"{gA}\u2212{gB}")
    ax[1].plot(range(len(order)), s.A_bar, "o-", ms=3, lw=1)
ax[0].axhline(0, color="k", lw=.8)
for a_ in ax:
    a_.set_xticks(range(len(order)))
    a_.set_xticklabels(order, rotation=55, ha="right", fontsize=6)
ax[0].set_ylabel(r"$\overline{\mathcal{B}}$")
ax[1].set_ylabel(r"$\overline{\mathcal{A}}$")
ax[0].legend(fontsize=6)
ax[0].set_title("a  Signed balance", fontsize=9)
ax[1].set_title("b  Asymmetry magnitude", fontsize=9)
fig.savefig(FIG / "fig6_robustness.pdf")
fig.savefig(FIG / "fig6_robustness.png"); plt.close(fig)

print("figures written:", sorted(p.name for p in FIG.iterdir()))
