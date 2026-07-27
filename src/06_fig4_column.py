#!/usr/bin/env python3
"""Regenerate fig4 at native two-column width (3.26in) with legible fonts.

Journals that typeset in two columns scale a 6.2in figure down by roughly half,
which renders the 7pt tick labels of fig4_bootstrap_forest illegible. This
script re-emits the same forest plot at final printed width so no downscaling
occurs.

Portability note (2026-07-27 refactor): absolute paths replaced with
project-relative pathlib paths plus DISPROT_OUT override. Figure size, fonts,
DPI and plotted values are unchanged.
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ------------------------------------------------------------------- paths
ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = Path(os.environ.get("DISPROT_OUT", ROOT / "results"))
TABLES_DIR = RESULTS_DIR / "tables"
FIG = RESULTS_DIR / "figures"
FIG.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({"font.size": 7, "axes.grid": True, "grid.alpha": .25,
                     "figure.dpi": 300, "savefig.bbox": "tight"})
p = pd.read_csv(TABLES_DIR / "thales_curated_sample_results.csv")
p["lab"] = (p.feature.map({"D_exp": "$D$", "af_very_low_content": "$V$"}) + ": "
            + p.group_A.str[:4] + "\u2212" + p.group_B.str[:4])
p = p.sort_values(["feature", "B_bar"])
fig, ax = plt.subplots(figsize=(3.26, 3.1))
y = np.arange(len(p))
ax.errorbar(p.B_bar, y, xerr=[p.B_bar - p.B_bar_lo, p.B_bar_hi - p.B_bar],
            fmt="o", ms=3, lw=.9, capsize=1.8, color="#1f4e79")
ax.axvline(0, color="k", lw=.7)
ax.set_yticks(y); ax.set_yticklabels(p.lab, fontsize=6)
ax.set_xlabel(r"$\overline{\mathcal{B}}$ (95% bootstrap CI)", fontsize=7)
ax.set_xlim(-1, 1); ax.tick_params(labelsize=6)
for ext in ("pdf", "png"):
    fig.savefig(FIG / f"fig4_bootstrap_forest_1col.{ext}")
print("written: fig4_bootstrap_forest_1col.pdf / .png")
