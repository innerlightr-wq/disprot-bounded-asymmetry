#!/usr/bin/env python3
"""Good-faith test: does a DeltaPhi-style convex triadic operator add
information beyond the component features, on the DisProt sample?

Instantiation (three normalized deviations in [0,1] per protein):
  dS = curated disorder fraction D_exp            (structural)
  dI = 1 - normalized AA-composition entropy      (informational)
  dC = AlphaFold very-low-confidence fraction     (coherence)
  DeltaPhi = a*dS + b*dI + g*dC,  a+b+g=1

This is a diagnostic reported in the repository record, not a manuscript
result: the conclusion was that the composite does not add value (taxon
ordering is weight-dependent, no regime boundary appears at the claimed
thresholds, and Human/Bacterial discrimination degrades from AUC 0.837 for the
best single component to 0.692 under the composite).

Portability note (2026-07-27 refactor): absolute paths replaced with
project-relative pathlib paths plus DISPROT_DATA / DISPROT_OUT overrides, and
input discovery by glob replaced with the same explicit ordered manifest used
by 01_build_tables.py. No computation was changed.
"""
import itertools
import json
import os
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# ------------------------------------------------------------------- paths
ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = Path(os.environ.get("DISPROT_DATA", ROOT / "data" / "raw"))
RESULTS_DIR = Path(os.environ.get("DISPROT_OUT", ROOT / "results"))
DIAG_DIR = RESULTS_DIR / "diagnostics"
DIAG_DIR.mkdir(parents=True, exist_ok=True)

SOURCE_FILES = [
    "disprot_taxon_9606_homo_sapiens.json",
    "disprot_taxon_559292_saccharomyces_cerevisiae.json",
    "disprot_taxon_83333_escherichia_coli_k12.json",
    "disprot_taxon_6239_caenorhabditis_elegans.json",
]

rows = []
for fname in SOURCE_FILES:
    path = DATA_DIR / fname
    if not path.exists():
        raise FileNotFoundError(
            f"Required input not found: {path}\n"
            f"Place the four DisProt JSON exports in {DATA_DIR}, or point the "
            f"DISPROT_DATA environment variable at the directory holding them."
        )
    with open(path) as fh:
        for r in json.load(fh)['data']:
            seq = r.get('sequence') or ''
            if not seq:
                continue
            c = Counter(seq); n = sum(c.values())
            p = np.array([v / n for v in c.values()])
            H = -(p * np.log(p)).sum() / np.log(20)      # normalized to [0,1]
            L = r['length']
            cov = set()
            for s in r['disprot_consensus'].get('Structural state', []):
                if s.get('type') == 'D':
                    cov.update(range(s['start'], s['end'] + 1))
            af = r.get('alphafold_very_low_content')
            rows.append(dict(acc=r['acc'], tax=r['ncbi_taxon_id'],
                             dS=len(cov) / L, dI=min(max(1 - H, 0), 1),
                             dC=af if af is not None else np.nan))
d = pd.DataFrame(rows)
d['group'] = d.tax.map({9606: 'Human', 559292: 'Fungal',
                        83333: 'Bacterial', 6239: 'Nematode'})
d['dC'] = pd.to_numeric(d['dC'], errors='coerce')
d = d.dropna(subset=['dC'])
print('n complete triads:', len(d))
print(d.groupby('group')[['dS', 'dI', 'dC']].median().round(4).to_string())

# component correlations: does the composite average redundant or unrelated axes?
print('\nSpearman between the three axes (all proteins):')
for a, b in itertools.combinations(['dS', 'dI', 'dC'], 2):
    rho, p = stats.spearmanr(d[a], d[b])
    print(f'  {a}-{b}: rho={rho:+.3f} p={p:.2g}')

# 1) is anything special at 1/3 or 0.40 under symmetric weights?
d['dphi'] = (d.dS + d.dI + d.dC) / 3
print('\nSymmetric-weight DeltaPhi distribution:')
print('  median %.4f  IQR %.4f-%.4f  min %.4f  max %.4f' % (
    d.dphi.median(), d.dphi.quantile(.25), d.dphi.quantile(.75),
    d.dphi.min(), d.dphi.max()))
print('  frac above 1/3: %.4f   frac above 0.40: %.4f' % (
    (d.dphi > 1 / 3).mean(), (d.dphi > 0.40).mean()))
# test for a density dip (regime boundary) near candidate thresholds
kde = stats.gaussian_kde(d.dphi)
gr = np.linspace(d.dphi.min(), d.dphi.max(), 400); dens = kde(gr)
loc = gr[np.r_[False, (dens[1:-1] < dens[:-2]) & (dens[1:-1] < dens[2:]), False]]
print('  interior density minima (candidate regime boundaries):',
      np.round(loc, 3) if len(loc) else 'NONE -> unimodal, no boundary')

# 2) does the group ranking depend on the weights?
print('\nGroup ordering by median DeltaPhi under different convex weights:')
for w in [(1 / 3, 1 / 3, 1 / 3), (.8, .1, .1), (.1, .8, .1), (.1, .1, .8),
          (.5, .5, 0), (0, .5, .5)]:
    a, b, g = w
    v = a * d.dS + b * d.dI + g * d.dC
    m = v.groupby(d.group).median().sort_values(ascending=False)
    print('  w=(%.2f,%.2f,%.2f) -> ' % w
          + ' > '.join(f'{k}({m[k]:.3f})' for k in m.index))

# 3) does DeltaPhi separate taxa better than its best single component?
print('\nHuman vs Bacterial separation (AUC = |rank-biserial| mapped to 0.5-1):')
h = d[d.group == 'Human']; bc = d[d.group == 'Bacterial']
for name in ['dS', 'dI', 'dC', 'dphi']:
    U = stats.mannwhitneyu(h[name], bc[name]).statistic
    auc = U / (len(h) * len(bc)); auc = max(auc, 1 - auc)
    print(f'  {name}: AUC={auc:.3f}')
d.to_csv(DIAG_DIR / 'dphi_operator_test.csv', index=False)
