# Bounded asymmetry profiles for protein disorder distributions across taxa

Reproducibility package for the methodological report *"Bounded asymmetry
profiles for comparing protein disorder distributions across taxa: a
methodological report with a negative concordance result for an
AlphaFold-derived proxy"* (E. De Jesús, 2026).

One command regenerates every table, figure and diagnostic in the manuscript
from the four deposited input files.

```bash
python run_all.py
```

---

## Project overview

Two-sample tests tell you *whether* two populations differ but not *where*
along the feature axis they differ. This study defines a bounded local asymmetry
profile over a continuous feature coordinate `x`,

```
S(x) = ( p_A(x) − p_B(x) ) / ( p_A(x) + p_B(x) + 2ε )
```

which is bounded in (−1, 1), antisymmetric under exchange of populations, and
zero exactly where the two density estimates agree. Three scalar summaries are
derived from it: the asymmetry magnitude `Ā`, the signed balance `B̄`, and the
directional consistency `C = |B̄|/Ā`.

The profile is evaluated on experimentally curated intrinsic disorder
annotations for **1,766 proteins** across four taxa (DisProt release 2026_06):
1,339 human, 223 *Saccharomyces cerevisiae*, 146 *Escherichia coli* K12, and 58
*Caenorhabditis elegans*.

The report has three findings, two of them negative:

1. **The profile adds information beyond global tests.** For human vs.
   *E. coli* curated disorder it resolves localised structure (`Ā = 0.143`,
   `C = 0.57`, two zero crossings) in a comparison that global tests read as
   null (KS *p* = 0.39, Mann–Whitney *p* = 0.85).
2. **The AlphaFold-derived very-low-confidence fraction is a poor
   protein-level proxy for curated disorder** in every taxon examined
   (Spearman ρ = 0.07–0.24), and its offset relative to curated disorder is
   taxon-dependent (median signed gap +0.001 human vs. −0.083 *E. coli*,
   Mann–Whitney *p* < 10⁻⁴). This is the study's gating result: it is reported
   even though it blocks the intended downstream application.
3. **The scalar magnitude index `Ā` is estimator-dependent.** It inflates from
   0.181 to 0.333 with histogram bin count because empty bins drive |S| → 1,
   and the signed index changes sign in 5 of 12 comparisons under histogram
   estimation while remaining stable across every kernel bandwidth and
   regularisation constant tested. Kernel estimation with a reported bandwidth
   is therefore required for the framework to be reproducible.

No network retrieval succeeded during the original analysis. Every external
database required by the protocol (AlphaFold DB, UniProt REST, MobiDB, the
DisProt API) returned HTTP 403 at the egress proxy. Those failures are logged
verbatim in `results/diagnostics/alphafold_retrieval_manifest.csv` rather than
being silently omitted, and they bound what the study can claim — see
§3.3 and §3.5 of the manuscript.

---

## Repository layout

```
.
├── README.md                   this file
├── CITATION.cff                citation metadata (software + preferred article citation)
├── LICENSE                     MIT — applies to src/ and run_all.py only
├── DATA_LICENSE.md             licensing for data and manuscript; DisProt attribution
├── requirements.txt            bounded runtime dependencies
├── requirements-lock.txt       exact versions used for the verification run
├── .gitignore
├── run_all.py                  runs the whole pipeline in dependency order
│
├── data/
│   ├── README.md               data provenance, taxonomy, known data-quality caveats
│   ├── raw/                    the four DisProt JSON exports, unmodified
│   │   ├── disprot_taxon_9606_homo_sapiens.json
│   │   ├── disprot_taxon_559292_saccharomyces_cerevisiae.json
│   │   ├── disprot_taxon_83333_escherichia_coli_k12.json
│   │   └── disprot_taxon_6239_caenorhabditis_elegans.json
│   └── metadata/
│       └── source_queries.csv  per-file provenance: taxon, counts, MD5, SHA-256
│
├── src/
│   ├── 01_build_tables.py      audit inputs; build protein- and region-level tables
│   ├── 02_asymmetry.py         profiles, scalar indices, bootstrap, robustness sweep
│   ├── 03_figures.py           Figures 1–6
│   ├── 04_retrieval_log.py     retrieval-failure manifest; unavailable accessions
│   ├── 05_dphi_test.py         composite-operator diagnostic (not a manuscript result)
│   └── 06_fig4_column.py       Figure 4 re-emitted at single-column width
│
├── results/
│   ├── README.md               output-to-manuscript map; verification record
│   ├── tables/                 protein/region tables and manuscript table sources
│   ├── diagnostics/            audit, robustness, calibration and retrieval CSVs
│   ├── figures/                Figures 1–6 as PDF and PNG
│   └── deposited_originals/    archival copies from the earlier deposit — DO NOT EDIT
│
└── paper/
    └── manuscript.pdf
```

Two notes on this layout:

- **`results/deposited_originals/` is not pipeline output.** It holds three
  files from the earlier archival deposit that the archived scripts cannot
  reproduce, preserved unchanged as historical record. `results/README.md`
  documents the differences in full. Nothing in the pipeline writes to that
  directory.
- `alphafold_retrieval_manifest.csv` lives in `results/diagnostics/` rather
  than `data/metadata/`, because it is generated by `04_retrieval_log.py` on
  every run. Keeping a second copy under `data/` would let the two drift apart.

---

## Installation

```bash
git clone https://github.com/innerlightr-wq/disprot-bounded-asymmetry.git
cd disprot-bounded-asymmetry

python -m venv venv
# macOS / Linux
source venv/bin/activate
# Windows (PowerShell)
venv\Scripts\Activate.ps1

pip install -r requirements.txt
```

To reproduce the verified run exactly, use the pinned versions instead:

```bash
pip install -r requirements-lock.txt
```

**NumPy 2.x is required.** `src/02_asymmetry.py` calls `numpy.trapezoid`, added
in NumPy 2.0. The scientific code is deliberately not rewritten for NumPy 1.x.

## Running

```bash
python run_all.py
```

Expected output:

```
Running 01_build_tables.py...
Running 02_asymmetry.py...
Running 03_figures.py...
Running 04_retrieval_log.py...
Running 05_dphi_test.py...
Running 06_fig4_column.py...
Pipeline completed successfully.
```

Runtime is a few minutes; the dominant cost is 2,000-replicate protein-level
bootstraps over 12 comparisons in `02_asymmetry.py`.

Scripts may also be run individually from anywhere, but **order matters**: `02`
reads what `01` writes, and `03` and `06` read what `02` writes. `04` needs only
`01`. Paths are resolved relative to the repository root via
`Path(__file__).resolve().parents[1]`, so the working directory does not matter.

### Environment overrides

Both variables are honoured by every script. Defaults are used when unset.

| Variable | Default | Purpose |
|---|---|---|
| `DISPROT_DATA` | `<repo>/data/raw` | directory holding the four JSON exports |
| `DISPROT_OUT` | `<repo>/results` | directory to write results into |

```bash
# macOS / Linux
DISPROT_OUT=/tmp/run1 python run_all.py

# Windows (PowerShell)
$env:DISPROT_OUT="C:\temp\run1"; python run_all.py
```

Output directories are created automatically. There are no absolute paths and
no POSIX assumptions anywhere in `src/`; all path handling is `pathlib`.

---

## Data provenance

Full detail is in [`data/README.md`](data/README.md). In summary:

- **Source.** Four per-organism JSON exports from **DisProt release 2026_06**,
  "with ambiguous evidences" variant. Included unmodified; MD5 and SHA-256 for
  each are in `data/metadata/source_queries.csv` and reproduced on every run in
  `results/diagnostics/source_file_audit.csv`.
- **Taxonomy.** Files are named by the **NCBI taxonomy ID actually carried by
  the records**: 9606 (*Homo sapiens*), 559292 (*S. cerevisiae* S288c), 83333
  (*E. coli* K12), 6239 (*C. elegans*). Two of these are strain-level
  identifiers, not the kingdom-level IDs 4751 (Fungi) and 2 (Bacteria) — using
  the kingdom IDs in the filenames would misstate what the files contain.
- **Group assignment.** Group labels are derived from `ncbi_taxon_id` inside
  each record, **never** from the filename, and `01_build_tables.py` raises if a
  file's records disagree with the taxon its name claims. The four exports are
  mutually disjoint: 1,766 records, 1,766 unique UniProt accessions, 1,766
  unique DisProt identifiers.
- **AlphaFold retrieval.** No residue-level pLDDT was ever retrieved. The
  `alphafold_very_low_content` field is taken as deposited in the DisProt
  export and is **unverified** against AlphaFold DB. Every attempt and its
  failure reason is in `results/diagnostics/alphafold_retrieval_manifest.csv`;
  the 103 proteins with no value are listed in
  `results/diagnostics/failed_accessions.csv`.
- **Mixed data types.** 34 `alphafold_very_low_content` values are stored as
  **strings** rather than floats in the export. `01_build_tables.py` coerces
  with `pd.to_numeric(..., errors="coerce")`. Reading this field without
  coercion will silently produce an object-dtype column and wrong statistics.

## Reproducibility

Every manuscript table, figure and diagnostic is regenerated by `run_all.py`
from `data/raw`; nothing is hand-edited or carried forward from a previous run.
Determinism comes from three things:

1. A fixed RNG seed (`np.random.default_rng(20260726)`) consumed in a fixed
   loop order, so the 2,000-replicate bootstrap CIs are reproducible to the last
   digit.
2. An explicit ordered input manifest in `src/01_build_tables.py` and
   `src/05_dphi_test.py`. The original scripts discovered inputs with
   `sorted(glob.glob("*.json"))`, which made table row order depend on
   arbitrary browser-assigned download suffixes; that is why the archived
   deposit's tables are ordered Bacterial-first. Row order is now fixed to the
   group order of manuscript Table 1.
3. Pinned versions in `requirements-lock.txt`.

The refactor that produced this repository was verified against the
pre-refactor run: **all 13 CSV outputs and all 13 stored profile arrays are
value-identical, and all 7 figure PNGs are byte-identical.** The full record,
including the three archived files the scripts cannot reproduce, is in
[`results/README.md`](results/README.md).

## Results map

| Manuscript item | Regenerated file(s) | Written by |
|---|---|---|
| Figure 1 | `results/figures/fig1_composition.{pdf,png}` | `03_figures.py` |
| Figure 2 | `results/figures/fig3_asymmetry_profiles.{pdf,png}` | `03_figures.py` |
| Figure 3 | `results/figures/fig4_bootstrap_forest.{pdf,png}` and `fig4_bootstrap_forest_1col.{pdf,png}` | `03_figures.py`, `06_fig4_column.py` |
| Figure 4 | `results/figures/fig5_concordance.{pdf,png}` | `03_figures.py` |
| Figure 5 | `results/figures/fig6_robustness.{pdf,png}` | `03_figures.py` |
| Table 1 | `results/tables/group_feature_summary.csv`, `results/tables/disprot_protein_level.csv` | `02_asymmetry.py`, `01_build_tables.py` |
| Table 2 | `results/tables/thales_curated_sample_results.csv` | `02_asymmetry.py` |
| §3.2 concordance (ρ, gaps) | `results/tables/cross_taxon_calibration_summary.csv`, `results/diagnostics/gap_by_taxon_tests.csv` | `02_asymmetry.py` |
| §3.4 estimator sensitivity | `results/diagnostics/robustness_sensitivity.csv`, `results/diagnostics/sign_stability.csv` | `02_asymmetry.py` |
| §2.3 audit / integrity claims | `results/diagnostics/source_file_audit.csv`, `results/diagnostics/phase1_integrity_checks.csv` | `01_build_tables.py` |
| §2.4 retrieval failures | `results/diagnostics/alphafold_retrieval_manifest.csv`, `results/diagnostics/failed_accessions.csv` | `04_retrieval_log.py` |

**Note on figure numbering.** Script filenames do not match manuscript figure
numbers. The manuscript has five figures in its body; the scripts emit six
(`fig1`–`fig6`), because `fig2_distributions` — the ECDF panel of feature
distributions — is not among them. `results/README.md` gives the full mapping.

## Environment

Tested on **Python 3.12.3**, Linux x86_64 (glibc 2.39), with NumPy 2.4.4,
pandas 3.0.2, SciPy 1.17.1 and matplotlib 3.10.8. Python 3.10 or newer should
work; the pipeline is pure Python plus the four scientific packages and has no
compiled extensions or OS-specific calls. matplotlib is used with the `Agg`
backend, so no display is required.

## Licensing

Three different sets of terms apply — see [`DATA_LICENSE.md`](DATA_LICENSE.md).
In short: the code is MIT; the DisProt input data in `data/raw/` is **CC BY 4.0**
(verified at https://disprot.org/about on 2026-07-27, which permits
redistribution with attribution); the manuscript is covered by neither.

## Related paper

The manuscript is in [`paper/manuscript.pdf`](paper/manuscript.pdf).

Preprint archived on Zenodo: https://doi.org/10.5281/zenodo.21628406
(version v8, published 2026-07-27, CC BY 4.0)

*TODO: replace the line above with the Zenodo **concept DOI**, which always
resolves to the latest version. Read it off the record page under
"Versions -> Cite all versions?". The version DOI above will go stale the next
time a new version is uploaded.*

**The Zenodo record and this repository are not the same artifact.** The
Zenodo deposit carries the preprint PDF together with the original,
pre-refactor scripts and their outputs. This repository carries the refactored,
portable pipeline. The two produce identical numbers -- that is verified in
`results/README.md` -- but the files differ, and the Zenodo copies of
`source_file_audit.csv`, `sign_stability.csv` and `gap_by_taxon_tests.csv` are
the richer variants that the archived scripts cannot regenerate (see
`results/deposited_originals/`).

The manuscript is not yet accepted for publication, so this repository carries
no journal name, volume, issue, pages or article DOI. `CITATION.cff` records
only what is presently true — title, author, version, release date and licence
— and journal metadata should be added there only after formal acceptance.

Note that `paper/manuscript.pdf` contains its own Data Availability Statement
and "How to Cite" block with their own identifiers. Those are part of the
manuscript and are deliberately left untouched by this repository; if they
disagree with anything here, the manuscript is the item that needs updating.

## Contact

Elias De Jesús — Independent Researcher, Virginia, USA
ORCID [0009-0007-0190-9143](https://orcid.org/0009-0007-0190-9143)
