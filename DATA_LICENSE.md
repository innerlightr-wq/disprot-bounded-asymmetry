# Licensing and attribution

This repository contains three kinds of material under three different sets of
terms. **The software licence does not cover the data or the manuscript.**

All three sets of terms are now confirmed; nothing in this file is provisional.

| Component | Paths | Terms |
|---|---|---|
| Analysis software | `src/`, `run_all.py` | MIT — see [`LICENSE`](LICENSE) |
| Derived results | `results/tables/`, `results/diagnostics/`, `results/figures/` | MIT, as software output; the underlying facts remain attributable to DisProt (see below) |
| Input data | `data/raw/*.json` | **CC BY 4.0** (DisProt) — attribution required; see §3 |
| Manuscript | `paper/manuscript.pdf` | © 2026 Elias De Jesús. All rights reserved unless superseded by the publishing agreement of the journal in which it appears. Not covered by the MIT licence. |

---

## 1. Software — MIT

The six analysis scripts in `src/` and the `run_all.py` orchestrator are
released under the MIT License, reproduced in full in [`LICENSE`](LICENSE).

## 2. Manuscript

`paper/manuscript.pdf` is a research manuscript, not software. It is **not**
covered by the MIT licence. Copyright is retained by the author, subject to
whatever terms the eventual publication agreement imposes. If the article is
published under a Creative Commons licence by the journal, that licence governs
the published version; this file should then be updated to name it.

## 3. Input data — DisProt

`data/raw/` contains four unmodified per-organism JSON exports from **DisProt
release 2026_06** ("with ambiguous evidences" variant). These are third-party
database records. They are included here so the pipeline is self-contained and
independently re-runnable, and they are byte-identical to the files as
downloaded (MD5 checksums are recorded in
`results/diagnostics/source_file_audit.csv`).

### Licence: CC BY 4.0 — confirmed

DisProt is distributed under the **Creative Commons Attribution 4.0
International License** (CC BY 4.0), https://creativecommons.org/licenses/by/4.0/

- Source of this statement: https://disprot.org/about ("License" section)
- Verified: 2026-07-27

CC BY 4.0 **permits redistribution**, including in a repository such as this
one, provided attribution is given. The four exports in `data/raw/` are
therefore included lawfully, unmodified, with the attribution below and with
checksums recorded so their integrity is verifiable.

Note that CC BY 4.0 applies to the DisProt data. It does **not** extend to this
repository's software (MIT, see §1) or to the manuscript (§2), and this
repository's author is not the rights holder for the DisProt data.

### Attribution

Any use of these data, or of results derived from them, should cite DisProt.
Primary citation for the release used here (2026_06):

> Nugnes, M. V., Bouhraoua, K. E. A., Zoubiri, M., Pancsa, R., Fichó, E.,
> DisProt Consortium, Tompa, P., Piovesan, D., Tosatto, S. C. E., &
> Aspromonte, M. C. (2026). DisProt in 2026: enhancing intrinsically disordered
> proteins accessibility, deposition, and annotation. *Nucleic Acids Research*,
> 54(D1), D383–D392. https://doi.org/10.1093/nar/gkaf1175
> (PubMed: 41249866)

Curation methodology, cited by DisProt alongside the release paper:

> Quaglia, F., Chasapi, A., Nugnes, M. V., Aspromonte, M. C., Leonardi, E.,
> Piovesan, D., & Tosatto, S. C. E. (2024). Best practices for the manual
> curation of intrinsically disordered proteins in DisProt. *Database*,
> baae009. https://doi.org/10.1093/database/baae009
> (PubMed: 38507044)

The previous release paper remains appropriate where earlier versions are
discussed:

> Aspromonte, M. C., Nugnes, M. V., Quaglia, F., Bouharoua, A., DisProt
> Consortium, Tosatto, S. C. E., & Piovesan, D. (2024). DisProt in 2024:
> improving function annotation of intrinsically disordered proteins.
> *Nucleic Acids Research*, 52(D1), D434–D441.
> https://doi.org/10.1093/nar/gkad928

The records additionally carry a precomputed `alphafold_very_low_content` field
derived from AlphaFold DB. Work relying on that field should also cite Jumper
et al. (2021), https://doi.org/10.1038/s41586-021-03819-2, and Varadi et al.
(2024), https://doi.org/10.1093/nar/gkad1011.

## 4. Third-party software

The pipeline depends on NumPy, pandas, SciPy and matplotlib, each distributed
under its own permissive licence (BSD-3-Clause or equivalent). Those licences
apply to those packages, not to this repository.
