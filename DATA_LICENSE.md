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

---

## AlphaFold validation data (added 2026-09-18)

The September 2026 validation audit downloaded per-residue confidence files from the AlphaFold
Protein Structure Database (`https://alphafold.ebi.ac.uk/files/AF-<acc>-F1-confidence_v6.json`) for
1,661 accessions, to verify what DisProt's `alphafold_very_low_content` field contains.

**Those files are not redistributed in this repository.** What is committed is:

* `data/alphafold_validation_manifest.csv` — one row per accession: source URL, retrieval timestamp,
  model version, sequence and model lengths, the **SHA-256 of the downloaded file**, and mapping
  status;
* `results/alphafold_validation/proxy_validation.csv` — the derived per-protein scalars (fractions
  below the confidence thresholds, mean and median pLDDT, coverage);
* `analysis/validate_alphafold_proxy.py` — the script that regenerates both from the accession list.

That combination reproduces the analysis without republishing third-party structure data. The
downloads are cached outside the repository (`--cache-dir`, default `/tmp/af_cache`, about 16 MB for
the full set).

This choice is deliberately conservative and does not depend on how AlphaFold DB is licensed. **If
you intend to redistribute the confidence files themselves, check the current licence terms on the
AlphaFold DB site first** — they were not verified here, because nothing in this repository requires
redistribution. Note also that AlphaFold DB serves only the current model version (`v6` at the time
of writing; `v1`–`v5` return HTTP 404), so a reconstruction run later may differ slightly from the
committed values; the manifest records the version and hash actually used.
