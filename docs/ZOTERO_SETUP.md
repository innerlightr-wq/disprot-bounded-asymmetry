# Zotero organisation for this project

Companion to [`NOVELTY_AND_PROVENANCE.md`](NOVELTY_AND_PROVENANCE.md). The authoritative
bibliography is [`references.bib`](../references.bib) (38 entries); Zotero is a convenience layer over
it, and nothing in the repository depends on Zotero being present.

## 1. Status

Applied 2026-09-18. The library now contains what this document specifies.

| Item | State |
|---|---|
| `references.bib` | **complete** — 38 entries, every DOI resolved against Crossref/DataCite, entries generated from registrar metadata rather than retyped |
| Collection tree | **created** — parent `DisProt — Bounded Asymmetry` plus all 17 subcollections |
| Items | **38 imported**, 0 failed, **0 duplicates**, each carrying its BibTeX citation key |
| Filing | **104 memberships** across the 38 items; all 17 subcollection counts match the specification exactly |
| Tags | **13 controlled tags applied**, counts matching the specification exactly; no tag outside the controlled list |
| Notes | **9 child notes**, one per prioritised source, all six fields present in each |
| Metadata agreement | **0 discrepancies in 38 items** against `references.bib` (citation key, title, year, DOI, first author, author count, journal, volume, pages) |

Library after the import: **180 top-level items, 106 collections, empty trash, 37 distinct tags.**
The six other project trees (`EOC — Collatz Conjecture Research`,
`OU Threshold — Heterogeneity and Coupling`, `Pell Spine — Optimal Separator`,
`Signed Context Decomposition — Cosmology`, `Strain–Vorticity Comparator Dynamics`,
`Structural Addresses Methodology`) were verified unchanged, with their subcollection counts intact.

Two notes on how this was done, for anyone repeating it:

* Writes to the Zotero local API need both an API key from `POST /api/local/authorize` — which waits
  on a confirmation dialog in the desktop application — and a `Zotero-Server-ID` header, which is
  returned on any read. Reads need neither. The key expires, and re-authorisation needs the dialog
  again.
* The import ran twice by accident, because the module that performs the writes was later imported
  for its BibTeX parser and had no `__main__` guard, so its top-level code executed a second time.
  That produced exactly one duplicate of each of the 38 items. The duplicates were identified by
  citation key, the earlier copy of each pair kept, and the 38 later copies deleted; the counts and
  the metadata check above were then re-verified from the live library, and the trash is empty. If
  you reuse the tooling, keep the write path behind a `__main__` guard.

Three of the 38 items are not journal articles: `Scott1992`, `SugiyamaEtAl2012` and `LeCam1986` are
books, so their DOIs live in the Extra field, which Zotero has no DOI field for on that item type.
`DeJesus2026disprot` is filed as a preprint.

## 2. Collection tree

Parent: **`DisProt — Bounded Asymmetry`**, with seventeen subcollections. The count in the last column
is how many of the 38 references file into each; **no subcollection is empty**, which is why this
list has no placeholders.

| Subcollection | Items |
|---|---|
| `00 — Reviews & Orientation` | 2 |
| `01 — DisProt & Intrinsic Disorder` | 4 |
| `02 — AlphaFold Confidence & pLDDT` | 6 |
| `03 — Disorder Prediction Benchmarks` | 6 |
| `04 — CAID & MobiDB` | 4 |
| `05 — Density Estimation / KDE` | 3 |
| `06 — Two-Sample Testing` | 4 |
| `07 — Density Difference & Relative Density Ratio` | 5 |
| `08 — Total Variation & Divergences` | 7 |
| `09 — Local Two-Sample / Witness Methods` | 2 |
| `10 — Permutation Calibration` | 2 |
| `11 — Unequal Sample Size / Bandwidth Effects` | 3 |
| `12 — Taxonomic Disorder Variation` | 3 |
| `13 — Proxy Validation & Calibration` | 3 |
| `14 — Negative Results / Method Evaluation` | 2 |
| `15 — Directly Cited in Repository` | **38** |
| `16 — Closest Prior Art / Novelty Checks` | 10 |

104 memberships across 38 distinct items: an item appears in several subcollections because the
subcollections record the *role* a source plays in the argument, not a partition of the bibliography.

`15 — Directly Cited in Repository` holds all 38 and is kept in exact agreement with
`references.bib` — that is its only job.

`16 — Closest Prior Art / Novelty Checks` is deliberately small and covers exactly the six
novelty questions the audit had to answer:

| Question | Sources |
|---|---|
| exact normalised density contrast | `LanceWilliams1966`, `BrayCurtis1957` |
| total variation and metric relations | `BrayCurtis1957`, `GibbsSu2002` |
| relative density-ratio / density-difference methods | `SugiyamaEtAl2012`, `SugiyamaEtAl2013` |
| KDE two-sample discrepancy statistics | `AndersonHallTitterington1994` |
| localised two-sample / witness methods | `Duong2013`, `Gretton2012` |
| AlphaFold–disorder benchmarking | `PiovesanEtAl2022`, `Necci2021` |

## 3. Controlled tags

Thirteen of the fifteen controlled tags are used. Counts:

| Tag | Items |
|---|---|
| `directly-cited` | 38 |
| `disorder-biology` | 14 |
| `classical-statistics` | 11 |
| `closest-prior-art` | 10 |
| `alphafold-context` | 6 |
| `proxy-validation` | 6 |
| `kde-methodology` | 5 |
| `known-reparameterized` | 5 |
| `two-sample-testing` | 5 |
| `total-variation` | 3 |
| `negative-result` | 2 |
| `permutation-calibration` | 2 |
| `computational-result` | 1 |

Items may carry several tags. No tag outside the controlled list is used.

**Two controlled tags are deliberately not created**, because creating them empty would be
misleading:

* `apparently-distinct` — this classification belongs to two *repository results* (the
  taxon-dependent offset, and the range compression of the proxy), not to any *source*. No source in
  this bibliography is itself apparently distinct from prior art. The classification lives in §18 of
  `NOVELTY_AND_PROVENANCE.md`, where it applies to claims.
* `uncertain-more-search` — nothing in the bibliography is in that state. The audit's open questions
  (§21) are about unretrieved *data* — residue-level pLDDT, the exact DisProt query strings — not
  about unresolved literature.

## 4. Notes

A child note in the following fixed format is intended for the nine sources that carry
interpretive weight:

```
REPOSITORY CLAIM
SOURCE RESULT
STATISTICAL / BIOLOGICAL RELATIONSHIP
ASSUMPTIONS
DOES SOURCE IMPLY REPO RESULT?
RECOMMENDED WORDING
```

Priority order, with the claim each note addresses:

| Source | Note addresses |
|---|---|
| `Quaglia2022` | what DisProt annotation is, and why it is not ground truth |
| `PiovesanEtAl2022` | the actual pLDDT/disorder relationship, and conditional folding as the competing explanation for §12.1 |
| `Necci2021` | CAID performance, i.e. why residue-level conclusions must not be drawn here |
| `LanceWilliams1966` | Canberra distance as the direct ancestor of `Ā` |
| `BrayCurtis1957` | Bray–Curtis as the `A_w`/total-variation ancestor |
| `GibbsSu2002` | total variation's place among the metrics, and what `A_w` therefore is |
| `AndersonHallTitterington1994` | plug-in KDE two-sample bias — why the nonzero null baseline is expected |
| `HemerikGoeman2018` | validity conditions for the permutation procedure actually used |
| `Duong2013` | calibrated localisation for KDE differences — the closest method to what the profile does informally |

In every one of the nine, the answer to **DOES SOURCE IMPLY REPO RESULT?** is *no*: each supplies an
ingredient, a standard equivalent, or a benchmark context, and none states the repository's
cross-taxon protein-level calibration result.

## 5. Reapplying this

The filing and tagging specification is machine-readable, one entry per citation key, and validated
against `references.bib` (the plan and the bibliography must contain exactly the same 38 keys). Two
invariants are checked before any write: no subcollection is empty, and no tag falls outside the
controlled list.
