# Manuscript status

Which document says what, and which of its claims still stand. Written 2026-09-18 during the
AlphaFold-proxy validation audit; see [`NOVELTY_AND_PROVENANCE.md`](NOVELTY_AND_PROVENANCE.md).

**No PDF was edited, and no Zenodo version was created.** This file exists because the repository
previously described one of the PDFs as a corrected manuscript when it is not one.

## 1. The three documents

| File | Bytes | SHA-256 |
|---|---|---|
| `paper/manuscript.pdf` | 356,597 | `89aca2675fadc76a3f5399a43d71e77663ccbd36d436b742cb46dc26bbdf70da` |
| `revision_addendum/bounded-asymmetry-profiles-cross-taxon-protein-disorder-2026.pdf` | 361,875 | `b6b67c5e9af92819964b540f3a74178318fd152bb7c3ace5900c8c6bde5a5087` |
| `revision_addendum/ADDENDUM_permutation_calibration.pdf` | 180,808 | `f229db0f91674a896bf888f2dde54d2b5d38996165046ea08b2055212b94ec2f` |

## 2. The two paper PDFs are the same paper

Extracted with `pdftotext -layout` and compared as word-token sequences:

* the **first 4,521 word tokens are byte-identical** (same MD5 of the token stream) — that is the
  whole of the title, abstract, introduction, methods, results and discussion;
* divergence begins at the archival-deposit DOI at the end of the body
  (`zenodo.21613281` → `zenodo.21627908`);
* the remaining 235 differing word tokens are entirely inside the reference list and
  acknowledgements, and are hyphenation and line-reflow differences (`dis-tance` vs `distance`,
  `https:` split across a line, and so on) plus one added author name in a reference;
* the string **"permut" appears 0 times in either paper PDF**, and 25 times in the addendum.

So `bounded-asymmetry-profiles-cross-taxon-protein-disorder-2026.pdf` is **a re-deposit of the same
manuscript, not a revision**. It contains no permutation calibration and no corrected claim.

## 3. Current status of each document

| Document | Status |
|---|---|
| `paper/manuscript.pdf` | **Superseded in part.** Its methodological headline claim is withdrawn (see §4). Its data, tables and figures reproduce exactly. |
| `…cross-taxon-protein-disorder-2026.pdf` | **Superseded in the same way and to the same extent** — it is the same text. Retained as the deposited record of that version. Not a correction. |
| `ADDENDUM_permutation_calibration.pdf` | **Current for the statistical claims.** Its §5 disposition table is authoritative for what was withdrawn, reversed and retained. |
| `docs/NOVELTY_AND_PROVENANCE.md` | **Current for prior art, and for the biological interpretation** after the September 2026 audits. |

## 4. Claims in the paper PDFs that are superseded

| Claim in both paper PDFs | Superseded by | Current status |
|---|---|---|
| "the profile … resolves localised structure (Ā = 0.143, C = 0.57) in a comparison that global tests read as null" (abstract) | addendum §5–§6 | **Withdrawn.** The observed Ā lies inside its own permutation null (median 0.111, 95th pct 0.174, raw p = 0.18). |
| Two zero crossings indicate real sign-change regions | addendum §5 | **Withdrawn as inferential**; survives only as description. The simultaneous band excludes zero nowhere in that comparison. |
| The profile adds sensitivity beyond scalar two-sample testing | addendum §5; re-derived in `NOVELTY_AND_PROVENANCE.md` §8.2 | **Reversed.** Under a matched per-feature Holm correction the profile rejects 5 comparisons and KS rejects 6, the profile's set being a strict subset. |
| "the AlphaFold-derived very-low-confidence fraction is a poor protein-level surrogate for curated disorder content in every taxon examined (Spearman ρ = 0.07–0.24)" | `NOVELTY_AND_PROVENANCE.md` §23 | **Materially qualified.** The correlation is suppressed by protein length, which is negatively associated with the curated fraction and positively with the proxy. Controlling for length raises the association substantially. The numbers as printed are correct; the interpretation "poor surrogate" is too strong. |
| The statistic is presented as newly defined, without prior art | `NOVELTY_AND_PROVENANCE.md` §3–§6 | **Attribution missing.** S(x) is a bounded density-ratio transform, Ā a regularised continuous Canberra distance, A_w exactly total variation at ε = 0. |

Claims that **stand unchanged**: every number, table and figure in the pipeline (reproduced exactly);
the estimator-dependence results; the fungal enrichment results (now supported under permutation);
and the observation that taxon-specific offsets exist between the proxy and curated disorder.

## 5. What a corrected version would need

1. An abstract that does not assert the withdrawn localisation claim.
2. Prior-art attribution for the statistic (`references.bib` supplies it).
3. The AlphaFold-proxy result restated with protein length controlled, and with the DisProt
   annotation-coverage caveat stated explicitly.
4. No change to any computed value.

Sequencing is the author's decision. Nothing here should be deposited before the abstract is fixed,
because the current abstract is the one carrying the withdrawn claim.
