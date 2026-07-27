# Input data

## Source database

**DisProt** — the manually curated database of intrinsically disordered
proteins. Release **2026_06**, "with ambiguous evidences" variant.

- Web: https://disprot.org
- Primary citation: Nugnes, M. V., Bouhraoua, K. E. A., Zoubiri, M., Pancsa, R.,
  Fichó, E., DisProt Consortium, Tompa, P., Piovesan, D., Tosatto, S. C. E., &
  Aspromonte, M. C. (2026). DisProt in 2026: enhancing intrinsically disordered
  proteins accessibility, deposition, and annotation. *Nucleic Acids Research*,
  54(D1), D383–D392. https://doi.org/10.1093/nar/gkaf1175
- **Licence: CC BY 4.0** (https://creativecommons.org/licenses/by/4.0/), per
  https://disprot.org/about, verified 2026-07-27. Redistribution is permitted
  with attribution, which is why `raw/` is included here. Full terms and
  attribution text: [`../DATA_LICENSE.md`](../DATA_LICENSE.md).

## Provenance table

| File | Organism (as recorded) | Taxon ID | Records | Retrieval date | Source |
|---|---|---|---|---|---|
| `raw/disprot_taxon_9606_homo_sapiens.json` | *Homo sapiens* | 9606 | 1,339 | *TODO — see below* | DisProt 2026_06, with ambiguous evidences |
| `raw/disprot_taxon_559292_saccharomyces_cerevisiae.json` | *Saccharomyces cerevisiae* (strain ATCC 204508 / S288c) | 559292 | 223 | *TODO — see below* | DisProt 2026_06, with ambiguous evidences |
| `raw/disprot_taxon_83333_escherichia_coli_k12.json` | *Escherichia coli* (strain K12) | 83333 | 146 | *TODO — see below* | DisProt 2026_06, with ambiguous evidences |
| `raw/disprot_taxon_6239_caenorhabditis_elegans.json` | *Caenorhabditis elegans* | 6239 | 58 | *TODO — see below* | DisProt 2026_06, with ambiguous evidences |

Machine-readable equivalent, with byte counts, MD5 and SHA-256 per file:
[`metadata/source_queries.csv`](metadata/source_queries.csv).

### Checksums

| File | Bytes | MD5 |
|---|---|---|
| `disprot_taxon_9606_homo_sapiens.json` | 12,221,224 | `0562d49c7a55bb2bbf018b49106c437a` |
| `disprot_taxon_559292_saccharomyces_cerevisiae.json` | 1,774,393 | `7f3444caf3ec364eedbf3a4e60801531` |
| `disprot_taxon_83333_escherichia_coli_k12.json` | 948,678 | `df0b07151610e12043d687362cb935c2` |
| `disprot_taxon_6239_caenorhabditis_elegans.json` | 377,380 | `8f652e6124b72fc87dd77a2521900697` |

These are regenerated on every run into
`../results/diagnostics/source_file_audit.csv`, so a reader can confirm their
copies match the files the reported numbers came from.

## Exact query and retrieval date — outstanding

> **TODO: Insert the exact DisProt query strings used for each taxonomic
> download from the original download history.**
>
> **TODO: Insert retrieval date(s) from the original download history.**
>
> Both fields are carried as these TODO placeholders in
> `metadata/source_queries.csv` (`exact_query_string`, `retrieval_date`).
> Neither is recoverable from the repository: the DisProt export format does
> not embed the query that produced it or the time it was run, and
> `disprot.org` was unreachable from the analysis environment (HTTP 403).
> Neither has been guessed, estimated or reconstructed.
>
> What *is* established from the file contents, and stated as fact rather than
> inference: each export is single-taxon, the four are mutually disjoint, and
> each is drawn from release 2026_06 in the with-ambiguous-evidences variant.
> Nothing further about the query parameters should be assumed.

## Filenames and taxonomy assignment

### Why these filenames

The files were downloaded as `DisProt_release_2026_06_with_ambiguous_evidences.json`,
`…__1_.json`, `…__2_.json` and `…__3_.json` — DisProt names exports by release,
not by organism, so four separate per-taxon queries all arrive with the same
name and the browser appends numeric suffixes. Those suffixes carry no
provenance, and they are **not stable across download sessions**: the archived
Zenodo deposit was produced in a session where the *E. coli* export was the
unsuffixed file, whereas in the session that built this repository it was the
*human* export. Renaming by taxon ID makes provenance legible and removes the
dependence on download accidents.

Mapping from original download names to repository names is recorded in
`metadata/source_queries.csv` (`original_download_filename` column).

### Strain-level, not kingdom-level, taxon IDs

Filenames use the **NCBI taxonomy ID that the records themselves carry**:

| Group label | Taxon ID | What it identifies |
|---|---|---|
| Human | 9606 | *Homo sapiens* |
| Fungal | 559292 | *S. cerevisiae* strain S288c — **not** 4751 (kingdom Fungi) |
| Bacterial | 83333 | *E. coli* strain K12 — **not** 2 (superkingdom Bacteria) |
| Nematode | 6239 | *C. elegans* |

The group labels "Fungal" and "Bacterial" are the manuscript's shorthand for
single-species samples, not claims about the kingdoms. Naming the files
`taxon_4751_fungi` and `taxon_2_bacteria` would overstate their scope.

### How group assignment works

`01_build_tables.py` maps `ncbi_taxon_id` → group label from inside each
record. **Filenames are never used to assign taxonomy.** As a guard against a
mis-rename, the script compares each file's recorded taxon against the taxon
its filename claims and raises `ValueError` on disagreement rather than
producing mislabelled output.

## AlphaFold retrieval notes

Each DisProt record carries a precomputed scalar `alphafold_very_low_content`,
denoted *V* in the manuscript and treated as the fraction of residues with
pLDDT < 50 (the AlphaFold DB "very low" confidence band).

**This interpretation is unverified.** It follows the field name and the
published confidence-band convention, but could not be checked against
AlphaFold DB: every required endpoint returned HTTP 403 with proxy reason
`host_not_allowed`. See
`../results/diagnostics/alphafold_retrieval_manifest.csv` for the timestamped
log of all five endpoints attempted.

Three consequences constrain what the data support:

1. **No residue-level pLDDT exists in this repository.** No per-residue
   sensitivity, specificity, ROC or MCC analysis was or can be performed here.
2. **The AlphaFold model version behind *V* is unknown.** This matters:
   AlphaFold DB has been rebuilt across releases, and AlphaFold3 produces
   markedly different pLDDT distributions in disordered regions than
   AlphaFold2.
3. ***V* is absent for 103 proteins** (92 human, 4 fungal, 3 bacterial, 4
   nematode). Only 18 of those exceed 2,700 residues, so model length limits do
   not explain most of the missingness, and the cause is unresolved. They are
   listed in `../results/diagnostics/failed_accessions.csv` and excluded
   pairwise, not imputed.

## Mixed data types requiring coercion

**This is the most likely way to get wrong numbers from these files.**

`alphafold_very_low_content` is not consistently typed in the export. Across the
1,766 records it appears as float, int, JSON `null`, and — for **34 records** —
as a **quoted string**. Loading it without coercion yields an object-dtype
pandas column; medians and correlations computed on that column are wrong or
raise.

The pipeline handles this explicitly:

```python
prot["af_very_low_content"] = pd.to_numeric(prot["af_very_low_content"],
                                            errors="coerce")
```

Observed types across the 1,766 records: 1,562 float, 67 int, 34 string, 103
JSON `null`.

One of the 34 strings is a trap worth naming. Accession **`P09651-2`** carries
the literal string **`"NaN"`**. It passes an `is not None` check but fails
numeric conversion, so the count of proteins "with a value" differs by one
depending on which test is used. This is the origin of the manuscript's 92-vs-93
/ 1,247-vs-1,246 off-by-one; see §6 of `../results/README.md` for the full
explanation.

Any reanalysis of these files must coerce, and must decide explicitly whether
"has a value" means non-null or numerically usable.

## Other data-quality caveats carried into the manuscript

1. **The exported `regions` array is a filtered subset of the curated set.**
   For 545 proteins it holds fewer entries than the record's own
   `regions_counter` field: 7,535 exported regions against 9,610 declared. The
   discrepancy could not be resolved against the live database (API blocked).
   The integrity check `regions_counter_matches_region_rows` reports 1,221 of
   1,766 proteins in agreement, and that number is expected to be below 1,766
   for this reason.
2. **`disorder_content` includes ambiguous residues.** For all 1,766 proteins
   the record-level `disorder_content` field equals the residue fraction covered
   by structural-state consensus segments of type **D or S** — i.e. it includes
   ambiguous molten-globule and pre-molten-globule residues. The manuscript's
   Eq. (5) therefore recomputes disorder fraction from **type-D segments only**
   (`D_exp`), and `D_exp` is what every reported statistic uses. Only 4 proteins
   carry type-S segments, so the choice matters for provenance more than for the
   numbers — but the two definitions are not interchangeable, and
   `disorder_content_source` is retained in the output table for comparison.
3. **Consensus type flags require explicit filtering.** Structural-state
   segments must be filtered on `type` before any coverage computation;
   unfiltered coverage silently mixes definitions. See `coverage()` in
   `01_build_tables.py`.
4. **Unannotated residues are not evidence of order.** `D_exp` treats every
   residue outside a curated disorder segment as contributing to the
   denominator. DisProt curation is region-targeted, so unannotated residues are
   largely uncharacterised rather than known-ordered. This is a stated
   limitation of the study, not a defect of the data.

## Ascertainment

DisProt proteins are present because their disorder was experimentally
investigated. They are not a random sample of any proteome, and the four groups
are not matched on length, function, subcellular localisation or study
intensity. No proteome-wide statement follows from any number in this
repository. Because UniProt was also unreachable, the planned comparison against
reference-proteome metadata was not performed, so the magnitude of the
ascertainment bias is **undocumented rather than small**.

## File descriptions

| Path | Description |
|---|---|
| `raw/*.json` | Four unmodified DisProt per-organism exports. Each is a JSON object with a `data` array of per-protein records carrying `acc`, `disprot_id`, `organism`, `ncbi_taxon_id`, `name`, `genes`, `length`, `sequence`, `regions`, `regions_counter`, `disorder_content`, `disprot_consensus`, `alphafold_very_low_content`, `UniParc`, `uniref50`, `released`. |
| `metadata/source_queries.csv` | Per-file provenance: group, repository and original filenames, release, variant, taxon ID, organism string, record and accession counts, byte count, MD5, SHA-256, and two fields (`exact_query_string`, `retrieval_date`) carried as TODO placeholders. |
