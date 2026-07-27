#!/usr/bin/env python3
"""
Phase 1: audit the four DisProt JSON exports and build protein-level and
region-level tables.

Provenance note: input files are DisProt release 2026_06 "with ambiguous
evidences" per-organism JSON exports (4 files, mutually disjoint).
No network access is used by this script.

Portability note (2026-07-27 refactor). Two changes only, neither of which
touches a computation:

  1. Absolute paths (/mnt/user-data/...) replaced with project-relative
     pathlib paths and DISPROT_DATA / DISPROT_OUT environment overrides.
  2. Input discovery by sorted(glob.glob("*.json")) replaced with the
     explicit SOURCE_FILES mapping below.

Change 2 fixes a real reproducibility defect. Under glob discovery, table row
order depended on arbitrary browser-assigned download suffixes ("__1_",
"__2_", ...), which differed between sessions; the archived deposit is ordered
Bacterial-first for exactly this reason. No computed value depends on load
order, because group membership is derived from ncbi_taxon_id (never from the
filename) and within-group record order comes from each file's internal
ordering. Load order is now fixed to the group order of manuscript Table 1.
"""
import hashlib
import json
import os
from pathlib import Path

import pandas as pd

# ------------------------------------------------------------------- paths
ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = Path(os.environ.get("DISPROT_DATA", ROOT / "data" / "raw"))
RESULTS_DIR = Path(os.environ.get("DISPROT_OUT", ROOT / "results"))
TABLES_DIR = RESULTS_DIR / "tables"
DIAG_DIR = RESULTS_DIR / "diagnostics"
for _d in (TABLES_DIR, DIAG_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# Explicit, ordered input manifest: (group label, NCBI taxon ID, filename).
# Taxon IDs are the strain-level identifiers carried by the records
# themselves, not kingdom-level approximations.
SOURCE_FILES = [
    ("Human",     9606,   "disprot_taxon_9606_homo_sapiens.json"),
    ("Fungal",    559292, "disprot_taxon_559292_saccharomyces_cerevisiae.json"),
    ("Bacterial", 83333,  "disprot_taxon_83333_escherichia_coli_k12.json"),
    ("Nematode",  6239,   "disprot_taxon_6239_caenorhabditis_elegans.json"),
]

GROUP_BY_TAXON = {
    9606: "Human",
    559292: "Fungal",
    83333: "Bacterial",
    6239: "Nematode",
}


def coverage(segs, want_types, length):
    """Residue count covered by consensus segments of the requested types."""
    cov = set()
    for s in segs:
        if s.get("type") in want_types:
            cov.update(range(s["start"], s["end"] + 1))
    cov = {r for r in cov if 1 <= r <= length}
    return len(cov)


prot_rows, region_rows, file_rows = [], [], []

for _group, _taxon, _fname in SOURCE_FILES:
    path = DATA_DIR / _fname
    if not path.exists():
        raise FileNotFoundError(
            f"Required input not found: {path}\n"
            f"Place the four DisProt JSON exports in {DATA_DIR}, or point "
            f"the DISPROT_DATA environment variable at the directory that "
            f"holds them."
        )
    raw = path.read_bytes()
    data = json.loads(raw)["data"]

    observed_taxon = data[0]["ncbi_taxon_id"]
    if observed_taxon != _taxon:
        raise ValueError(
            f"{path.name}: expected NCBI taxon {_taxon} ({_group}) but the "
            f"records carry taxon {observed_taxon}. The input files may have "
            f"been renamed incorrectly; refusing to proceed rather than "
            f"mislabel provenance."
        )

    file_rows.append({
        "file": path.name,
        "bytes": len(raw),
        "md5": hashlib.md5(raw).hexdigest(),
        "n_records": len(data),
        "n_unique_acc": len({r["acc"] for r in data}),
        "organism": data[0]["organism"],
        "ncbi_taxon_id": data[0]["ncbi_taxon_id"],
    })

    for r in data:
        L = r["length"]
        ss = r.get("disprot_consensus", {}).get("Structural state", [])
        n_d = coverage(ss, {"D"}, L)
        n_ds = coverage(ss, {"D", "S"}, L)
        af = r.get("alphafold_very_low_content")

        prot_rows.append({
            "acc": r["acc"],
            "disprot_id": r["disprot_id"],
            "organism": r["organism"],
            "ncbi_taxon_id": r["ncbi_taxon_id"],
            "group": GROUP_BY_TAXON.get(r["ncbi_taxon_id"], "Other"),
            "protein_name": r["name"],
            "gene": (r["genes"][0]["name"]["value"]
                     if r.get("genes") and r["genes"][0].get("name") else None),
            "length": L,
            "seq_len_check": len(r.get("sequence") or ""),
            "n_regions": r["regions_counter"],
            "n_disorder_res_D": n_d,
            "n_disorder_res_D_or_S": n_ds,
            "D_exp": n_d / L,
            "D_exp_ambiguous": n_ds / L,
            "disorder_content_source": r["disorder_content"],
            "af_very_low_content": af,
            "af_available": af is not None,
            "has_ambiguous_S": n_ds != n_d,
            "uniparc": r.get("UniParc"),
            "uniref50": r.get("uniref50"),
            "released": r.get("released"),
            "source_file": path.name,
        })

        for g in r.get("regions", []):
            region_rows.append({
                "acc": r["acc"],
                "disprot_id": r["disprot_id"],
                "region_id": g.get("region_id"),
                "organism": r["organism"],
                "ncbi_taxon_id": r["ncbi_taxon_id"],
                "group": GROUP_BY_TAXON.get(r["ncbi_taxon_id"], "Other"),
                "protein_length": L,
                "start": g.get("start"),
                "end": g.get("end"),
                "region_len": (g["end"] - g["start"] + 1)
                              if g.get("start") and g.get("end") else None,
                "term_namespace": g.get("term_namespace"),
                "term_id": g.get("term_id"),
                "term_ontology": g.get("term_ontology"),
                "ec_id": g.get("ec_id"),
                "ec_name": g.get("ec_name"),
                "ec_ontology": g.get("ec_ontology"),
                "reference_id": g.get("reference_id"),
                "curator_name": g.get("curator_name"),
                "released": g.get("released"),
                "version": g.get("version"),
            })

prot = pd.DataFrame(prot_rows)
prot["af_very_low_content"] = pd.to_numeric(prot["af_very_low_content"],
                                            errors="coerce")
reg = pd.DataFrame(region_rows)
files = pd.DataFrame(file_rows)

# ---- integrity checks -------------------------------------------------
checks = []


def chk(name, value):
    checks.append({"check": name, "value": value})


chk("n_protein_records", len(prot))
chk("n_unique_acc", prot["acc"].nunique())
chk("n_unique_disprot_id", prot["disprot_id"].nunique())
chk("duplicate_acc_across_files", int(prot["acc"].duplicated().sum()))
chk("length_vs_sequence_mismatch",
    int((prot["length"] != prot["seq_len_check"]).sum()))
chk("missing_accession", int(prot["acc"].isna().sum()))
chk("D_exp_out_of_unit_interval",
    int(((prot["D_exp"] < 0) | (prot["D_exp"] > 1)).sum()))
chk("proteins_with_ambiguous_S_segments", int(prot["has_ambiguous_S"].sum()))
chk("disorder_content_matches_D_or_S_coverage",
    int((abs(prot["D_exp_ambiguous"] - prot["disorder_content_source"])
         < 1e-6).sum()))
chk("af_value_missing", int((~prot["af_available"]).sum()))
chk("n_regions_total", len(reg))
chk("regions_missing_boundaries",
    int(reg[["start", "end"]].isna().any(axis=1).sum()))
chk("regions_with_end_beyond_length",
    int((reg["end"] > reg["protein_length"]).sum()))
chk("regions_counter_matches_region_rows",
    int((prot.set_index("acc")["n_regions"]
         == reg.groupby("acc").size().reindex(prot["acc"]).fillna(0).values
         ).sum()))
chk("duplicate_region_ids", int(reg["region_id"].duplicated().sum()))

prot.to_csv(TABLES_DIR / "disprot_protein_level.csv", index=False)
reg.to_csv(TABLES_DIR / "disprot_region_level.csv", index=False)
files.to_csv(DIAG_DIR / "source_file_audit.csv", index=False)
pd.DataFrame(checks).to_csv(DIAG_DIR / "phase1_integrity_checks.csv",
                            index=False)

print(files.to_string(index=False))
print()
print(pd.DataFrame(checks).to_string(index=False))
print()
print(prot.groupby("group").agg(
    n=("acc", "size"),
    af_ok=("af_available", "sum"),
    med_len=("length", "median"),
    med_Dexp=("D_exp", "median"),
    med_af=("af_very_low_content", "median"),
).to_string())
