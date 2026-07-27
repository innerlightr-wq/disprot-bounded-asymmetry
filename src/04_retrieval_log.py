#!/usr/bin/env python3
"""Retrieval provenance log and list of accessions with no AlphaFold value.

Every external database required by the protocol was unreachable from the
generating environment (HTTP 403, proxy reason host_not_allowed). This script
records those attempts verbatim so the failure is part of the archived record
rather than an undocumented gap, and lists the 103 accessions for which the
DisProt export carries no alphafold_very_low_content value.

Portability note (2026-07-27 refactor): absolute paths replaced with
project-relative pathlib paths plus DISPROT_OUT override, and the module was
reformatted for legibility. The logged rows and the derived accession list are
unchanged.
"""
import os
from pathlib import Path

import pandas as pd

# ------------------------------------------------------------------- paths
ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = Path(os.environ.get("DISPROT_OUT", ROOT / "results"))
TABLES_DIR = RESULTS_DIR / "tables"
DIAG_DIR = RESULTS_DIR / "diagnostics"
DIAG_DIR.mkdir(parents=True, exist_ok=True)

rows = [
    {"source": "AlphaFold Protein Structure Database",
     "endpoint": "https://alphafold.ebi.ac.uk/api/prediction/{acc}",
     "purpose": "residue-level pLDDT, model id/version",
     "attempt_date": "2026-07-26",
     "http_status": 403,
     "proxy_deny_reason": "host_not_allowed",
     "outcome": "FAILED",
     "note": "Sandbox egress allowlist excludes this host; no accessions retrieved."},
    {"source": "AlphaFold Protein Structure Database",
     "endpoint": "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-confidence_v4.json",
     "purpose": "confidenceScore array (pLDDT)",
     "attempt_date": "2026-07-26",
     "http_status": None,
     "proxy_deny_reason": "host_not_allowed",
     "outcome": "FAILED",
     "note": "Not reachable from container; also not fetchable via the assistant web tool (URL not in prior results)."},
    {"source": "UniProt REST",
     "endpoint": "https://rest.uniprot.org/uniprotkb/{acc}.json",
     "purpose": "sequence, length, review status, isoform status",
     "attempt_date": "2026-07-26",
     "http_status": 403,
     "proxy_deny_reason": "host_not_allowed",
     "outcome": "FAILED",
     "note": "Sequence and length were instead taken from the DisProt export and internally cross-checked."},
    {"source": "MobiDB",
     "endpoint": "https://mobidb.org/api/download?acc={acc}&format=json",
     "purpose": "consensus disorder annotation",
     "attempt_date": "2026-07-26",
     "http_status": 403,
     "proxy_deny_reason": "host_not_allowed",
     "outcome": "FAILED",
     "note": "No MobiDB comparison possible."},
    {"source": "DisProt API",
     "endpoint": "https://disprot.org/api/{disprot_id}",
     "purpose": "verify export completeness vs regions_counter",
     "attempt_date": "2026-07-26",
     "http_status": 403,
     "proxy_deny_reason": "host_not_allowed",
     "outcome": "FAILED",
     "note": "regions_counter discrepancy could not be resolved against the live database."},
]
pd.DataFrame(rows).to_csv(DIAG_DIR / "alphafold_retrieval_manifest.csv",
                          index=False)

p = pd.read_csv(TABLES_DIR / "disprot_protein_level.csv")
f = p[~p.af_available][["acc", "disprot_id", "organism", "group", "length",
                        "protein_name"]].copy()
f["reason"] = ("alphafold_very_low_content absent from DisProt export; "
               "AlphaFold DB not reachable to determine cause (candidate "
               "causes: no AFDB model, length outside AFDB limits, or "
               "non-canonical isoform)")
f.to_csv(DIAG_DIR / "failed_accessions.csv", index=False)

print("retrieval manifest rows:", len(rows))
print("failed/unavailable accessions:", len(f))
print(f.groupby("group").size().to_string())
print("\nlength of unavailable, human:", sorted(f[f.group == 'Human'].length)[-5:])
print("n unavailable with length > 2700:", (f.length > 2700).sum())
