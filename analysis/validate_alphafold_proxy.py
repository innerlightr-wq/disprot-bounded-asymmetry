#!/usr/bin/env python3
"""Validate what DisProt's `alphafold_very_low_content` field actually is.

The repository never computes this field: `src/01_build_tables.py` copies it
verbatim out of the DisProt JSON export. Its meaning is therefore entirely
DisProt's, and the field name is the only documentation available -- DisProt's
public pages and API do not define it. This script settles the question
empirically, by reconstructing the fraction of very-low-confidence residues
directly from AlphaFold DB and comparing.

What it does, for a list of UniProt accessions:

  1. fetches AF-<acc>-F1-confidence_v<N>.json from AlphaFold DB (cached on disk);
  2. computes, from the per-residue pLDDT array,
         f50  = count(pLDDT <  50) / n_model
         f70  = count(pLDDT <  70) / n_model
         fcatD = count(confidenceCategory == "D") / n_model
     plus mean/median pLDDT;
  3. joins the stored DisProt field and DisProt's own sequence length;
  4. writes one row per accession, and a provenance manifest with the SHA-256
     of every downloaded file.

Nothing here trains a model, tunes a threshold, or selects a statistic: the
thresholds 50 and 70 are AlphaFold DB's own published confidence bands, and
category "D" is the provider's own label for the very-low band.

Network use is confined to AlphaFold DB. Downloads are cached under
--cache-dir, which is NOT part of the repository: the manifest plus this script
reconstruct them.

Usage
-----
    python analysis/validate_alphafold_proxy.py \
        --protein-table results/tables/disprot_protein_level.csv \
        --out-dir results/alphafold_validation \
        --cache-dir /tmp/af_cache \
        --strata          # deterministic stratified sample (default)

    python analysis/validate_alphafold_proxy.py --accessions P04637 P38398 ...
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
AF_URL = "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-confidence_v{ver}.json"
# AlphaFold DB serves only the current model version; older ones return 404.
# Versions are tried in this order and the first that resolves is recorded.
AF_VERSIONS = (6, 5, 4)
UA = "disprot-bounded-asymmetry/validation (mailto:dejesuselias10@gmail.com)"


def fetch_confidence(acc: str, cache: Path, pause: float = 0.2):
    """Return (payload, model_version, sha256, source_url, from_cache)."""
    for ver in AF_VERSIONS:
        dest = cache / f"AF-{acc}-F1-confidence_v{ver}.json"
        if dest.exists():
            raw = dest.read_bytes()
            return (json.loads(raw), ver, hashlib.sha256(raw).hexdigest(),
                    AF_URL.format(acc=acc, ver=ver), True)
    last = None
    for ver in AF_VERSIONS:
        url = AF_URL.format(acc=acc, ver=ver)
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                raw = resp.read()
        except urllib.error.HTTPError as exc:
            last = f"HTTP {exc.code}"
            continue
        except Exception as exc:                        # network-level failure
            last = repr(exc)
            continue
        cache.mkdir(parents=True, exist_ok=True)
        (cache / f"AF-{acc}-F1-confidence_v{ver}.json").write_bytes(raw)
        time.sleep(pause)
        return (json.loads(raw), ver, hashlib.sha256(raw).hexdigest(), url, False)
    return None, None, None, last, False


def summarise(payload) -> dict:
    scores = np.asarray(payload["confidenceScore"], dtype=float)
    cats = payload.get("confidenceCategory") or []
    n = scores.size
    if n == 0:
        # An empty score array must fail loudly: silently returning 0.0 would
        # enter the comparison as a perfectly-confident protein.
        raise ValueError("confidenceScore is empty; refusing to report a fraction")
    if not np.isfinite(scores).all():
        raise ValueError("confidenceScore contains non-finite values")
    out = {
        "n_model": int(n),
        "f50": float((scores < 50).sum() / n),
        "f70": float((scores < 70).sum() / n),
        "n_lt50": int((scores < 50).sum()),
        "mean_plddt": float(scores.mean()),
        "median_plddt": float(np.median(scores)),
    }
    if len(cats) == n:
        out["fcatD"] = float(sum(1 for c in cats if c == "D") / n)
        # the provider's own band edge, read off its own data
        d = [s for s, c in zip(scores, cats) if c == "D"]
        out["catD_max_score"] = float(max(d)) if d else float("nan")
    return out


def pick_strata(prot: pd.DataFrame, per_group: int) -> pd.DataFrame:
    """Deterministic stratified sample: sorted by accession within each stratum."""
    p = prot.dropna(subset=["af_very_low_content"]).copy()
    strata = {
        "full_disorder": p[p.D_exp >= 1.0],
        "near_zero_disorder": p[p.D_exp < 0.05],
        "intermediate": p[(p.D_exp >= 0.20) & (p.D_exp <= 0.60)],
    }
    frames = []
    for name, sub in strata.items():
        # spread across taxa: round-robin over groups, accession-sorted
        sub = sub.sort_values("acc")
        picks, per_taxon = [], {}
        for _, row in sub.iterrows():
            per_taxon.setdefault(row.group, []).append(row)
        order = sorted(per_taxon)
        i = 0
        while len(picks) < per_group and any(per_taxon[g] for g in order):
            g = order[i % len(order)]
            if per_taxon[g]:
                picks.append(per_taxon[g].pop(0))
            i += 1
        frames.append(pd.DataFrame(picks).assign(stratum=name))
    return pd.concat(frames, ignore_index=True)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--protein-table",
                    default=os.environ.get("DISPROT_PROTEIN_TABLE",
                                           ROOT / "results/tables/disprot_protein_level.csv"))
    ap.add_argument("--out-dir", default=ROOT / "results/alphafold_validation")
    ap.add_argument("--cache-dir", default=os.environ.get("AF_CACHE", "/tmp/af_cache"))
    ap.add_argument("--accessions", nargs="*", default=None,
                    help="explicit accession list; overrides the stratified sample")
    ap.add_argument("--per-group", type=int, default=20)
    ap.add_argument("--all", action="store_true",
                    help="use every accession in the protein table")
    args = ap.parse_args(argv)

    prot = pd.read_csv(args.protein_table)
    out_dir, cache = Path(args.out_dir), Path(args.cache_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.accessions:
        sel = prot[prot.acc.isin(args.accessions)].copy(); sel["stratum"] = "explicit"
    elif args.all:
        sel = prot.dropna(subset=["af_very_low_content"]).copy(); sel["stratum"] = "all"
    else:
        sel = pick_strata(prot, args.per_group)
    sel = sel.sort_values(["stratum", "acc"]).reset_index(drop=True)
    print(f"{len(sel)} accessions selected "
          f"({sel.stratum.value_counts().to_dict()})", flush=True)

    rows, manifest = [], []
    for i, r in sel.iterrows():
        payload, ver, sha, url, cached = fetch_confidence(r.acc, cache)
        status = "ok" if payload else "model_unavailable"
        rec = {"acc": r.acc, "disprot_id": r.disprot_id, "group": r.group,
               "stratum": r.stratum, "length_disprot": int(r.length),
               "D_exp": float(r.D_exp), "stored": float(r.af_very_low_content),
               "model_version": ver, "mapping_status": status}
        if payload:
            s = summarise(payload)
            rec.update(s)
            rec["coverage"] = s["n_model"] / int(r.length)
            if s["n_model"] != int(r.length):
                rec["mapping_status"] = "length_mismatch"
            rec["stored_minus_f50"] = rec["stored"] - s["f50"]
            rec["stored_times_len"] = rec["stored"] * int(r.length)
        rows.append(rec)
        manifest.append({"acc": r.acc,
                         "source_url": url,
                         "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                         "from_cache": cached,
                         "model_version": ver,
                         "length_disprot": int(r.length),
                         "n_model": rec.get("n_model"),
                         "sha256": sha,
                         "mapping_status": rec["mapping_status"]})
        if (i + 1) % 10 == 0:
            print(f"  {i+1}/{len(sel)}", flush=True)

    val = pd.DataFrame(rows)
    val.to_csv(out_dir / "proxy_validation.csv", index=False)
    pd.DataFrame(manifest).to_csv(out_dir / "alphafold_validation_manifest.csv", index=False)
    print(f"\nwrote {out_dir/'proxy_validation.csv'} and the provenance manifest")

    ok = val[val.mapping_status.isin(["ok", "length_mismatch"])].dropna(subset=["f50"])
    if not len(ok):
        print("no models retrieved; nothing to compare")
        return 1
    d = ok["stored"] - ok["f50"]
    print(f"\nstored vs reconstructed f50, n = {len(ok)}")
    print(f"  max |stored - f50|    : {d.abs().max():.6f}")
    print(f"  median |stored - f50| : {d.abs().median():.6f}")
    print(f"  mean  (stored - f50)  : {d.mean():+.6f}")
    print(f"  Spearman              : {ok['stored'].corr(ok['f50'], method='spearman'):.4f}")
    print(f"  Pearson               : {ok['stored'].corr(ok['f50']):.4f}")
    exact = (d.abs() < 1e-9).sum()
    print(f"  exact matches         : {exact}/{len(ok)}")
    if "fcatD" in ok:
        print(f"  max |f50 - fcatD|     : {(ok['f50']-ok['fcatD']).abs().max():.2e} "
              f"(provider category D vs pLDDT<50)")
    print("\nby stratum (median):")
    print(ok.groupby("stratum")[["D_exp", "stored", "f50", "f70", "median_plddt", "coverage"]]
            .median().round(4).to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
