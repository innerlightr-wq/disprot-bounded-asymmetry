"""Unit tests for analysis/validate_alphafold_proxy.py.

No network access: every test builds its own tiny confidence payload. The module
is loaded by path with importlib, matching the convention already used by
tests/test_null_calibration_grid.py, because `analysis/` is not a package.
"""
import importlib.util
import json
import math
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "analysis" / "validate_alphafold_proxy.py"

_spec = importlib.util.spec_from_file_location("validate_alphafold_proxy", MODULE_PATH)
vap = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(vap)


def payload(scores, cats=None):
    out = {"residueNumber": list(range(1, len(scores) + 1)), "confidenceScore": list(scores)}
    if cats is not None:
        out["confidenceCategory"] = list(cats)
    return out


# --------------------------------------------------------------- f50 / f70
def test_f50_counts_strictly_below_fifty():
    # 50.0 itself must NOT count: AlphaFold DB's "very low" band is pLDDT < 50
    s = summarise = vap.summarise(payload([10.0, 49.99, 50.0, 50.01, 90.0]))
    assert s["n_lt50"] == 2
    assert s["f50"] == pytest.approx(2 / 5)


def test_f70_counts_strictly_below_seventy():
    s = vap.summarise(payload([10.0, 69.99, 70.0, 95.0]))
    assert s["f70"] == pytest.approx(2 / 4)


def test_all_low_and_all_high_are_the_endpoints():
    assert vap.summarise(payload([1.0, 2.0, 3.0]))["f50"] == 1.0
    assert vap.summarise(payload([99.0, 98.0]))["f50"] == 0.0


# ------------------------------------------------------------- denominator
def test_denominator_is_the_number_of_modelled_residues():
    """f50 divides by len(confidenceScore), never by a stored sequence length."""
    s = vap.summarise(payload([10.0] * 3 + [90.0] * 7))
    assert s["n_model"] == 10
    assert s["f50"] == pytest.approx(0.3)


def test_fraction_not_percentage():
    """The field must be a fraction in [0,1], not a percentage."""
    s = vap.summarise(payload([10.0] * 50 + [90.0] * 50))
    assert 0.0 <= s["f50"] <= 1.0
    assert s["f50"] == pytest.approx(0.5)
    assert s["f50"] != pytest.approx(50.0)


# ------------------------------------------------- provider category check
def test_provider_category_D_matches_below_fifty():
    scores = [30.0, 49.9, 50.1, 80.0, 95.0]
    cats = ["D", "D", "L", "M", "H"]
    s = vap.summarise(payload(scores, cats))
    assert s["fcatD"] == pytest.approx(s["f50"])
    assert s["catD_max_score"] == pytest.approx(49.9)


def test_category_array_of_wrong_length_is_ignored():
    s = vap.summarise(payload([30.0, 90.0], ["D"]))       # mismatched length
    assert "fcatD" not in s


# ------------------------------------------------------- degenerate inputs
def test_single_residue_model():
    s = vap.summarise(payload([20.0]))
    assert s["n_model"] == 1 and s["f50"] == 1.0


def test_empty_model_raises_rather_than_returning_zero():
    """An empty score array must not silently yield 0.0 or NaN."""
    with pytest.raises(Exception):
        vap.summarise(payload([]))


def test_median_and_mean_are_reported():
    s = vap.summarise(payload([10.0, 20.0, 90.0]))
    assert s["median_plddt"] == pytest.approx(20.0)
    assert s["mean_plddt"] == pytest.approx(40.0)


# ---------------------------------------------------------- cache / fetch
def test_fetch_prefers_the_cache_and_reports_the_version(tmp_path):
    """A cached file must be used without any network call, and its version kept."""
    body = json.dumps(payload([10.0, 90.0], ["D", "H"])).encode()
    (tmp_path / "AF-P12345-F1-confidence_v6.json").write_bytes(body)
    data, ver, sha, url, cached = vap.fetch_confidence("P12345", tmp_path)
    assert cached is True
    assert ver == 6
    assert data["confidenceScore"] == [10.0, 90.0]
    assert len(sha) == 64
    assert "AF-P12345-F1-confidence_v6.json" in url


def test_missing_accession_is_reported_not_raised(tmp_path, monkeypatch):
    """An unavailable model yields a None payload so the caller can record it."""
    def boom(*a, **k):
        raise OSError("no network in tests")
    monkeypatch.setattr(vap.urllib.request, "urlopen", boom)
    data, ver, sha, info, cached = vap.fetch_confidence("P00000", tmp_path)
    assert data is None and ver is None and sha is None


# -------------------------------------------------------- strata selection
def test_strata_selection_is_deterministic_and_taxon_spread():
    import pandas as pd
    rows = []
    for i in range(40):
        rows.append({"acc": f"A{i:03d}", "disprot_id": f"DP{i:05d}",
                     "group": ["Human", "Fungal", "Bacterial", "Nematode"][i % 4],
                     "length": 300, "D_exp": 1.0 if i < 20 else 0.01,
                     "af_very_low_content": 0.2})
    prot = pd.DataFrame(rows)
    a = vap.pick_strata(prot, per_group=8)
    b = vap.pick_strata(prot, per_group=8)
    assert list(a.acc) == list(b.acc)                      # deterministic
    full = a[a.stratum == "full_disorder"]
    assert len(full) == 8
    assert full.group.nunique() == 4                       # spread across taxa
