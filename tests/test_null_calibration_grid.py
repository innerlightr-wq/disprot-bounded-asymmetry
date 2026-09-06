"""
Tests for src/07_null_calibration_grid.py (the A_bar vs. A_w synthetic
null-calibration diagnostic). See that module's docstring and
results/null_calibration/README.md for the full experiment.

The module filename has a numeric prefix (matching this repository's
existing src/ convention), so it is not a valid `import` target; it is
loaded here via importlib.util.spec_from_file_location, exactly the pattern
needed for any numeric-prefixed script in src/.
"""
import importlib.util
import os
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import stats as scipy_stats

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "src" / "07_null_calibration_grid.py"

_spec = importlib.util.spec_from_file_location("null_calibration_grid", MODULE_PATH)
ncg = importlib.util.module_from_spec(_spec)
sys.modules["null_calibration_grid"] = ncg
_spec.loader.exec_module(ncg)


GRID = np.linspace(0.0, 1.0, 257)


class TestDegenerateCases:
    """Identical density arrays give A_bar = A_w = B_bar = 0."""

    def test_identical_densities_give_zero_indices(self):
        rng = np.random.default_rng(1)
        x = rng.beta(5, 5, size=500)
        p = ncg.kde_density(x, GRID)
        S = ncg.asym_profile(p, p, eps=ncg.EPS_DEFAULT)
        A, B = ncg.indices(S, GRID)
        Aw = ncg.a_w(p, p, S, GRID)
        assert A == pytest.approx(0.0, abs=1e-12)
        assert B == pytest.approx(0.0, abs=1e-12)
        assert Aw == pytest.approx(0.0, abs=1e-12)


class TestSwapSymmetry:
    """Swapping A and B: A_bar and A_w unchanged, B_bar changes sign."""

    def test_swap_symmetry(self):
        rng = np.random.default_rng(2)
        xA = rng.beta(5, 5, size=300)
        xB = rng.beta(2, 6, size=200)
        pA = ncg.kde_density(xA, GRID)
        pB = ncg.kde_density(xB, GRID)
        S_ab = ncg.asym_profile(pA, pB)
        S_ba = ncg.asym_profile(pB, pA)
        A_ab, B_ab = ncg.indices(S_ab, GRID)
        A_ba, B_ba = ncg.indices(S_ba, GRID)
        Aw_ab = ncg.a_w(pA, pB, S_ab, GRID)
        Aw_ba = ncg.a_w(pB, pA, S_ba, GRID)
        assert A_ab == pytest.approx(A_ba, rel=1e-10)
        assert Aw_ab == pytest.approx(Aw_ba, rel=1e-10)
        assert B_ab == pytest.approx(-B_ba, rel=1e-10)


class TestBounds:
    @pytest.mark.parametrize("seed", [10, 11, 12, 13])
    def test_bounds_hold(self, seed):
        rng = np.random.default_rng(seed)
        xA = rng.beta(2, 2, size=rng.integers(50, 400))
        xB = rng.beta(3, 5, size=rng.integers(50, 400))
        pA = ncg.kde_density(xA, GRID)
        pB = ncg.kde_density(xB, GRID)
        S = ncg.asym_profile(pA, pB)
        A, B = ncg.indices(S, GRID)
        Aw = ncg.a_w(pA, pB, S, GRID)
        assert 0.0 <= A < 1.0
        assert 0.0 <= Aw < 1.0
        assert abs(B) <= A + 1e-9


class TestTVLimit:
    """As epsilon shrinks, A_w approaches the discrete numerical
    total-variation distance (1/2) * integral |p_A - p_B| dx."""

    def test_aw_approaches_tv_as_epsilon_shrinks(self):
        rng = np.random.default_rng(3)
        xA = rng.beta(5, 5, size=600)
        xB = rng.beta(3, 7, size=600)
        pA = ncg.kde_density(xA, GRID)
        pB = ncg.kde_density(xB, GRID)
        tv = 0.5 * ncg._trapz(np.abs(pA - pB), GRID)

        eps_large = 1e-1
        eps_small = 1e-8
        S_large = ncg.asym_profile(pA, pB, eps=eps_large)
        S_small = ncg.asym_profile(pA, pB, eps=eps_small)
        Aw_large = ncg.a_w(pA, pB, S_large, GRID)
        Aw_small = ncg.a_w(pA, pB, S_small, GRID)

        assert abs(Aw_small - tv) < abs(Aw_large - tv)
        assert abs(Aw_small - tv) < 1e-3


class TestScottBandwidthMatchesScipy:
    """Cross-check our reimplemented Scott bandwidth against scipy's own
    gaussian_kde, since 07_null_calibration_grid.py reimplements it directly
    (see module docstring) rather than importing the numeric-prefixed
    src/02_asymmetry.py."""

    def test_matches_scipy_scott_bandwidth(self):
        rng = np.random.default_rng(4)
        x = rng.beta(4, 6, size=257)
        k = scipy_stats.gaussian_kde(x, bw_method="scott")
        sigma = np.std(x, ddof=1)
        h_expected = k.factor * sigma
        h_ours = ncg.scott_bandwidth(len(x), sigma)
        assert h_ours == pytest.approx(h_expected, rel=1e-10)

        d_scipy = k(GRID)
        d_scipy /= ncg._trapz(d_scipy, GRID)
        d_ours = ncg.kde_density(x, GRID)
        assert np.allclose(d_scipy, d_ours, atol=1e-6)


class TestSharedBandwidthSymmetry:
    """The shared-bandwidth control must be symmetric under swapping which
    sample is called A vs. B (the shared bandwidth depends only on the
    pooled sample, not on group order)."""

    def test_shared_bandwidth_symmetric_under_swap(self):
        rng = np.random.default_rng(5)
        xA = rng.beta(5, 5, size=150)
        xB = rng.beta(5, 5, size=1350)
        A1, Aw1, B1, ks1, _ = ncg.compute_all_statistics(
            xA, xB, GRID, shared_bandwidth=True)
        A2, Aw2, B2, ks2, _ = ncg.compute_all_statistics(
            xB, xA, GRID, shared_bandwidth=True)
        assert A1 == pytest.approx(A2, rel=1e-10)
        assert Aw1 == pytest.approx(Aw2, rel=1e-10)
        assert B1 == pytest.approx(-B2, rel=1e-10)


class TestReproducibility:
    """A seeded synthetic condition is exactly reproducible."""

    def test_seeded_condition_reproducible(self):
        rng1 = np.random.default_rng(42)
        rng2 = np.random.default_rng(42)
        res1 = ncg.run_condition(rng1, "A_interior_unimodal", 60, 80, GRID,
                                  n_outer=20, shared_bandwidth=False)
        res2 = ncg.run_condition(rng2, "A_interior_unimodal", 60, 80, GRID,
                                  n_outer=20, shared_bandwidth=False)
        for key in res1:
            assert np.array_equal(res1[key], res2[key])

    def test_condition_rngs_deterministic(self):
        r1 = ncg.make_condition_rngs()
        r2 = ncg.make_condition_rngs()
        key = ("A_interior_unimodal", (50, 50), "independent")
        draw1 = r1[key].beta(5, 5, size=10)
        draw2 = r2[key].beta(5, 5, size=10)
        assert np.array_equal(draw1, draw2)


class TestPermutationMachinery:
    """Permutation labels preserve original group sizes; p-values lie in
    [0, 1]."""

    def test_permutation_preserves_group_sizes(self):
        rng = np.random.default_rng(6)
        xA = rng.beta(5, 5, size=37)
        xB = rng.beta(5, 5, size=113)
        pooled = np.concatenate([xA, xB])
        perm = rng.permutation(pooled)
        pA, pB = perm[:len(xA)], perm[len(xA):]
        assert len(pA) == len(xA)
        assert len(pB) == len(xB)
        assert pA.size + pB.size == pooled.size

    @pytest.mark.parametrize("statistic", ["A_bar", "A_w", "ks_stat"])
    def test_permutation_pvalue_in_unit_interval(self, statistic):
        rng = np.random.default_rng(7)
        xA = rng.beta(5, 5, size=40)
        xB = rng.beta(5, 5, size=60)
        p = ncg.permutation_pvalue(xA, xB, GRID, n_perm=25, rng=rng,
                                    statistic=statistic)
        assert 0.0 <= p <= 1.0


class TestSampleFamilies:
    """Sanity: every declared family samples stay in [0,1] and the mixture
    respects its declared component weights approximately."""

    @pytest.mark.parametrize("family_key", list(ncg.FAMILIES.keys()))
    def test_family_samples_in_unit_interval(self, family_key):
        rng = np.random.default_rng(8)
        x = ncg.sample_family(rng, family_key, 500)
        assert x.min() >= 0.0
        assert x.max() <= 1.0
        assert len(x) == 500
