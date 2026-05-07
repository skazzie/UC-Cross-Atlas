"""Tests for concordance metrics on synthetic data.

Headline convention is larger-is-stronger (cell-type-level Z-scores or
regression coefficients). Tests synthesize Z-scores accordingly.
"""

import math

import numpy as np
import pytest

from metrics import (
    DEFAULT_SEED,
    bootstrap_spearman_ci,
    concordance,
    fdr_concordance,
    spearman,
    top_k_jaccard,
)


CELL_TYPES = ["T", "B", "Plasma", "Mono", "DC", "Mast", "Neut",
              "NK", "ILC", "Enterocyte", "Goblet", "Stem", "Fib", "Endo", "SMC"]


def perfect_zscores():
    """Z-scores from 1.0 (weakest) to 15.0 (strongest)."""
    return {ct: float(i + 1) for i, ct in enumerate(CELL_TYPES)}


def perfect_qvals(n_sig=8):
    """First n_sig cell types significant at FDR<0.05; rest non-significant."""
    return {ct: 0.001 if i < n_sig else 0.5 for i, ct in enumerate(CELL_TYPES)}


def big_cell_counts():
    """All cell types pass the min-cell-count threshold."""
    return {ct: 1000 for ct in CELL_TYPES}


def test_spearman_perfect_agreement():
    z = perfect_zscores()
    rho, _, n, _ = spearman(z, z)
    assert n == len(CELL_TYPES)
    assert rho == pytest.approx(1.0)


def test_spearman_perfect_disagreement():
    z_a = perfect_zscores()
    z_b = {ct: z_a[CELL_TYPES[len(CELL_TYPES) - 1 - i]] for i, ct in enumerate(CELL_TYPES)}
    rho, _, _, _ = spearman(z_a, z_b)
    assert rho == pytest.approx(-1.0)


def test_spearman_drops_non_shared_cell_types():
    z_a = perfect_zscores()
    z_b = {ct: z_a[ct] for ct in CELL_TYPES[:10]}
    z_b["Extra"] = 0.5
    rho, _, n, _ = spearman(z_a, z_b)
    assert n == 10
    assert rho == pytest.approx(1.0)


def test_min_cells_filter_excludes_low_count():
    z_a = perfect_zscores()
    z_b = perfect_zscores()
    counts_a = big_cell_counts()
    counts_b = {ct: 1000 for ct in CELL_TYPES}
    counts_b["T"] = 10  # below threshold
    counts_b["B"] = 30  # below threshold
    rho, _, n, dropped = spearman(z_a, z_b, counts_a, counts_b, min_cells=50)
    assert n == len(CELL_TYPES) - 2
    assert dropped == 2
    assert rho == pytest.approx(1.0)


def test_jaccard_full_agreement_zscores():
    z = perfect_zscores()
    assert top_k_jaccard(z, z, k=5, larger_is_stronger=True) == pytest.approx(1.0)


def test_jaccard_partial_zscores():
    # Top-5 in atlas_a is f,g,h,i,j (Z = 5..1 reversed: actually big values)
    # Use explicit small example.
    z_a = {"a": 10, "b": 9, "c": 8, "d": 7, "e": 6,
           "f": 1, "g": 2, "h": 3, "i": 4, "j": 5}
    z_b = {"a": 10, "b": 9, "c": 8, "f": 7, "g": 6,
           "d": 1, "e": 2, "h": 3, "i": 4, "j": 5}
    # Top-5 a: a,b,c,d,e ; top-5 b: a,b,c,f,g ; intersection {a,b,c}, union 7
    assert top_k_jaccard(z_a, z_b, k=5, larger_is_stronger=True) == pytest.approx(3 / 7)


def test_jaccard_disjoint():
    z_a = {f"x{i}": float(i) for i in range(10)}
    z_b = {f"x{i}": float(10 - i) for i in range(10)}
    assert top_k_jaccard(z_a, z_b, k=3, larger_is_stronger=True) == pytest.approx(0.0)


def test_jaccard_smaller_is_stronger_path():
    # If passing p-values instead, set larger_is_stronger=False
    p_a = {"a": 0.001, "b": 0.002, "c": 0.003, "d": 0.004, "e": 0.005,
           "f": 0.10, "g": 0.20, "h": 0.30, "i": 0.40, "j": 0.50}
    p_b = {"a": 0.001, "b": 0.002, "c": 0.003, "f": 0.004, "g": 0.005,
           "d": 0.10, "e": 0.20, "h": 0.30, "i": 0.40, "j": 0.50}
    assert top_k_jaccard(p_a, p_b, k=5, larger_is_stronger=False) == pytest.approx(3 / 7)


def test_kappa_perfect_agreement_mixed_marginals():
    q = {"a": 0.001, "b": 0.001, "c": 0.5, "d": 0.5, "e": 0.5}
    kappa, n_a, n_b, _ = fdr_concordance(q, q)
    assert kappa == pytest.approx(1.0)
    assert n_a == 2 and n_b == 2


def test_kappa_chance_agreement():
    rng = np.random.default_rng(42)
    cells = [f"ct{i}" for i in range(200)]
    q_a = {c: float(rng.random()) for c in cells}
    q_b = {c: float(rng.random()) for c in cells}
    kappa, _, _, _ = fdr_concordance(q_a, q_b, threshold=0.5)
    assert abs(kappa) < 0.2


def test_kappa_all_significant_returns_nan():
    q = {ct: 0.001 for ct in CELL_TYPES}
    kappa, n_a, n_b, _ = fdr_concordance(q, q)
    assert math.isnan(kappa)
    assert n_a == len(CELL_TYPES) and n_b == len(CELL_TYPES)


def test_kappa_all_nonsignificant_returns_nan():
    q = {ct: 0.9 for ct in CELL_TYPES}
    kappa, n_a, n_b, _ = fdr_concordance(q, q)
    assert math.isnan(kappa)
    assert n_a == 0 and n_b == 0


def test_concordance_dataclass_aggregates_broad_tier():
    z = perfect_zscores()
    q = perfect_qvals(n_sig=8)
    res = concordance(z, z, q, q, is_fine_tier=False)
    assert res.spearman_rho == pytest.approx(1.0)
    assert res.jaccard_top5 == pytest.approx(1.0)
    assert res.jaccard_top10 == pytest.approx(1.0)
    assert math.isnan(res.jaccard_top20)  # not computed at broad tier
    assert res.kappa == pytest.approx(1.0)
    assert res.kappa_threshold == 0.05
    assert not res.kappa_threshold_used_due_to_saturation
    assert res.n_sig_a == 8 and res.n_sig_b == 8
    assert res.n_common == len(CELL_TYPES)


def test_concordance_dataclass_includes_jaccard20_at_fine_tier():
    keys = [f"ct{i}" for i in range(30)]
    z = {k: float(30 - i) for i, k in enumerate(keys)}
    q = {k: 0.001 if i < 10 else 0.5 for i, k in enumerate(keys)}
    res = concordance(z, z, q, q, is_fine_tier=True)
    assert res.jaccard_top20 == pytest.approx(1.0)


def test_kappa_saturation_contingency_promotes_to_fdr_001():
    # Construct a scenario where >=80% pass FDR<0.05 in both atlases
    keys = [f"ct{i}" for i in range(20)]
    z = {k: float(i) for i, k in enumerate(keys)}
    # 18/20 = 90% saturation at FDR<0.05; but only 5 each pass FDR<0.01
    q_a = {k: 0.001 if i < 5 else (0.04 if i < 18 else 0.5) for i, k in enumerate(keys)}
    q_b = {k: 0.001 if i < 5 else (0.04 if i < 18 else 0.5) for i, k in enumerate(keys)}
    res = concordance(z, z, q_a, q_b)
    assert res.kappa_threshold_used_due_to_saturation
    assert res.kappa_threshold == 0.01
    assert res.n_sig_a == 5 and res.n_sig_b == 5


def test_bootstrap_ci_brackets_point_estimate_and_uses_locked_seed():
    rng = np.random.default_rng(0)
    n = 30
    keys = [f"ct{i}" for i in range(n)]
    base = rng.normal(size=n)
    z_a = dict(zip(keys, base))
    z_b = dict(zip(keys, base + rng.normal(scale=0.3, size=n)))
    rho, _, _, _ = spearman(z_a, z_b)
    lo, hi = bootstrap_spearman_ci(z_a, z_b, n_iter=500)
    # Default seed = 42 (locked)
    assert lo <= rho <= hi
    assert hi - lo < 1.0


def test_bootstrap_ci_is_deterministic_under_locked_seed():
    rng = np.random.default_rng(7)
    n = 25
    keys = [f"ct{i}" for i in range(n)]
    base = rng.normal(size=n)
    z_a = dict(zip(keys, base))
    z_b = dict(zip(keys, base + rng.normal(scale=0.3, size=n)))
    lo1, hi1 = bootstrap_spearman_ci(z_a, z_b, n_iter=200, seed=DEFAULT_SEED)
    lo2, hi2 = bootstrap_spearman_ci(z_a, z_b, n_iter=200, seed=DEFAULT_SEED)
    assert lo1 == lo2 and hi1 == hi2


def test_too_few_common_cell_types():
    z_a = {"a": 1.0}
    z_b = {"a": 1.0}
    rho, _, n, _ = spearman(z_a, z_b)
    assert math.isnan(rho)
    assert n == 1
