"""Tests for concordance metrics on synthetic data."""

import math

import numpy as np
import pytest

from metrics import (
    bootstrap_spearman_ci,
    concordance,
    fdr_concordance,
    spearman,
    top_k_jaccard,
)


CELL_TYPES = ["T", "B", "Plasma", "Mono", "DC", "Mast", "Neut",
              "NK", "ILC", "Enterocyte", "Goblet", "Stem", "Fib", "Endo", "SMC"]


def perfect_pvals():
    return {ct: 10 ** -(len(CELL_TYPES) - i) for i, ct in enumerate(CELL_TYPES)}


def test_spearman_perfect_agreement():
    p = perfect_pvals()
    rho, _, n = spearman(p, p)
    assert n == len(CELL_TYPES)
    assert rho == pytest.approx(1.0)


def test_spearman_perfect_disagreement():
    p_a = perfect_pvals()
    p_b = {ct: p_a[CELL_TYPES[len(CELL_TYPES) - 1 - i]] for i, ct in enumerate(CELL_TYPES)}
    rho, _, _ = spearman(p_a, p_b)
    assert rho == pytest.approx(-1.0)


def test_spearman_drops_non_shared_cell_types():
    p_a = perfect_pvals()
    p_b = {ct: p_a[ct] for ct in CELL_TYPES[:10]}
    p_b["Extra"] = 0.5
    rho, _, n = spearman(p_a, p_b)
    assert n == 10
    assert rho == pytest.approx(1.0)


def test_jaccard_full_agreement():
    p = perfect_pvals()
    assert top_k_jaccard(p, p, k=5) == pytest.approx(1.0)


def test_jaccard_partial():
    p_a = {"a": 0.001, "b": 0.002, "c": 0.003, "d": 0.004, "e": 0.005,
           "f": 0.10, "g": 0.20, "h": 0.30, "i": 0.40, "j": 0.50}
    p_b = {"a": 0.001, "b": 0.002, "c": 0.003, "f": 0.004, "g": 0.005,
           "d": 0.10, "e": 0.20, "h": 0.30, "i": 0.40, "j": 0.50}
    assert top_k_jaccard(p_a, p_b, k=5) == pytest.approx(3 / 7)


def test_jaccard_disjoint():
    p_a = {f"x{i}": i / 100 for i in range(10)}
    p_b = {f"x{i}": (10 - i) / 100 for i in range(10)}
    assert top_k_jaccard(p_a, p_b, k=3) == pytest.approx(0.0)


def test_kappa_perfect_agreement_mixed_marginals():
    q_a = {"a": 0.001, "b": 0.001, "c": 0.5, "d": 0.5, "e": 0.5}
    q_b = q_a
    kappa, n_a, n_b, _ = fdr_concordance(q_a, q_b)
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
    q_a = {ct: 0.001 for ct in CELL_TYPES}
    q_b = {ct: 0.001 for ct in CELL_TYPES}
    kappa, n_a, n_b, _ = fdr_concordance(q_a, q_b)
    assert math.isnan(kappa)
    assert n_a == len(CELL_TYPES) and n_b == len(CELL_TYPES)


def test_kappa_all_nonsignificant_returns_nan():
    q_a = {ct: 0.9 for ct in CELL_TYPES}
    q_b = {ct: 0.9 for ct in CELL_TYPES}
    kappa, n_a, n_b, _ = fdr_concordance(q_a, q_b)
    assert math.isnan(kappa)
    assert n_a == 0 and n_b == 0


def test_concordance_dataclass_aggregates():
    p = perfect_pvals()
    q = {ct: 0.001 if i < 8 else 0.5 for i, ct in enumerate(CELL_TYPES)}
    res = concordance(p, p, q, q)
    assert res.spearman_rho == pytest.approx(1.0)
    assert res.jaccard_top5 == pytest.approx(1.0)
    assert res.jaccard_top10 == pytest.approx(1.0)
    assert res.kappa == pytest.approx(1.0)
    assert res.n_sig_a == 8 and res.n_sig_b == 8
    assert res.n_common == len(CELL_TYPES)


def test_bootstrap_ci_brackets_point_estimate():
    rng = np.random.default_rng(0)
    n = 30
    keys = [f"ct{i}" for i in range(n)]
    base = rng.normal(size=n)
    p_a = dict(zip(keys, base))
    p_b = dict(zip(keys, base + rng.normal(scale=0.3, size=n)))
    rho, _, _ = spearman(p_a, p_b)
    lo, hi = bootstrap_spearman_ci(p_a, p_b, n_iter=500, seed=1)
    assert lo <= rho <= hi
    assert hi - lo < 1.0


def test_too_few_common_cell_types():
    p_a = {"a": 0.01}
    p_b = {"a": 0.01}
    rho, _, n = spearman(p_a, p_b)
    assert math.isnan(rho)
    assert n == 1
