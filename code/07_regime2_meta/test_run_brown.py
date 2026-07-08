"""Tests for run_brown.py — Brown's method math + E2E with fixtures.

run_brown.py's helpers are defined inside main(), so we re-derive them
here for the math tests (kept in sync by docstring reference + a
sanity check that the driver still works against fixtures). The E2E
test exercises the full pipeline against synthetic per-atlas results
+ null tensors.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import chi2


_HERE = Path(__file__).resolve().parent
_REPO = _HERE.parents[1]
_DRIVER = _HERE / "run_brown.py"


def _brown_cov_from_corr(r: float) -> float:
    """Mirror of the cubic in run_brown.py (Brown 1975 polynomial)."""
    return 3.263 * r + 0.710 * r * r + 0.027 * r ** 3


def _brown_combine(pvals: np.ndarray, corr: np.ndarray) -> float:
    pvals = np.clip(pvals, 1e-300, 1.0)
    k = len(pvals)
    if k == 1:
        return float(pvals[0])
    T = -2.0 * np.sum(np.log(pvals))
    E = 2.0 * k
    var = 4.0 * k
    for i in range(k):
        for j in range(i + 1, k):
            var += 2.0 * _brown_cov_from_corr(float(corr[i, j]))
    c = var / (2.0 * E)
    f = 2.0 * E * E / var
    return float(chi2.sf(T / c, df=f))


def test_brown_cov_from_corr_zero_correlation_is_zero() -> None:
    """At r=0, the polynomial returns 0 — independent test statistics
    contribute no extra covariance term."""
    assert _brown_cov_from_corr(0.0) == pytest.approx(0.0)


def test_brown_cov_from_corr_perfect_correlation_matches_brown_1975() -> None:
    """Brown 1975: at r=1.0 the cov term should approach 4 — the
    maximum possible covariance of -2 ln p_i values when the two
    statistics are perfectly correlated. Polynomial gives 4.0."""
    # 3.263 + 0.710 + 0.027 = 4.000
    assert _brown_cov_from_corr(1.0) == pytest.approx(4.0, abs=1e-9)


def test_brown_combine_single_pvalue_passthrough() -> None:
    """For k=1 the combined p-value is the input p-value itself —
    Brown's method has no meaning."""
    combined = _brown_combine(np.array([0.04]), np.array([[1.0]]))
    assert combined == pytest.approx(0.04)


def test_brown_combine_independent_matches_fisher() -> None:
    """With zero correlation, Brown reduces to Fisher's combined p:
    T = -2 sum log p ~ chi^2(2k). Fisher's combined p = chi2.sf(T, 2k)."""
    pvals = np.array([0.01, 0.05, 0.10])
    corr = np.eye(3)
    combined = _brown_combine(pvals, corr)
    # Fisher's exact form
    T = -2.0 * np.sum(np.log(pvals))
    expected = float(chi2.sf(T, df=2 * len(pvals)))
    assert combined == pytest.approx(expected, abs=1e-9)


def test_brown_combine_with_correlation_inflates_p_relative_to_fisher() -> None:
    """When the underlying statistics are correlated, the effective k
    shrinks and the combined p-value should be LARGER than Fisher's
    independence-assuming p-value."""
    pvals = np.array([0.01, 0.05, 0.10])
    corr_zero = np.eye(3)
    corr_strong = np.ones((3, 3)) * 0.7
    np.fill_diagonal(corr_strong, 1.0)
    fisher = _brown_combine(pvals, corr_zero)
    brown_correlated = _brown_combine(pvals, corr_strong)
    assert brown_correlated > fisher


def test_brown_combine_handles_tiny_pvalue_without_overflow() -> None:
    """A clipped 1e-300 input should not throw — clip protects log(0)."""
    pvals = np.array([1e-308, 0.5])
    corr = np.eye(2)
    combined = _brown_combine(pvals, corr)
    assert np.isfinite(combined)
    assert combined < 0.5  # very small p_1 should still drive combined low


# ---- E2E with synthetic null tensors --------------------------------


def _make_regime1_dir_scdrs(scdrs_dir: Path, atlas: str, gwas: str,
                            tier: str, cell_types: list[str],
                            pvals: list[float]) -> None:
    df = pd.DataFrame({
        "group": cell_types,
        "assoc_mcz": [0.0] * len(cell_types),
        "assoc_mcp": pvals,
        "assoc_mcq": pvals,
        "n_cell": [1000] * len(cell_types),
    })
    path = scdrs_dir / f"{atlas}_{gwas}" / f"cell_type_{tier}" \
        / f"{atlas}_{gwas}.scdrs_group"
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def _make_null_npz(null_dir: Path, atlas: str, tier: str,
                   cell_types: list[str], n_draws: int = 200,
                   seed: int = 0) -> None:
    rng = np.random.default_rng(seed=seed)
    nulls = rng.normal(size=(len(cell_types), n_draws)).astype(np.float32)
    null_dir.mkdir(parents=True, exist_ok=True)
    np.savez(
        null_dir / f"{atlas}_{tier}_nulls.npz",
        cell_types=np.array(cell_types, dtype=object),
        nulls=nulls,
    )


def test_run_brown_e2e_writes_combined_pvalues(tmp_path: Path) -> None:
    """Full driver run against synthetic regime-1 + null fixtures."""
    scdrs_dir = tmp_path / "results" / "scdrs"
    null_dir = tmp_path / "results" / "null_draws"
    cell_types = ["T cell", "B cell", "fibroblast"]
    # Per-atlas p-values: vary a bit so the combination is non-trivial.
    _make_regime1_dir_scdrs(scdrs_dir, "smillie", "delange", "broad",
                            cell_types, [0.01, 0.5, 0.05])
    _make_regime1_dir_scdrs(scdrs_dir, "garrido_trigo", "delange", "broad",
                            cell_types, [0.02, 0.4, 0.07])
    _make_regime1_dir_scdrs(scdrs_dir, "taurus", "delange", "broad",
                            cell_types, [0.005, 0.6, 0.04])
    for atlas, seed in [("smillie", 1), ("garrido_trigo", 2), ("taurus", 3)]:
        _make_null_npz(null_dir, atlas, "broad", cell_types, seed=seed)

    out = tmp_path / "brown_results.tsv"
    rc = subprocess.call([
        sys.executable, str(_DRIVER),
        "--method", "scdrs",
        "--tier", "broad",
        "--gwas", "delange",
        "--regime1-dir", str(scdrs_dir),
        "--null-draws-dir", str(null_dir),
        "--atlases", "smillie", "garrido_trigo", "taurus",
        "--out", str(out),
    ], cwd=str(_REPO))
    assert rc == 0, f"driver exited {rc}"
    assert out.exists()
    df = pd.read_csv(out, sep="\t")
    assert set(df["cell_type"]) == set(cell_types)
    for col in ("combined_pval", "n_atlases_combined",
                "correlation_fallback", "atlases"):
        assert col in df.columns
    # All three atlases combined for every cell type in this fixture.
    assert (df["n_atlases_combined"] == 3).all()
    # P-values should be in (0, 1].
    assert (df["combined_pval"] > 0).all()
    assert (df["combined_pval"] <= 1).all()
    # T cell had the smallest per-atlas p-values; combined should be smaller
    # than for B cell which had p ~ 0.5 everywhere.
    t_p = float(df.loc[df["cell_type"] == "T cell", "combined_pval"].iloc[0])
    b_p = float(df.loc[df["cell_type"] == "B cell", "combined_pval"].iloc[0])
    assert t_p < b_p


def test_run_brown_handles_missing_null_with_fallback(tmp_path: Path) -> None:
    """If one atlas's null tensor is absent, the driver falls back to
    the median cross-atlas correlation and flags the row."""
    scdrs_dir = tmp_path / "results" / "scdrs"
    null_dir = tmp_path / "results" / "null_draws"
    cell_types = ["T cell"]
    for atlas in ("smillie", "garrido_trigo", "taurus"):
        _make_regime1_dir_scdrs(scdrs_dir, atlas, "delange", "broad",
                                cell_types, [0.05])
    # Provide nulls for only two of three atlases.
    _make_null_npz(null_dir, "smillie", "broad", cell_types, seed=1)
    _make_null_npz(null_dir, "garrido_trigo", "broad", cell_types, seed=2)

    out = tmp_path / "brown_fallback.tsv"
    rc = subprocess.call([
        sys.executable, str(_DRIVER),
        "--method", "scdrs", "--tier", "broad", "--gwas", "delange",
        "--regime1-dir", str(scdrs_dir),
        "--null-draws-dir", str(null_dir),
        "--atlases", "smillie", "garrido_trigo", "taurus",
        "--out", str(out),
    ], cwd=str(_REPO))
    assert rc == 0
    df = pd.read_csv(out, sep="\t")
    assert len(df) == 1
    # The fallback should be flagged.
    assert bool(df.iloc[0]["correlation_fallback"]) is True
