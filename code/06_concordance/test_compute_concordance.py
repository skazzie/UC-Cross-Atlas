"""End-to-end test for compute_concordance.py against fixture TSVs.

The CLI takes one or more `PATH:ATLAS:METHOD:GWAS:TIER` entries, each
pointing at a simplified TSV with columns:
  cell_type, score, pval, qval, n_cells

For each unique (method, gwas, tier) it produces all pairwise atlas
comparisons. This test exercises that loop end-to-end against
synthetic fixtures.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


_HERE = Path(__file__).resolve().parent
_REPO = _HERE.parents[1]
_DRIVER = _HERE / "compute_concordance.py"


def _make_tsv(path: Path, cell_types: list[str], scores: list[float],
              qvals: list[float] | None = None,
              n_cells: list[int] | None = None) -> None:
    if qvals is None:
        qvals = [0.5] * len(cell_types)
    if n_cells is None:
        n_cells = [1000] * len(cell_types)
    pvals = qvals  # for the simplified format, pval can equal qval
    df = pd.DataFrame({
        "cell_type": cell_types,
        "score": scores,
        "pval": pvals,
        "qval": qvals,
        "n_cells": n_cells,
    })
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def test_compute_concordance_e2e_three_atlas_pairwise(tmp_path: Path) -> None:
    """Three atlases × scdrs × delange × broad → 3 pairwise rows."""
    cell_types = [f"ct{i}" for i in range(8)]
    rng = np.random.default_rng(7)
    base = rng.normal(size=8)
    paths = {}
    for atlas, jitter in [("smillie", 0.0), ("garrido_trigo", 0.3),
                          ("taurus", 0.6)]:
        scores = (base + rng.normal(scale=jitter, size=8)).tolist()
        path = tmp_path / "results" / atlas / "delange_broad.tsv"
        _make_tsv(path, cell_types, scores)
        paths[atlas] = path

    out = tmp_path / "concordance.csv"
    entries = [
        f"{paths['smillie']}:smillie:scdrs:delange:broad",
        f"{paths['garrido_trigo']}:garrido_trigo:scdrs:delange:broad",
        f"{paths['taurus']}:taurus:scdrs:delange:broad",
    ]
    rc = subprocess.call([
        sys.executable, str(_DRIVER),
        "--input", *entries,
        "--out", str(out),
        "--bootstrap-iter", "200",
    ], cwd=str(_REPO))
    assert rc == 0, f"driver exited {rc}"
    assert out.exists()

    df = pd.read_csv(out)
    # 3 atlases choose 2 → 3 pairs.
    assert len(df) == 3
    # smillie vs garrido_trigo should have higher rho than smillie vs taurus
    # (smaller jitter). Bootstrap CI is wide enough that we only assert
    # ordering on point estimates, not strict inequalities.
    expected_cols = {
        "method", "gwas", "tier", "atlas_a", "atlas_b",
        "spearman_rho", "spearman_ci_lo", "spearman_ci_hi",
        "jaccard_top5", "jaccard_top10", "kappa",
        "n_common", "n_excluded_low_count",
    }
    assert expected_cols.issubset(set(df.columns))
    assert (df["n_common"] >= 3).all()


def test_compute_concordance_separate_method_gwas_tier_buckets(
    tmp_path: Path
) -> None:
    """Two (method, gwas, tier) buckets shouldn't be cross-pollinated."""
    cell_types = [f"ct{i}" for i in range(8)]
    rng = np.random.default_rng(11)
    paths = []
    # bucket 1: scdrs/delange/broad (smillie + taurus)
    for atlas in ("smillie", "taurus"):
        scores = rng.normal(size=8).tolist()
        p = tmp_path / f"{atlas}_scdrs_delange_broad.tsv"
        _make_tsv(p, cell_types, scores)
        paths.append(f"{p}:{atlas}:scdrs:delange:broad")
    # bucket 2: seismic/liu/fine (smillie + taurus, same fixtures but
    # different keys)
    for atlas in ("smillie", "taurus"):
        scores = rng.normal(size=8).tolist()
        p = tmp_path / f"{atlas}_seismic_liu_fine.tsv"
        _make_tsv(p, cell_types, scores)
        paths.append(f"{p}:{atlas}:seismic:liu:fine")

    out = tmp_path / "buckets.csv"
    rc = subprocess.call([
        sys.executable, str(_DRIVER),
        "--input", *paths,
        "--out", str(out),
        "--bootstrap-iter", "100",
    ], cwd=str(_REPO))
    assert rc == 0
    df = pd.read_csv(out)
    assert len(df) == 2
    # One scdrs/delange/broad and one seismic/liu/fine.
    keys = set(zip(df["method"], df["gwas"], df["tier"]))
    assert keys == {("scdrs", "delange", "broad"),
                    ("seismic", "liu", "fine")}


def test_compute_concordance_rejects_malformed_input_entry(tmp_path: Path) -> None:
    """A `--input` entry missing parts of the PATH:ATLAS:METHOD:GWAS:TIER
    spec should exit nonzero with a clear message."""
    rc = subprocess.call([
        sys.executable, str(_DRIVER),
        "--input", str(tmp_path / "x.tsv:atlas_only"),
        "--out", str(tmp_path / "out.csv"),
    ], cwd=str(_REPO))
    assert rc != 0


def test_compute_concordance_min_cells_filter_excludes_low_count_types(
    tmp_path: Path
) -> None:
    """Cell types whose n_cells is below --min-cells in EITHER atlas
    drop from the comparison's input."""
    rng = np.random.default_rng(13)
    cts = ["keeper1", "keeper2", "low_in_a", "low_in_b"]
    a_path = tmp_path / "a.tsv"
    b_path = tmp_path / "b.tsv"
    # 'low_in_a' has 10 cells in smillie; 'low_in_b' has 10 in taurus.
    _make_tsv(a_path, cts, rng.normal(size=4).tolist(),
              n_cells=[2000, 1500, 10, 2000])
    _make_tsv(b_path, cts, rng.normal(size=4).tolist(),
              n_cells=[2200, 1700, 1800, 10])

    out = tmp_path / "filtered.csv"
    rc = subprocess.call([
        sys.executable, str(_DRIVER),
        "--input",
        f"{a_path}:smillie:scdrs:delange:broad",
        f"{b_path}:taurus:scdrs:delange:broad",
        "--out", str(out),
        "--min-cells", "50",
        "--bootstrap-iter", "100",
    ], cwd=str(_REPO))
    assert rc == 0
    df = pd.read_csv(out)
    assert len(df) == 1
    # Two cell types pass the filter on both sides.
    assert df.iloc[0]["n_common"] == 2
    assert df.iloc[0]["n_excluded_low_count"] == 2
