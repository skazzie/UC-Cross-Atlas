"""End-to-end test for run_cross_gwas.py against synthetic fixtures.

For each (atlas, method, tier), the driver compares cell-type-level
results under de Lange vs Liu. Same metrics as 08 but the axis is
GWAS instead of method.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


_HERE = Path(__file__).resolve().parent
_REPO = _HERE.parents[1]
_DRIVER = _HERE / "run_cross_gwas.py"


def _make_scdrs(scdrs_dir: Path, atlas: str, gwas: str, tier: str,
                rho_target: float, n_types: int = 8) -> None:
    rng = np.random.default_rng(seed=hash((atlas, gwas, tier)) % (2**31))
    truth = rng.normal(size=n_types)
    noise = rng.normal(scale=np.sqrt(max(1e-9, 1 - rho_target**2)), size=n_types)
    df = pd.DataFrame({
        "group": [f"ct{i:02d}" for i in range(n_types)],
        "assoc_mcz": rho_target * truth + noise,
        "assoc_mcp": rng.uniform(low=1e-9, high=1.0, size=n_types),
        "assoc_mcq": rng.uniform(low=1e-9, high=1.0, size=n_types),
        "n_cell": rng.integers(low=80, high=5000, size=n_types),
    })
    path = scdrs_dir / f"{atlas}_{gwas}" / f"cell_type_{tier}" \
        / f"{atlas}_{gwas}.scdrs_group"
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def _make_seismic(seismic_dir: Path, atlas: str, gwas: str, tier: str,
                  rho_target: float, n_types: int = 8) -> None:
    rng = np.random.default_rng(seed=(hash((atlas, gwas, tier)) + 7) % (2**31))
    truth = rng.normal(size=n_types)
    noise = rng.normal(scale=np.sqrt(max(1e-9, 1 - rho_target**2)), size=n_types)
    df = pd.DataFrame({
        "cell_type": [f"ct{i:02d}" for i in range(n_types)],
        "coefficient": rho_target * truth + noise,
        "se": rng.uniform(low=0.05, high=0.20, size=n_types),
        "pvalue": rng.uniform(low=1e-9, high=1.0, size=n_types),
        "n_cells": rng.integers(low=80, high=5000, size=n_types),
    })
    path = seismic_dir / f"{atlas}_{gwas}_{tier}.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def test_run_cross_gwas_e2e_writes_expected_schema(tmp_path: Path) -> None:
    scdrs_dir = tmp_path / "results" / "scdrs"
    seis_dir = tmp_path / "results" / "seismic"
    # Need both delange + liu present for one (atlas, method, tier).
    for gwas in ("delange", "liu"):
        _make_scdrs(scdrs_dir, "smillie", gwas, "broad", 0.6)
        _make_seismic(seis_dir, "smillie", gwas, "broad", 0.6)

    out = tmp_path / "cross_gwas.tsv"
    rc = subprocess.call([
        sys.executable, str(_DRIVER),
        "--atlases", "smillie",
        "--methods", "scdrs",
        "--tiers", "broad",
        "--scdrs-dir", str(scdrs_dir),
        "--seismic-dir", str(seis_dir),
        "--out", str(out),
        "--bootstrap-n", "200",
    ], cwd=str(_REPO))
    assert rc == 0
    df = pd.read_csv(out, sep="\t")
    assert len(df) == 1
    row = df.iloc[0]
    assert row["atlas"] == "smillie"
    assert row["method"] == "scdrs"
    assert row["tier"] == "broad"
    assert -1.0 <= row["spearman_rho"] <= 1.0
    assert "n_sig_delange" in df.columns
    assert "n_sig_liu" in df.columns


def test_run_cross_gwas_seismic_path_covers_fdr_bh(tmp_path: Path) -> None:
    """seismic outputs raw pvalue; driver must BH-adjust before
    calling the metrics with qvalues. Smoke test: it runs."""
    scdrs_dir = tmp_path / "scdrs"  # unused but path arg is required
    seis_dir = tmp_path / "seismic"
    for gwas in ("delange", "liu"):
        _make_seismic(seis_dir, "smillie", gwas, "broad", 0.5)

    out = tmp_path / "cross_gwas_seismic.tsv"
    rc = subprocess.call([
        sys.executable, str(_DRIVER),
        "--atlases", "smillie",
        "--methods", "seismic",
        "--tiers", "broad",
        "--scdrs-dir", str(scdrs_dir),
        "--seismic-dir", str(seis_dir),
        "--out", str(out),
        "--bootstrap-n", "200",
    ], cwd=str(_REPO))
    assert rc == 0
    df = pd.read_csv(out, sep="\t")
    assert len(df) == 1
    assert df.iloc[0]["method"] == "seismic"


def test_run_cross_gwas_exits_clean_when_no_combos_match(tmp_path: Path) -> None:
    """If nothing's on disk, the driver should exit nonzero with a
    clear message — never write an empty TSV."""
    out = tmp_path / "empty.tsv"
    rc = subprocess.call([
        sys.executable, str(_DRIVER),
        "--atlases", "smillie",
        "--methods", "scdrs",
        "--tiers", "broad",
        "--scdrs-dir", str(tmp_path / "missing_scdrs"),
        "--seismic-dir", str(tmp_path / "missing_seismic"),
        "--out", str(out),
    ], cwd=str(_REPO))
    assert rc != 0
    assert not out.exists()
