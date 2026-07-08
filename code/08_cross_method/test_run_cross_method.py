"""End-to-end test for run_cross_method.py against synthetic fixtures.

The driver loads scDRS group + seismicGWAS results for the UC trio
of atlases × 2 GWAS × 2 tiers, computes Spearman / Jaccard / kappa,
and writes a long-format TSV. This test:

1. Builds fixtures in tmp_path matching `_shared/result_loading.py`'s
   expected layout.
2. Spawns the driver as a subprocess (so its own sys.path setup runs).
3. Reads the output TSV and asserts the schema, the bootstrap CI
   determinism (seed locked at 42 per DECISIONS), and at least one
   sensible value.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


_HERE = Path(__file__).resolve().parent
_REPO = _HERE.parents[1]
_DRIVER = _HERE / "run_cross_method.py"


def _make_scdrs_fixture(scdrs_dir: Path, atlas: str, gwas: str, tier: str,
                        rho_target: float, n_types: int = 8) -> None:
    """Build a scDRS group file where z_mean correlates with a target
    coefficient pattern, so the downstream concordance has a known
    direction (not exact rho — bootstrap CI bands will still vary)."""
    rng = np.random.default_rng(seed=hash((atlas, gwas, tier)) % (2**31))
    cell_types = [f"ct{i:02d}" for i in range(n_types)]
    truth = rng.normal(size=n_types)
    noise = rng.normal(scale=np.sqrt(max(1e-9, 1 - rho_target**2)), size=n_types)
    z = rho_target * truth + noise
    df = pd.DataFrame({
        "group": cell_types,
        "assoc_mcz": z,
        "assoc_mcp": rng.uniform(low=1e-9, high=1.0, size=n_types),
        "assoc_mcq": rng.uniform(low=1e-9, high=1.0, size=n_types),
        "n_cell": rng.integers(low=80, high=5000, size=n_types),
    })
    path = scdrs_dir / f"{atlas}_{gwas}" / f"cell_type_{tier}" \
        / f"{atlas}_{gwas}.scdrs_group"
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def _make_seismic_fixture(seismic_dir: Path, atlas: str, gwas: str, tier: str,
                          rho_target: float, n_types: int = 8) -> None:
    rng = np.random.default_rng(seed=(hash((atlas, gwas, tier)) + 1) % (2**31))
    cell_types = [f"ct{i:02d}" for i in range(n_types)]
    truth = rng.normal(size=n_types)
    noise = rng.normal(scale=np.sqrt(max(1e-9, 1 - rho_target**2)), size=n_types)
    coef = rho_target * truth + noise
    df = pd.DataFrame({
        "cell_type": cell_types,
        "coefficient": coef,
        "se": rng.uniform(low=0.05, high=0.20, size=n_types),
        "pvalue": rng.uniform(low=1e-9, high=1.0, size=n_types),
        "n_cells": rng.integers(low=80, high=5000, size=n_types),
    })
    path = seismic_dir / f"{atlas}_{gwas}_{tier}.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def test_run_cross_method_e2e_writes_expected_schema(tmp_path: Path) -> None:
    scdrs_dir = tmp_path / "results" / "scdrs"
    seis_dir = tmp_path / "results" / "seismic"
    # One atlas × one gwas × one tier — enough for the driver to write
    # one row and exercise every code path.
    _make_scdrs_fixture(scdrs_dir, "smillie", "delange", "broad",
                        rho_target=0.6)
    _make_seismic_fixture(seis_dir, "smillie", "delange", "broad",
                          rho_target=0.6)

    out = tmp_path / "cross_method.tsv"
    rc = subprocess.call([
        sys.executable, str(_DRIVER),
        "--atlases", "smillie",
        "--gwas", "delange",
        "--tiers", "broad",
        "--scdrs-dir", str(scdrs_dir),
        "--seismic-dir", str(seis_dir),
        "--out", str(out),
        "--bootstrap-n", "200",  # cheap; correctness test, not power
    ], cwd=str(_REPO))
    assert rc == 0, f"driver exited {rc}"
    assert out.exists()

    df = pd.read_csv(out, sep="\t")
    assert len(df) == 1
    expected_cols = {
        "atlas", "gwas", "tier", "spearman_rho", "ci_lo", "ci_hi",
        "jaccard_top5", "jaccard_top10",
        "kappa", "kappa_threshold", "kappa_saturated",
        "n_sig_scdrs", "n_sig_seismic", "n_common",
        "tool_version", "git_sha",
    }
    assert expected_cols.issubset(set(df.columns))
    row = df.iloc[0]
    assert row["atlas"] == "smillie"
    assert row["gwas"] == "delange"
    assert row["tier"] == "broad"
    assert -1.0 <= row["spearman_rho"] <= 1.0
    assert row["ci_lo"] <= row["spearman_rho"] <= row["ci_hi"]
    assert row["n_common"] >= 3


def test_run_cross_method_skips_when_only_some_inputs_present(tmp_path: Path) -> None:
    """Missing input for a (atlas, gwas, tier) should warn-and-skip,
    not crash, as long as at least one combo produces output."""
    scdrs_dir = tmp_path / "results" / "scdrs"
    seis_dir = tmp_path / "results" / "seismic"
    # smillie/delange/broad has both
    _make_scdrs_fixture(scdrs_dir, "smillie", "delange", "broad", 0.5)
    _make_seismic_fixture(seis_dir, "smillie", "delange", "broad", 0.5)
    # garrido_trigo/delange/broad has only scDRS — seismic missing
    _make_scdrs_fixture(scdrs_dir, "garrido_trigo", "delange", "broad", 0.5)

    out = tmp_path / "cross_method_skip.tsv"
    rc = subprocess.call([
        sys.executable, str(_DRIVER),
        "--atlases", "smillie", "garrido_trigo",
        "--gwas", "delange",
        "--tiers", "broad",
        "--scdrs-dir", str(scdrs_dir),
        "--seismic-dir", str(seis_dir),
        "--out", str(out),
        "--bootstrap-n", "100",
    ], cwd=str(_REPO))
    assert rc == 0
    df = pd.read_csv(out, sep="\t")
    # Only smillie completed.
    assert len(df) == 1
    assert df.iloc[0]["atlas"] == "smillie"


def test_run_cross_method_seed_determinism(tmp_path: Path) -> None:
    """Two runs with the same fixtures + locked seed must produce
    identical CI bands (the bootstrap is the only stochastic step)."""
    scdrs_dir = tmp_path / "results" / "scdrs"
    seis_dir = tmp_path / "results" / "seismic"
    _make_scdrs_fixture(scdrs_dir, "smillie", "delange", "broad", 0.5)
    _make_seismic_fixture(seis_dir, "smillie", "delange", "broad", 0.5)

    out1 = tmp_path / "run1.tsv"
    out2 = tmp_path / "run2.tsv"
    for out in (out1, out2):
        rc = subprocess.call([
            sys.executable, str(_DRIVER),
            "--atlases", "smillie", "--gwas", "delange",
            "--tiers", "broad",
            "--scdrs-dir", str(scdrs_dir),
            "--seismic-dir", str(seis_dir),
            "--out", str(out),
            "--bootstrap-n", "200",
        ], cwd=str(_REPO))
        assert rc == 0

    df1 = pd.read_csv(out1, sep="\t")
    df2 = pd.read_csv(out2, sep="\t")
    for col in ("spearman_rho", "ci_lo", "ci_hi"):
        assert df1.iloc[0][col] == df2.iloc[0][col], \
            f"{col} not deterministic under locked seed"
