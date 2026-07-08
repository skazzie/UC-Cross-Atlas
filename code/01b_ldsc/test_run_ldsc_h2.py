"""Tests for run_ldsc_h2.py — command construction + intercept TSV writer."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_ldsc_h2 import build_ldsc_command, write_intercept_tsv  # noqa: E402


def test_build_ldsc_command_has_all_required_flags(tmp_path: Path) -> None:
    cmd = build_ldsc_command(
        ldsc_dir=tmp_path / "ldsc",
        sumstats=tmp_path / "uc_delange.sumstats.gz",
        ld_dir=tmp_path / "eur_w_ld_chr",
        out_prefix=tmp_path / "results" / "uc_delange",
        python_bin="python2.7",
    )
    assert cmd[0] == "python2.7"
    assert cmd[1].endswith("ldsc.py")
    assert "--h2" in cmd
    assert "--ref-ld-chr" in cmd
    assert "--w-ld-chr" in cmd
    assert "--out" in cmd
    out_idx = cmd.index("--out")
    assert cmd[out_idx + 1].endswith("uc_delange")


def test_build_ldsc_command_ld_dir_has_trailing_slash(tmp_path: Path) -> None:
    """LDSC requires the LD-score directory to end in `/` so it can
    construct per-chromosome paths like `eur_w_ld_chr/1.l2.ldscore.gz`."""
    cmd = build_ldsc_command(
        ldsc_dir=tmp_path,
        sumstats=tmp_path / "x.gz",
        ld_dir=tmp_path / "eur_w_ld_chr",
        out_prefix=tmp_path / "y",
    )
    for flag in ("--ref-ld-chr", "--w-ld-chr"):
        idx = cmd.index(flag)
        assert cmd[idx + 1].endswith("/"), f"{flag} value must end with '/'"


def test_write_intercept_tsv_ok_status_with_full_fields(tmp_path: Path) -> None:
    parsed = {
        "h2": 0.07, "h2_se": 0.003,
        "intercept": 1.02, "intercept_se": 0.01,
        "ratio": 0.08, "ratio_se": 0.03,
        "lambda_gc": 1.17, "mean_chisq": 1.27,
        "status": "ok",
    }
    path = write_intercept_tsv(tmp_path / "uc_delange", "uc_delange", parsed)
    assert path.exists()
    df = pd.read_csv(path, sep="\t")
    assert len(df) == 1
    row = df.iloc[0]
    assert row["gwas"] == "uc_delange"
    assert row["intercept"] == pytest.approx(1.02)
    assert row["status"] == "ok"
    assert bool(row["in_band"]) is True  # intercept 1.02 <= 1.10


def test_write_intercept_tsv_partial_status_writes_nulls(tmp_path: Path) -> None:
    parsed = {
        "h2": None, "h2_se": None,
        "intercept": 1.05, "intercept_se": 0.01,
        "ratio": None, "ratio_se": None,
        "lambda_gc": None, "mean_chisq": None,
        "status": "partial",
    }
    path = write_intercept_tsv(tmp_path / "scz", "scz_trubetskoy", parsed)
    df = pd.read_csv(path, sep="\t")
    assert df.iloc[0]["status"] == "partial"
    # in_band still computable from the intercept alone.
    assert bool(df.iloc[0]["in_band"]) is True


def test_write_intercept_tsv_missing_status_in_band_false(tmp_path: Path) -> None:
    parsed = {k: None for k in (
        "h2", "h2_se", "intercept", "intercept_se",
        "ratio", "ratio_se", "lambda_gc", "mean_chisq",
    )}
    parsed["status"] = "missing"
    path = write_intercept_tsv(tmp_path / "fail", "broken_gwas", parsed)
    df = pd.read_csv(path, sep="\t")
    assert df.iloc[0]["status"] == "missing"
    assert bool(df.iloc[0]["in_band"]) is False
