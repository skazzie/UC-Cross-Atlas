"""Unit tests for parse_log.py — LDSC `.log` parser + F10 acceptance."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from parse_log import in_band_for_f10, parse  # noqa: E402


_TYPICAL_LOG = """\
*********************************************************************
* LD Score Regression (LDSC)
* Version 1.0.1
*********************************************************************
Beginning analysis at Sat Jun  7 19:30:00 2026

Reading summary statistics from /path/to/uc_delange.sumstats.gz ...
Read summary statistics for 1156421 SNPs.
After merging with regression SNP LD, 1146421 SNPs remain.

Total Observed scale h2: 0.0742 (0.0034)
Lambda GC: 1.1721
Mean Chi^2: 1.2746
Intercept: 1.0231 (0.0091)
Ratio: 0.084 (0.033)

Analysis finished at Sat Jun  7 19:30:25 2026
Total time elapsed: 25.3s
"""


def test_parse_extracts_all_fields(tmp_path: Path) -> None:
    log = tmp_path / "uc_delange.log"
    log.write_text(_TYPICAL_LOG)
    result = parse(log)
    assert result["h2"] == pytest.approx(0.0742)
    assert result["h2_se"] == pytest.approx(0.0034)
    assert result["intercept"] == pytest.approx(1.0231)
    assert result["intercept_se"] == pytest.approx(0.0091)
    assert result["ratio"] == pytest.approx(0.084)
    assert result["ratio_se"] == pytest.approx(0.033)
    assert result["lambda_gc"] == pytest.approx(1.1721)
    assert result["mean_chisq"] == pytest.approx(1.2746)
    assert result["status"] == "ok"


def test_parse_missing_file_returns_missing_file_status(tmp_path: Path) -> None:
    result = parse(tmp_path / "does_not_exist.log")
    assert result["status"] == "missing-file"
    assert result["intercept"] is None


def test_parse_partial_log_returns_partial_status(tmp_path: Path) -> None:
    log = tmp_path / "partial.log"
    log.write_text(
        "Total Observed scale h2: 0.05 (0.01)\n"
        "Intercept: 1.05 (0.02)\n"
        # No Ratio line.
    )
    result = parse(log)
    assert result["intercept"] == pytest.approx(1.05)
    assert result["ratio"] is None
    assert result["status"] == "partial"


def test_parse_empty_log_returns_missing(tmp_path: Path) -> None:
    log = tmp_path / "empty.log"
    log.write_text("LDSC crashed before computing anything.\n")
    result = parse(log)
    assert result["status"] == "missing"


def test_parse_handles_liability_scale_h2(tmp_path: Path) -> None:
    log = tmp_path / "liability.log"
    log.write_text(
        "Total Liability scale h2: 0.123 (0.045)\n"
        "Intercept: 1.0 (0.01)\n"
        "Ratio: 0.0 (0.01)\n"
    )
    result = parse(log)
    assert result["h2"] == pytest.approx(0.123)
    assert result["h2_se"] == pytest.approx(0.045)
    assert result["status"] == "ok"


def test_in_band_intercept_below_threshold_passes() -> None:
    assert in_band_for_f10(intercept=1.05, ratio=None) is True


def test_in_band_intercept_at_threshold_passes() -> None:
    assert in_band_for_f10(intercept=1.10, ratio=None) is True


def test_in_band_intercept_above_threshold_fails_without_ratio() -> None:
    assert in_band_for_f10(intercept=1.20, ratio=None) is False


def test_in_band_ratio_rescues_above_threshold_intercept() -> None:
    # Elevated intercept but ratio low → still in band (the ratio gate
    # passes when most of the inflation IS polygenic, regardless of
    # the raw intercept value).
    assert in_band_for_f10(intercept=1.30, ratio=0.10) is True


def test_in_band_both_none_returns_false() -> None:
    # F10 cannot be claimed without evidence.
    assert in_band_for_f10(intercept=None, ratio=None) is False


def test_parse_handles_realistic_scientific_notation(tmp_path: Path) -> None:
    log = tmp_path / "sci.log"
    log.write_text(
        "Intercept: 1.0231 (0.0091)\n"
        "Ratio: 8.4e-2 (3.3e-2)\n"
    )
    result = parse(log)
    assert result["ratio"] == pytest.approx(0.084)
    assert result["ratio_se"] == pytest.approx(0.033)
