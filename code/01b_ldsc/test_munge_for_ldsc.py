"""Tests for munge_for_ldsc.py — merge + command construction."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from munge_for_ldsc import build_ldsc_command, build_merged_sumstats  # noqa: E402


def _write_snp_loc(path: Path, rows: list[tuple]) -> None:
    pd.DataFrame(rows, columns=["SNP", "CHR", "BP"]).to_csv(
        path, sep="\t", header=False, index=False
    )


def _write_pval(path: Path, rows: list[tuple]) -> None:
    pd.DataFrame(rows, columns=["SNP", "P", "N"]).to_csv(
        path, sep="\t", index=False
    )


def test_merge_without_original_sumstats_drops_alleles(tmp_path: Path) -> None:
    snp_loc = tmp_path / "uc.snp.loc"
    pval = tmp_path / "uc.pval"
    _write_snp_loc(snp_loc, [
        ("rs1", 1, 100),
        ("rs2", 1, 200),
        ("rs3", 2, 300),
    ])
    _write_pval(pval, [
        ("rs1", 0.5, 36160),
        ("rs2", 1e-8, 36160),
        ("rs3", 0.9, 36160),
    ])
    merged = build_merged_sumstats(
        snp_loc, pval, original_sumstats=None,
        original_a1_col=None, original_a2_col=None, original_snp_col=None,
    )
    assert list(merged.columns) == ["SNP", "CHR", "BP", "P", "N"]
    assert len(merged) == 3
    assert merged["N"].iloc[0] == 36160


def test_merge_with_original_sumstats_joins_alleles(tmp_path: Path) -> None:
    snp_loc = tmp_path / "uc.snp.loc"
    pval = tmp_path / "uc.pval"
    orig = tmp_path / "uc_original.tsv"
    _write_snp_loc(snp_loc, [("rs1", 1, 100), ("rs2", 1, 200)])
    _write_pval(pval, [("rs1", 0.5, 36160), ("rs2", 1e-8, 36160)])
    pd.DataFrame({
        "hm_rsid": ["rs1", "rs2", "rs3"],
        "hm_other_allele": ["A", "G", "C"],
        "hm_effect_allele": ["T", "C", "A"],
    }).to_csv(orig, sep="\t", index=False)

    merged = build_merged_sumstats(
        snp_loc, pval, original_sumstats=orig,
        original_snp_col="hm_rsid",
        original_a1_col="hm_effect_allele",
        original_a2_col="hm_other_allele",
    )
    assert set(merged.columns) == {"SNP", "CHR", "BP", "P", "N", "A1", "A2"}
    assert len(merged) == 2  # rs3 was in orig but not in our pval; inner join drops it
    a1_for_rs1 = merged.loc[merged["SNP"] == "rs1", "A1"].iloc[0]
    assert a1_for_rs1 == "T"


def test_merge_pval_missing_n_column_errors(tmp_path: Path) -> None:
    snp_loc = tmp_path / "uc.snp.loc"
    pval = tmp_path / "uc.pval"
    _write_snp_loc(snp_loc, [("rs1", 1, 100)])
    pd.DataFrame({"SNP": ["rs1"], "P": [0.5]}).to_csv(
        pval, sep="\t", index=False
    )
    with pytest.raises(SystemExit, match="no N column"):
        build_merged_sumstats(
            snp_loc, pval, original_sumstats=None,
            original_a1_col=None, original_a2_col=None, original_snp_col=None,
        )


def test_build_ldsc_command_has_alleles_flag_when_a1a2_present(tmp_path: Path) -> None:
    cmd = build_ldsc_command(
        ldsc_dir=tmp_path / "ldsc",
        snplist=tmp_path / "w_hm3.snplist",
        in_tsv=tmp_path / "merged.tsv",
        out_prefix=tmp_path / "out" / "uc_delange",
        has_alleles=True,
    )
    assert "--a1" in cmd and "--a2" in cmd
    assert "--no-alleles" not in cmd


def test_build_ldsc_command_no_alleles_when_a1a2_absent(tmp_path: Path) -> None:
    cmd = build_ldsc_command(
        ldsc_dir=tmp_path / "ldsc",
        snplist=tmp_path / "w_hm3.snplist",
        in_tsv=tmp_path / "merged.tsv",
        out_prefix=tmp_path / "out" / "uc_delange",
        has_alleles=False,
    )
    assert "--no-alleles" in cmd
    assert "--a1" not in cmd


def test_build_ldsc_command_signed_sumstats_default(tmp_path: Path) -> None:
    cmd = build_ldsc_command(
        ldsc_dir=tmp_path, snplist=tmp_path / "snp",
        in_tsv=tmp_path / "m.tsv", out_prefix=tmp_path / "o",
        has_alleles=False,
    )
    idx = cmd.index("--signed-sumstats")
    assert cmd[idx + 1] == "P,0"


def test_build_ldsc_command_required_flags_present(tmp_path: Path) -> None:
    cmd = build_ldsc_command(
        ldsc_dir=tmp_path / "ldsc",
        snplist=tmp_path / "snplist",
        in_tsv=tmp_path / "in.tsv",
        out_prefix=tmp_path / "out",
        has_alleles=True,
    )
    for flag in ("--sumstats", "--out", "--merge-alleles",
                 "--signed-sumstats", "--snp", "--p", "--N-col"):
        assert flag in cmd, f"missing flag: {flag}"
