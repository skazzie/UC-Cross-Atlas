"""Unit tests for result_loading.py — schema normalization + filters."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

# Allow `from result_loading import ...` regardless of pytest CWD.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from result_loading import (  # noqa: E402
    load_scdrs_group,
    load_seismic_results,
    require_files,
    scdrs_group_path,
    seismic_path,
    shared_cell_types,
    to_lookup,
)


def _write_scdrs(scdrs_dir: Path, atlas: str, gwas: str, tier: str,
                 df: pd.DataFrame) -> Path:
    """Write a scDRS group file at the canonical layout."""
    path = scdrs_group_path(scdrs_dir, atlas, gwas, tier)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    return path


def _write_seismic(seismic_dir: Path, atlas: str, gwas: str, tier: str,
                   df: pd.DataFrame) -> Path:
    path = seismic_path(seismic_dir, atlas, gwas, tier)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    return path


def test_load_scdrs_group_normalizes_short_column_names(tmp_path: Path) -> None:
    """scDRS ships group / assoc_mcz / assoc_mcp / assoc_mcq / n_cell;
    the loader normalizes to cell_type / z_mean / pvalue / fdr / n_cells."""
    raw = pd.DataFrame({
        "group": ["T cell", "B cell", "fibroblast"],
        "assoc_mcz": [2.5, -0.3, 1.1],
        "assoc_mcp": [0.005, 0.62, 0.13],
        "assoc_mcq": [0.02, 0.62, 0.20],
        "n_cell": [2300, 1500, 800],
    })
    _write_scdrs(tmp_path, "smillie", "delange", "broad", raw)
    df = load_scdrs_group(tmp_path, "smillie", "delange", "broad")
    assert list(df.columns) == ["cell_type", "z_mean", "pvalue", "fdr", "n_cells"]
    assert df.iloc[0]["cell_type"] == "T cell"
    assert df.iloc[0]["z_mean"] == 2.5
    assert df.iloc[2]["n_cells"] == 800


def test_load_scdrs_group_accepts_canonical_names(tmp_path: Path) -> None:
    """If a file already uses canonical names, the loader leaves them."""
    raw = pd.DataFrame({
        "cell_type": ["T cell"],
        "z_mean": [2.5],
        "pvalue": [0.005],
        "fdr": [0.02],
        "n_cells": [2300],
    })
    _write_scdrs(tmp_path, "garrido_trigo", "delange", "broad", raw)
    df = load_scdrs_group(tmp_path, "garrido_trigo", "delange", "broad")
    assert df.iloc[0]["cell_type"] == "T cell"


def test_load_scdrs_group_raises_on_missing_columns(tmp_path: Path) -> None:
    raw = pd.DataFrame({
        "group": ["T cell"],
        "assoc_mcz": [2.5],
        # missing pvalue, fdr, n_cells
    })
    _write_scdrs(tmp_path, "smillie", "delange", "broad", raw)
    with pytest.raises(KeyError, match="missing columns"):
        load_scdrs_group(tmp_path, "smillie", "delange", "broad")


def test_load_seismic_results_validates_columns(tmp_path: Path) -> None:
    raw = pd.DataFrame({
        "cell_type": ["T cell"],
        "coefficient": [0.15],
        "se": [0.04],
        "pvalue": [0.001],
        "n_cells": [2300],
    })
    _write_seismic(tmp_path, "smillie", "delange", "broad", raw)
    df = load_seismic_results(tmp_path, "smillie", "delange", "broad")
    assert df.iloc[0]["coefficient"] == 0.15


def test_load_seismic_results_raises_on_missing_columns(tmp_path: Path) -> None:
    raw = pd.DataFrame({
        "cell_type": ["T cell"],
        "coefficient": [0.15],
        # missing se, pvalue, n_cells
    })
    _write_seismic(tmp_path, "smillie", "delange", "broad", raw)
    with pytest.raises(KeyError, match="missing columns"):
        load_seismic_results(tmp_path, "smillie", "delange", "broad")


def test_shared_cell_types_filters_by_min_cells() -> None:
    df_a = pd.DataFrame({
        "cell_type": ["T cell", "B cell", "fibroblast"],
        "n_cells":   [2300,     1500,     30],
    })
    df_b = pd.DataFrame({
        "cell_type": ["T cell", "fibroblast", "endothelium"],
        "n_cells":   [3000,     800,         200],
    })
    shared = shared_cell_types(df_a, df_b, min_cells=50)
    # fibroblast fails min_cells in df_a (30); endothelium absent in df_a;
    # B cell absent in df_b. Only T cell qualifies.
    assert shared == ["T cell"]


def test_shared_cell_types_default_min_cells_is_50() -> None:
    df_a = pd.DataFrame({"cell_type": ["x"], "n_cells": [49]})
    df_b = pd.DataFrame({"cell_type": ["x"], "n_cells": [100]})
    assert shared_cell_types(df_a, df_b) == []
    df_a = pd.DataFrame({"cell_type": ["x"], "n_cells": [50]})
    assert shared_cell_types(df_a, df_b) == ["x"]


def test_to_lookup_builds_dict_from_columns() -> None:
    df = pd.DataFrame({
        "cell_type": ["T cell", "B cell"],
        "z_mean":    [2.5,      -0.3],
    })
    d = to_lookup(df, "cell_type", "z_mean")
    assert d == {"T cell": 2.5, "B cell": -0.3}


def test_require_files_passes_when_all_exist(tmp_path: Path) -> None:
    files = [tmp_path / "a.tsv", tmp_path / "b.tsv"]
    for p in files:
        p.write_text("ok")
    require_files(files)  # should not raise


def test_require_files_raises_with_listed_paths(tmp_path: Path) -> None:
    present = tmp_path / "present.tsv"
    present.write_text("ok")
    absent = tmp_path / "absent.tsv"
    with pytest.raises(SystemExit, match="Missing input"):
        require_files([present, absent])
