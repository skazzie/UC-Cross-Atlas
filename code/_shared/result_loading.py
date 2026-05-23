"""Loaders for scDRS + seismicGWAS result files.

Paths match the conventions documented in code/03_scdrs/README.md and
code/04_seismic/README.md:

  results/scdrs/<atlas>_<gwas>/cell_type_<tier>/<atlas>_<gwas>.scdrs_group
  results/seismic/<atlas>_<gwas>_<tier>.tsv

Both 08_cross_method and 09_cross_gwas drivers need these; do not duplicate.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pandas as pd


def scdrs_group_path(scdrs_dir: Path, atlas: str, gwas: str, tier: str) -> Path:
    return (
        scdrs_dir
        / f"{atlas}_{gwas}"
        / f"cell_type_{tier}"
        / f"{atlas}_{gwas}.scdrs_group"
    )


def seismic_path(seismic_dir: Path, atlas: str, gwas: str, tier: str) -> Path:
    return seismic_dir / f"{atlas}_{gwas}_{tier}.tsv"


def load_scdrs_group(
    scdrs_dir: Path, atlas: str, gwas: str, tier: str
) -> pd.DataFrame:
    """Return df with columns: cell_type, z_mean, pvalue, fdr, n_cells.

    scDRS' ``perform-downstream --group-analysis`` writes a tab-separated
    table with one row per cell-type group. Column names vary slightly
    between scDRS versions; this loader normalizes to the canonical set
    used by downstream metrics and raises a clear error if the expected
    columns are absent.
    """
    path = scdrs_group_path(scdrs_dir, atlas, gwas, tier)
    df = pd.read_csv(path, sep="\t")
    rename_map = {
        "group": "cell_type",
        "assoc_mcz": "z_mean",
        "assoc_mcp": "pvalue",
        "assoc_mcq": "fdr",
        "n_cell": "n_cells",
    }
    for src, dst in rename_map.items():
        if src in df.columns and dst not in df.columns:
            df = df.rename(columns={src: dst})
    required = {"cell_type", "z_mean", "pvalue", "fdr", "n_cells"}
    missing = required - set(df.columns)
    if missing:
        raise KeyError(
            f"{path} is missing columns {missing}; got {list(df.columns)}"
        )
    return df[["cell_type", "z_mean", "pvalue", "fdr", "n_cells"]].copy()


def load_seismic_results(
    seismic_dir: Path, atlas: str, gwas: str, tier: str
) -> pd.DataFrame:
    """Return df with columns: cell_type, coefficient, se, pvalue, n_cells.

    Written by ``code/04_seismic/run_seismic.R``; column names follow the
    seismicGWAS convention.
    """
    path = seismic_path(seismic_dir, atlas, gwas, tier)
    df = pd.read_csv(path, sep="\t")
    required = {"cell_type", "coefficient", "se", "pvalue", "n_cells"}
    missing = required - set(df.columns)
    if missing:
        raise KeyError(
            f"{path} is missing columns {missing}; got {list(df.columns)}"
        )
    return df[["cell_type", "coefficient", "se", "pvalue", "n_cells"]].copy()


def shared_cell_types(
    df_a: pd.DataFrame, df_b: pd.DataFrame, min_cells: int = 50
) -> list[str]:
    """Cell-type intersection of two result tables, each side >= min_cells.

    Both inputs must have ``cell_type`` and ``n_cells`` columns.
    """
    a_ok = set(df_a.loc[df_a["n_cells"] >= min_cells, "cell_type"])
    b_ok = set(df_b.loc[df_b["n_cells"] >= min_cells, "cell_type"])
    return sorted(a_ok & b_ok)


def to_lookup(
    df: pd.DataFrame, key: str, value: str
) -> dict[str, float]:
    """Build a {cell_type: value} dict from a result table."""
    return dict(zip(df[key].astype(str), df[value].astype(float)))


def require_files(paths: Iterable[Path]) -> None:
    """Raise SystemExit(2) with a clear message for any missing path."""
    missing = [str(p) for p in paths if not Path(p).exists()]
    if missing:
        msg = "Missing input file(s):\n  " + "\n  ".join(missing)
        raise SystemExit(msg)
