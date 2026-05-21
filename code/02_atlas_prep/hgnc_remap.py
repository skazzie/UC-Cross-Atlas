"""HGNC remap utility.

CELLxGENE deposits store gene identifiers as Ensembl IDs in ``var_names``
and the corresponding HGNC symbol in ``var['feature_name']``. Downstream
MAGMA / scDRS / seismicGWAS code matches on HGNC symbols, so every
loader's final step is to remap ``var_names`` to symbols and drop
duplicates / non-canonical entries.

Spec: ``code/02_atlas_prep/atlas_schemas.md`` and DECISIONS.md correction
2026-05-20 (4/7), (5/7).
"""

from __future__ import annotations

import gzip
import logging
import os
import urllib.request
from pathlib import Path

from anndata import AnnData
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

NCBI_GENE_INFO_URL = (
    "https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
)


def _gene_info_cache_path() -> Path:
    base = os.environ.get("UCC_DATA")
    if base:
        cache_dir = Path(base) / "reference"
    else:
        cache_dir = Path("data") / "reference"
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir / "gene_info.tsv.gz"


def _load_ncbi_symbol_set() -> set[str]:
    """Authoritative HGNC-approved symbol set from NCBI gene_info.

    Cached at ``$UCC_DATA/reference/gene_info.tsv.gz`` (or
    ``./data/reference/gene_info.tsv.gz``); refetched only if missing.
    """
    cache = _gene_info_cache_path()
    if not cache.exists():
        logger.info("Downloading NCBI gene_info to %s", cache)
        urllib.request.urlretrieve(NCBI_GENE_INFO_URL, cache)
    symbols: set[str] = set()
    with gzip.open(cache, "rt", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").lstrip("#").split("\t")
        sym_idx = header.index("Symbol")
        syn_idx = header.index("Synonyms")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(sym_idx, syn_idx):
                continue
            sym = parts[sym_idx].strip()
            if sym and sym != "-":
                symbols.add(sym)
            syn = parts[syn_idx].strip()
            if syn and syn != "-":
                for s in syn.split("|"):
                    s = s.strip()
                    if s:
                        symbols.add(s)
    return symbols


def ensembl_to_hgnc(adata: AnnData) -> AnnData:
    """Set ``adata.var_names`` to HGNC symbols.

    Strategy:
      1. Prefer ``var['feature_name']`` (HGNC symbols shipped by CELLxGENE).
      2. Otherwise fall back to the NCBI ``gene_info`` symbol set: keep
         var_names that are already authoritative symbols.
      3. Drop duplicate symbols, keeping the row with highest summed
         expression (proxy for highest-expressed isoform/loci).
      4. Drop symbols not in the NCBI authoritative list.
      5. Log dropped counts.

    Returns the updated AnnData (a copy; original is not modified).
    """
    if "feature_name" in adata.var.columns:
        new_names = adata.var["feature_name"].astype(str).to_numpy()
        source = "var['feature_name']"
    else:
        new_names = adata.var_names.astype(str).to_numpy()
        source = "var_names (no feature_name column)"
    logger.info("HGNC source: %s", source)

    n_in = adata.n_vars
    adata = adata.copy()
    adata.var_names = pd.Index(new_names)
    adata.var_names_make_unique = False  # we resolve duplicates explicitly below

    # 3) Resolve duplicate symbols by keeping the row with highest summed expression.
    if not adata.var_names.is_unique:
        n_dup = int((adata.var_names.value_counts() > 1).sum())
        logger.info("Resolving %d duplicate HGNC symbol(s) by max-expression", n_dup)
        X = adata.X
        if hasattr(X, "toarray"):
            sums = np.asarray(X.sum(axis=0)).ravel()
        else:
            sums = np.asarray(X).sum(axis=0).ravel()
        order = np.argsort(-sums, kind="stable")
        sorted_names = adata.var_names.to_numpy()[order]
        _, first_idx = np.unique(sorted_names, return_index=True)
        keep_positions = order[first_idx]
        mask = np.zeros(adata.n_vars, dtype=bool)
        mask[keep_positions] = True
        adata = adata[:, mask].copy()
        logger.info(
            "Kept %d unique symbols after deduplication (was %d)",
            adata.n_vars,
            n_in,
        )

    # 4) Drop symbols not in NCBI authoritative list.
    try:
        approved = _load_ncbi_symbol_set()
    except Exception as exc:  # pragma: no cover
        logger.warning(
            "Could not load NCBI gene_info (%s); skipping symbol-validity filter.",
            exc,
        )
        approved = None

    if approved is not None:
        in_approved = adata.var_names.isin(approved)
        n_drop = int((~in_approved).sum())
        if n_drop:
            logger.info("Dropping %d symbols not in NCBI authoritative list", n_drop)
            adata = adata[:, in_approved].copy()

    logger.info("HGNC remap: %d -> %d genes", n_in, adata.n_vars)
    return adata
