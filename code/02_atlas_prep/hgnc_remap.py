"""HGNC remap utility — pinned, reproducible.

CELLxGENE deposits store gene identifiers as Ensembl IDs in ``var_names``
and the corresponding HGNC symbol in ``var['feature_name']``. Loaders
that ship raw symbols (Smillie SCP259, Mennillo GEO) skip the
``feature_name`` step. Downstream MAGMA / scDRS / seismicGWAS code
matches on HGNC symbols, so every loader's final step is to remap
``var_names`` to symbols, drop duplicates, and filter to the
NCBI-authoritative symbol set.

Reproducibility pin
-------------------
The authoritative symbol set is read from a **committed, dated NCBI
``gene_info`` snapshot** at ``data/reference/gene_info.tsv.gz``. NCBI
updates monthly; the live-fetch version of this module produced runs
that were not byte-reproducible. The committed snapshot is the single
source of truth; no live URL fetch happens. To refresh the pin,
download a new ``Homo_sapiens.gene_info.gz``, update
``GENE_INFO_PIN_DATE`` below, replace the committed file in one commit,
and document the date bump in DECISIONS.md.

The approved set is built from the ``Symbol`` column **only**. The
previous version also pulled in the ``Synonyms`` column, which made the
membership filter near-permissive and let deprecated aliases pass
through unchanged. If alias resolution is ever needed, it belongs as an
explicit alias → approved remap, not as a membership test.

Canonical-hit survival gate
---------------------------
After remap + dedup + symbol-validity filter, the five canonical UC
GWAS hits in ``CANONICAL_UC_HITS`` must survive (≥95% of the list).
This is a hard gate: if a fundamental immune/IBD locus disappears from
an atlas after harmonization, the pipeline is silently wrong and
downstream scores are misleading.

Spec: ``code/02_atlas_prep/README.md`` (HGNC pin section);
DECISIONS.md corrections 2026-05-20 (5/7) and 2026-06-04 (11).
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

from anndata import AnnData
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# NCBI gene_info dump committed to data/reference/. Update both fields
# together if you re-pin to a fresher snapshot.
GENE_INFO_PIN_DATE = "2026-05-21"
_PINNED_GENE_INFO_REL = Path("data") / "reference" / "gene_info.tsv.gz"

# Canonical UC GWAS hits; ≥95% must survive remap to pass the gate.
# These are well-established IBD/UC loci (IL23R, JAK2, TYK2 = JAK-STAT
# axis; NKX2-3, ATG16L1 = IBD susceptibility). If even one drops out
# silently, the harmonization step is broken in a way that biases every
# downstream score.
CANONICAL_UC_HITS: tuple[str, ...] = (
    "IL23R", "JAK2", "TYK2", "NKX2-3", "ATG16L1",
)
CANONICAL_SURVIVAL_THRESHOLD: float = 0.95


def _pinned_gene_info_path() -> Path:
    """Resolve the committed gene_info snapshot relative to the repo root."""
    # __file__ -> code/02_atlas_prep/hgnc_remap.py; up three is the repo root.
    return Path(__file__).resolve().parents[2] / _PINNED_GENE_INFO_REL


def _load_ncbi_symbol_set() -> set[str]:
    """Read the committed NCBI ``gene_info`` snapshot and return its
    ``Symbol`` column as a set.

    No synonyms. Raises ``FileNotFoundError`` if the snapshot is missing
    (the pin is required — no silent fallback).
    """
    path = _pinned_gene_info_path()
    if not path.exists():
        raise FileNotFoundError(
            f"HGNC remap: pinned gene_info snapshot not found at {path}. "
            f"This file is committed to the repo at data/reference/; "
            f"check your working copy (it should not be gitignored). "
            f"Live-fetch is intentionally disabled — see module docstring."
        )
    symbols: set[str] = set()
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").lstrip("#").split("\t")
        sym_idx = header.index("Symbol")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= sym_idx:
                continue
            sym = parts[sym_idx].strip()
            if sym and sym != "-":
                symbols.add(sym)
    logger.info(
        "Loaded %d NCBI-approved Symbols (pin date %s)",
        len(symbols), GENE_INFO_PIN_DATE,
    )
    return symbols


def ensembl_to_hgnc(adata: AnnData) -> AnnData:
    """Set ``adata.var_names`` to HGNC symbols.

    Strategy:
      1. Prefer ``var['feature_name']`` (HGNC symbols shipped by
         CELLxGENE).
      2. Otherwise fall back to the existing ``var_names`` (loaders that
         ship raw symbols, e.g. Smillie SCP259).
      3. Resolve duplicate symbols by keeping the row with highest
         summed expression (proxy for highest-expressed locus / isoform).
      4. Drop symbols not in the pinned NCBI ``Symbol`` set.
      5. Assert ≥95% of ``CANONICAL_UC_HITS`` survive; raise otherwise.
      6. Log dropped counts.

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

    # 3) Resolve duplicate symbols by max-summed-expression. Summing the
    # duplicate rows (rather than keeping the highest-expressed one) would
    # also be defensible, but max-keep preserves a row's full sparsity
    # pattern and avoids merging genes that share a symbol but have
    # genuinely distinct Ensembl IDs.
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
            adata.n_vars, n_in,
        )

    # 4) Drop symbols not in the pinned NCBI Symbol set.
    approved = _load_ncbi_symbol_set()
    in_approved = adata.var_names.isin(approved)
    n_drop = int((~in_approved).sum())
    if n_drop:
        logger.info(
            "Dropping %d symbols not in pinned NCBI Symbol set (pin date %s)",
            n_drop, GENE_INFO_PIN_DATE,
        )
        adata = adata[:, in_approved].copy()

    # 5) Canonical-hit survival gate. Loud raise — a missing IBD GWAS hit
    # post-remap means the harmonization is silently wrong.
    final_names = set(adata.var_names)
    present = [g for g in CANONICAL_UC_HITS if g in final_names]
    survival = len(present) / len(CANONICAL_UC_HITS)
    if survival < CANONICAL_SURVIVAL_THRESHOLD:
        missing = [g for g in CANONICAL_UC_HITS if g not in final_names]
        raise ValueError(
            f"HGNC remap: canonical UC-hit survival {survival:.0%} < "
            f"{CANONICAL_SURVIVAL_THRESHOLD:.0%}. Missing: {missing} "
            f"(checked: {list(CANONICAL_UC_HITS)}). The input atlas does "
            f"not carry one of the load-bearing IBD/UC GWAS loci after "
            f"symbol harmonization — verify (a) the input gene set, "
            f"(b) the gene_info pin (date {GENE_INFO_PIN_DATE}), and "
            f"(c) whether the missing symbol is shipped as an alias and "
            f"needs an explicit alias→approved remap."
        )
    logger.info(
        "Canonical-hit survival: %d/%d (%.0f%%) — %s",
        len(present), len(CANONICAL_UC_HITS), survival * 100, present,
    )

    logger.info("HGNC remap: %d -> %d genes", n_in, adata.n_vars)
    return adata
