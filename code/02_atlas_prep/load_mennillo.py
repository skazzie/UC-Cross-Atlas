"""Loader for Mennillo 2024 UC core atlas (SKELETON — deferred to next session).

Source: GEO GSE229072 (verify accession on download).
Paper: Mennillo 2024, Nat Commun 15:1493.

Deferred pending GEO download to Hummingbird scratch.

Per DECISIONS.md correction 2026-05-20 (5/7), after download the loader
must apply log1p(CP10k) normalization so the input distribution is
comparable to the log-normalized CELLxGENE atlases.
"""

from __future__ import annotations

from anndata import AnnData


def load(
    h5ad_path: str,
    apply_v1_filter: bool = True,
    raw_count_mode: bool = False,
) -> AnnData:
    """Not yet implemented; GEO download required first.

    TODO once Mennillo .h5ad is on disk:
      - Verify GEO accession (expected GSE229072)
      - Discover obs schema
      - Subset to pre-treatment baseline samples only
        (anti-integrin therapy time course)
      - Identify broad and fine tier columns
      - Apply log1p(CP10k) normalization (raw counts -> log-normalized)
        to match the CELLxGENE atlases (Correction 5/7)
    """
    raise NotImplementedError(
        "Mennillo loader is deferred. Download GSE229072 from GEO "
        "(verify accession) to Hummingbird scratch before this loader "
        "can be completed. See DECISIONS.md correction 2026-05-20 (5/7) "
        "for the normalization requirement."
    )
