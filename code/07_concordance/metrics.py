"""
Cross-atlas concordance metrics for cell-type prioritization rankings.

Three metrics from spec section 2.7:
  - Spearman rank correlation between cell-type rank vectors
  - Top-k Jaccard overlap of nominated cell types
  - Cohen's kappa for FDR-significance agreement (with marginals)

Plus a bootstrap CI helper that resamples cell types (cheap; for resampling
cells, see compute_concordance.py).

All functions take dicts keyed by cell-type label so callers don't have to
worry about ordering or about which atlases share which cell types.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np
from scipy.stats import spearmanr
from sklearn.metrics import cohen_kappa_score


@dataclass
class ConcordanceResult:
    spearman_rho: float
    spearman_p: float
    n_common: int
    jaccard_top5: float
    jaccard_top10: float
    kappa: float
    n_sig_a: int
    n_sig_b: int


def _common_keys(a: Mapping, b: Mapping) -> list:
    return sorted(set(a) & set(b))


def spearman(scores_a: Mapping[str, float],
             scores_b: Mapping[str, float]) -> tuple[float, float, int]:
    """Spearman ρ on shared cell types. `scores` are values where smaller = stronger
    (e.g. p-values or q-values); the correlation is invariant to that, but stay
    consistent with `top_k_jaccard`."""
    common = _common_keys(scores_a, scores_b)
    if len(common) < 3:
        return float("nan"), float("nan"), len(common)
    a = [scores_a[c] for c in common]
    b = [scores_b[c] for c in common]
    rho, p = spearmanr(a, b)
    return float(rho), float(p), len(common)


def top_k_jaccard(scores_a: Mapping[str, float],
                  scores_b: Mapping[str, float],
                  k: int = 5) -> float:
    """Jaccard overlap of the top-k cell types in each atlas, ranked by smallest
    score (so pass p-values or q-values, not z-statistics)."""
    common = _common_keys(scores_a, scores_b)
    if len(common) < k:
        return float("nan")
    top_a = set(sorted(common, key=lambda c: scores_a[c])[:k])
    top_b = set(sorted(common, key=lambda c: scores_b[c])[:k])
    union = top_a | top_b
    return len(top_a & top_b) / len(union) if union else float("nan")


def fdr_concordance(qvals_a: Mapping[str, float],
                    qvals_b: Mapping[str, float],
                    threshold: float = 0.05) -> tuple[float, int, int, int]:
    """Cohen's κ on binary 'q < threshold' calls, plus marginals.

    Returns (kappa, n_sig_a, n_sig_b, n_common). Always report the marginals
    alongside κ — when nearly every cell type is significant in both atlases,
    κ collapses even though raw agreement is high (spec §2.5)."""
    common = _common_keys(qvals_a, qvals_b)
    if len(common) < 2:
        return float("nan"), 0, 0, len(common)
    sig_a = [qvals_a[c] < threshold for c in common]
    sig_b = [qvals_b[c] < threshold for c in common]
    if all(sig_a) and all(sig_b):
        # sklearn returns nan with a warning; surface it explicitly.
        return float("nan"), sum(sig_a), sum(sig_b), len(common)
    if not any(sig_a) and not any(sig_b):
        return float("nan"), 0, 0, len(common)
    kappa = cohen_kappa_score(sig_a, sig_b)
    return float(kappa), sum(sig_a), sum(sig_b), len(common)


def concordance(pvals_a: Mapping[str, float],
                pvals_b: Mapping[str, float],
                qvals_a: Mapping[str, float],
                qvals_b: Mapping[str, float],
                fdr_threshold: float = 0.05) -> ConcordanceResult:
    rho, rho_p, n = spearman(pvals_a, pvals_b)
    j5 = top_k_jaccard(pvals_a, pvals_b, k=5)
    j10 = top_k_jaccard(pvals_a, pvals_b, k=10)
    kappa, ns_a, ns_b, _ = fdr_concordance(qvals_a, qvals_b, threshold=fdr_threshold)
    return ConcordanceResult(
        spearman_rho=rho, spearman_p=rho_p, n_common=n,
        jaccard_top5=j5, jaccard_top10=j10,
        kappa=kappa, n_sig_a=ns_a, n_sig_b=ns_b,
    )


def bootstrap_spearman_ci(scores_a: Mapping[str, float],
                          scores_b: Mapping[str, float],
                          n_iter: int = 1000,
                          alpha: float = 0.05,
                          seed: int = 0) -> tuple[float, float]:
    """Percentile bootstrap CI for Spearman ρ by resampling cell types with
    replacement. This is the *cheap* CI — it does NOT account for cell-level
    sampling noise upstream. For the headline numbers in the paper, also run
    cell-level bootstraps (resample cells within each atlas, re-run scDRS,
    re-compute ρ); 100 iterations is enough for a usable 95% CI per spec §2.7."""
    common = _common_keys(scores_a, scores_b)
    if len(common) < 5:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    a = np.array([scores_a[c] for c in common])
    b = np.array([scores_b[c] for c in common])
    n = len(common)
    rhos = np.empty(n_iter)
    for i in range(n_iter):
        idx = rng.integers(0, n, size=n)
        if len(set(idx)) < 3:
            rhos[i] = np.nan
            continue
        rho, _ = spearmanr(a[idx], b[idx])
        rhos[i] = rho
    rhos = rhos[~np.isnan(rhos)]
    if len(rhos) < 10:
        return float("nan"), float("nan")
    lo = float(np.percentile(rhos, 100 * alpha / 2))
    hi = float(np.percentile(rhos, 100 * (1 - alpha / 2)))
    return lo, hi
