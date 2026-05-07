"""
Cross-atlas concordance metrics for cell-type prioritization rankings.

Locked v1 commitments (PLAN.md, DECISIONS.md):

- **Headline metric: Spearman rho on cell-type-level Z-scores (scDRS)
  and regression coefficients (seismicGWAS) — NOT p-values.** scDRS
  p-values use an atlas-specific Monte Carlo null; same biology in two
  atlases produces different p-values purely from compositional differences.
- **Tied-rank handling: average-rank tie-breaking (scipy default).**
- **Concordance computed on shared cell-type intersection per pair, with
  minimum cell-count threshold (>=50 cells in BOTH atlases per cell type).**
  Atlas-specific cell types are reported separately, not entered into
  concordance metrics in either direction.
- **Bootstrap 95% CIs on every reported Spearman rho:** 1000 iterations,
  resampling over cell types within shared intersection, percentile method,
  seed = 42 (locked random seed).
- **Top-k Jaccard:** k = 5, 10 at broad tier; k = 5, 10, 20 at fine tier.
- **Cohen's kappa on FDR-significance** with marginal-saturation contingency:
  if >=80% of cell types pass FDR < 0.05 in both atlases, report
  kappa @ FDR < 0.01 as the headline kappa instead.

All metrics take dicts keyed by cell-type label so callers do not have to
worry about ordering or about which atlases share which cell types.

Score-direction convention
--------------------------
This module's API treats "scores" as values where LARGER = STRONGER
(cell-type-level Z-scores, regression coefficients). This is the locked
v1 headline convention. Spearman rho is invariant under monotone
transforms, but `top_k_jaccard` is not — pass `larger_is_stronger=False`
when ranking by p-values instead of Z-scores.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Optional

import numpy as np
from scipy.stats import spearmanr
from sklearn.metrics import cohen_kappa_score

# Locked random seed (DECISIONS.md)
DEFAULT_SEED = 42

# Locked min cell-count threshold (DECISIONS.md)
DEFAULT_MIN_CELLS = 50

# Locked saturation threshold for kappa contingency (PLAN.md, DECISIONS.md)
KAPPA_SATURATION_THRESHOLD = 0.80


@dataclass
class ConcordanceResult:
    """Aggregated concordance for a single (atlas_a, atlas_b) pair under
    a fixed (method, GWAS, granularity)."""
    spearman_rho: float
    spearman_p: float
    spearman_ci_lo: float
    spearman_ci_hi: float
    n_common: int
    jaccard_top5: float
    jaccard_top10: float
    jaccard_top20: float  # only meaningful at fine tier; NaN at broad
    kappa: float
    kappa_threshold: float  # 0.05 by default; 0.01 under saturation contingency
    kappa_threshold_used_due_to_saturation: bool
    n_sig_a: int
    n_sig_b: int
    excluded_low_count: int  # cell types dropped by min-cell-count filter


def _intersect_with_min_cells(
    scores_a: Mapping[str, float],
    scores_b: Mapping[str, float],
    cell_counts_a: Optional[Mapping[str, int]],
    cell_counts_b: Optional[Mapping[str, int]],
    min_cells: int,
) -> tuple[list, int]:
    """Return the shared cell-type set with each side meeting min_cells.

    If cell counts are not provided, no min-cell filter is applied
    (concordance is computed on raw shared keys). Returns (kept, n_dropped).
    """
    common = sorted(set(scores_a) & set(scores_b))
    if cell_counts_a is None or cell_counts_b is None:
        return common, 0
    kept = [
        c for c in common
        if cell_counts_a.get(c, 0) >= min_cells and cell_counts_b.get(c, 0) >= min_cells
    ]
    return kept, len(common) - len(kept)


def spearman(
    scores_a: Mapping[str, float],
    scores_b: Mapping[str, float],
    cell_counts_a: Optional[Mapping[str, int]] = None,
    cell_counts_b: Optional[Mapping[str, int]] = None,
    min_cells: int = DEFAULT_MIN_CELLS,
) -> tuple[float, float, int, int]:
    """Spearman rho on shared cell types meeting the min-cell-count threshold.

    Tied ranks are broken by average-rank (scipy.stats.spearmanr default,
    locked in DECISIONS.md).

    Returns (rho, p, n_kept, n_dropped_low_count).
    """
    kept, dropped = _intersect_with_min_cells(
        scores_a, scores_b, cell_counts_a, cell_counts_b, min_cells
    )
    if len(kept) < 3:
        return float("nan"), float("nan"), len(kept), dropped
    a = [scores_a[c] for c in kept]
    b = [scores_b[c] for c in kept]
    rho, p = spearmanr(a, b)
    return float(rho), float(p), len(kept), dropped


def top_k_jaccard(
    scores_a: Mapping[str, float],
    scores_b: Mapping[str, float],
    k: int = 5,
    cell_counts_a: Optional[Mapping[str, int]] = None,
    cell_counts_b: Optional[Mapping[str, int]] = None,
    min_cells: int = DEFAULT_MIN_CELLS,
    larger_is_stronger: bool = True,
) -> float:
    """Jaccard overlap of the top-k cell types in each atlas.

    `larger_is_stronger=True` (locked default): scores are Z-scores or
    regression coefficients; rank descending.

    `larger_is_stronger=False`: scores are p-values; rank ascending.
    """
    kept, _ = _intersect_with_min_cells(
        scores_a, scores_b, cell_counts_a, cell_counts_b, min_cells
    )
    if len(kept) < k:
        return float("nan")
    if larger_is_stronger:
        top_a = set(sorted(kept, key=lambda c: -scores_a[c])[:k])
        top_b = set(sorted(kept, key=lambda c: -scores_b[c])[:k])
    else:
        top_a = set(sorted(kept, key=lambda c: scores_a[c])[:k])
        top_b = set(sorted(kept, key=lambda c: scores_b[c])[:k])
    union = top_a | top_b
    return len(top_a & top_b) / len(union) if union else float("nan")


def fdr_concordance(
    qvals_a: Mapping[str, float],
    qvals_b: Mapping[str, float],
    threshold: float = 0.05,
    cell_counts_a: Optional[Mapping[str, int]] = None,
    cell_counts_b: Optional[Mapping[str, int]] = None,
    min_cells: int = DEFAULT_MIN_CELLS,
) -> tuple[float, int, int, int]:
    """Cohen's kappa on binary 'q < threshold' calls, plus marginals.

    Returns (kappa, n_sig_a, n_sig_b, n_kept). Always report the marginals
    alongside kappa — when nearly every cell type is significant in both
    atlases, kappa collapses even though raw agreement is high. Use
    `concordance(...)` to apply the saturation contingency automatically.
    """
    kept, _ = _intersect_with_min_cells(
        qvals_a, qvals_b, cell_counts_a, cell_counts_b, min_cells
    )
    if len(kept) < 2:
        return float("nan"), 0, 0, len(kept)
    sig_a = [qvals_a[c] < threshold for c in kept]
    sig_b = [qvals_b[c] < threshold for c in kept]
    if all(sig_a) and all(sig_b):
        return float("nan"), sum(sig_a), sum(sig_b), len(kept)
    if not any(sig_a) and not any(sig_b):
        return float("nan"), 0, 0, len(kept)
    kappa = cohen_kappa_score(sig_a, sig_b)
    return float(kappa), sum(sig_a), sum(sig_b), len(kept)


def bootstrap_spearman_ci(
    scores_a: Mapping[str, float],
    scores_b: Mapping[str, float],
    n_iter: int = 1000,
    alpha: float = 0.05,
    seed: int = DEFAULT_SEED,
    cell_counts_a: Optional[Mapping[str, int]] = None,
    cell_counts_b: Optional[Mapping[str, int]] = None,
    min_cells: int = DEFAULT_MIN_CELLS,
) -> tuple[float, float]:
    """Percentile bootstrap 95% CI for Spearman rho.

    Locked v1 default: 1000 iterations, percentile method, seed=42,
    resampling over cell types within the shared cell-type intersection
    (post min-cell-count filter). Cheap (~seconds per comparison). BCa
    deferred to revision if reviewers request.

    Note: this resamples cell types — it does NOT capture cell-level
    sampling noise. Donor-level uncertainty is reported separately as
    LOO jackknife ranges (see PLAN.md §"Donor-LOO uncertainty intervals").
    """
    kept, _ = _intersect_with_min_cells(
        scores_a, scores_b, cell_counts_a, cell_counts_b, min_cells
    )
    if len(kept) < 5:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    a = np.array([scores_a[c] for c in kept])
    b = np.array([scores_b[c] for c in kept])
    n = len(kept)
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


def concordance(
    scores_a: Mapping[str, float],
    scores_b: Mapping[str, float],
    qvals_a: Mapping[str, float],
    qvals_b: Mapping[str, float],
    fdr_threshold: float = 0.05,
    fdr_threshold_strict: float = 0.01,
    saturation_threshold: float = KAPPA_SATURATION_THRESHOLD,
    cell_counts_a: Optional[Mapping[str, int]] = None,
    cell_counts_b: Optional[Mapping[str, int]] = None,
    min_cells: int = DEFAULT_MIN_CELLS,
    n_bootstrap: int = 1000,
    seed: int = DEFAULT_SEED,
    is_fine_tier: bool = False,
    larger_is_stronger: bool = True,
) -> ConcordanceResult:
    """Compute the locked-v1 concordance bundle for one (atlas_a, atlas_b) pair.

    Headline scores must be larger-is-stronger (cell-type-level Z-scores
    for scDRS, regression coefficients for seismicGWAS). Pass p/q-values via
    `qvals_a` / `qvals_b` for the kappa metric and within-atlas saturation
    contingency.
    """
    rho, rho_p, n_kept, dropped = spearman(
        scores_a, scores_b, cell_counts_a, cell_counts_b, min_cells
    )
    ci_lo, ci_hi = bootstrap_spearman_ci(
        scores_a, scores_b, n_iter=n_bootstrap, seed=seed,
        cell_counts_a=cell_counts_a, cell_counts_b=cell_counts_b,
        min_cells=min_cells,
    )

    j5 = top_k_jaccard(
        scores_a, scores_b, k=5,
        cell_counts_a=cell_counts_a, cell_counts_b=cell_counts_b,
        min_cells=min_cells, larger_is_stronger=larger_is_stronger,
    )
    j10 = top_k_jaccard(
        scores_a, scores_b, k=10,
        cell_counts_a=cell_counts_a, cell_counts_b=cell_counts_b,
        min_cells=min_cells, larger_is_stronger=larger_is_stronger,
    )
    if is_fine_tier:
        j20 = top_k_jaccard(
            scores_a, scores_b, k=20,
            cell_counts_a=cell_counts_a, cell_counts_b=cell_counts_b,
            min_cells=min_cells, larger_is_stronger=larger_is_stronger,
        )
    else:
        j20 = float("nan")

    # Cohen's kappa with marginal-saturation contingency (PLAN.md, DECISIONS.md)
    kappa_05, n_sig_a, n_sig_b, n_kept_q = fdr_concordance(
        qvals_a, qvals_b, threshold=fdr_threshold,
        cell_counts_a=cell_counts_a, cell_counts_b=cell_counts_b,
        min_cells=min_cells,
    )
    kappa_used = kappa_05
    threshold_used = fdr_threshold
    saturated = False
    if n_kept_q > 0:
        frac_a = n_sig_a / n_kept_q
        frac_b = n_sig_b / n_kept_q
        if frac_a >= saturation_threshold and frac_b >= saturation_threshold:
            kappa_strict, n_sig_a_strict, n_sig_b_strict, _ = fdr_concordance(
                qvals_a, qvals_b, threshold=fdr_threshold_strict,
                cell_counts_a=cell_counts_a, cell_counts_b=cell_counts_b,
                min_cells=min_cells,
            )
            kappa_used = kappa_strict
            threshold_used = fdr_threshold_strict
            saturated = True
            n_sig_a = n_sig_a_strict
            n_sig_b = n_sig_b_strict

    return ConcordanceResult(
        spearman_rho=rho,
        spearman_p=rho_p,
        spearman_ci_lo=ci_lo,
        spearman_ci_hi=ci_hi,
        n_common=n_kept,
        jaccard_top5=j5,
        jaccard_top10=j10,
        jaccard_top20=j20,
        kappa=kappa_used,
        kappa_threshold=threshold_used,
        kappa_threshold_used_due_to_saturation=saturated,
        n_sig_a=n_sig_a,
        n_sig_b=n_sig_b,
        excluded_low_count=dropped,
    )
