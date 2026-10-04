"""Driver: load HCA Gut (Elmentaite 2021), remap cell_type_broad to the
canonical vocab, precompute neighbors, and write h5ad + partial covariate
file.

Source pin: CELLxGENE deposit ``f34d2b82-9265-4a73-bda4-852933bf2a8d.h5ad``;
see ``code/02_atlas_prep/atlas_schemas.md`` and the HCA Gut entry in
DECISIONS.md (broad-comparator role, overlap-with-Pan-GI caveat).

Covariate structure — DECISIONS-locked (open-items section of the
Hummingbird rebuild handoff): ``const``, ``log_n_genes``, ``log_n_counts``,
then pre-expanded float 0.0/1.0 dummies for ``assay``, ``batch``,
``Fraction``, and ``sex``. HCA Gut has no biopsy-level sample identifier
separate from ``batch`` (per-library 10x run), and donor-level variation
is absorbed by ``batch`` under the sample-in-donor nesting invariant
checked below. There is NO donor dummy block in the cov file: batch
nests in donor, so donor dummies would be a linear combination of batch
dummies and the design would be rank-deficient (same failure mode as
TAURUS's donor+sample+intercept crash).

Why pre-expanded, not raw categoricals: scDRS's ``category2dummy`` path
emits pandas bool dummies, then ``df_cov.values`` upcasts to ``object``
and np.linalg.solve dies with ``Cannot cast dtype('O') to float64``.
Pre-expanding to explicit float dummies here hands scDRS a fully numeric
matrix. An all-numeric assertion runs before write.

NO disease / health covariate. HCA Gut's ``disease`` column MUST NOT
enter the cov — that would scrub the very signal scDRS tests
(no-disease-covariate rule).

Mirrors ``run_pangi_load.py``:
- scDRS-replica pre-filter (``min_genes=250`` / ``min_cells=50``) then
  ``sc.pp.pca(n_comps=20)`` + ``sc.pp.neighbors(n_neighbors=15,
  n_pcs=20)``, so ``obsp['connectivities']`` ships with the h5ad and
  ``scdrs perform-downstream --flag-filter-data False`` doesn't rebuild
  the kNN (~45 min saved per run). ``uns['scdrs_prefilter']`` stamped.
- Object-dtype obs columns coerced to str with NaN -> "NA" before
  ``write_h5ad`` — h5py cannot serialise mixed/NaN object columns; this
  killed the first Pan-GI load on ``obs['batch']``.
- Explicit ``--out-h5ad`` / ``--out-cov`` (no ``../../`` defaults that
  resolve outside the repo and fail after a successful full load).
- Fail-loud gate on unmapped ``cell_type_broad`` labels against
  ``_broad_vocab._BROAD_VOCAB``, with the full label list in the error
  so the HCA_TO_BROAD map (empty v0) can be populated iteratively.
- Depth source is logged so a future CELLxGENE schema change is visible
  in the driver log, not a silent miscounting of depth.
"""

from __future__ import annotations

import argparse
import logging

import numpy as np
import pandas as pd
import scanpy as sc

from load_hca_gut import load
from _broad_vocab import _BROAD_VOCAB

logger = logging.getLogger(__name__)

# scDRS 1.0.2 load-time filter defaults (scdrs/util.py:82-88). Match
# here so the on-disk cell set == what --flag-filter-data True would
# produce.
SCDRS_MIN_GENES_PER_CELL = 250
SCDRS_MIN_CELLS_PER_GENE = 50

# HCA Gut donor identifier (per load_hca_gut.load() → obs['donor']).
DONOR_COL = "donor"

# Categorical cov columns pre-expanded to float 0/1 dummies.
# Locked by DECISIONS (open-items section of 2026-10-03 handoff): the
# HCA covariate block is {assay, batch, Fraction, sex}. These are built
# as a single hstack of pd.get_dummies(...).astype(float) frames below.
HCA_COV_CATEGORICALS: tuple[str, ...] = ("assay", "batch", "Fraction", "sex")

# Map from HCA Gut `category` (= cell_type_broad after load_hca_gut.load()
# standardises obs) into the 15-term canonical broad vocab. v0 ships
# EMPTY — the gate-2 assertion below raises on first run with the full
# label list so Muskaan can populate it in one biology-review commit,
# same iterative pattern as PANGI_TO_BROAD in run_pangi_load.py.
#
# If `category` already uses canonical labels verbatim (as expected for
# HCA Gut's broad category axis), the gate-2 check becomes a no-op guard
# and the map stays empty. If the check raises, add entries here.
HCA_TO_BROAD: dict[str, str] = {}

# Gate (1): every value the map ships must be in the canonical vocab.
# Vacuously true while the map is empty; protects the day we start
# filling it in.
_unmapped_broad_out = set(HCA_TO_BROAD.values()) - _BROAD_VOCAB
if _unmapped_broad_out:
    raise ValueError(
        f"run_hca_load.HCA_TO_BROAD ships broad values outside "
        f"_BROAD_VOCAB: {sorted(_unmapped_broad_out)}. Typo on the "
        f"value side of the map; see canonical_broad_DRAFT.md."
    )
del _unmapped_broad_out


def _raw_counts_matrix(adata):
    """Depth covariates (log_n_counts, log_n_genes) MUST be computed
    from raw UMIs. Prefer ``layers['counts']``, then ``adata.raw.X``,
    then ``adata.X`` with a loud warning (would miscount depth if X is
    log-normalised). HCA Gut on CELLxGENE typically ships raw counts in
    ``.raw.X`` and log1p(CP10k) in ``.X``.
    """
    if "counts" in adata.layers:
        logger.info("[driver] depth from layers['counts']")
        return adata.layers["counts"]
    if getattr(adata, "raw", None) is not None:
        logger.info("[driver] depth from adata.raw.X (no layers['counts'])")
        return adata.raw.X
    logger.warning(
        "[driver] no layers['counts'] and no adata.raw; using adata.X "
        "for depth. If X is log-normalized this MISCOUNTS depth — fix "
        "upstream in load_hca_gut.load()."
    )
    return adata.X


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--h5ad", required=True,
                   help="HCA Gut CELLxGENE h5ad "
                        "(f34d2b82-9265-4a73-bda4-852933bf2a8d.h5ad)")
    # No ../../ defaults — those resolve outside the repo and fail at
    # write time after a successful full load. Caller must specify.
    p.add_argument("--out-h5ad", required=True)
    p.add_argument("--out-cov", required=True)
    a = p.parse_args()

    adata = load(a.h5ad)  # apply_v1_filter=True; Adult x colon
    logger.info(
        "[driver] loaded: %d cells x %d genes", adata.n_obs, adata.n_vars,
    )

    # ---- Cell-type broad gate (fail-loud on unmapped labels) ----
    raw_broad = adata.obs["cell_type_broad"].astype(str)
    unique_broad = sorted(raw_broad.unique())
    if HCA_TO_BROAD:
        unmapped = sorted(set(unique_broad) - set(HCA_TO_BROAD))
        if unmapped:
            raise KeyError(
                f"HCA Gut driver: {len(unmapped)} cell_type_broad labels "
                f"have no HCA_TO_BROAD entry. Extend the map "
                f"(run_hca_load.HCA_TO_BROAD). Unmapped labels (full "
                f"list): {unmapped}"
            )
        broad = raw_broad.map(HCA_TO_BROAD)
    else:
        broad = raw_broad
    # Gate 2: emitted values must be in the canonical vocab, regardless
    # of whether we routed through HCA_TO_BROAD or passed labels verbatim.
    outside = sorted(set(broad.dropna().unique()) - _BROAD_VOCAB)
    if outside:
        raise KeyError(
            f"HCA Gut driver: {len(outside)} cell_type_broad values are "
            f"not in _BROAD_VOCAB. Either HCA_TO_BROAD needs these keys "
            f"(current map: {sorted(HCA_TO_BROAD) if HCA_TO_BROAD else 'empty'}), "
            f"or _BROAD_VOCAB needs extending (see canonical_broad_DRAFT.md). "
            f"Unmapped values (full list): {outside}"
        )
    adata.obs["cell_type_broad"] = broad.astype("category")
    logger.info(
        "[driver] cell_type_broad gate passed: %d unique canonical labels",
        int(adata.obs["cell_type_broad"].nunique()),
    )

    # ---- scDRS-replica pre-filter (matches add_neighbors.py so
    #      perform-downstream can be invoked with --flag-filter-data
    #      False without a cell-set mismatch). filter_cells counts
    #      nonzero-X genes per cell; log1p(0)==0 so the cell set is
    #      identical to what scDRS would produce on raw counts.
    n0_cells, n0_genes = adata.n_obs, adata.n_vars
    sc.pp.filter_cells(adata, min_genes=SCDRS_MIN_GENES_PER_CELL)
    sc.pp.filter_genes(adata, min_cells=SCDRS_MIN_CELLS_PER_GENE)
    print(
        f"[driver] scDRS-replica filter: "
        f"{n0_cells}->{adata.n_obs} cells "
        f"(min_genes={SCDRS_MIN_GENES_PER_CELL}), "
        f"{n0_genes}->{adata.n_vars} genes "
        f"(min_cells={SCDRS_MIN_CELLS_PER_GENE})"
    )

    # ---- PCA + kNN so obsp['connectivities'] ships with the h5ad and
    #      scdrs perform-downstream skips its ~45-min rebuild. Params
    #      match scDRS group-analysis defaults (knn_n_pcs=20,
    #      knn_n_neighbors=15).
    sc.pp.pca(adata, n_comps=20)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=20)
    assert "connectivities" in adata.obsp, (
        "sc.pp.neighbors did not populate obsp['connectivities']"
    )

    adata.uns["scdrs_prefilter"] = {
        "min_genes_per_cell": SCDRS_MIN_GENES_PER_CELL,
        "min_cells_per_gene": SCDRS_MIN_CELLS_PER_GENE,
        "n_pcs": 20,
        "n_neighbors": 15,
        "invoke_perform_downstream_with":
            "--flag-filter-data False --flag-raw-count False",
    }

    # ---- Sanitize obs for h5py serialization ----
    # CELLxGENE deposits ship several object-dtype obs columns with
    # mixed str/NaN that h5py can't write (TypeError: "Can't implicitly
    # convert non-string objects to strings"). Coerce every object-dtype
    # obs column to str; NaN becomes "NA". Applied generically — any
    # future schema change on HCA Gut would trip the same crash
    # otherwise, and whack-a-mole across atlases is the pattern we're
    # avoiding (Pan-GI already paid for this on obs['batch']).
    coerced = []
    for col in adata.obs.columns:
        if pd.api.types.is_object_dtype(adata.obs[col]):
            adata.obs[col] = adata.obs[col].fillna("NA").astype(str)
            coerced.append(col)
    if coerced:
        logger.info(
            "[driver] coerced %d object-dtype obs cols to str for h5py: %s",
            len(coerced), coerced,
        )

    adata.write_h5ad(a.out_h5ad)
    print(
        f"[driver] wrote {a.out_h5ad}: {adata.n_obs} cells x "
        f"{adata.n_vars} genes; obsp['connectivities'] populated "
        f"(n_pcs=20, n_neighbors=15)."
    )

    # ---- Cov file build ----
    # Verify every DECISIONS-locked cov categorical is present before we
    # fail halfway through building dummies.
    missing_cov_cols = [c for c in HCA_COV_CATEGORICALS
                        if c not in adata.obs.columns]
    if missing_cov_cols:
        raise KeyError(
            f"[driver] HCA Gut obs is missing cov columns {missing_cov_cols}. "
            f"DECISIONS locks cov = {list(HCA_COV_CATEGORICALS)}; the loader "
            f"(load_hca_gut.load) must carry these through."
        )

    # Sample-in-donor nesting gate. Treat `batch` as the biopsy/library
    # proxy (one-10x-run granularity): the gate requires each batch to
    # map to exactly one donor so that batch dummies fully absorb donor-
    # level variation and donor can be safely omitted from the cov.
    # Same logic as TAURUS/Pan-GI's sample-in-donor check — HCA just
    # uses `batch` because it ships no finer sample id.
    nest_df = pd.DataFrame({
        "batch": adata.obs["batch"].astype(str).values,
        "donor": adata.obs[DONOR_COL].astype(str).values,
    })
    donors_per_batch = nest_df.groupby("batch")["donor"].nunique()
    multi_donor_batches = donors_per_batch[donors_per_batch > 1]
    if len(multi_donor_batches) > 0:
        raise SystemExit(
            f"[driver] Batch-in-donor nesting violated: "
            f"{len(multi_donor_batches)} batch(es) span multiple donors. "
            f"Donor-level variation would NOT be absorbed by batch dummies "
            f"and the cov would under-correct. Top 10 offenders "
            f"(batch -> n_donors): {multi_donor_batches.head(10).to_dict()}. "
            f"Fix: either add a biopsy-level column that nests in donor, "
            f"or restore an explicit donor dummy block (and accept the "
            f"rank-deficiency hit) in HCA_COV_CATEGORICALS."
        )
    print(
        f"[driver] batch-in-donor nesting gate passed: "
        f"{donors_per_batch.size} batches across "
        f"{nest_df['donor'].nunique()} donors, each batch -> one donor "
        f"(donor dummies REDUNDANT with batch dummies; donor DROPPED)."
    )

    X = _raw_counts_matrix(adata)
    n_counts = np.asarray(X.sum(axis=1)).ravel()
    n_genes  = np.asarray((X > 0).sum(axis=1)).ravel()

    # Pre-expand each cov categorical into explicit float 0.0/1.0
    # dummies with a column-name prefix so no two cov columns can
    # collide. ``pd.get_dummies(..., prefix=<name>).astype(float)``
    # bypasses scDRS's category2dummy bool-dummy path that upcasts
    # ``df_cov.values`` to object dtype.
    dummy_frames = []
    dummy_counts = {}
    for col in HCA_COV_CATEGORICALS:
        vals = adata.obs[col].astype(str).values
        dummies = pd.get_dummies(
            pd.Series(vals, index=adata.obs_names, name=col),
            prefix=col,
        ).astype(float)
        dummy_frames.append(dummies)
        dummy_counts[col] = dummies.shape[1]
    logger.info(
        "[driver] cov dummy block widths: %s (total %d columns)",
        dummy_counts, sum(dummy_counts.values()),
    )

    depth = pd.DataFrame({
        "const":        1.0,
        "log_n_genes":  np.log1p(n_genes),
        "log_n_counts": np.log1p(n_counts),
    }, index=adata.obs_names)
    cov = pd.concat([depth, *dummy_frames], axis=1)
    cov.index.name = "cell"

    # Belt-and-suspenders: every column must be a numeric dtype so
    # scDRS's np.asarray(df_cov.values, dtype=float) doesn't hit
    # dtype('O'). Regression here would be caught in the driver rather
    # than deep inside scDRS.
    non_numeric = [c for c in cov.columns
                   if not pd.api.types.is_numeric_dtype(cov[c])]
    if non_numeric:
        raise SystemExit(
            f"[driver] cov has non-numeric columns after pre-expansion: "
            f"{non_numeric}. scDRS will crash on df_cov.values upcast."
        )

    cov.to_csv(a.out_cov, sep="\t")
    print(
        f"[driver] wrote {a.out_cov}: {len(cov.columns)} numeric cols "
        f"(const + 2 depth + "
        f"{sum(dummy_counts.values())} float dummies from "
        f"{list(HCA_COV_CATEGORICALS)}) — donor DROPPED (batch nests in "
        f"donor); no disease/health per no-disease-covariate rule."
    )


if __name__ == "__main__":
    main()
