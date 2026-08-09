"""Driver: load TAURUS-IBD, apply UC×colonic×baseline filter, precompute
neighbors, and write h5ad + partial covariate file.

Source pin (DECISIONS 16): Zenodo v3 ``10.5281/zenodo.14007626``,
``TAURUS_raw_counts_annotated_final.h5ad``, md5
``c1bd13b92cacb164a401c6c4a4e7912c``.

DATA VERSION CAVEAT — disregard any TAURUS metadata dated before
2024-10-30. Earlier Zenodo revisions ship incorrect donor/sample
metadata that will fail the loader's Supp-Table-1B donor invariant.
Use v3 (posted 2024-10-30 or later).

Covariate structure: depth (log_n_genes, log_n_counts) + SAMPLE dummies
only. NO donor column and NO disease/health covariate. Reason to drop
donor: in TAURUS each sample_id nests within exactly one donor (donor
UC1 → sample CID003356-1, ...), so donor + sample dummies + const are
linearly dependent (~23 redundant columns → rank 54 of 77 →
np.linalg.solve singular in scDRS's covariate regression). Sample-
level dummies fully absorb donor-level variation, consistent with
sample-level batch correction as used for Smillie. NO disease/health
covariate per the "no-disease-covariate" rule (that would scrub the
signal we test).

Sample dummy prefers ``obs['sample_id']`` (true biopsy-level unit that
nests within donor); if the loader didn't carry it through, we fall
back to ``obs['region']`` (=Site: Ascending_Colon / Descending_Colon /
Rectum / Sigmoid) but the nesting assertion below will fail loud in
that case — Site crosses donors, so dropping donor while using Site as
the sample dummy would under-correct.

Neighbors precomputed inline (min_genes=250 / min_cells=50 scDRS-replica
filter, then ``sc.pp.pca(n_comps=20)`` + ``sc.pp.neighbors(n_neighbors=15,
n_pcs=20)``) so the shipped h5ad works with ``scdrs perform-downstream
--flag-filter-data False --flag-raw-count False`` — same recipe as
add_neighbors.py. ``uns['scdrs_prefilter']`` stamped for provenance.
"""

import argparse
import numpy as np
import pandas as pd
import scanpy as sc

from load_taurus import load

# scDRS 1.0.2 load-time filter defaults (scdrs/util.py:82-88). Match here
# so on-disk cell set == what --flag-filter-data True would produce.
SCDRS_MIN_GENES_PER_CELL = 250
SCDRS_MIN_CELLS_PER_GENE = 50

# Extra gate on top of the loader's 22-donor Supp-Table-1B invariant:
# never ship a TAURUS h5ad with < 8 UC donors. This is a defense against
# an unnoticed loader change or a wrong Zenodo version silently thinning
# the cohort — 8 is well below the 22 the loader itself asserts, so
# tripping this means something has gone badly wrong upstream.
MIN_UC_DONORS_POST_FILTER = 8


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--h5ad", required=True,
                   help="TAURUS pooled h5ad (TAURUS_raw_counts_annotated_final.h5ad)")
    # VM layout: outputs live under scratch/data/atlases/ (not
    # ../../data/atlases/ — that's the code-tree repo path, no
    # room for a 12 GB TAURUS write). Defaults resolve relative to
    # `code/02_atlas_prep/`, which is where run_*.py is invoked from.
    p.add_argument("--out-h5ad",
                   default="../../scratch/data/atlases/taurus_uc_colonic_baseline.h5ad")
    p.add_argument("--out-cov",
                   default="../../scratch/data/atlases/taurus_uc_colonic_baseline.cov.tsv")
    a = p.parse_args()

    adata = load(a.h5ad)  # apply_v1_filter=True default; runs 22-donor gate

    n_uc_donors = int(adata.obs["donor"].nunique())
    if n_uc_donors < MIN_UC_DONORS_POST_FILTER:
        raise SystemExit(
            f"[driver] TAURUS UC-donor gate failed: got {n_uc_donors}, "
            f"required >= {MIN_UC_DONORS_POST_FILTER}. Confirm the Zenodo "
            f"h5ad is v3 or later (2024-10-30+); pre-Oct-30-2024 revisions "
            f"have wrong donor/sample metadata."
        )
    print(f"[driver] UC donors post-filter: {n_uc_donors} "
          f"(>= {MIN_UC_DONORS_POST_FILTER} gate)")

    # scDRS-replica pre-filter — matches add_neighbors.py so
    # perform-downstream can be invoked with --flag-filter-data False
    # without a cell-set mismatch. filter_cells counts nonzero-X genes
    # per cell; log1p(0)==0, so the cell set is identical to what
    # scDRS would produce on raw counts.
    n0_cells, n0_genes = adata.n_obs, adata.n_vars
    sc.pp.filter_cells(adata, min_genes=SCDRS_MIN_GENES_PER_CELL)
    sc.pp.filter_genes(adata, min_cells=SCDRS_MIN_CELLS_PER_GENE)
    print(
        f"[driver] scDRS-replica filter: "
        f"{n0_cells}->{adata.n_obs} cells (min_genes={SCDRS_MIN_GENES_PER_CELL}), "
        f"{n0_genes}->{adata.n_vars} genes (min_cells={SCDRS_MIN_CELLS_PER_GENE})"
    )

    # Precompute PCA + kNN so obsp['connectivities'] ships with the h5ad
    # and scdrs perform-downstream skips its ~45-min rebuild. Params
    # match scDRS group-analysis defaults (knn_n_pcs=20, knn_n_neighbors=15).
    sc.pp.pca(adata, n_comps=20)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=20)
    assert "connectivities" in adata.obsp, "sc.pp.neighbors did not populate obsp"

    # Provenance stamp — mirrors add_neighbors.py so a future scDRS
    # default drift is catchable by comparing this against
    # scdrs.util.load_h5ad's filter constants.
    adata.uns["scdrs_prefilter"] = {
        "min_genes_per_cell": SCDRS_MIN_GENES_PER_CELL,
        "min_cells_per_gene": SCDRS_MIN_CELLS_PER_GENE,
        "n_pcs": 20,
        "n_neighbors": 15,
        "invoke_perform_downstream_with":
            "--flag-filter-data False --flag-raw-count False",
    }

    adata.write_h5ad(a.out_h5ad)
    print(
        f"[driver] wrote {a.out_h5ad}: {adata.n_obs} cells x {adata.n_vars} genes, "
        f"obsp['connectivities'] populated (n_pcs=20, n_neighbors=15)."
    )

    # Depth + SAMPLE dummies only (no donor, no disease/health). See
    # module docstring — donor is dropped because sample nests in donor,
    # so donor dummies are a linear combination of sample dummies and
    # the joint design is singular (scDRS's covariate regression uses
    # np.linalg.solve, which crashes on the resulting rank-deficient
    # normal-equations matrix).
    X = adata.layers["counts"]
    n_counts = np.asarray(X.sum(axis=1)).ravel()
    n_genes  = np.asarray((X > 0).sum(axis=1)).ravel()
    if "sample_id" in adata.obs.columns:
        sample_col_name = "sample_id"
        sample_values = adata.obs["sample_id"].astype(str).values
    else:
        sample_col_name = "region"  # =Site; site-level, not biopsy-level
        sample_values = adata.obs["region"].astype(str).values

    # Nesting gate: each sample maps to exactly one donor. This must
    # hold for the "drop donor, keep sample" strategy to be valid —
    # otherwise sample dummies leave donor variation uncorrected and
    # the covariate matrix under-corrects rather than over-parametrises.
    nest_df = pd.DataFrame({
        "sample": sample_values,
        "donor":  adata.obs["donor"].astype(str).values,
    })
    donors_per_sample = nest_df.groupby("sample")["donor"].nunique()
    multi_donor_samples = donors_per_sample[donors_per_sample > 1]
    if len(multi_donor_samples) > 0:
        raise SystemExit(
            f"[driver] Sample-in-donor nesting violated: "
            f"{len(multi_donor_samples)} sample(s) span multiple donors — "
            f"cannot safely drop donor from cov. sample_col='{sample_col_name}'. "
            f"First 10 offenders (sample -> n_donors): "
            f"{multi_donor_samples.head(10).to_dict()}. "
            f"Fix: use a truly biopsy-level column for sample, or revert "
            f"to keeping donor and use a rank-tolerant regression."
        )
    print(
        f"[driver] Sample nesting check passed: {donors_per_sample.size} "
        f"unique samples across {nest_df['donor'].nunique()} donors, "
        f"each sample → exactly one donor (sample_col='{sample_col_name}')."
    )

    cov = pd.DataFrame({
        "const": 1,
        "log_n_genes":  np.log1p(n_genes),
        "log_n_counts": np.log1p(n_counts),
        "sample": sample_values,
    }, index=adata.obs_names)
    cov.index.name = "cell"
    cov.to_csv(a.out_cov, sep="\t")
    print(f"[driver] wrote {a.out_cov}: {list(cov.columns)} "
          f"(sample <- obs['{sample_col_name}']; donor DROPPED — absorbed "
          f"by sample dummies; no disease/health per no-disease-covariate "
          f"rule)")


if __name__ == "__main__":
    main()
