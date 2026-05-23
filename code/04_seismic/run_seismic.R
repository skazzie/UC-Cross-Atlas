#!/usr/bin/env Rscript
# run_seismic.R — seismicGWAS regime-1 pipeline for one (atlas, GWAS, tier).
#
# Spec: code/04_seismic/README.md and DECISIONS.md.
#
# CLI:
#   Rscript code/04_seismic/run_seismic.R \
#     --atlas garrido_trigo --gwas delange --tier broad \
#     --h5ad-path data/atlases/garrido_trigo.h5ad \
#     --magma-z results/magma/delange_gene_z.tsv \
#     --out-dir results/seismic \
#     --permutations 1000 \
#     --seed 42
#
# Behaviour:
#   1. Load h5ad as SCE via zellkonverter; validate cell_type_<tier> column.
#   2. Compute (and cache) per-atlas specificity in long format.
#   3. Regression with explicit confounders (gene_length_log, ld_score,
#      transcript_count) — overrides package defaults if they differ.
#   4. Write headline TSV: cell_type, coefficient, se, pvalue, n_genes, n_cells.
#   5. M=permutations of the gene-Z vector; per-permutation coefficients
#      saved as long-format feather.
#   6. Test-retest gate: re-run regression with identical inputs; assert
#      spearman correlation of coefficients >= 0.999.

suppressPackageStartupMessages({
  library(optparse)
  library(seismicGWAS)
  library(SingleCellExperiment)
  library(zellkonverter)
  library(arrow)
})

# ---- CLI ------------------------------------------------------------------

option_list <- list(
  make_option("--atlas",        type = "character"),
  make_option("--gwas",         type = "character"),
  make_option("--tier",         type = "character", default = "broad",
              help = "broad or fine"),
  make_option("--h5ad-path",    type = "character"),
  make_option("--magma-z",      type = "character"),
  make_option("--out-dir",      type = "character", default = "results/seismic"),
  make_option("--permutations", type = "integer",   default = 1000L),
  make_option("--seed",         type = "integer",   default = 42L),
  make_option("--recompute-spec", action = "store_true", default = FALSE,
              help = "Force recomputation of the cached specificity table."),
  make_option("--skip-permutations", action = "store_true", default = FALSE,
              help = "Skip the M=permutations Brown null draws (debug only).")
)

opt <- parse_args(OptionParser(option_list = option_list))

required <- c("atlas", "gwas", "tier", "h5ad-path", "magma-z")
missing <- required[vapply(required, function(x) is.null(opt[[x]]), logical(1))]
if (length(missing) > 0) {
  stop("Missing required arguments: ", paste(missing, collapse = ", "))
}

require_path <- function(p, descr) {
  if (!file.exists(p)) {
    message(sprintf("ERROR: missing %s: %s", descr, p))
    quit(status = 2)
  }
}
require_path(opt[["h5ad-path"]], "h5ad")
require_path(opt[["magma-z"]],   "MAGMA gene-Z TSV")

dir.create(opt[["out-dir"]], recursive = TRUE, showWarnings = FALSE)
spec_dir <- file.path(opt[["out-dir"]], "..", "..",
                      "code", "04_seismic", "specificity_long")
dir.create(spec_dir, recursive = TRUE, showWarnings = FALSE)
perm_dir <- file.path(opt[["out-dir"]], "permutations")
dir.create(perm_dir, recursive = TRUE, showWarnings = FALSE)

t0 <- Sys.time()
message(sprintf("[seismic] atlas=%s gwas=%s tier=%s perm=%d seed=%d",
                opt$atlas, opt$gwas, opt$tier,
                opt$permutations, opt$seed))

# ---- 1. Load h5ad as SCE --------------------------------------------------

sce <- readH5AD(opt[["h5ad-path"]])
ct_col <- paste0("cell_type_", opt$tier)
if (!(ct_col %in% colnames(colData(sce)))) {
  stop(sprintf("colData(sce) has no '%s' column. Available: %s",
               ct_col, paste(colnames(colData(sce)), collapse = ", ")))
}
message(sprintf("[seismic] loaded %d cells x %d genes",
                ncol(sce), nrow(sce)))

# ---- 2. Specificity (cached per atlas x tier) -----------------------------

spec_path <- file.path(
  spec_dir,
  paste0(opt$atlas, "_", opt$tier, "_specificity.feather")
)

if (file.exists(spec_path) && !opt[["recompute-spec"]]) {
  message(sprintf("[seismic] specificity cache hit: %s", spec_path))
  spec_long <- arrow::read_feather(spec_path)
  spec_obj <- spec_long  # downstream API may differ; see note below
} else {
  message("[seismic] computing specificity (this can take several minutes)")
  spec_obj <- calc_specificity(sce, ct_label_col = ct_col)
  # Reshape to long format: atlas, granularity, cell_type, gene, specificity.
  spec_mat <- as.matrix(spec_obj)
  spec_long <- data.frame(
    atlas       = opt$atlas,
    granularity = opt$tier,
    cell_type   = rep(colnames(spec_mat), each = nrow(spec_mat)),
    gene        = rep(rownames(spec_mat), times = ncol(spec_mat)),
    specificity = as.numeric(spec_mat),
    stringsAsFactors = FALSE
  )
  spec_long <- spec_long[spec_long$specificity != 0, ]
  arrow::write_feather(spec_long, spec_path)
}

# ---- 3. Regression --------------------------------------------------------

# Confounders verification — DECISIONS.md requires the explicit triple
# (gene length, gene-gene LD, transcript count). If package defaults
# differ we override. See code/04_seismic/README.md §Confounders.
CONFOUNDERS <- c("gene_length_log", "ld_score", "transcript_count")

magma_z <- read.table(opt[["magma-z"]], header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)

run_regression <- function(spec, mz) {
  get_ct_trait_associations(spec, mz, confounders = CONFOUNDERS)
}

res <- run_regression(spec_obj, magma_z)
n_cells_per_ct <- as.integer(table(colData(sce)[[ct_col]]))
names(n_cells_per_ct) <- names(table(colData(sce)[[ct_col]]))

headline <- data.frame(
  cell_type   = res$cell_type,
  coefficient = res$coefficient,
  se          = res$se,
  pvalue      = res$pvalue,
  n_genes     = res$n_genes,
  n_cells     = n_cells_per_ct[as.character(res$cell_type)],
  stringsAsFactors = FALSE
)
headline_path <- file.path(
  opt[["out-dir"]],
  paste0(opt$atlas, "_", opt$gwas, "_", opt$tier, ".tsv")
)
write.table(headline, headline_path, sep = "\t",
            row.names = FALSE, quote = FALSE)
message(sprintf("[seismic] wrote %s (n_cell_types=%d)",
                headline_path, nrow(headline)))

# ---- 4. M permutations of MAGMA gene-Z ------------------------------------

if (!opt[["skip-permutations"]]) {
  set.seed(opt$seed)
  message(sprintf("[seismic] running %d permutations", opt$permutations))
  perm_rows <- vector("list", opt$permutations)
  for (i in seq_len(opt$permutations)) {
    mz_perm <- magma_z
    mz_perm$z <- sample(magma_z$z)
    r <- run_regression(spec_obj, mz_perm)
    perm_rows[[i]] <- data.frame(
      permutation_idx = i,
      cell_type       = r$cell_type,
      coefficient     = r$coefficient,
      stringsAsFactors = FALSE
    )
    if (i %% 100 == 0) {
      message(sprintf("[seismic]   ... permutation %d/%d", i, opt$permutations))
    }
  }
  perm_long <- do.call(rbind, perm_rows)
  perm_path <- file.path(
    perm_dir,
    paste0(opt$atlas, "_", opt$gwas, "_", opt$tier, "_permnulls.feather")
  )
  arrow::write_feather(perm_long, perm_path)
  message(sprintf("[seismic] wrote %s (%d rows)", perm_path, nrow(perm_long)))
}

# ---- 5. Test-retest gate --------------------------------------------------

res_rerun <- run_regression(spec_obj, magma_z)
ord <- order(res$cell_type)
ord2 <- order(res_rerun$cell_type)
rho <- cor(res$coefficient[ord], res_rerun$coefficient[ord2],
           method = "spearman")
if (is.na(rho) || rho < 0.999) {
  message(sprintf("[seismic] FAIL: test-retest spearman=%.6f < 0.999", rho))
  quit(status = 1)
}
message(sprintf("[seismic] test-retest PASS: spearman=%.6f", rho))

elapsed <- format(round(Sys.time() - t0, 1))
message(sprintf("[seismic] DONE in %s", elapsed))
