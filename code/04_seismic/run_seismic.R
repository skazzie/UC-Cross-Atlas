#!/usr/bin/env Rscript
# run_seismic.R — seismicGWAS pipeline for one (atlas, GWAS, tier).
#
# Writes <out-dir>/<atlas>_<gwas>_<tier>.tsv with columns:
#   cell_type, pvalue, FDR, n_cells
#
# seismicGWAS 1.0.0's get_ct_trait_associations returns only cell_type,
# pvalue, FDR — no coefficient/se/n_genes. This driver targets that API.
#
# Spec: code/04_seismic/README.md and DECISIONS.md.
#
# CLI (either --sce-rds OR --h5ad-path is required):
#   # Preferred on the GCP VM (bypasses zellkonverter/basilisk):
#   Rscript code/04_seismic/run_seismic.R \
#     --atlas garrido_trigo --gwas delange --tier broad \
#     --sce-rds data/atlases/garrido_trigo.sce.rds \
#     --magma-z results/magma/delange_gene_z.tsv \
#     --out-dir results/seismic \
#     --seed 42
#
#   # Legacy (needs a working zellkonverter/basilisk install):
#   Rscript code/04_seismic/run_seismic.R \
#     --atlas garrido_trigo --gwas delange --tier broad \
#     --h5ad-path data/atlases/garrido_trigo.h5ad \
#     ...
#
#   Optional (off by default so the headline TSV writes reliably):
#     --run-retest         re-run regression, spearman(-log10 p) >= 0.999
#     --run-permutations   M shuffles of gene-Z, per-perm -log10(p) feather
#     --permutations INT   M value when --run-permutations is set (default 1000)
#
# Behaviour:
#   1. Load SCE — from --sce-rds via readRDS (pure R, no basilisk) OR
#      from --h5ad-path via zellkonverter::readH5AD. Validate
#      cell_type_<tier> column.
#   2. calc_specificity on the SCE (recomputed every run — a few minutes;
#      no on-disk cache, which lost gene rownames through feather).
#   3. Regression via get_ct_trait_associations with magma_gene_col=SYMBOL,
#      magma_z_col=ZSTAT (our MAGMA tables use those column names). The
#      1.0.0 API takes no `confounders` argument — passing one errors.
#   4. Write headline TSV: cell_type, pvalue, FDR, n_cells.
#   5. Optional (--run-permutations): shuffle gene-Z M times, save
#      per-permutation -log10(pvalue) as long-format feather.
#   6. Optional (--run-retest): re-run regression with identical inputs;
#      assert spearman correlation of -log10(pvalue) >= 0.999.

suppressPackageStartupMessages({
  library(optparse)
  library(seismicGWAS)
  library(SingleCellExperiment)
  library(arrow)
})

# ---- CLI ------------------------------------------------------------------

option_list <- list(
  make_option("--atlas",        type = "character"),
  make_option("--gwas",         type = "character"),
  make_option("--tier",         type = "character", default = "broad",
              help = "broad or fine"),
  make_option("--sce-rds",      type = "character",
              help = "Pre-built SingleCellExperiment .rds (preferred; pure R load)."),
  make_option("--h5ad-path",    type = "character",
              help = "h5ad file (requires working zellkonverter/basilisk)."),
  make_option("--magma-z",      type = "character"),
  make_option("--out-dir",      type = "character", default = "results/seismic"),
  make_option("--permutations", type = "integer",   default = 1000L,
              help = "M value used only when --run-permutations is set."),
  make_option("--seed",         type = "integer",   default = 42L),
  make_option("--run-permutations", action = "store_true", default = FALSE,
              help = "Run M gene-Z shuffles (default off; headline TSV only)."),
  make_option("--run-retest", action = "store_true", default = FALSE,
              help = "Run test-retest reproducibility gate (default off)."),
  # Kept for backward compat with existing slurm callers; no effect.
  make_option("--skip-permutations", action = "store_true", default = FALSE,
              help = "Deprecated no-op; permutations are opt-in via --run-permutations."),
  make_option("--recompute-spec", action = "store_true", default = FALSE,
              help = "Deprecated no-op; specificity is always recomputed.")
)

opt <- parse_args(OptionParser(option_list = option_list))

required <- c("atlas", "gwas", "tier", "magma-z")
missing <- required[vapply(required, function(x) is.null(opt[[x]]), logical(1))]
if (length(missing) > 0) {
  stop("Missing required arguments: ", paste(missing, collapse = ", "))
}
if (is.null(opt[["sce-rds"]]) && is.null(opt[["h5ad-path"]])) {
  stop("Must provide either --sce-rds (preferred) or --h5ad-path.")
}
if (!is.null(opt[["sce-rds"]]) && !is.null(opt[["h5ad-path"]])) {
  stop("Provide exactly one of --sce-rds or --h5ad-path, not both.")
}

require_path <- function(p, descr) {
  if (!file.exists(p)) {
    message(sprintf("ERROR: missing %s: %s", descr, p))
    quit(status = 2)
  }
}
if (!is.null(opt[["sce-rds"]]))  require_path(opt[["sce-rds"]],  "SCE .rds")
if (!is.null(opt[["h5ad-path"]])) require_path(opt[["h5ad-path"]], "h5ad")
require_path(opt[["magma-z"]],   "MAGMA gene-Z TSV")

dir.create(opt[["out-dir"]], recursive = TRUE, showWarnings = FALSE)
perm_dir <- file.path(opt[["out-dir"]], "permutations")
if (isTRUE(opt[["run-permutations"]])) {
  dir.create(perm_dir, recursive = TRUE, showWarnings = FALSE)
}

t0 <- Sys.time()
message(sprintf("[seismic] atlas=%s gwas=%s tier=%s seed=%d run_perm=%s run_retest=%s",
                opt$atlas, opt$gwas, opt$tier, opt$seed,
                isTRUE(opt[["run-permutations"]]),
                isTRUE(opt[["run-retest"]])))

# ---- 1. Load SCE ---------------------------------------------------------

if (!is.null(opt[["sce-rds"]])) {
  message(sprintf("[seismic] readRDS %s", opt[["sce-rds"]]))
  sce <- readRDS(opt[["sce-rds"]])
} else {
  # zellkonverter drags in basilisk, which on our sudo-less VM tries to
  # source-build Python 3.14 and dies. Only load it on the path where we
  # actually need it.
  message(sprintf("[seismic] zellkonverter::readH5AD %s", opt[["h5ad-path"]]))
  suppressPackageStartupMessages(library(zellkonverter))
  sce <- readH5AD(opt[["h5ad-path"]])
}
ct_col <- paste0("cell_type_", opt$tier)
if (!(ct_col %in% colnames(colData(sce)))) {
  stop(sprintf("colData(sce) has no '%s' column. Available: %s",
               ct_col, paste(colnames(colData(sce)), collapse = ", ")))
}
message(sprintf("[seismic] loaded %d cells x %d genes",
                ncol(sce), nrow(sce)))

# ---- 2. Specificity -------------------------------------------------------

# No on-disk cache: the prior feather round-trip stored specificity as a
# long-format data.frame with a `gene` column, but on cache hit the code
# assigned that data.frame straight to spec_obj — losing the wide-matrix
# shape and rownames that get_ct_trait_associations needs for gene
# overlap. Result was a 1-gene overlap and check_overlap failure.
# Recomputing is a few minutes and cheap relative to a whole failed run.
message("[seismic] calc_specificity ...")
t_spec <- Sys.time()
spec_obj <- calc_specificity(sce, ct_label_col = ct_col)
message(sprintf("[seismic] specificity done in %s",
                format(round(Sys.time() - t_spec, 1))))

# ---- 3. Regression --------------------------------------------------------

# seismicGWAS 1.0.0 signature:
#   get_ct_trait_associations(sscore, magma,
#     magma_gene_col = "GENE", magma_z_col = "ZSTAT", model = "linear")
# There is NO `confounders` argument in 1.0.0 — passing one errors out
# with "unused argument". Our MAGMA gene-Z tables use SYMBOL/ZSTAT
# columns (see code/01_magma/make_scdrs_gs.py), so we override
# magma_gene_col to match. Confounder adjustment tracked separately in
# DECISIONS.md; not enforceable through this API version.

magma_z <- read.table(opt[["magma-z"]], header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)

run_regression <- function(spec, mz) {
  get_ct_trait_associations(spec, mz,
                            magma_gene_col = "SYMBOL",
                            magma_z_col    = "ZSTAT")
}

res <- run_regression(spec_obj, magma_z)

# seismicGWAS 1.0.0 returns only cell_type, pvalue, FDR.
n_cells_per_ct <- as.integer(table(colData(sce)[[ct_col]]))
names(n_cells_per_ct) <- names(table(colData(sce)[[ct_col]]))

headline <- data.frame(
  cell_type = as.character(res$cell_type),
  pvalue    = as.numeric(res$pvalue),
  FDR       = as.numeric(res$FDR),
  n_cells   = n_cells_per_ct[as.character(res$cell_type)],
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

# ---- 4. Optional: M permutations of MAGMA gene-Z --------------------------

# seismic 1.0.0 has no coefficient in its output, so the null we track
# per permutation is -log10(pvalue) (mirrors the score column the
# laptop driver synthesizes for cross-method comparison).
if (isTRUE(opt[["run-permutations"]])) {
  set.seed(opt$seed)
  M <- opt$permutations
  message(sprintf("[seismic] running %d permutations", M))
  perm_rows <- vector("list", M)
  for (i in seq_len(M)) {
    mz_perm <- magma_z
    mz_perm$z <- sample(magma_z$z)
    r <- run_regression(spec_obj, mz_perm)
    perm_rows[[i]] <- data.frame(
      permutation_idx = i,
      cell_type       = as.character(r$cell_type),
      neg_log10_p     = -log10(as.numeric(r$pvalue)),
      stringsAsFactors = FALSE
    )
    if (i %% 100 == 0) {
      message(sprintf("[seismic]   ... permutation %d/%d", i, M))
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

# ---- 5. Optional: test-retest gate ---------------------------------------

if (isTRUE(opt[["run-retest"]])) {
  res_rerun <- run_regression(spec_obj, magma_z)
  ord  <- order(res$cell_type)
  ord2 <- order(res_rerun$cell_type)
  # 1.0.0 exposes pvalue/FDR only — use -log10(pvalue) as the
  # reproducibility signal (the coefficient-based gate isn't available).
  score1 <- -log10(as.numeric(res$pvalue))[ord]
  score2 <- -log10(as.numeric(res_rerun$pvalue))[ord2]
  rho <- suppressWarnings(cor(score1, score2, method = "spearman"))
  if (is.na(rho) || rho < 0.999) {
    message(sprintf("[seismic] FAIL: test-retest spearman(-log10 p)=%.6f < 0.999", rho))
    quit(status = 1)
  }
  message(sprintf("[seismic] test-retest PASS: spearman(-log10 p)=%.6f", rho))
}

elapsed <- format(round(Sys.time() - t0, 1))
message(sprintf("[seismic] DONE in %s", elapsed))
