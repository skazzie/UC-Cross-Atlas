# Laptop-only seismicGWAS runner — bypasses zellkonverter via flat-file
# SCE reconstruction (Windows can't source-build zellkonverter without
# Rtools 4.5). For HB use code/04_seismic/run_seismic.R, which uses
# zellkonverter::readH5AD directly.
#
# Differences vs run_seismic.R:
#   - SCE comes from scripts/h5ad_to_flat.py output via sce_from_flat.R,
#     not from readH5AD.
#   - Accepts MULTIPLE magma-z files in one invocation so the
#     SCE-build + specificity step is amortized across both MHC variants
#     (the laptop dry run needs MHC-excluded + MHC-included for the
#     trade-off panel).
#   - No M=permutations Brown null draws (HB-only — would take hours on
#     laptop and isn't needed for the dry-run cross-method comparison).
#   - No test-retest gate (the dry run already validated the mechanic;
#     test-retest is a Phase 5 deliverable, not a Phase-1 validation).
#
# CLI:
#   conda run -n uc-cross-atlas-r Rscript code/04_seismic/run_seismic_laptop.R \
#     --flat-dir data/atlases/garrido_flat/ \
#     --atlas garrido \
#     --tier broad \
#     --magma-z results/magma/uc_delange_geneZ.tsv:delange \
#               results/magma/uc_delange_mhc_geneZ.tsv:delange_mhc \
#     --out-dir results/seismic/

suppressPackageStartupMessages({
  library(optparse)
  library(seismicGWAS)
})

option_list <- list(
  make_option("--flat-dir",  type = "character"),
  make_option("--atlas",     type = "character"),
  make_option("--tier",      type = "character", default = "broad",
              help = "broad or fine"),
  make_option("--magma-z",   type = "character", action = "callback",
              callback = function(parser, opt, value, ...) {
                # optparse doesn't support nargs='+' natively; accept comma-separated.
                strsplit(value, ",")[[1]]
              }),
  make_option("--out-dir",   type = "character", default = "results/seismic")
)

parser <- OptionParser(option_list = option_list)
parsed <- parse_args(parser, positional_arguments = TRUE)
opt <- parsed$options
pos <- parsed$args

# Accept --magma-z as comma-separated OR as repeated positional args after the flag.
# We also accept the positional convention: every positional arg ending in .tsv
# is treated as a magma-z spec.
magma_specs <- if (!is.null(opt[["magma-z"]])) opt[["magma-z"]] else pos[grepl("\\.tsv", pos)]
if (length(magma_specs) == 0) {
  stop("No magma-z specs given. Use --magma-z PATH:TAG[,PATH:TAG,...] or positional PATH:TAG args.")
}

for (req in c("flat-dir", "atlas", "tier")) {
  if (is.null(opt[[req]])) stop(sprintf("Missing required arg: --%s", req))
}

# Resolve paths.
here <- normalizePath(getwd())
flat_dir <- normalizePath(opt[["flat-dir"]], mustWork = TRUE)
out_dir  <- opt[["out-dir"]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Build SCE from flat files -----------------------------------------

source(file.path("code", "04_seismic", "sce_from_flat.R"))
sce <- build_sce_from_flat(flat_dir)
ct_col <- paste0("cell_type_", opt[["tier"]])
if (!(ct_col %in% colnames(SummarizedExperiment::colData(sce)))) {
  stop(sprintf("colData has no '%s'. Available: %s", ct_col,
               paste(colnames(SummarizedExperiment::colData(sce)), collapse = ", ")))
}

# ---- 2. Specificity (per atlas x tier; cheap to recompute on Garrido) -----

t0 <- Sys.time()
message(sprintf("[seismic] calc_specificity on %s/%s ...", opt[["atlas"]], opt[["tier"]]))
spec_obj <- calc_specificity(sce, ct_label_col = ct_col)
message(sprintf("[seismic] specificity done in %s",
                format(round(Sys.time() - t0, 1))))

# ---- 3. Loop over magma-z specs -------------------------------------------

# Each spec is PATH:TAG (e.g. "results/magma/uc_delange_geneZ.tsv:delange").
parse_spec <- function(s) {
  # rsplit on the LAST colon — Windows path drive letters use colons too.
  pos <- max(gregexpr(":", s)[[1]])
  list(path = substr(s, 1, pos - 1),
       tag  = substr(s, pos + 1, nchar(s)))
}

for (s in magma_specs) {
  spec <- parse_spec(s)
  magma_path <- normalizePath(spec$path, mustWork = TRUE)
  tag <- spec$tag
  message(sprintf("[seismic] === %s : %s ===", tag, magma_path))

  magma_z <- read.table(magma_path, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE)
  # seismicGWAS::translate_gene_ids defaults assume column names; check what we have.
  # Our prepare_gwas pipeline writes SYMBOL + ZSTAT (see code/01_magma/make_scdrs_gs.py).
  if (!all(c("SYMBOL", "ZSTAT") %in% colnames(magma_z))) {
    stop(sprintf("magma-z TSV %s missing SYMBOL or ZSTAT columns; got: %s",
                 magma_path, paste(colnames(magma_z), collapse = ", ")))
  }

  t1 <- Sys.time()
  message(sprintf("[seismic] get_ct_trait_associations (%s) ...", tag))
  # seismicGWAS::get_ct_trait_associations defaults to magma_gene_col='GENE'
  # and magma_z_col='ZSTAT'. The GENE column is matched against the SCE
  # rownames, which for Garrido are HGNC symbols (the loader's
  # ensembl_to_hgnc step strips Ensembl). So GENE in the magma data
  # must be the SYMBOL column from our gene-Z TSV, not the ENTREZ id.
  mz <- data.table::data.table(GENE = magma_z$SYMBOL, ZSTAT = magma_z$ZSTAT)
  res <- get_ct_trait_associations(spec_obj, mz)
  message(sprintf("[seismic] regression done in %s",
                  format(round(Sys.time() - t1, 1))))

  # seismicGWAS::get_ct_trait_associations returns ONLY cell_type, pvalue,
  # FDR (one-sided p from the linear-model slope). No coefficient/se. We
  # synthesize 'score' = -log10(pvalue) so the cross-method TSV layer
  # (result_loading.shared_cell_types + 06_concordance) has a comparable
  # signed magnitude column matching scDRS's z_mean. Larger-is-stronger
  # convention preserved (mirrors DEFAULT in 08_cross_method test fixtures).
  n_cells_per_ct <- as.integer(table(SummarizedExperiment::colData(sce)[[ct_col]]))
  names(n_cells_per_ct) <- names(table(SummarizedExperiment::colData(sce)[[ct_col]]))

  headline <- data.frame(
    cell_type   = as.character(res$cell_type),
    pvalue      = as.numeric(res$pvalue),
    fdr         = as.numeric(res$FDR),
    score       = -log10(as.numeric(res$pvalue)),
    n_cells     = n_cells_per_ct[as.character(res$cell_type)],
    stringsAsFactors = FALSE
  )
  out_path <- file.path(out_dir,
                        paste0(opt[["atlas"]], "_", tag, "_",
                               opt[["tier"]], ".tsv"))
  write.table(headline, out_path, sep = "\t",
              row.names = FALSE, quote = FALSE)
  message(sprintf("[seismic] wrote %s (n=%d cell types)",
                  out_path, nrow(headline)))
}

message(sprintf("[seismic] DONE total %s", format(round(Sys.time() - t0, 1))))
