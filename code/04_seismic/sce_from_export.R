#!/usr/bin/env Rscript
# sce_from_export.R — assemble a SingleCellExperiment from the flat-file
# export produced by h5ad_to_sce_export.py, then serialize to .rds.
#
# Pure R: reads counts.mtx via Matrix::readMM, no reticulate, no basilisk,
# no zellkonverter. This is the second stage of the seismic pipeline's
# h5ad-bypass path (Python exports flat → R builds SCE + saves .rds →
# run_seismic reads .rds).
#
# The primary and only assay is `counts` (raw). Seismic's calc_specificity
# expects raw counts; we intentionally do not synthesize logcounts here.
#
# CLI:
#   Rscript code/04_seismic/sce_from_export.R \
#     --export-dir data/atlases/garrido_trigo_export/ \
#     --out-rds    data/atlases/garrido_trigo.sce.rds

suppressPackageStartupMessages({
  library(optparse)
  library(Matrix)
  library(SingleCellExperiment)
})

build_sce_from_export <- function(export_dir) {
  export_dir <- normalizePath(export_dir, mustWork = TRUE)
  message(sprintf("[sce_from_export] reading %s", export_dir))

  barcodes <- readLines(file.path(export_dir, "barcodes.tsv"))
  genes    <- readLines(file.path(export_dir, "genes.tsv"))
  obs      <- read.table(file.path(export_dir, "obs.tsv"),
                         sep = "\t", header = TRUE, check.names = FALSE,
                         row.names = 1, stringsAsFactors = FALSE,
                         comment.char = "", quote = "")

  if (nrow(obs) != length(barcodes)) {
    stop(sprintf("obs.tsv rows (%d) != barcodes.tsv lines (%d)",
                 nrow(obs), length(barcodes)))
  }
  if (!identical(rownames(obs), barcodes)) {
    stop("obs.tsv row order does not match barcodes.tsv")
  }

  required_cols <- c("cell_type_broad", "cell_type_fine", "donor")
  missing_cols <- setdiff(required_cols, colnames(obs))
  if (length(missing_cols) > 0) {
    stop(sprintf("obs.tsv missing required columns: %s (have: %s)",
                 paste(missing_cols, collapse = ", "),
                 paste(colnames(obs), collapse = ", ")))
  }

  # counts.mtx is cells x genes (Python/AnnData orientation). SCE expects
  # genes x cells — transpose after reading.
  counts <- Matrix::readMM(file.path(export_dir, "counts.mtx"))
  if (nrow(counts) != length(barcodes) || ncol(counts) != length(genes)) {
    stop(sprintf("counts.mtx shape (%d x %d) doesn't match cells x genes (%d x %d)",
                 nrow(counts), ncol(counts),
                 length(barcodes), length(genes)))
  }
  counts <- Matrix::t(counts)
  rownames(counts) <- genes
  colnames(counts) <- barcodes

  sce <- SingleCellExperiment(
    assays  = list(counts = counts),
    colData = S4Vectors::DataFrame(obs),
    rowData = S4Vectors::DataFrame(gene = genes, row.names = genes)
  )
  message(sprintf("[sce_from_export] built SCE: %d cells x %d genes (assay: counts)",
                  ncol(sce), nrow(sce)))
  sce
}

if (sys.nframe() == 0L) {
  option_list <- list(
    make_option("--export-dir", type = "character",
                help = "Directory produced by h5ad_to_sce_export.py"),
    make_option("--out-rds", type = "character",
                help = "Destination .rds path for the SingleCellExperiment")
  )
  opt <- parse_args(OptionParser(option_list = option_list))

  for (req in c("export-dir", "out-rds")) {
    if (is.null(opt[[req]])) stop(sprintf("Missing required arg: --%s", req))
  }

  sce <- build_sce_from_export(opt[["export-dir"]])

  out_rds <- opt[["out-rds"]]
  dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)
  saveRDS(sce, out_rds)
  message(sprintf("[sce_from_export] wrote %s", normalizePath(out_rds)))
}
