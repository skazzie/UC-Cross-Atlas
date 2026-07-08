# Build a SingleCellExperiment from the flat-file dump produced by
# scripts/h5ad_to_flat.py. Used by run_seismic_laptop.R on the Windows
# laptop where zellkonverter can't be installed without Rtools.
#
# Usage (sourced, not run directly):
#   source("code/04_seismic/sce_from_flat.R")
#   sce <- build_sce_from_flat("data/atlases/garrido_flat/")

suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
})

build_sce_from_flat <- function(flat_dir) {
  flat_dir <- normalizePath(flat_dir, mustWork = TRUE)
  message(sprintf("[sce_from_flat] reading %s", flat_dir))

  # cells / genes metadata.
  cells <- read.table(file.path(flat_dir, "cells.tsv"),
                      sep = "\t", header = TRUE, check.names = FALSE,
                      row.names = 1, stringsAsFactors = FALSE,
                      comment.char = "")
  genes <- read.table(file.path(flat_dir, "genes.tsv"),
                      sep = "\t", header = TRUE, check.names = FALSE,
                      row.names = 1, stringsAsFactors = FALSE,
                      comment.char = "")
  message(sprintf("[sce_from_flat] cells=%d, genes=%d",
                  nrow(cells), nrow(genes)))

  # X is cells x genes in the mtx file (Python anndata convention).
  # SCE expects assays as genes x cells — transpose.
  X <- Matrix::readMM(file.path(flat_dir, "X.mtx"))
  if (nrow(X) != nrow(cells) || ncol(X) != nrow(genes)) {
    stop(sprintf("X.mtx shape (%d x %d) doesn't match cells x genes (%d x %d)",
                 nrow(X), ncol(X), nrow(cells), nrow(genes)))
  }
  X <- Matrix::t(X)
  rownames(X) <- rownames(genes)
  colnames(X) <- rownames(cells)

  assays_list <- list(logcounts = X)

  counts_path <- file.path(flat_dir, "counts.mtx")
  if (file.exists(counts_path)) {
    counts <- Matrix::readMM(counts_path)
    counts <- Matrix::t(counts)
    rownames(counts) <- rownames(genes)
    colnames(counts) <- rownames(cells)
    assays_list[["counts"]] <- counts
    message("[sce_from_flat] loaded counts layer")
  }

  sce <- SingleCellExperiment(
    assays = assays_list,
    colData = S4Vectors::DataFrame(cells),
    rowData = S4Vectors::DataFrame(genes)
  )
  message(sprintf("[sce_from_flat] built SCE: %d cells x %d genes (assays: %s)",
                  ncol(sce), nrow(sce),
                  paste(SummarizedExperiment::assayNames(sce), collapse = ", ")))
  sce
}
