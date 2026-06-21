# Isolated seismicGWAS install — skips scOntoMatch (Mingw-w64 crash trigger)
# and zellkonverter (Rtools-source-build needed). seismicGWAS itself does NOT
# depend on zellkonverter; that was only in code/04_seismic/run_seismic.R as
# the h5ad loader. We bypass it via Python-side flat-file dump (see
# scripts/garrido_to_flat.py and code/04_seismic/run_seismic_laptop.R).
#
# Usage:
#   conda run -n uc-cross-atlas-r Rscript scripts/install_seismic_only.R <SEISMIC_LOCAL_PATH>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: install_seismic_only.R <path-to-seismicGWAS-clone>")
seismic_path <- args[[1]]
if (!dir.exists(seismic_path)) stop(sprintf("Not a directory: %s", seismic_path))

repos <- "https://cloud.r-project.org"

cat("=== Phase 1: ensure CRAN deps ===\n")
cran_pkgs <- c("data.table", "ggrepel", "magrittr", "speedglm")
need <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(need)) {
  cat(sprintf("  installing: %s\n", paste(need, collapse = ", ")))
  install.packages(need, repos = repos)
} else {
  cat("  all present\n")
}

cat("=== Phase 2: ensure Bioc deps ===\n")
bioc_pkgs <- c("SingleCellExperiment", "SummarizedExperiment",
               "AnnotationDbi", "org.Hs.eg.db")
need <- setdiff(bioc_pkgs, rownames(installed.packages()))
if (length(need)) {
  cat(sprintf("  installing via BiocManager: %s\n", paste(need, collapse = ", ")))
  BiocManager::install(need, update = FALSE, ask = FALSE)
} else {
  cat("  all present\n")
}

cat("=== Phase 3: install seismicGWAS from local clone ===\n")
if (!requireNamespace("seismicGWAS", quietly = TRUE)) {
  remotes::install_local(seismic_path, dependencies = FALSE, upgrade = "never")
} else {
  cat("  seismicGWAS already installed; reinstalling\n")
  remove.packages("seismicGWAS")
  remotes::install_local(seismic_path, dependencies = FALSE, upgrade = "never")
}

cat("\n=== Final smoke test ===\n")
for (p in c("data.table", "ggrepel", "speedglm",
            "SingleCellExperiment", "SummarizedExperiment",
            "AnnotationDbi", "org.Hs.eg.db",
            "Matrix", "seismicGWAS")) {
  ok <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("  %-22s %s\n", p, ok))
}
cat("\n[install_seismic_only] done\n")
