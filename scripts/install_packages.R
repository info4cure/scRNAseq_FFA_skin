#!/usr/bin/env Rscript

# Installs R/Bioconductor packages not available via conda
# Run this script after creating and activating your conda env.

message("🔧 Installing missing CRAN & Bioconductor packages...")

# Install BiocManager if not present
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# ---- Bioconductor core ----
BiocManager::install(c(
  "SingleCellExperiment",
  "batchelor",
  "limma",
  "scran",
  "scater",
  "org.Hs.eg.db"
), ask = FALSE, update = TRUE)

# ---- scRNA-seq specific ----
BiocManager::install(c(
  "monocle3",
  "GSVA",
  "clusterProfiler",
  "CellChat"
), ask = FALSE, update = TRUE)

# ---- Extra CRAN packages ----
cran_pkgs <- c("patchwork", "cowplot", "RColorBrewer", "irlba")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

message("✅ All packages installed or already present.")
