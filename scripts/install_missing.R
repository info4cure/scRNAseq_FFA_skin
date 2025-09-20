#!/usr/bin/env Rscript

# Script to install missing R packages for scRNAseq_FFA_skin project
# Run with: Rscript scripts/install_missing.R

message("🔧 Installing missing R packages...")

# Ensure BiocManager and remotes are available
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")

# Core dependencies
pkgs_cran <- c("expm", "ggpubr")
pkgs_bioc <- c("ComplexHeatmap", "BiocNeighbors")

# Install from CRAN
for (pkg in pkgs_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Install from Bioconductor
BiocManager::install(pkgs_bioc, ask = FALSE, update = TRUE)

# Install CellChat from GitHub
if (!requireNamespace("CellChat", quietly = TRUE)) {
  remotes::install_github("sqjin/CellChat")
}

message("✅ All missing packages should now be installed.")
