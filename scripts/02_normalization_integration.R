############################################################
# 02_normalization_integration.R
# Normalization and integration of scRNA-seq samples
# Project: scRNAseq_FFA_skin
# Author: Juan Ruano, MD, PhD, MSc
# Date: 2025-09-18
############################################################

# 📌 Qué hace este script
# 
# Lee los objetos filtrados de 01_qc.R.
# Aplica SCTransform (corrige efectos técnicos y mito%).
# Integra las muestras con 3000 genes comunes.
# Guarda el objeto integrado (integrated_seurat.rds).
# Exporta un PDF rápido con distribuciones (integration_diagnostics.pdf).

# ---- Load libraries ----
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})

# ---- Define input and output ----
input_dir <- "results/01_qc_filtered/"
output_dir <- "results/02_integration/"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---- Define samples ----
samples <- c("FFA1", "FFA2", "FFA3", "FFA4", "CTRL1", "CTRL2", "CTRL3", "CTRL4")

# ---- Load filtered objects ----
seurat_list <- list()
for (s in samples) {
  message("Loading: ", s)
  seurat_list[[s]] <- readRDS(file.path(input_dir, paste0(s, "_filtered.rds")))
}

# ---- Normalize with SCTransform ----
for (s in names(seurat_list)) {
  message("Running SCTransform: ", s)
  seurat_list[[s]] <- SCTransform(seurat_list[[s]],
                                  vars.to.regress = "percent.mt",
                                  verbose = FALSE)
}

# ---- Integration ----
# Select features for integration
features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)

# Prep SCT integration
seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = features)

# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = seurat_list,
                                  normalization.method = "SCT",
                                  anchor.features = features)

# Integrate data
combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# ---- Save integrated object ----
saveRDS(combined, file.path(output_dir, "integrated_seurat.rds"))

# ---- Quick diagnostics ----
pdf(file.path(output_dir, "integration_diagnostics.pdf"), width = 10, height = 6)
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA"), group.by = "group")
dev.off()

message("Integration completed. Output saved in: ", output_dir)
