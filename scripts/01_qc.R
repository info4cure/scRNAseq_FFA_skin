############################################################
# 01_qc.R
# Quality control for scRNA-seq data (10X Genomics outputs)
# Project: scRNAseq_FFA_skin
# Author: Juan Ruano, MD, PhD, MSc
# Date: 2025-09-18
############################################################

# 📌 Qué hace este script
# 
# Carga las 8 muestras desde las carpetas data/FFA1/filtered_feature_bc_matrix/, etc.
# Crea objetos Seurat con metadatos (sample_id, group).
# Calcula % de mitocondriales (percent.mt).
# Aplica filtros típicos:
#   nFeature_RNA > 500 & < 6000
#   nCount_RNA > 1000
#   percent.mt < 15%
# Guarda .rds individuales y un objeto combinado (combined_filtered.rds).
# Exporta un PDF con gráficas QC (violines de genes, UMIs y %mt).

# Lo importante es que la estructura del script ya está preparada: solo tienes que decirle dónde están los archivos de CellRanger (matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz).
# 
# Por defecto suelen estar en el el subdirectorio filtered_feature_bc_matrix:
#   
#   FFA1/outs/filtered_feature_bc_matrix/
#   FFA1/outs/filtered_feature_bc_matrix/
#   ...
# 
# pero a veces el proveedor te los entrega ya renombrados (ej. data/FFA1/filtered_feature_bc_matrix/).

# ---- Load libraries ----
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})

# ---- Define input and output paths ----
# Adjust these paths depending on your structure
data_dir <- "data/"      # Directory containing CellRanger outputs
results_dir <- "results/01_qc_filtered/"

if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# ---- Define samples ----
# Example: each sample has a folder with matrix.mtx, barcodes.tsv, features.tsv
samples <- c("FFA1", "FFA2", "FFA3", "FFA4", "CTRL1", "CTRL2", "CTRL3", "CTRL4")
groups  <- c(rep("FFA", 4), rep("Control", 4))

# ---- Load 10X data ----
seurat_list <- list()

for (i in seq_along(samples)) {
  sample_id <- samples[i]
  group_id  <- groups[i]
  
  message("Loading sample: ", sample_id)
  
  # Path to CellRanger filtered_feature_bc_matrix
  sample_path <- file.path(data_dir, sample_id, "filtered_feature_bc_matrix")
  
  seurat_obj <- Read10X(sample_path) %>%
    CreateSeuratObject(project = "scRNAseq_FFA_skin",
                       min.cells = 3,
                       min.features = 200)
  
  # Add metadata
  seurat_obj$sample_id <- sample_id
  seurat_obj$group <- group_id
  
  # ---- QC metrics ----
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # Save raw object (optional)
  saveRDS(seurat_obj, file = file.path(results_dir, paste0(sample_id, "_raw.rds")))
  
  seurat_list[[sample_id]] <- seurat_obj
}

# ---- Apply QC filters ----
# Typical thresholds – adjust for your dataset
for (s in names(seurat_list)) {
  obj <- seurat_list[[s]]
  
  seurat_list[[s]] <- subset(
    obj,
    subset = nFeature_RNA > 500 & nFeature_RNA < 6000 &
      nCount_RNA > 1000 & percent.mt < 15
  )
  
  message(s, ": ", ncol(seurat_list[[s]]), " cells after filtering")
  
  # Save filtered object
  saveRDS(seurat_list[[s]], file = file.path(results_dir, paste0(s, "_filtered.rds")))
}

# ---- Summary QC plots ----
pdf(file.path(results_dir, "QC_metrics.pdf"), width = 10, height = 6)
for (s in names(seurat_list)) {
  VlnPlot(seurat_list[[s]],
          features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
          ncol = 3) + ggtitle(s)
}
dev.off()

# ---- Save combined object (optional) ----
combined <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = samples)
saveRDS(combined, file = file.path(results_dir, "combined_filtered.rds"))

message("QC step completed. Results saved in: ", results_dir)
