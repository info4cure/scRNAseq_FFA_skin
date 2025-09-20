#!/usr/bin/env Rscript

# ============================
# Script 04: Marker detection & cell type annotation
# ============================

# Aquí vamos a usar dos métodos complementarios:
# - Manual annotation (usando tu librería markers_skin_extended.csv).
# - Automated annotation con SingleR (usando celulas de piel/public datasets o el HumanPrimaryCellAtlas como referencia).
# Qué tendrás al final:
# -cluster_markers.csv: tabla con marcadores por cluster.
# -Dos anotaciones:
# *celltype_manual (tu librería de marcadores de piel).
# *celltype_singleR (HumanPrimaryCellAtlas).
# -Figura umap_annotation.pdf con ambos UMAP comparativos.
# -seurat_annotated.rds con toda la metadata lista.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(SingleR)
  library(celldex)   # Para referencias estandarizadas
  library(SingleCellExperiment)
})

# ----------------------------
# Load clustered object
# ----------------------------
seurat_obj <- readRDS("results/seurat_clusters.rds")

# ----------------------------
# Find cluster markers
# ----------------------------
markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Save markers table
write_csv(markers, "results/tables/cluster_markers.csv")

# ----------------------------
# Manual annotation (skin-specific markers)
# ----------------------------
skin_markers <- read_csv("resources/markers_skin_extended.csv")
print(table(skin_markers$CellType))

seurat_obj$celltype_manual <- "Unassigned"

for (cell_type in unique(skin_markers$CellType)) {
  genes <- skin_markers %>%
    filter(CellType == cell_type) %>%
    pull(Gene)
  
  expressed_clusters <- markers %>%
    filter(gene %in% genes) %>%
    pull(cluster) %>%
    unique()
  
  if (length(expressed_clusters) > 0) {
    seurat_obj$celltype_manual[seurat_obj$seurat_clusters %in% expressed_clusters] <- cell_type
  }
}

# ----------------------------
# Automated annotation with SingleR
# ----------------------------
# Convert Seurat -> SingleCellExperiment
sce <- as.SingleCellExperiment(seurat_obj)

# Use Human Primary Cell Atlas reference
ref <- celldex::HumanPrimaryCellAtlasData()

pred <- SingleR(test = sce, ref = ref, labels = ref$label.main)

# Add SingleR labels
seurat_obj$celltype_singleR <- pred$labels

# ----------------------------
# Save annotated plots
# ----------------------------
pdf("results/plots/umap_annotation.pdf", width = 10, height = 6)

# Manual annotation
print(DimPlot(seurat_obj, reduction = "umap", group.by = "celltype_manual", label = TRUE) +
        ggtitle("UMAP: Manual annotation (skin markers)"))

# Automated annotation
print(DimPlot(seurat_obj, reduction = "umap", group.by = "celltype_singleR", label = TRUE) +
        ggtitle("UMAP: Automated annotation (SingleR/HumanPrimaryCellAtlas)"))

dev.off()

# ----------------------------
# Save annotated object
# ----------------------------
saveRDS(seurat_obj, "results/seurat_annotated.rds")

message("✅ Script 04_marker_annotation.R completed successfully.")
