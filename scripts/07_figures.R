#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(EnhancedVolcano)
  library(pheatmap)
})

# Genera las figuras finales para el manuscrito o un preprint.
# Aquí la idea es consolidar:
# -UMAPs con clusters y condiciones (FFA vs control).
# -Heatmaps de genes marcadores.
# -DotPlots de pathways.
# -Volcano plots (si usamos pseudobulk o DEG).
# 
# 📌 Output esperado:
# *UMAP_clusters_condition.tiff → UMAP doble (clusters y condición).
# *Heatmap_top10_markers.pdf → Heatmap con los top10 genes.
# *DotPlot_custom_pathways.pdf → Dotplot/heatmap con pathways de apoptosis, stress, senescencia, etc.
# *Volcano_FFA_vs_Control.pdf → Volcano plot DEG.
# Guardar todo en TIFF/PNG a 300 dpi.

# ----------------------------
# Load integrated object
# ----------------------------
seurat_file <- "results/seurat/integrated_clustering.rds"
if (!file.exists(seurat_file)) stop("❌ Run 03_dimreduction_clustering.R first.")
seurat_obj <- readRDS(seurat_file)

# Create results directory
fig_dir <- "results/figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# ----------------------------
# 1. UMAPs
# ----------------------------
p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Clusters")
p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "condition") +
  ggtitle("FFA vs Control")

ggsave(file.path(fig_dir, "UMAP_clusters_condition.tiff"),
       p1 + p2, width = 10, height = 5, dpi = 300)

# ----------------------------
# 2. Heatmap of top markers
# ----------------------------
# Load markers (from script 04_marker_annotation.R)
markers_file <- "results/markers/top_markers.csv"
if (file.exists(markers_file)) {
  markers <- read.csv(markers_file)
  top10 <- markers %>%
    group_by(cluster) %>%
    top_n(10, wt = avg_log2FC)
  heatmap_genes <- unique(top10$gene)
  
  pdf(file.path(fig_dir, "Heatmap_top10_markers.pdf"), width = 8, height = 10)
  DoHeatmap(seurat_obj, features = heatmap_genes, size = 3) +
    ggtitle("Top 10 marker genes per cluster")
  dev.off()
}

# ----------------------------
# 3. DotPlot of pathways
# ----------------------------
if (file.exists("results/pathways/GSVA_custom_scores.csv")) {
  gsva_scores <- read.csv("results/pathways/GSVA_custom_scores.csv", row.names = 1)
  gsva_subset <- gsva_scores[c("Apoptosis","Cellular_Stress","Senescence","Tcell_Exhaustion","Mitochondrial_Stress"),]
  
  pdf(file.path(fig_dir, "DotPlot_custom_pathways.pdf"), width = 8, height = 6)
  pheatmap(gsva_subset, cluster_rows = TRUE, cluster_cols = TRUE,
           main = "Custom pathway activity")
  dev.off()
}

# ----------------------------
# 4. Volcano plots
# ----------------------------
deg_file <- "results/deg/FFA_vs_Control_DEG.csv"
if (file.exists(deg_file)) {
  degs <- read.csv(deg_file)
  pdf(file.path(fig_dir, "Volcano_FFA_vs_Control.pdf"), width = 7, height = 6)
  EnhancedVolcano(degs,
                  lab = degs$gene,
                  x = "log2FoldChange",
                  y = "p_val_adj",
                  pCutoff = 0.05,
                  FCcutoff = 0.5,
                  title = "FFA vs Control")
  dev.off()
}

message("✅ Final figures generated in results/figures/")
