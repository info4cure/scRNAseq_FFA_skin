############################################################
# 03_dimreduction_clustering.R
# Dimensionality reduction and clustering
# Project: scRNAseq_FFA_skin
# Author: Juan Ruano, MD, PhD, MSc
# Date: 2025-09-18
############################################################

# Qué hace este script:
# 
# Cargar el objeto integrado de 02_normalization_integration.R.
# Realizar reducción de dimensionalidad (PCA).
# Calcular embeddings no lineales (UMAP/t-SNE).
# Identificar clusters con FindNeighbors y FindClusters.
# Guardar el objeto con la metadata de clustering.
# Exportar gráficos (UMAP coloreado por grupo y por cluster).

# ---- Load libraries ----
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})

# ---- Define input and output ----
input_file <- "results/02_integration/integrated_seurat.rds"
output_dir <- "results/03_clustering/"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---- Load integrated object ----
combined <- readRDS(input_file)

# ---- PCA ----
combined <- RunPCA(combined, npcs = 50, verbose = FALSE)

# Elbow plot para elegir PCs
pdf(file.path(output_dir, "PCA_elbow_plot.pdf"))
ElbowPlot(combined, ndims = 50)
dev.off()

# ---- UMAP ----
combined <- RunUMAP(combined, dims = 1:30)  # ajusta el número de PCs según el elbow plot

# ---- t-SNE (opcional) ----
# combined <- RunTSNE(combined, dims = 1:30)

# ---- Clustering ----
combined <- FindNeighbors(combined, dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5) # prueba 0.3–1.0

# ---- Save object ----
saveRDS(combined, file.path(output_dir, "clustered_seurat.rds"))

# ---- Plots ----
pdf(file.path(output_dir, "UMAP_plots.pdf"), width = 10, height = 5)
p1 <- DimPlot(combined, reduction = "umap", group.by = "group") + ggtitle("UMAP by group")
p2 <- DimPlot(combined, reduction = "umap", label = TRUE) + ggtitle("UMAP by clusters")
print(p1 + p2)
dev.off()

message("Dimensionality reduction and clustering completed. Results saved in: ", output_dir)
