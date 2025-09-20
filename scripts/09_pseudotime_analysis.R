#!/usr/bin/env Rscript


# Ahora vamos con el script 09_pseudotime_analysis.R, donde trabajaremos 
# con Monocle3 (o Slingshot si prefieres) para reconstruir trayectorias celulares
# y ver cómo cambian genes/pathways dinámicamente en las muestras de FFA vs controles.
# 
# 📌 Qué hace este script:
# -Convierte tu objeto Seurat integrado en un objeto Monocle3 cds.
# -Reconstruye trayectorias celulares en espacio UMAP.
# -Ordena células en pseudotime (puedes definir raíz = keratinocitos basales en controles).
# -Exporta valores de pseudotime y genera plots:
#   .Trayectoria coloreada por pseudotime.
#   .Trayectoria coloreada por condición (FFA vs Control).
#   -Corre graph_test para identificar genes regulados dinámicamente en pseudotime.
# 
# 📊 Resultados guardados en results/pseudotime/:
# *pseudotime_values.csv → pseudotime por célula.
# *trajectory_pseudotime.tiff → pseudotime visual.
# *trajectory_condition.tiff → comparando FFA vs controles.
# *DEG_pseudotime.csv → genes con expresión dinámica.

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(dplyr)
  library(ggplot2)
})


# ----------------------------
# Load Seurat object
# ----------------------------
seurat_file <- "results/seurat/integrated_clustering.rds"
if (!file.exists(seurat_file)) stop("❌ Run 03_dimreduction_clustering.R first.")
seurat_obj <- readRDS(seurat_file)

# ----------------------------
# Convert Seurat -> Monocle3
# ----------------------------
cds <- as.cell_data_set(seurat_obj)

# Transfer cluster and condition metadata
cds@clusters$UMAP <- seurat_obj@meta.data$seurat_clusters
colData(cds)$condition <- seurat_obj@meta.data$condition

# ----------------------------
# Preprocess and reduce dimension
# ----------------------------
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, reduction_method = "UMAP")

# ----------------------------
# Learn trajectory graph
# ----------------------------
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds)

# ----------------------------
# Order cells along pseudotime
# ----------------------------
# Root cells: can be chosen manually (e.g., basal keratinocytes in controls)
cds <- order_cells(cds)

# ----------------------------
# Save pseudotime values
# ----------------------------
pseudotime_vals <- pseudotime(cds)
out_file <- "results/pseudotime/pseudotime_values.csv"
if (!dir.exists("results/pseudotime")) dir.create("results/pseudotime", recursive = TRUE)
write.csv(pseudotime_vals, out_file)

# ----------------------------
# Plot pseudotime trajectories
# ----------------------------
p1 <- plot_cells(cds,
                 color_cells_by = "pseudotime",
                 label_groups_by_cluster = FALSE,
                 label_leaves = TRUE,
                 label_branch_points = TRUE) +
  ggtitle("Pseudotime trajectory (Monocle3)")

p2 <- plot_cells(cds,
                 color_cells_by = "condition",
                 label_groups_by_cluster = TRUE) +
  ggtitle("Trajectory colored by condition (FFA vs Control)")

ggsave("results/pseudotime/trajectory_pseudotime.tiff", plot = p1, width = 7, height = 6, dpi = 300)
ggsave("results/pseudotime/trajectory_condition.tiff", plot = p2, width = 7, height = 6, dpi = 300)

# ----------------------------
# Differential gene expression along pseudotime
# ----------------------------
deg_pseudotime <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
deg_pseudotime <- deg_pseudotime %>%
  arrange(q_value) %>%
  filter(q_value < 0.05)

write.csv(deg_pseudotime, "results/pseudotime/DEG_pseudotime.csv")

message("✅ Pseudotime analysis complete. Results in results/pseudotime/")
