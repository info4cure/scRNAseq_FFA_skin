#!/usr/bin/env Rscript

# ============================
# Script 05: Differential Expression
# ============================

# 📌 Qué tendremos al final:
# *DEG_global_FFA_vs_Control.csv: todos los genes diferencialmente expresados entre condiciones.
# *DEG_cluster_X_FFA_vs_Control.csv: DEGs por cluster.
# *DEG_summary.csv: resumen con nº de genes detectados por cluster.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

# ----------------------------
# Load annotated object
# ----------------------------
seurat_obj <- readRDS("results/seurat_annotated.rds")

# ----------------------------
# Check metadata
# ----------------------------
print(head(seurat_obj@meta.data))

# Se espera que exista una columna "condition" con valores "FFA" o "Control"
if (!"condition" %in% colnames(seurat_obj@meta.data)) {
  stop("❌ Metadata must include a 'condition' column (e.g., FFA vs Control).")
}

# ----------------------------
# Differential expression (global)
# ----------------------------
deg_global <- FindMarkers(
  seurat_obj,
  ident.1 = "FFA",
  ident.2 = "Control",
  group.by = "condition",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

# Save results
write_csv(deg_global, "results/tables/DEG_global_FFA_vs_Control.csv")

# ----------------------------
# Differential expression (per cluster)
# ----------------------------
clusters <- levels(seurat_obj$seurat_clusters)

for (cl in clusters) {
  deg_cluster <- FindMarkers(
    subset(seurat_obj, subset = seurat_clusters == cl),
    ident.1 = "FFA",
    ident.2 = "Control",
    group.by = "condition",
    logfc.threshold = 0.25,
    min.pct = 0.1
  )
  
  out_file <- paste0("results/tables/DEG_cluster_", cl, "_FFA_vs_Control.csv")
  write_csv(deg_cluster, out_file)
}

# ----------------------------
# Save DEG summary
# ----------------------------
summary_table <- data.frame(
  Cluster = clusters,
  nDEG = sapply(clusters, function(cl) {
    fn <- paste0("results/tables/DEG_cluster_", cl, "_FFA_vs_Control.csv")
    if (file.exists(fn)) nrow(read_csv(fn)) else 0
  })
)

write_csv(summary_table, "results/tables/DEG_summary.csv")

message("✅ Script 05_differential_expression.R completed successfully.")
