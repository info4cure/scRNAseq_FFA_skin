#!/usr/bin/env Rscript

# 📌 Qué hace este script:
# Carga el objeto integrado de Seurat.
# -Ejecuta GSVA con tus gene sets personalizados (definidos en resources/custom_genesets.gmt).
# -Guarda una tabla con scores de GSVA por célula.
# -Añade esos scores al metadata del objeto Seurat.
# -Genera ridge plots para los pathways clave (Apoptosis, Stress, Senescence, Exhaustion, Mitochondria).
# 
# -📊 Output esperado en results/gsva/:
# *GSVA_scores.csv → tabla con los scores por célula.
# *RidgePlot_Apoptosis.tiff, RidgePlot_Senescence.tiff, etc. → gráficos de densidad por cluster.
# *seurat_with_gsva.rds → objeto Seurat con los GSVA scores añadidos.

suppressPackageStartupMessages({
  library(Seurat)
  library(GSVA)
  library(ggplot2)
  library(ggridges)
  library(dplyr)
})

# ----------------------------
# Load Seurat object
# ----------------------------
seurat_file <- "results/seurat/integrated_clustering.rds"
if (!file.exists(seurat_file)) stop("❌ Run 03_dimreduction_clustering.R first.")
seurat_obj <- readRDS(seurat_file)

# Create results directory
gsva_dir <- "results/gsva"
if (!dir.exists(gsva_dir)) dir.create(gsva_dir, recursive = TRUE)

# ----------------------------
# Load custom gene sets
# ----------------------------
geneset_file <- "resources/custom_genesets.gmt"
if (!file.exists(geneset_file)) stop("❌ Missing custom gene sets file.")
gene_sets <- GSEABase::getGmt(geneset_file)

# ----------------------------
# Prepare expression matrix (log-normalized)
# ----------------------------
expr_mat <- as.matrix(seurat_obj@assays$RNA@data)

# Run GSVA
gsva_res <- gsva(expr_mat, gene_sets, method = "gsva", kcdf = "Poisson", parallel.sz = 4)

# Save GSVA results
write.csv(gsva_res, file.path(gsva_dir, "GSVA_scores.csv"))

# ----------------------------
# Add GSVA scores back to Seurat metadata (per cell)
# ----------------------------
seurat_obj <- AddMetaData(seurat_obj, t(gsva_res))

# ----------------------------
# Ridge plots for selected pathways
# ----------------------------
selected_pathways <- c("Apoptosis", "Cellular_Stress", "Senescence",
                       "Mitochondrial_Stress", "Tcell_Exhaustion")

for (pathway in selected_pathways) {
  if (!(pathway %in% colnames(seurat_obj@meta.data))) next
  p <- ggplot(seurat_obj@meta.data, aes(x = !!sym(pathway), y = seurat_clusters, fill = condition)) +
    geom_density_ridges(alpha = 0.6, scale = 1.2) +
    theme_minimal(base_size = 12) +
    ggtitle(paste("GSVA scores -", pathway)) +
    xlab("GSVA enrichment score") +
    ylab("Clusters")
  
  ggsave(file.path(gsva_dir, paste0("RidgePlot_", pathway, ".tiff")),
         plot = p, width = 7, height = 5, dpi = 300)
}

saveRDS(seurat_obj, file.path(gsva_dir, "seurat_with_gsva.rds"))

message("✅ GSVA and ridge plots generated in results/gsva/")
