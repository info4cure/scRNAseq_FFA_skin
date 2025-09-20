#!/usr/bin/env Rscript

# Este script explora comunicación intercelular con la librería CellChat. 
# Esto es especialmente interesante para tus datos de piel (FFA vs controles) 
# porque permite ver cómo se modifican las interacciones entre queratinocitos, 
# fibroblastos, células inmunes, sebocitos, etc.
# 
# 📌 Qué hace este script:
# -Carga tu objeto Seurat anotado con celltypes.
# -Usa la base de datos de ligando–receptor de CellChat para humanos.
# -Estima probabilidades de comunicación entre tipos celulares.
# -Genera:
#   *Redes globales (número e intensidad de interacciones).
#   *Redes específicas de vías (ej. IFN, TGF-β, TNF, CXCL).
#   *Bubble plots de ligando–receptor.
#   *Diagramas tipo chord para pathways seleccionados.
# -Exporta tablas con pesos de interacción y pathways detectados.
# 
# 📊 Resultados en results/cellchat/:
# *network_weights.csv, network_pathways.csv
# *cellchat_overview.pdf
# *pathway_IFN-II.tiff, etc.
# *bubble_plots.pdf
# *chord_IFN_TGF.pdf


suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(patchwork)
})

# ----------------------------
# Load Seurat object with annotations
# ----------------------------
seurat_file <- "results/seurat/annotated_clusters.rds"
if (!file.exists(seurat_file)) stop("❌ Run 04_marker_annotation.R first.")
seurat_obj <- readRDS(seurat_file)

# ----------------------------
# Subset metadata
# ----------------------------
meta <- seurat_obj@meta.data
if (!"celltype" %in% colnames(meta)) {
  stop("❌ celltype annotations missing in metadata. Please run annotation step.")
}

# ----------------------------
# Prepare data for CellChat
# ----------------------------
data_input <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")

cellchat <- createCellChat(object = data_input, meta = meta, group.by = "celltype")
cellchat <- addMeta(cellchat, meta = meta)

# Set ligand-receptor database (human skin context)
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB

# ----------------------------
# Preprocess
# ----------------------------
cellchat <- subsetData(cellchat)       # subset to signaling genes
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

# Compute communication probability
cellchat <- computeCommunProb(cellchat, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Aggregate network
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# ----------------------------
# Save communication tables
# ----------------------------
if (!dir.exists("results/cellchat")) dir.create("results/cellchat", recursive = TRUE)

write.csv(cellchat@net$weight, "results/cellchat/network_weights.csv")
write.csv(cellchat@netP$pathways, "results/cellchat/network_pathways.csv")

# ----------------------------
# Visualization
# ----------------------------
pdf("results/cellchat/cellchat_overview.pdf", width = 10, height = 8)
netVisual_circle(cellchat@net$count, vertex.weight = as.numeric(table(cellchat@idents)),
                 weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = as.numeric(table(cellchat@idents)),
                 weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction strength")
dev.off()

# Pathway-specific networks
pathways.show <- c("IFN-II", "TGFB", "TNF", "EGF", "CXCL") # customize as relevant
for (pathway in pathways.show) {
  p <- netVisual_aggregate(cellchat, signaling = pathway, layout = "circle")
  ggsave(paste0("results/cellchat/pathway_", pathway, ".tiff"),
         plot = p, width = 7, height = 6, dpi = 300)
}

# Bubble plot for ligand-receptor pairs
pdf("results/cellchat/bubble_plots.pdf", width = 10, height = 6)
netVisual_bubble(cellchat, sources.use = NULL, targets.use = NULL, remove.isolate = TRUE)
dev.off()

# Chord diagram for selected pathways
pdf("results/cellchat/chord_IFN_TGF.pdf", width = 10, height = 8)
netVisual_chord_gene(cellchat, signaling = c("IFN-II", "TGFB"))
dev.off()

# Save object
saveRDS(cellchat, "results/cellchat/cellchat_obj.rds")

message("✅ CellChat analysis complete. Results in results/cellchat/")
