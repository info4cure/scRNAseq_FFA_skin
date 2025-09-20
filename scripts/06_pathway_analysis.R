#!/usr/bin/env Rscript

# ============================
# Script 06: Pathway Analysis
# ============================


# 📌 Qué genera este script:
#   
# *Archivos CSV con resultados de GO, KEGG, Reactome y GSEA.
# *Figuras (dotplots, ridgeplots) en la carpeta results/figures/.



suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(readr)
  library(dplyr)
})

# ----------------------------
# Load DEGs
# ----------------------------
deg_file <- "results/tables/DEG_global_FFA_vs_Control.csv"
if (!file.exists(deg_file)) stop("❌ DEG file not found. Run 05_differential_expression.R first.")

deg_global <- read_csv(deg_file)

# Ensure gene symbols column
if (!"gene" %in% colnames(deg_global)) {
  stop("❌ DEG file must contain a column named 'gene'.")
}

# ----------------------------
# Define gene lists
# ----------------------------
# Upregulated
genes_up <- deg_global %>%
  filter(avg_log2FC > 0.25 & p_val_adj < 0.05) %>%
  pull(gene)

# Downregulated
genes_down <- deg_global %>%
  filter(avg_log2FC < -0.25 & p_val_adj < 0.05) %>%
  pull(gene)

# Ranked gene list for GSEA
gene_ranks <- deg_global %>%
  arrange(desc(avg_log2FC)) %>%
  select(gene, avg_log2FC) %>%
  deframe()

# ----------------------------
# GO Enrichment (Biological Process)
# ----------------------------
ego_up <- enrichGO(gene = genes_up,
                   OrgDb = org.Hs.eg.db,
                   keyType = "SYMBOL",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)

ego_down <- enrichGO(gene = genes_down,
                     OrgDb = org.Hs.eg.db,
                     keyType = "SYMBOL",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05)

# ----------------------------
# KEGG and Reactome
# ----------------------------
ekegg <- enrichKEGG(gene = genes_up,
                    organism = 'hsa',
                    pvalueCutoff = 0.05)

ereactome <- enrichPathway(gene = genes_up,
                           organism = "human",
                           pvalueCutoff = 0.05,
                           readable = TRUE)

# ----------------------------
# GSEA
# ----------------------------
egsea <- GSEA(geneList = gene_ranks,
              OrgDb = org.Hs.eg.db,
              keyType = "SYMBOL",
              ont = "BP",
              pvalueCutoff = 0.05)

# ----------------------------
# Save results
# ----------------------------
write_csv(as.data.frame(ego_up), "results/pathways/GO_upregulated.csv")
write_csv(as.data.frame(ego_down), "results/pathways/GO_downregulated.csv")
write_csv(as.data.frame(ekegg), "results/pathways/KEGG_upregulated.csv")
write_csv(as.data.frame(ereactome), "results/pathways/Reactome_upregulated.csv")
write_csv(as.data.frame(egsea), "results/pathways/GSEA_global.csv")

# ----------------------------
# Plots
# ----------------------------
pdf("results/figures/pathway_dotplots.pdf", width = 8, height = 6)
dotplot(ego_up, showCategory = 20) + ggtitle("GO BP - Upregulated genes")
dotplot(ego_down, showCategory = 20) + ggtitle("GO BP - Downregulated genes")
dotplot(ekegg, showCategory = 20) + ggtitle("KEGG - Upregulated genes")
dotplot(ereactome, showCategory = 20) + ggtitle("Reactome - Upregulated genes")
dev.off()

pdf("results/figures/GSEA_global.pdf", width = 8, height = 6)
ridgeplot(egsea, showCategory = 20)
dev.off()

message("✅ Script 06_pathway_analysis.R completed successfully.")
