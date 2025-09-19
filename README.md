# scRNAseq_FFA_skin
Reproducible pipeline for single-cell RNA-seq analysis of frontal fibrosing alopecia (FFA) and healthy controls.

## 📌 Project description
This repository hosts a reproducible pipeline for the analysis of single-cell RNA sequencing (scRNA-seq) data from **frontal fibrosing alopecia (FFA) patients** and **healthy controls**.  
The project aims to provide a transparent and shareable workflow that can be reused, adapted, and cited by other researchers.

---

## 🎯 Objectives
- Preprocess and quality-control scRNA-seq data.  
- Perform clustering, cell type annotation, and trajectory inference.  
- Identify differentially expressed genes (DEGs) between FFA and controls.  
- Explore pathway enrichment and regulatory networks.  
- Provide a reproducible and publishable workflow following FAIR principles.  

---

## 🗂 Repository structure
```plaintext
scRNAseq_FFA_skin/
├── data/               # Symbolic links to raw data or metadata (not uploaded)
├── scripts/            # Analysis scripts (R/Seurat pipeline)
│   ├── 01_preprocessing.R
│   ├── 02_qc_filtering.R
│   ├── 03_clustering.R
│   ├── 04_diff_expr.R
│   ├── 05_pathway_analysis.R
│   ├── 06_celltype_annotation.R
│   ├── 07_ridgeplots_gsva.R
│   ├── 08_pseudotime_analysis.R
│   └── 09_cellchat_analysis.R
├── results/            # Intermediate outputs (RDS objects, plots, tables)
├── docs/               # Documentation and protocol notes
├── config/             # Configuration files (YAML, metadata)
├── environment.yml     # Conda environment for reproducibility
└── README.md           # Project description
```
# Analysis pipeline for scRNA-seq in FFA skin biopsies

This directory contains R scripts that implement a reproducible workflow for single-cell RNA-seq data analysis using Seurat and complementary packages.

### Workflow overview

1. **01_qc.R**
   - Load CellRanger outputs (`matrix.mtx`, `barcodes.tsv`, `features.tsv`)
   - Create Seurat object(s)
   - Perform initial quality control (mitochondrial %, number of genes, number of UMIs)
   - Save filtered objects for downstream steps

2. **02_normalization_integration.R**
   - Normalize data with `SCTransform`
   - Integrate multiple samples (CCA or RPCA-based integration)
   - Save integrated Seurat object

3. **03_dimreduction_clustering.R**
   - Run PCA and UMAP/t-SNE
   - Identify clusters (`FindNeighbors`, `FindClusters`)
   - Save Seurat object with clustering metadata

4. **04_marker_annotation.R**
   - Find cluster markers (`FindAllMarkers`)
   - Annotate cell types using canonical markers and/or automated tools (e.g., `SingleR`)
   - Save annotation results

5. **05_differential_expression.R**
   - Compare expression profiles between FFA vs Controls
   - Export DEG tables for downstream validation

6. **06_pathway_analysis.R**
   - Run enrichment analysis (e.g., GO, Reactome, KEGG)
   - Apply GSVA or GSEA to cluster-level or pseudobulk data
   - Save results

7. **07_figures.R**
   - Generate final publication-ready figures:
     - UMAP plots
     - Heatmaps
     - Dot plots
     - Volcano plots
   - Export figures in 300 dpi TIFF/PNG

### Advanced modules

8. **08_ridgeplots_gsva.R**
   - Performs GSVA using curated gene sets (immune pathways, keratinocyte biology, fibrosis, etc.).
   - Generates ridge plots to visualize pathway activity across cell subtypes or clusters.
   - Outputs tables of GSVA enrichment scores and publication-ready plots.

9. **09_pseudotime_analysis.R**
   - Constructs pseudotime trajectories using Monocle3 or Slingshot.
   - Identifies genes/pathways dynamically regulated along differentiation axes.
   - Exports pseudotime plots and gene expression trends.

10. **10_cellchat_analysis.R**
   - Infers intercellular communication networks using `CellChat`.
   - Identifies enriched ligand–receptor interactions between FFA vs controls.
   - Generates chord diagrams, bubble plots, and network graphs.


### Notes
- All scripts should be run sequentially.
- Parameters are stored in `config/config.yaml`.
- Intermediate objects are saved in `results/` for reproducibility.
- Use `set.seed()` for consistent results across runs.

