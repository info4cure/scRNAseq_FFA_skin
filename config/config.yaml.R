# Configuration file for scRNAseq_FFA_skin project
# Author: Juan Ruano
# Repository: https://github.com/info4cure/scRNAseq_FFA_skin

project:
  name: "scRNAseq_FFA_skin"
description: "Single-cell RNA sequencing analysis of frontal fibrosing alopecia (FFA) and healthy control scalp skin"
organism: "Homo sapiens"
genome: "GRCh38"

data:
  input_format: "10x"
samples:
  - id: FFA_1
group: FFA
path: "data/FFA_1/"
- id: FFA_2
group: FFA
path: "data/FFA_2/"
- id: FFA_3
group: FFA
path: "data/FFA_3/"
- id: FFA_4
group: FFA
path: "data/FFA_4/"
- id: HC_1
group: Healthy
path: "data/HC_1/"
- id: HC_2
group: Healthy
path: "data/HC_2/"
- id: HC_3
group: Healthy
path: "data/HC_3/"
- id: HC_4
group: Healthy
path: "data/HC_4/"

qc:
  min_genes_per_cell: 200
max_genes_per_cell: 6000
max_mito_percent: 10
min_cells_per_gene: 3

normalization:
  method: "lognormalize"
scale_factor: 10000
variable_features: 2000

dimensionality_reduction:
  pca_dims: 30
umap_neighbors: 15
umap_min_dist: 0.3

clustering:
  method: "leiden"
resolution: 0.6

differential_expression:
  test: "wilcox"
logfc_threshold: 0.25
adj_pval_threshold: 0.05

visualization:
  color_by: ["group", "cluster", "sample"]
highlight_genes: ["KRT14", "CD3E", "COL1A1", "IFNG", "CXCL9", "CXCL10"]
