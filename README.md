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
scRNAseq_FFA_skin/
├── data/ # Scripts or symbolic links to raw data (not uploaded)
├── scripts/ # Analysis scripts (R/Seurat pipeline)
│ ├── 01_preprocessing.R
│ ├── 02_qc_filtering.R
│ ├── 03_clustering.R
│ ├── 04_diff_expr.R
│ └── 05_pathway_analysis.R
├── results/ # Intermediate outputs (RDS objects, plots, tables)
├── docs/ # Documentation and protocol notes
├── environment.yml # Conda environment for reproducibility
└── README.md # Project description

---

## ⚙️ Installation & requirements
We recommend using **Conda** for reproducibility.

```bash
# Create environment
conda env create -f environment.yml

# Activate environment
conda activate scRNAseq_FFA
Dependencies will include R (Seurat, tidyverse) or Python (Scanpy, anndata), depending on the chosen workflow.

🚀 Usage

Example (R/Seurat-based):

# Preprocessing
Rscript scripts/01_preprocessing.R

# QC and filtering
Rscript scripts/02_qc_filtering.R


Example (Python/Scanpy-based):

python scripts/01_preprocessing.py
python scripts/02_qc_filtering.py

📊 Results

Figures and tables will be stored in results/.

Documentation of each step will be provided in docs/.

📖 Citation

Once the project is finalized, all releases will be archived on Zenodo, with a DOI that can be cited in publications.

📜 License

This project is released under the MIT License. You are free to use, modify, and distribute it, provided that proper credit is given.

👨‍🔬 Author

Juan Ruano, MD, PhD, MSc
Immune-mediated Inflammatory Skin Diseases Laboratory (GC29),
IMIBIC / University of Córdoba,
Reina Sofía University Hospital, Spain.
