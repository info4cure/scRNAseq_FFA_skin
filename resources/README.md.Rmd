---
title: "cell type–specific markers for human skin"
author: "Juan Ruano"
date: "2025-09-19"
output: html_document
---
# 📂 Resources

This folder contains reference files used throughout the **scRNAseq_FFA_skin** pipeline.  
These files are not raw data but curated resources that support analysis and annotation.

## Contents

- **markers_skin_extended.csv**  
  Curated list of cell type–specific markers for human skin, compiled from the literature and previous studies.  
  Includes keratinocyte subtypes, fibroblast subsets, melanocytes, sebocytes, adipocytes, Schwann cells, dendritic cell subtypes, MAIT cells, ILCs, and other immune and stromal populations relevant for skin biology and disease.  

## Usage

These resources are typically loaded in scripts such as:
- `04_marker_annotation.R` → for automated/semi-automated cell type annotation.  
- `06_pathway_analysis.R` → as background gene sets for enrichment analysis.  

## Notes

- The list is extendable: add new markers or gene sets as knowledge evolves.  
- Each file should include a short header with **source/reference** for transparency.  
- Do not include raw sequencing data in this folder — only supporting metadata and curated references.  
