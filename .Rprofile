# .Rprofile for scRNAseq_FFA_skin

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Load commonly used packages (optional, comment out if too slow)
pkgs <- c(
  "Seurat",
  "tidyverse",
  "data.table",
  "patchwork",
  "cowplot",
  "RColorBrewer"
)

for (p in pkgs) {
  if (requireNamespace(p, quietly = TRUE)) {
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
}

# Message on startup
packageStartupMessage("✅ scRNAseq_FFA_skin project environment loaded")

