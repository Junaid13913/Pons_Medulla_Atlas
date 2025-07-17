# Comprehensive Single-Cell Transcriptomic Atlas of Mouse Pons & Medulla

This repository provides metadata, scripts, and documentation for the large-scale transcriptomic atlas of the pons and medulla, compiled from 8 independent single-cell/single-nucleus RNA sequencing (sc/snRNA-seq) datasets. 

The integration was performed using a standardized bioinformatic workflow, clustering, and marker-based annotation. The dataset comprises 317,985 quality-passed cells across 45 cell types.

This integrated large-scale dataset serves as a valuable single-cell neuroscience dataset to explore region-specific molecular insights and cellular diversity in the brainstem.

---

## Dataset Overview
- Cells: 317,985 quality passed cells
- Cell Types: 45 distinct cell types
- Methodology: scRNA-seq/snRNA-seq integration
- Data Format:`.rds` files (Seurat objects)
- Integration Tool: [Seurat v4.3.0](https://satijalab.org/seurat/)

## Download the Dataset
The `.rds` files containing normalized expression data can be accessed via Figshare:
1. `Final_Integrated_Pons_Medulla.rds` | Full dataset with all identified cell types |Download](https://doi.org/10.6084/m9.figshare.28342025.v1) |
2. `Pons_Medulla_Neurons_level_3.rds` | Neuronal subtypes dataset | [Download](https://doi.org/10.6084/m9.figshare.28342025.v1)

## Reproducing the Analysis

### Prerequisites
- R version: ≥ 4.2.1
- Key R packages:
  - `devtools`
  - `reticulate`
  - `Seurat`
  - `Matrix`
  - `dplyr`
  - `ggplot2`
  - `tidyr`
  - `tidyverse`
  - `readxl`
  - `ComplexHeatmap`
  - `RColorBrewer`
  - `clusterProfiler`

Install dependencies in R:
```R
install.packages(c("devtools", "reticulate", "Seurat", "Matrix", "dplyr", 
                   "ggplot2", "tidyr", "tidyverse", "readxl", "ComplexHeatmap", 
                   "RColorBrewer"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("clusterProfiler")

## How to Use the Data
1. Download the `.rds` files from Figshare.
2. Load the dataset into R:
   ```r
   library(Seurat)
   pons_medulla <- readRDS("Final_Integrated_Pons_Medulla.rds")

