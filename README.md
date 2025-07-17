# Comprehensive Single-Cell Transcriptomic Atlas of Mouse Pons & Medulla

This repository provides metadata, scripts, and documentation for the large-scale transcriptomic atlas of the mouse pons and medulla subregions, compiled from 8 independent single-cell/single-nucleus RNA sequencing (sc/snRNA-seq) datasets. 

The integration was performed using a standardized bioinformatic workflow, clustering, and marker-based annotation. The dataset comprises 317,985 quality-passed cells across 45 cell types.

This integrated large-scale dataset serves as a valuable single-cell neuroscience dataset to explore region-specific molecular insights and cellular diversity in the brainstem.

---

## Dataset Overview
- Cells: 317,985 quality passed cells
- Species: Mouse (Mus musculus)
- Cell Types: 45 distinct cell types
- Methodology: scRNA-seq/snRNA-seq integration
- Data Format:`.rds` files (Seurat objects)
- Integration Tool: [Seurat v4.4.2](https://satijalab.org/seurat/)

## Download the Dataset
The `.rds` files containing normalized expression data can be accessed via Figshare:
1. `Final_Integrated_Pons_Medulla.rds` | Full dataset with all identified cell types |Download](https://doi.org/10.6084/m9.figshare.28342025.v1) |
2. `Pons_Medulla_Neurons_level_3.rds` | Neuronal subtypes dataset | [Download](https://doi.org/10.6084/m9.figshare.28342025.v1)
3. Differentially Expressed Genes (DEGs)

File: DEG_Level1_Pons_Medulla.csv

Description: Cluster-level transcriptional differences across all cells (Level 1 annotations).

File: DEG_Level1_Neurons.csv

Description: Cluster-level transcriptional differences within neuronal populations only.

## Reproducing the Analysis
## Code and Technical Details
All technical details for:

Preprocessing

Integration

Clustering

Cell type annotation

Subclustering and subtyping

are provided in the code section of this repository. 

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

# Check metadata
head(pons_medulla@meta.data)

# View the annotated UMAP
DimPlot(pons_medulla, reduction = "umap", group.by = "cell_type")

# View all annotated cell types
unique(pons_medulla$cell_type)

# Cell counts per dataset
table(pons_medulla$Study)

# Cell counts per cell type
table(pons_medulla$cell_type)

# Final UMAP with all major cell types and subtypes
custom_colors <- c(
  "GABA1" = "#1b9e77", "GABA2" = "#d95f02", "GABA3" = "#7570b3", 
  "GABA4" = "#e7298a", "GABA5" = "#66a61e", "GABA6" = "#e6ab02",
  "Glut1" = "#a6761d", "Glut2" = "#666666", "Glut3" = "#1f78b4", 
  "Glut4" = "#b2df8a", "Glut5" = "#33a02c", "Glut6" = "#fb9a99", 
  "Glut7" = "#e31a1c",
  "Autonomic projection neurons" = "#d53e4f", "Cholinergic neurons" = "#ff7f00", 
  "Haemoglobin" = "#6a3d9a", "Ion_trans gilal cells" = "#8dd3c7", "Myelinating neurons" = "#b3de69", 
  "Projection neurons" = "#bc80bd", "Synaptic neurons" = "#ccebc5", 
  "Sensory neurons" = "#ffed6f",
  "Neu_1" = "#8c564b", "Neu_2" = "#9467bd", "Neu_3" = "#17becf", 
  "Neu_4" = "#d62728", "Neu_5" = "#2ca02c", "Neu_6" = "#ff7f0e",
  "Olig1" = "#FF6F61", "Olig2" = "#6A0572", "Olig3" = "#2A9D8F", 
  "Olig4" = "#264653", "Olig5" = "#F4A261", "Olig6" = "#E63946",
  "Olig7" = "#2166AC",
  "Neurons" = "#E63946",                
  "Astro1" = "#F4A261",                
  "Astro2" = "#2A9D8F",                 
  "Oligodendrocytes" = "#CC99FF",        
  "Polydendrocytes" = "#E76F51",         
  "Mature oligodendrocytes" = "#b15928",
  "Microglia" = "#66CC99",               
  "Endothelial" = "#457B9D",             
  "Noradrenergic_neurons" = "#A8DADC",   
  "Vascular Leptomeningeal cells" = "#E9C46A", 
  "Newly formed oligodendrocytes" = "#6A4C93", 
  "Purkinje cells" = "#FF6F61",          
  "Metabolic Glio-Pericyte" = "#FF9966"  
)

# custom color palette to UMAP plot
DimPlot(pons_medulla, reduction = "umap", group.by = "cell_type", raster = FALSE, pt.size = 0.2, label = FALSE) + 
  scale_color_manual(values = custom_colors) 

# Cell type proportion group by each study
# data frame for cell type counts
cell_type_counts <- as.data.frame(table(pons_medulla$cell_type))
colnames(cell_type_counts) <- c("Cell_Type", "Count")

# Bar plot
ggplot(cell_type_counts, aes(x = reorder(Cell_Type, -Count), y = Count, fill = Cell_Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = custom_colors) +  # Use your defined colors
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Adjust text
        axis.title = element_text(size = 12, face = "bold"),
        legend.position = "none") +
  labs(title = "Cell Type Distribution", x = "Cell Type", y = "Count")

