#Seurat Tutorial Go-through 1: https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html

install.packages("Seurat")
library(dplyr)
library(Seurat)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "hg19") #Data imported from the local directory 

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
##?? Interpretation of this violin plot? What is nFeature_RNA? 

