# ============================================================
# Understanding Seurat
# 01 - Loading and Inspecting Single-Cell RNA-seq Data
# ============================================================

# Load Seurat - Loads the Seurat package into the current R session 
library(Seurat) 

# Check the installed Seurat version - this is important for reproducibilty
packageVersion("Seurat")

# For this first learning exercise, we'll use the classic PBMC 3k dataset.
# This is a small single-cell dataset containing roughly 2,700 peripheral blood mononuclear cells.
# The official Seurat introductory workflow uses this dataset.

# ============================================================
# Locate the raw 10X data
# ============================================================

# Here we are creating a variable containing the location of the data.
data_dir <- "data/raw/filtered_gene_bc_matrices/hg19"
data_dir

# Check whether the data directory exists
dir.exists(data_dir)
# List the files in the 10X data directory
list.files(data_dir)
# Read the 10X expression matrix
counts <- Read10X(data.dir = data_dir) 
# Examine the object returned by Read10X()
class(counts)
# Examine its dimensions: meaning number of rows and columns
dim(counts)
# Number of genes/features
nrow(counts)
# Number of cells: meaning number of cell barcodes
ncol(counts)
# Look at the first few genes: gives us the features
head(rownames(counts))
# Look at the first few cell barcodes: think of these as identifiers assigned to individual droplets/cell containing partitions during the 10X experiment.
head(colnames(counts))
# Look at a small section of the expression matrix
counts[1:10, 1:10]

# ============================================================
# Inspect the raw expression matrix
# ============================================================

# Number of genes
nrow(counts)
# Number of cells
ncol(counts)
# Number of detected gene-cell combinations
sum(counts > 0)

# ============================================================
# Create the Seurat object
# ============================================================

# The left-hand counts is the argument name expected by the function
# The right-hand counts is the R object containing the expression matrix
# Therefore one can say take the object called counts and provide it to the function's counts argument
seurat_obj <- CreateSeuratObject(counts = counts)

# ============================================================
# Inspect the Seurat object
# ============================================================

class(seurat_obj)
seurat_obj
head(seurat_obj@meta.data)
Assays(seurat_obj)
DefaultAssay(seurat_obj)

# ============================================================
# Save the initial Seurat object
# ============================================================

saveRDS(
  seurat_obj,
  file = "data/seurat_obj_raw.rds"
)