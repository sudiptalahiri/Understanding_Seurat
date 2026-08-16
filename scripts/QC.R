# ============================================================
# 02_QC.R
# Quality control of the PBMC3K scRNA-seq dataset
# ============================================================

# ============================================================
# QC ANALYSIS
# ============================================================
#
# This script examines and filters cells based on:
#
# 1. Number of detected genes (nFeature_RNA)
# 2. Total RNA counts (nCount_RNA)
# 3. Mitochondrial RNA percentage (percent.mt)
#
# The QC checkpoint is saved before filtering so that
# alternative thresholds can be tested without rebuilding
# the Seurat object.
#
# ============================================================

# ============================================================
# OPTIONAL: Reload previously saved QC checkpoint
# ============================================================
#
# Activate this section if you want to restart QC analysis
# from the saved QC checkpoint rather than rebuilding the
# Seurat object from the raw 10X data.
#
# To activate it, remove the "#" from the beginning of
# each line below.
#
# seurat_obj <- readRDS(
#   "data/seurat_obj_qc_metrics.rds"
# )
# ============================================================

# ============================================================
# Load required package
# ============================================================

library(Seurat)
library(ggplot2)

# ============================================================
# Load the raw Seurat object
# ============================================================

# We don't want to rerun Read10X() and CreateSeuratObject() every time we work on QC
seurat_obj <- readRDS(
  file = "data/seurat_obj_raw.rds"
)

# ============================================================
# Verify the Seurat object
# ============================================================

seurat_obj

# ============================================================
# Examine basic QC metrics
# ============================================================
# These metrics will decides the quality of the assay.Check if they are too low or too high.

summary(seurat_obj$nFeature_RNA)
summary(seurat_obj$nCount_RNA)

# ============================================================
# Visualize basic QC metrics
# ============================================================
# Use the Seurat object as the source of the data with these variables (features)
# c() function in R means combine and ncol controls the number of columns in the plot arrangement

VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA"),
  ncol = 2
)

# ============================================================
# Calculate mitochondrial RNA percentage
# ============================================================

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
  seurat_obj,
  pattern = "^MT-"
)

# ============================================================
# Examine mitochondrial QC
# ============================================================

summary(seurat_obj$percent.mt)

VlnPlot(
  seurat_obj,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  ncol = 3
)

# ============================================================
# Examine relationships between mitochondrial
# percentage and other QC metrics
# ============================================================

FeatureScatter(
  seurat_obj,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
)

FeatureScatter(
  seurat_obj,
  feature1 = "nFeature_RNA",
  feature2 = "percent.mt"
)

# ============================================================
# Save QC checkpoint
# ============================================================

saveRDS(
  seurat_obj,
  file = "data/seurat_obj_qc_metrics.rds"
)

# ============================================================
# Step 9: Explore potential QC thresholds
# ============================================================
#
# At this stage we have NOT filtered any cells.
#
# We are going to examine how many cells would be affected
# by different possible QC thresholds.
#
# This allows us to choose thresholds based on the actual
# distribution of this dataset rather than blindly copying
# thresholds from a tutorial.
# ============================================================

# Number of cells below different nFeature_RNA thresholds
# e.g How many cells have fewer than 200 detected genes? 

sum(seurat_obj$nFeature_RNA < 200)
sum(seurat_obj$nFeature_RNA < 300)
sum(seurat_obj$nFeature_RNA < 400)
sum(seurat_obj$nFeature_RNA < 500)

# Number of cells above different nFeature_RNA thresholds
# We are trying to examine whether the upper tail represents a  
# handful of extreme cells or a meaningful population

sum(seurat_obj$nFeature_RNA > 1500)
sum(seurat_obj$nFeature_RNA > 2000)
sum(seurat_obj$nFeature_RNA > 2500)
sum(seurat_obj$nFeature_RNA > 3000)

# Number of cells above different mitochondrial thresholds

sum(seurat_obj$percent.mt > 3)
sum(seurat_obj$percent.mt > 5)
sum(seurat_obj$percent.mt > 10)
sum(seurat_obj$percent.mt > 15)
sum(seurat_obj$percent.mt > 20)

# Total number of cells

n_cells <- ncol(seurat_obj)
n_cells

# Percentage of cells above mitochondrial thresholds

sum(seurat_obj$percent.mt > 3) / n_cells * 100
sum(seurat_obj$percent.mt > 5) / n_cells * 100
sum(seurat_obj$percent.mt > 10) / n_cells * 100
sum(seurat_obj$percent.mt > 15) / n_cells * 100
sum(seurat_obj$percent.mt > 20) / n_cells * 100

qc_threshold_summary <- data.frame(
  
  metric = c(
    "nFeature_RNA",
    "nFeature_RNA",
    "nFeature_RNA",
    "nFeature_RNA",
    "nCount_RNA",
    "nCount_RNA",
    "percent.mt",
    "percent.mt",
    "percent.mt",
    "percent.mt",
    "percent.mt"
  ),
  
  direction = c(
    "below",
    "below",
    "below",
    "below",
    "above",
    "above",
    "above",
    "above",
    "above",
    "above",
    "above"
  ),
  
  threshold = c(
    200,
    300,
    400,
    500,
    10000,
    15000,
    3,
    5,
    10,
    15,
    20
  ),
  
  cells_affected = c(
    sum(seurat_obj$nFeature_RNA < 200),
    sum(seurat_obj$nFeature_RNA < 300),
    sum(seurat_obj$nFeature_RNA < 400),
    sum(seurat_obj$nFeature_RNA < 500),
    sum(seurat_obj$nCount_RNA > 10000),
    sum(seurat_obj$nCount_RNA > 15000),
    sum(seurat_obj$percent.mt > 3),
    sum(seurat_obj$percent.mt > 5),
    sum(seurat_obj$percent.mt > 10),
    sum(seurat_obj$percent.mt > 15),
    sum(seurat_obj$percent.mt > 20)
  )
)

qc_threshold_summary
qc_threshold_summary$percent_affected <-
  +     qc_threshold_summary$cells_affected /
  +     n_cells * 100
qc_threshold_summary

# Create a directory in the path results/qc relative to the project directory
# Make sure the results/qc directory exist and if it does not
# then create one and any missing parent directories. 
# If the directory exist then do not show warnings

dir.create(
  "results/qc",
  recursive = TRUE,
  showWarnings = FALSE
)

# Create a csv file of the qc_threshold_summary dataframe and write it to a file
# Save the file and get rid of the row numbers of the R dataframes

write.csv(
  qc_threshold_summary,
  file = "results/qc/qc_threshold_summary.csv",
  row.names = FALSE
)

saveRDS(
  qc_threshold_summary,
  file = "results/qc/qc_threshold_summary.rds"
)

# ============================================================
# QC filtering parameters
# ============================================================

min_features <- 200
max_percent_mt <- 5

# ============================================================
# Identify cells that will be removed
# ============================================================

qc_remove <- subset(
  seurat_obj@meta.data,
  nFeature_RNA <= min_features |
    percent.mt >= max_percent_mt
)

nrow(qc_remove)
nrow(qc_remove) / ncol(seurat_obj) * 100

# ============================================================
# Apply QC filters
# ============================================================

seurat_obj_qc <- subset(
  seurat_obj,
  subset =
    nFeature_RNA > min_features &
    percent.mt < max_percent_mt
)

# ============================================================
# Before vs. after QC summary
# ============================================================

qc_filter_summary <- data.frame(
  metric = c(
    "cells_before",
    "cells_removed",
    "cells_after",
    "percent_removed",
    "percent_retained"
  ),
  
  value = c(
    ncol(seurat_obj),
    ncol(seurat_obj) - ncol(seurat_obj_qc),
    ncol(seurat_obj_qc),
    (ncol(seurat_obj) - ncol(seurat_obj_qc)) /
      ncol(seurat_obj) * 100,
    ncol(seurat_obj_qc) /
      ncol(seurat_obj) * 100
  )
)

qc_filter_summary

write.csv(
  qc_filter_summary,
  file = "results/qc/qc_filter_summary.csv",
  row.names = FALSE
)

# ============================================================
# Save QC-filtered Seurat object
# ============================================================

saveRDS(
  seurat_obj_qc,
  file = "data/seurat_obj_qc.rds"
)

qc_vln_after <- VlnPlot(
  seurat_obj_qc,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  ncol = 3
)

ggplot2::ggsave(
  filename = "results/qc/qc_violin_after_filtering.png",
  plot = qc_vln_after,
  width = 12,
  height = 5,
  dpi = 300
)