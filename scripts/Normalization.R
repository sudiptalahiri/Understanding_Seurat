# ============================================================
# Understanding Seurat
# Step: Normalization
# ============================================================


# ------------------------------------------------------------
# Optional: reload QC-filtered Seurat object
# ------------------------------------------------------------

seurat_obj_qc <- readRDS(
  "data/seurat_obj_qc.rds"
)

# Check that the object loaded correctly

seurat_obj_qc

ncol(seurat_obj_qc)
nrow(seurat_obj_qc)
summary(seurat_obj_qc$percent.mt)
DefaultAssay(seurat_obj_qc)
seurat_obj_qc[["RNA"]]

# ------------------------------------------------------------
# Normalize RNA expression data
# ------------------------------------------------------------

seurat_obj_qc <- NormalizeData(
  seurat_obj_qc
)

# ------------------------------------------------------------
# Verify normalization
# ------------------------------------------------------------

seurat_obj_qc[["RNA"]]

dim(
  GetAssayData(
    seurat_obj_qc,
    assay = "RNA",
    layer = "data"
  )
)

# ------------------------------------------------------------
# Save normalized checkpoint
# ------------------------------------------------------------

saveRDS(
  seurat_obj_qc,
  file = "data/seurat_obj_normalized.rds"
)