## make_obj_with_RNA_assay_fullcells.R
## server3
## Enrico Barrozo, Ph.D.
## 4.17.26

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

# -------------------- USER CONFIG --------------------
setwd("/path/to/local_home/visiumhd_mnp/results/audra")

obj_path <- "/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/SCT_chunks_mapped/audra/mnp_obj_w_nn_predicted_celltype.rds"

# output names
out_name <- "mnp_obj_w_nn_predicted_celltype_fullcells_with_RNAassay.rds"
out_wd   <- file.path(getwd(), out_name)
out_raid <- file.path(dirname(obj_path), out_name)

# -------------------- LOAD --------------------
obj <- readRDS(obj_path)

# keep all cells
obj_full <- obj

# -------------------- BUILD RNA ASSAY FROM SCT --------------------
assay_names <- tryCatch(
  names(obj_full@assays),
  error = function(e) NULL
)

if (is.null(assay_names)) {
  stop("Could not retrieve assay names from obj_full. Check object class.")
}

print(assay_names)

if (!"SCT" %in% assay_names) {
  stop("SCT assay not found in obj_full.")
}

DefaultAssay(obj_full) <- "SCT"

# log1p-normalized SCT data
sct_data <- GetAssayData(obj_full, assay = "SCT", slot = "data")

# sparse-safe pseudo-count reconstruction from SCT data
# NOTE: these are pseudo-counts, not original raw RNA counts
sct_counts <- sct_data
sct_counts@x <- expm1(sct_counts@x)
sct_counts@x[sct_counts@x < 0] <- 0
sct_counts@x <- as.numeric(round(sct_counts@x))

# create RNA assay
obj_full[["RNA"]] <- CreateAssayObject(counts = sct_counts)

# populate RNA data slot with same log1p matrix used above
obj_full[["RNA"]]@data <- sct_data

# set default assay for downstream compatibility
DefaultAssay(obj_full) <- "RNA"

# -------------------- CHECKS --------------------
cat("Assays:\n")
print(Assays(obj_full))

cat("\nDefault assay:\n")
print(DefaultAssay(obj_full))

cat("\nRNA data dim:\n")
print(dim(GetAssayData(obj_full, assay = "RNA", slot = "data")))

cat("\nRNA counts dim:\n")
print(dim(GetAssayData(obj_full, assay = "RNA", slot = "counts")))

cat("\nCell count:\n")
print(ncol(obj_full))

# -------------------- SAVE --------------------
saveRDS(obj_full, out_wd)
saveRDS(obj_full, out_raid)

cat("\nSaved to:\n")
cat(out_wd, "\n")
cat(out_raid, "\n")





#### 
obj <- readRDS("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/SCT_chunks_mapped/audra/mnp_obj_w_nn_predicted_celltype_fullcells_with_RNAassay.rds")
DefaultAssay(obj) <- "RNA"
#obj <- NormalizeData(obj, verbose = FALSE)   # optional if she wants to overwrite RNA@data from counts
#obj <- FindVariableFeatures(obj, verbose = FALSE)
#obj <- ScaleData(obj, verbose = FALSE)
#obj <- RunPCA(obj, verbose = FALSE)