# RCTD_spatialmnp_visiumv1.R




# 8.02.25
  # server2

# Enrico Barrozo, Ph.D. 
# Assistant Professor
# Department of Obstetrics & Gynecology
# Baylor College of Medicine and Texas Children’s Hospital

# Load session information
sessionInfo()

## Set working directory and seed
setwd("/path/to/local_home/visiumhd_mnp/results/seurat_all_v1")
set.seed(1)
options(future.globals.maxSize = 300e9) ## 300 gb ram 
## free -h ## in terminal

# Integration with scRNA-seq data (deconvolution); Robust Cell Type Decomposition RCTD
library(spacexr)
# load the reference 
reference <- readRDS("/path/to/local_storage/raid2/project/spatial_MNP/placental_ref_Jul_2025.rds") ## this is  11.5 gb
DefaultAssay(reference) <- "RNA"
Idents(reference) <- "cell_type_expanded"
levels(reference)  # Should list only placenta + decidua cell types
t1 <- table(Idents(reference))  # Confirm cell counts
write.table(t1, "reference_cell_type_expanded.counts.txt", sep="\t")

counts <- GetAssayData(reference, assay = "RNA", layer = "counts")
counts <- round(counts)  # Ensure whole numbers
# counts <- as.matrix(counts)  # Convert to standard matrix format
# storage.mode(counts) <- "integer"  # Set as integer type

cluster <- as.factor(reference$cell_type_expanded)
names(cluster) <- colnames(reference)
nUMI <- reference$nCount_RNA
names(nUMI) <- colnames(reference)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
rm(reference)
gc()
reference.ref <- Reference(counts, cluster, nUMI)





## Set working directory and seed
setwd("/path/to/local_home/visiumhd_mnp/results/seurat_all_v1")
library(Seurat)
library(spacexr) # for SpatialRNA
library(purrr)   # for list operations
# --- LIBRARIES ---
library(tidyverse)
library(data.table)
seurat_all <- readRDS("/path/to/local_storage/raid2/project/visiumhd_mnp/all_spatial_seurat_objects.rds") ## this is  11.5 gb
# 1. Define your sample groups
#  "HP2","HP2", 
table(names(seurat_all))
# You'll see HP2 appears twice

# Let's manually fix:
# (assuming order is as shown in your post)
names(seurat_all)[which(names(seurat_all) == "HP2")[1]] <- "HP2_xen"
names(seurat_all)[which(names(seurat_all) == "HP2")[1]] <- "HP2_hd"
for (sample in names(seurat_all)) {
  obj <- seurat_all[[sample]]
  obj$orig.ident <- sample
  seurat_all[[sample]] <- obj
}

# Confirm
table(names(seurat_all))


query_list  <- list()
coords_list <- list()
umi_thresholds <- c()  # Save for record-keeping


# Helper: optionally visualize UMI distributions per sample
plot_umi_dist <- function(umis, sample_name) {
    png(sprintf("umi_dist_%s.png", sample_name), width=900, height=600)
    hist(umis, breaks=100, main=paste("UMI per pixel:", sample_name), xlab="UMI", log="x")
    abline(v=median(umis), col="red", lty=2)
    dev.off()
}
for (sample in names(seurat_all)) {
    obj <- seurat_all[[sample]]
    assay_name <- DefaultAssay(obj)
    counts <- GetAssayData(obj, assay = assay_name, layer = "counts")
    img_names <- names(obj@images)
    coords <- NULL
    if (length(img_names) > 0 && "centroids" %in% names(obj@images[[img_names[1]]]@boundaries)) {
        centroids_obj <- obj@images[[img_names[1]]]@boundaries[["centroids"]]
        barcodes <- centroids_obj@cells
        coords_matrix <- centroids_obj@coords
        rownames(coords_matrix) <- barcodes
        coords_df <- as.data.frame(coords_matrix)
        colnames(coords_df) <- c("x", "y")
        coords <- coords_df
    } else {
        coords <- GetTissueCoordinates(obj, which = "centroids")
        if ("cell" %in% colnames(coords)) {
            rownames(coords) <- coords$cell
            coords$cell <- NULL
        }
    }
    # Defensive: intersect and match order
    common_barcodes <- intersect(colnames(counts), rownames(coords))
    if (length(common_barcodes) == 0) {
        warning("No matching barcodes for ", sample)
        next
    }
    counts <- counts[, common_barcodes, drop = FALSE]
    coords <- coords[common_barcodes, c("x", "y"), drop = FALSE]
    stopifnot(all(colnames(counts) == rownames(coords)))

    # --- DYNAMIC UMI THRESHOLD ---
    umis <- Matrix::colSums(counts)
    plot_umi_dist(umis, sample)  # Optional, remove/comment out if not desired

    # Pick a dynamic threshold (example: 10th percentile, but you can use elbow, fixed value, etc.)
    umi_threshold <- quantile(umis, 0.10)  # Change 0.10 as needed for your cohort/platform
    umi_threshold <- max(10, round(umi_threshold)) # enforce a sensible minimum

    umi_thresholds[sample] <- umi_threshold  # for record keeping

    # Filter barcodes by threshold
    keep_barcodes <- names(umis)[umis >= umi_threshold]
    counts_filt   <- counts[, keep_barcodes, drop = FALSE]
    coords_filt   <- coords[keep_barcodes, , drop = FALSE]
    stopifnot(all(colnames(counts_filt) == rownames(coords_filt)))

    # Build query object for RCTD
    query <- SpatialRNA(coords_filt, counts_filt, Matrix::colSums(counts_filt))
    query_list[[sample]] <- query
    coords_list[[sample]] <- coords_filt

    # Optional: log counts/kept ratio per sample
    message(sample, ": kept ", length(keep_barcodes), " of ", length(umis), " pixels (UMI threshold: ", umi_threshold, ")")
}

# Save umi_thresholds for supplement/QC:
write.csv(data.frame(sample=names(umi_thresholds), umi_threshold=umi_thresholds), "umi_thresholds_per_sample.csv", row.names=FALSE)





#xenium_samples <- c("HP1","HP2_xen", "HT1_xen", "HT2a_xen", "HT2b", "HT3_xen", "LP1_xen", "LP2", "LT1_xen", "LT2", "HP13b_xen", "NC4_xen")
#visHD_samples  <- c("HT3", "HP2_hd", "HT2_hd", "LP1_hd")
visium_v1_samples <- c("HP13b_v1", "NC1a", "NC1b", 
    "NC1c", "ND5", "SP11", "SP8", "SP9", 
    "SP10", "ND6", "ND7", "SP12", "HP13a", 
    "NC2", "NC3", "NC4_v1")

visium_v1_samples <- c("NC1a", "NC1b", 
    "NC1c", "ND5", "SP11", "SP8", "SP9", 
    "SP10", "ND6", "ND7", "SP12", "HP13a", 
    "NC2", "NC3", "NC4_v1")
visium_v1_samples <- c("HP13b_v1", "NC1a", "NC1b", 
    "NC1c", "ND5", "SP11", "SP8", "SP9", 
    "SP10", "ND6", "ND7", "SP12", "HP13a", 
    "NC2", "NC3", "NC4_v1")

# [1] "LP1_v1" "LT1_v1" "HP14" missing in seurat_all? yes

visium_v1_samples <- c("HP13b_v1")

results_dir <- "/path/to/local_home/visiumhd_mnp/results/seurat_all_v1"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# Confirm groups
#cat("Xenium: ", paste(xenium_samples, collapse = ", "), "\n")
#cat("Visium HD: ", paste(visHD_samples, collapse = ", "), "\n")
cat("Visium v1/v2: ", paste(visium_v1_samples, collapse = ", "), "\n")


run_rctd_batch <- function(sample_names, doublet_mode, label) {
    for (sample in sample_names) {
        message(label, " - Processing: ", sample)
        query <- query_list[[sample]]

common_genes <- intersect(rownames(query@counts), rownames(reference.ref@counts))
ordered_genes <- sort(common_genes)
counts_query_ordered     <- query@counts[ordered_genes, ]
counts_reference_ordered <- reference.ref@counts[ordered_genes, ]
nUMI_reference_ordered   <- Matrix::colSums(counts_reference_ordered)
# Optional: coerce rownames to character
rownames(counts_query_ordered) <- as.character(rownames(counts_query_ordered))
rownames(counts_reference_ordered) <- as.character(rownames(counts_reference_ordered))
query_sub <- spacexr::SpatialRNA(coords = query@coords, counts = counts_query_ordered)
reference_sub <- spacexr::Reference(
    counts = counts_reference_ordered, 
    cell_types = reference.ref@cell_types, 
    nUMI = nUMI_reference_ordered
)
stopifnot(identical(rownames(query_sub@counts), rownames(reference_sub@counts)))
        rctd_obj <- create.RCTD(query_sub, reference_sub, max_cores = 1)
        saveRDS(rctd_obj, file.path(results_dir, paste0("rctd_", sample, "_pre.rds")))
        rctd_obj <- run.RCTD(rctd_obj, doublet_mode = doublet_mode) # <-- This is correct
        saveRDS(rctd_obj, file.path(results_dir, paste0("rctd_", sample, ".rds")))
    }
}
# Then call with:
run_rctd_batch(visium_v1_samples, doublet_mode = "doublet", label = "Visium v1/v2")





library(Seurat)
library(spacexr)

results_dir <- "/path/to/local_home/visiumhd_mnp/results/seurat_all_v1"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# samples <- names(query_list)
# Use only samples with existing RCTD results
rctd_result_files <- list.files(results_dir, pattern = "\\.rds$", full.names = FALSE)
samples <- gsub("^rctd_|\\.rds$", "", rctd_result_files[grep("^rctd_", rctd_result_files)])
# Or, if you want to stick to a specific group, set samples = visHD_samples, etc.
filter_log <- data.frame(Sample = character(), Total = integer(), Annotated = integer(),
                         Unannotated = integer(), Filtered = integer(), FilterType = character(),
                         stringsAsFactors = FALSE)

# Helper for universal barcode harmonization
clean_barcode <- function(x) gsub(".*_", "", gsub("-[12]$", "", x))

for (sample in samples) {
    message("Annotating and filtering: ", sample)
    rctd_obj <- readRDS(file.path(results_dir, paste0("rctd_", sample, ".rds")))
    seurat_obj <- seurat_all[[sample]]
    annotations.df <- rctd_obj@results$results_df

    # ----- CHECK: skip sample if annotation is invalid or empty -----
 if (!is.data.frame(annotations.df) || nrow(annotations.df) == 0) {
    warning(paste("No annotation data for sample", sample, "- skipping"))
    total_n <- suppressWarnings(tryCatch(ncol(seurat_obj), error = function(e) 0))
    if (length(total_n) == 0) total_n <- 0
    filter_log <- rbind(filter_log, data.frame(
        Sample = as.character(sample), 
        Total = as.integer(total_n), 
        Annotated = as.integer(0), 
        Unannotated = as.integer(total_n),
        Filtered = as.integer(0), 
        FilterType = as.character("no annotation"), 
        stringsAsFactors = FALSE
    ))
    next
}
    # ----- END CHECK -----

    # --- BARCODE HARMONIZATION ---
    annotation_barcodes_clean <- clean_barcode(rownames(annotations.df))
    seurat_barcodes_clean     <- clean_barcode(colnames(seurat_obj))
    rownames(annotations.df)  <- annotation_barcodes_clean
    names(seurat_barcodes_clean) <- colnames(seurat_obj) # Keep mapping to original barcodes

    # Build matched annotation dataframe, in the order of seurat_obj columns
    ann_aligned <- annotations.df[match(seurat_barcodes_clean, annotation_barcodes_clean), , drop=FALSE]
    rownames(ann_aligned) <- colnames(seurat_obj) # so AddMetaData matches directly

    # Add all RCTD metadata (will fill NA for unmatched barcodes)
    seurat_obj <- AddMetaData(seurat_obj, metadata = ann_aligned)

    # Assign predicted.celltype as character, regardless of column type
    if ("first_type" %in% colnames(seurat_obj@meta.data)) {
        seurat_obj$predicted.celltype <- as.character(seurat_obj$first_type)
    }

    # Assign confidence score if available (can support first_type_score or legacy names)
    if ("first_type_score" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$predicted.celltype.score <- seurat_obj$first_type_score
    } else if (!("predicted.celltype.score" %in% colnames(seurat_obj@meta.data))) {
        seurat_obj$predicted.celltype.score <- NA_real_
    }


    # Log counts
    total_n <- ncol(seurat_obj)
    annotated_n <- sum(!is.na(seurat_obj$predicted.celltype))
    unannotated_n <- sum(is.na(seurat_obj$predicted.celltype))

    # Save full annotated Seurat object
    saveRDS(seurat_obj, file.path(results_dir, paste0("seurat_", sample, "_annotated.rds")))

  # No filtering; export all cells
obj_filtered <- seurat_obj
filter_type <- "no filtering"
filtered_n <- ncol(obj_filtered)

# Handle case with no cells in object (should only happen if object is empty)
if (filtered_n == 0) {
    warning(paste("No cells present for sample", sample))
    filter_log <- rbind(filter_log, data.frame(
        Sample = as.character(sample),
        Total = as.integer(total_n),
        Annotated = as.integer(annotated_n),
        Unannotated = as.integer(unannotated_n),
        Filtered = as.integer(0),
        FilterType = as.character(filter_type),
        stringsAsFactors = FALSE
    ))
    next
}

# Output table: orig.ident x predicted.celltype (filtered)
t1 <- table(obj_filtered$orig.ident, obj_filtered$predicted.celltype)
outname <- paste0("rctd_", sample, "_celltype_counts_", gsub(" ", "_", filter_type), ".csv")
write.csv(t1, file.path(results_dir, outname))

}

# Save logs
write.csv(filter_log, file.path(results_dir, "rctd_filtering_log.csv"), row.names = FALSE)



cat(sprintf("%s: %d of %d cells annotated (%.1f%%)\n", sample, annotated_n, total_n, 100 * annotated_n / total_n))

pdf(file.path(results_dir, paste0("rctd_", sample, "_umap_predicted_celltype.pdf")))
print(DimPlot(seurat_obj, group.by = "predicted.celltype", label = TRUE))
dev.off()



        ## probably delete HP2_hd and HP2_xen 
library(Seurat)
seurat_files <- list.files(results_dir, pattern = "^seurat_.*_annotated.rds$", full.names = TRUE)
seurat_list <- lapply(seurat_files, readRDS)

obj_merged_annot <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = samples)
saveRDS(obj_merged_annot, file.path(results_dir, "seurat_all_RCTD_annotated_merged.rds"))

table(obj_merged_annot$predicted.celltype)
with(obj_merged_annot@meta.data, table(orig.ident, predicted.celltype))
annot_table <- obj_merged_annot@meta.data[, c("orig.ident", "predicted.celltype", "predicted.celltype.score")]
write.csv(annot_table, file.path(results_dir, "RCTD_annotations_per_cell.csv"))
not_annotated <- colnames(seurat_obj)[is.na(seurat_obj$predicted.celltype)]
writeLines(not_annotated, file.path(results_dir, paste0("unannotated_barcodes_", sample, ".txt")))






#Step-by-step Best Practice:
#Filter your Seurat object up front
#(before normalization, chunking, or SCTransform, e.g. obj <- subset(obj, nCount_Xenium.500 > 0))
#Run RCTD on the filtered data
#(RCTD will only analyze the barcodes/pixels present in your filtered data)
#Barcode harmonization
#(clean barcodes in both the RCTD annotation table and Seurat object as needed)
#Immediately before transferring annotations:
#Subset your annotation dataframe (from RCTD) to match only the barcodes present in your filtered Seurat object:
barcodes_to_keep <- colnames(filtered_seurat_obj)
ann_aligned <- ann_aligned[barcodes_to_keep, , drop=FALSE]
rownames(ann_aligned) <- barcodes_to_keep
filtered_seurat_obj <- AddMetaData(filtered_seurat_obj, ann_aligned)









cd /path/to/local_home/visiumhd_mnp/results/seurat_all_v1
ls
for f in *annotated.rds;
do
cp $f /path/to/local_storage/raid2/project/spatial_MNP/seurat_all_ref_merged
done




/path/to/local_storage/raid2/project/spatial_MNP/seurat_all_ref_merged


