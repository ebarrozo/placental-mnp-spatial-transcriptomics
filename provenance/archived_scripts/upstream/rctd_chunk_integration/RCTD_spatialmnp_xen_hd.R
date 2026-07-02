RCTD_spatialmnp_xen_hd.R
# 8.02.25
  # server1

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

library(Seurat)
library(spacexr) # for SpatialRNA
library(purrr)   # for list operations
# --- LIBRARIES ---
library(tidyverse)
library(data.table)


## Set working directory and seed
setwd("/path/to/local_home/visiumhd_mnp/results/seurat_all_v1")

reference <- readRDS("/path/to/local_storage/raid2/project/spatial_MNP/placental_ref_Jul_2025.rds") ## this is  11.5 gb
# Integration with scRNA-seq data (deconvolution); Robust Cell Type Decomposition RCTD
library(spacexr)
# load the reference 
#reference <- readRDS("/path/to/local_storage/raid2/project/spatial_MNP/placental_ref_Jul_2025.rds") ## this is  11.5 gb
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
reference.ref <- Reference(counts, cluster, nUMI)




# Helper: optionally visualize UMI distributions per sample
plot_umi_dist <- function(umis, sample_name) {
    png(sprintf("umi_dist_%s.png", sample_name), width=900, height=600)
    hist(umis, breaks=100, main=paste("UMI per pixel:", sample_name), xlab="UMI", log="x")
    abline(v=median(umis), col="red", lty=2)
    dev.off()
}
seurat_all <- readRDS("/path/to/local_storage/raid2/project/visiumhd_mnp/all_spatial_seurat_objects.rds") ## this is  11.5 gb

query_list  <- list()
coords_list <- list()
umi_thresholds <- c()  # Save for record-keeping

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


# 1. Define your sample groups
#  "HP2","HP2", 
table(names(seurat_all))
# You'll see HP2 appears twice

# Let's manually fix:
# (assuming order is as shown in your post)
#names(seurat_all)[which(names(seurat_all) == "HP2")[1]] <- "HP2_xen"
names(seurat_all)[which(names(seurat_all) == "HP2")[1]] <- "HP2_hd"
for (sample in names(seurat_all)) {
  obj <- seurat_all[[sample]]
  obj$orig.ident <- sample
  seurat_all[[sample]] <- obj
}

# Confirm
table(names(seurat_all))


# 1. Define your sample groups

xenium_samples <- c("HP1", "HP2", "HT1_xen", "HT2a_xen", "HT2b", "HT3_xen", "LP1_xen", "LP2", "LT1_xen", "LT2", "HP13b_xen", "NC4_xen")
visHD_samples  <- c("HP2_hd", "HT3", "HT2_hd", "LP1_hd")
visHD_samples  <- c("HP2_hd", "HT2_hd", "LP1_hd")

# [1] "LP1_v1" "LT1_v1" "HP14" missing in seurat_all? yes

results_dir <- "/path/to/local_home/visiumhd_mnp/results/seurat_all_v1"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# Confirm groups
#cat("Xenium: ", paste(xenium_samples, collapse = ", "), "\n")
cat("Visium HD: ", paste(visHD_samples, collapse = ", "), "\n")
#cat("Visium v1/v2: ", paste(visium_v1_samples, collapse = ", "), "\n")


# 2. RCTD loop function (modular)
run_rctd_batch <- function(sample_names, doublet_mode, label) {
    for (sample in sample_names) {
        message(label, " - Processing: ", sample)
        query <- query_list[[sample]]
        rctd_obj <- create.RCTD(query, reference.ref, max_cores = 1) 
        saveRDS(rctd_obj, file.path(results_dir, paste0("rctd_", sample, "_pre.rds")))
        rctd_obj <- run.RCTD(rctd_obj, doublet_mode = doublet_mode)
        saveRDS(rctd_obj, file.path(results_dir, paste0("rctd_", sample, ".rds")))
    }
}

# 3. Run by group (change doublet_mode as needed)
#run_rctd_batch(xenium_samples,     doublet_mode = "full",    label = "Xenium")
run_rctd_batch(visHD_samples,      doublet_mode = "full",    label = "Visium HD")
#run_rctd_batch(visium_v1_samples,  doublet_mode = "doublet", label = "Visium v1/v2")


HP1_pre <- readRDS("rctd_HP1_pre.rds") 
HP1_rctd <- readRDS("rctd_HP1.rds") 
n, Unannotated = unannotated_n,
+         Filtered = filtered_n, FilterType = filter_type, stringsAsFactors = FALSE
+     ))
+ }
Annotating and filtering: HP1_pre
Error in data.frame(Sample = sample, Total = ncol(seurat_obj), Annotated = 0,  : 
  arguments imply differing number of rows: 1, 0
In addition: Warning message:
No weights data for sample HP1_pre - skipping 


library(Seurat)
library(spacexr)

rctd_result_files <- list.files(results_dir, pattern = "\\.rds$", full.names = FALSE)
samples <- gsub("^rctd_|\\.rds$", "", rctd_result_files[grep("^rctd_", rctd_result_files)])

filter_log <- data.frame(Sample = character(), Total = integer(), Annotated = integer(),
                         Unannotated = integer(), Filtered = integer(), FilterType = character(),
                         stringsAsFactors = FALSE)

for (sample in samples) {
    message("Annotating and filtering: ", sample)
    rctd_obj <- readRDS(file.path(results_dir, paste0("rctd_", sample, ".rds")))
    seurat_obj <- seurat_all[[sample]]

    # -- Get barcodes from Seurat or from spatialRNA if empty --
    seurat_barcodes <- colnames(seurat_obj)
    if (is.null(seurat_barcodes) || length(seurat_barcodes) == 0) {
        # If Seurat object is empty, build a new object from spatialRNA slot
        if ("spatialRNA" %in% slotNames(seurat_obj)) {
            pixel_counts <- seurat_obj@spatialRNA@counts
            pixel_barcodes <- pixel_counts@Dimnames[[2]]
            pixel_features <- pixel_counts@Dimnames[[1]]
            # Convert to standard matrix if needed
            pixel_counts_mat <- as.matrix(pixel_counts)
            # Create new Seurat object
            seurat_obj <- CreateSeuratObject(counts = pixel_counts_mat, assay = "Spatial")
            seurat_barcodes <- colnames(seurat_obj)
        } else {
            warning(paste("No barcodes in Seurat object or spatialRNA slot for", sample, "- skipping"))
            next
        }
    }

    weights.df <- as.matrix(rctd_obj@results$weights)
    weights_barcodes <- rownames(weights.df)

    # --- Check and match barcodes directly ---
    if (length(seurat_barcodes) == 0 || length(weights_barcodes) == 0) {
        warning(paste("No barcodes for sample", sample, "- skipping"))
        next
    }

    common_barcodes <- intersect(seurat_barcodes, weights_barcodes)
    if (length(common_barcodes) == 0) {
        warning(paste("No matching barcodes for sample", sample, "- skipping"))
        next
    }

    weights.df <- weights.df[common_barcodes, , drop = FALSE]
    seurat_obj <- subset(seurat_obj, cells = common_barcodes)
    seurat_barcodes <- colnames(seurat_obj)

    if (nrow(weights.df) == 0 || length(seurat_barcodes) == 0) {
        warning(paste("Nothing left after subsetting for", sample, "- skipping"))
        next
    }

    main_cell_type_idx <- apply(weights.df, 1, which.max)
    main_cell_type     <- colnames(weights.df)[main_cell_type_idx]
    max_score          <- apply(weights.df, 1, max)

    if (!(length(main_cell_type) == length(max_score) && length(main_cell_type) == nrow(weights.df))) {
        cat("Diagnostic for sample:", sample, "\n")
        cat("  nrow(weights.df):", nrow(weights.df), "\n")
        cat("  length(main_cell_type):", length(main_cell_type), "\n")
        cat("  length(max_score):", length(max_score), "\n")
        warning(paste("Vector length mismatch for sample", sample, "- skipping"))
        next
    }

    weights.df <- as.data.frame(weights.df)
    weights_annot_df <- data.frame(
        barcode = rownames(weights.df),
        predicted.celltype = main_cell_type,
        predicted.celltype.score = max_score,
        stringsAsFactors = FALSE
    )
    weights_annot_df <- cbind(weights_annot_df, weights.df)
    rownames(weights_annot_df) <- weights_annot_df$barcode

    ann_aligned <- weights_annot_df[seurat_barcodes, , drop=FALSE]
    rownames(ann_aligned) <- seurat_barcodes

    # Add metadata to Seurat object
    seurat_obj <- AddMetaData(seurat_obj, ann_aligned)

    # Log counts
    total_n <- ncol(seurat_obj)
    annotated_n <- sum(!is.na(seurat_obj$predicted.celltype))
    unannotated_n <- sum(is.na(seurat_obj$predicted.celltype))

    # Save full annotated Seurat object
    saveRDS(seurat_obj, file.path(results_dir, paste0("seurat_", sample, "_annotated.rds")))

    # Filtering by celltype confidence score if present
    if ("predicted.celltype.score" %in% colnames(seurat_obj@meta.data)) {
        obj_filtered <- subset(seurat_obj, subset = predicted.celltype.score > 0.1)
        filter_type <- "score > 0.1"
    } else {
        obj_filtered <- seurat_obj
        filter_type <- "no score filter"
        message("No predicted.celltype.score column for ", sample, "; exporting all cells")
    }
    filtered_n <- ncol(obj_filtered)

    # Handle case with no passing cells
    if (filtered_n == 0) {
        warning(paste("No cells pass filter for sample", sample))
        filter_log <- rbind(filter_log, data.frame(
            Sample = sample, Total = total_n, Annotated = annotated_n, Unannotated = unannotated_n,
            Filtered = 0, FilterType = filter_type, stringsAsFactors = FALSE
        ))
        next
    }

    # Output table: orig.ident x predicted.celltype (filtered)
    t1 <- table(obj_filtered$orig.ident, obj_filtered$predicted.celltype)
    outname <- paste0("rctd_", sample, "_celltype_counts_", gsub(" ", "_", filter_type), ".csv")
    write.csv(t1, file.path(results_dir, outname))

    # Update filter log
    filter_log <- rbind(filter_log, data.frame(
        Sample = sample, Total = total_n, Annotated = annotated_n, Unannotated = unannotated_n,
        Filtered = filtered_n, FilterType = filter_type, stringsAsFactors = FALSE
    ))
}

# Save logs
write.csv(filter_log, file.path(results_dir, "rctd_filtering_log.csv"), row.names = FALSE)

all_annotated_files <- list.files(pattern = "_annotated\\.rds$")
write.csv(all_annotated_files, file="annotated.filenames.csv")






library(Seurat)
library(Matrix)

# Clean barcode helper
clean_barcode <- function(x) gsub(".*_", "", gsub("-[12]$", "", x))

# Load objects if not already loaded
# HP1_pre <- readRDS("/path/to/HP1_pre.rds")
# HP1_rctd <- readRDS("/path/to/rctd_HP1_pre.rds")

# 1. Extract and check weights
weights_mat <- as.matrix(HP1_rctd@results$weights)
if (is.null(weights_mat) || nrow(weights_mat) == 0 || ncol(weights_mat) == 0) {
    stop("RCTD weights matrix is empty for HP1_pre.")
}

# 2. Harmonize barcodes
weights_barcodes_clean <- clean_barcode(rownames(weights_mat))
seurat_barcodes_clean  <- clean_barcode(colnames(HP1_pre))

rownames(weights_mat) <- weights_barcodes_clean
names(seurat_barcodes_clean) <- colnames(HP1_pre)

# 3. Check for matching barcodes
matched_idx <- match(seurat_barcodes_clean, weights_barcodes_clean)
if (all(is.na(matched_idx))) stop("No matching barcodes after harmonization for HP1_pre.")

# 4. Assign main cell type and score
main_cell_type_idx <- apply(weights_mat, 1, which.max)
main_cell_type     <- colnames(weights_mat)[main_cell_type_idx]
max_score          <- apply(weights_mat, 1, max)

if (length(main_cell_type) == 0 || all(is.na(main_cell_type))) stop("No valid main cell type for HP1_pre.")

# 5. Build annotation data.frame
weights_annot_df <- data.frame(
    barcode = rownames(weights_mat),
    predicted.celltype = main_cell_type,
    predicted.celltype.score = max_score,
    weights_mat,
    stringsAsFactors = FALSE
)
rownames(weights_annot_df) <- weights_annot_df$barcode

# 6. Align annotation to Seurat object column order
ann_aligned <- weights_annot_df[match(seurat_barcodes_clean, weights_barcodes_clean), , drop=FALSE]
rownames(ann_aligned) <- colnames(HP1_pre)

# 7. Add to Seurat object (NAs for unmatched)
HP1_pre <- AddMetaData(HP1_pre, ann_aligned)

# 8. Log annotation stats
cat("Total barcodes:", ncol(HP1_pre), "\n")
cat("Annotated barcodes:", sum(!is.na(HP1_pre$predicted.celltype)), "\n")
cat("Unannotated barcodes:", sum(is.na(HP1_pre$predicted.celltype)), "\n")
cat("Annotation complete!\n")

# 9. Save if needed
# saveRDS(HP1_pre, "HP1_pre_annotated.rds")



if (all(is.na(matched_idx))) stop("No matching barcodes after harmonization for HP1_pre.")
Error: No matching barcodes after harmonization for HP1_pre.