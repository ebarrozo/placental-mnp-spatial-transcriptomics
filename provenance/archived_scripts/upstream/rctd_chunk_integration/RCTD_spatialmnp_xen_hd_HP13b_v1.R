# RCTD_spatialmnp_visiumv1.R

library(Seurat)
library(spacexr)
library(Matrix)

# Set your sample names and reference
visium_v1_samples <- c("HP13b_v1")
results_dir <- "/path/to/local_home/visiumhd_mnp/results/seurat_all_v1/HP13b_v1"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# Load objects
seurat_all <- readRDS("/path/to/local_storage/raid2/project/visiumhd_mnp/all_spatial_seurat_objects.rds")
query_obj <- seurat_all[["HP13b_v1"]]
reference <- readRDS("/path/to/local_storage/raid2/project/spatial_MNP/placental_ref_Jul_2025.rds") ## this is  11.5 gb
DefaultAssay(reference) <- "RNA"
Idents(reference) <- "cell_type_expanded"
levels(reference) <- gsub("/", "-", levels(reference))
t1 <- table(Idents(reference))
write.table(t1, "reference_cell_type_expanded.counts.txt", sep="\t")

counts <- GetAssayData(reference, assay = "RNA", layer = "counts")
counts <- round(counts)
cluster <- as.factor(reference$cell_type_expanded)
names(cluster) <- colnames(reference)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
reference.ref <- Reference(counts, cluster, nUMI)

# Helper function (unchanged)
get_query_for_sample <- function(obj, umi_pct = 0.10, min_umi = 10) {
    assay_name <- DefaultAssay(obj)
    counts <- GetAssayData(obj, assay = assay_name, layer = "counts")
    coords <- tryCatch({
        img_names <- names(obj@images)
        if (length(img_names) > 0 && "centroids" %in% names(obj@images[[img_names[1]]]@boundaries)) {
            centroids_obj <- obj@images[[img_names[1]]]@boundaries[["centroids"]]
            barcodes <- centroids_obj@cells
            coords_matrix <- centroids_obj@coords
            rownames(coords_matrix) <- barcodes
            coords_df <- as.data.frame(coords_matrix)
            colnames(coords_df) <- c("x", "y")
            coords_df
        } else {
            coords <- GetTissueCoordinates(obj, which = "centroids")
            if ("cell" %in% colnames(coords)) {
                rownames(coords) <- coords$cell
                coords$cell <- NULL
            }
            coords
        }
    }, error = function(e) { NULL })
    common_barcodes <- intersect(colnames(counts), rownames(coords))
    if (length(common_barcodes) == 0) return(NULL)
    counts <- counts[, common_barcodes, drop = FALSE]
    coords <- coords[common_barcodes, c("x", "y"), drop = FALSE]
    umis <- Matrix::colSums(counts)
    umi_threshold <- max(min_umi, round(quantile(umis, umi_pct)))
    keep_barcodes <- names(umis)[umis >= umi_threshold]
    if (length(keep_barcodes) == 0) return(NULL)
    counts_filt <- counts[, keep_barcodes, drop = FALSE]
    coords_filt <- coords[keep_barcodes, , drop = FALSE]
    if (!all(colnames(counts_filt) == rownames(coords_filt))) return(NULL)
    list(counts = counts_filt, coords = coords_filt)
}

# --- NEW: Remove all-zero genes before building RCTD objects ---
run_rctd_batch <- function(sample_names, reference.ref, doublet_mode = "doublet", label = "Visium v1/v2", results_dir = ".") {
    for (sample in sample_names) {
        message(label, " - Processing: ", sample)
        if (!sample %in% names(seurat_all)) {
            warning("Sample ", sample, " not found in seurat_all, skipping.")
            next
        }
        obj <- seurat_all[[sample]]
        query_info <- get_query_for_sample(obj)
        if (is.null(query_info)) {
            warning("Could not build query object for sample: ", sample)
            next
        }
        counts_query <- query_info$counts
        coords_query <- query_info$coords

        # Harmonize gene sets
        common_genes <- intersect(rownames(counts_query), rownames(reference.ref@counts))
        if (length(common_genes) < 2000) {
            warning("Low gene overlap for sample ", sample, " (n=", length(common_genes), "), skipping.")
            next
        }
        ordered_genes <- sort(common_genes)
        counts_query_ordered     <- counts_query[ordered_genes, , drop=FALSE]
        counts_reference_ordered <- reference.ref@counts[ordered_genes, , drop=FALSE]
        nUMI_reference_ordered   <- Matrix::colSums(counts_reference_ordered)

        # --- Remove all-zero genes ---
        genes_nonzero_query     <- rownames(counts_query_ordered)[Matrix::rowSums(counts_query_ordered) > 0]
        genes_nonzero_reference <- rownames(counts_reference_ordered)[Matrix::rowSums(counts_reference_ordered) > 0]
        good_genes <- intersect(genes_nonzero_query, genes_nonzero_reference)
        cat("Genes nonzero in both: ", length(good_genes), "\n")
        if (length(good_genes) < 2000) {
            warning("Low number of nonzero genes after filtering for sample ", sample, ". Skipping.")
            next
        }
        counts_query_final     <- counts_query_ordered[good_genes, , drop=FALSE]
        counts_reference_final <- counts_reference_ordered[good_genes, , drop=FALSE]
        nUMI_reference_final   <- Matrix::colSums(counts_reference_final)

        cat("Sample:", sample, "\n")
        cat("Query (final): ", nrow(counts_query_final), "genes x", ncol(counts_query_final), "barcodes\n")
        cat("Reference (final): ", nrow(counts_reference_final), "genes x", ncol(counts_reference_final), "cells\n")
        cat("First 5 genes (final): ", paste(head(rownames(counts_query_final), 5), collapse=", "), "\n")

        # Build RCTD query/reference
        query_sub <- spacexr::SpatialRNA(coords = coords_query, counts = counts_query_final)
        reference_sub <- spacexr::Reference(
            counts = counts_reference_final, 
            cell_types = reference.ref@cell_types, 
            nUMI = nUMI_reference_final
        )
        stopifnot(identical(rownames(query_sub@counts), rownames(reference_sub@counts)))

        # RCTD run
        rctd_obj <- create.RCTD(query_sub, reference_sub, max_cores = 1)
        saveRDS(rctd_obj, file.path(results_dir, paste0("rctd_", sample, "_pre.rds")))
        rctd_obj <- run.RCTD(rctd_obj, doublet_mode = doublet_mode)
        saveRDS(rctd_obj, file.path(results_dir, paste0("rctd_", sample, ".rds")))
        message("Done: ", sample)
    }
}

# Call the batch function
run_rctd_batch(
    sample_names = visium_v1_samples,
    reference.ref = reference.ref,
    doublet_mode = "doublet",
    label = "Visium v1/v2",
    results_dir = results_dir
)




library(Seurat)
library(spacexr)

results_dir <- "/path/to/local_home/visiumhd_mnp/results/seurat_all_v1/HP13b_v1"
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



annotated <- readRDS("/path/to/local_home/visiumhd_mnp/results/seurat_all_v1/HP13b_v1/seurat_HP13b_v1_annotated.rds") ## this is  11.5 gb



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







cd /path/to/local_home/visiumhd_mnp/results/seurat_all_v1/HP13b_v1
for f in *annotated.rds;
do
cp $f /path/to/local_storage/raid3/project/seurat_all_ref_merged
done



cd /path/to/local_home/visiumhd_mnp/results/seurat_all_v1
ls
for f in *annotated.rds;
do
cp $f /path/to/local_storage/raid2/project/spatial_MNP/seurat_all_ref_merged
done




/path/to/local_storage/raid2/project/spatial_MNP/seurat_all_ref_merged


