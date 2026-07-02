## HP2_hd.R
# 7.31.25
  # server2
## I noticed HP2 is the seurat orig.id for both the xenium and visiumHD object so going back and running through the workflow 
	# renamining both 



# Load session information
sessionInfo()

## Set working directory and seed
setwd("/path/to/local_home/visiumhd_mnp/results/seurat_all_v1")
set.seed(1)
options(future.globals.maxSize = 90e9) ## 300 gb ram 
## free -h ## in terminal

# --- LIBRARIES ---
library(Seurat)
library(tidyverse)
library(data.table)

# --- METADATA ---
metadata <- read.csv("/path/to/local_home/visiumhd_mnp/docs/spatial_cohort_annotation_spatialMedian_quickannot_HP2.csv", stringsAsFactors=FALSE)

# --- XENIUM (HARDCODED, NAMED WITH XenID_seurat) ---
xenium_objs <- list()

# Map: XeniumID_seq  → XenID_seurat
xenium_map <- setNames(metadata$XenID_seurat, metadata$XeniumID_seq)

# Hardcoded Xenium loads with object names set from metadata
xenium_objs[[xenium_map["S31542"]]]   <- LoadXenium("/path/to/local_storage/raid2/project/xenium_1_2025/Xenium/20250108__200739__Belfort_1/20250108__200739__Belfort_1/output-XETG00320__0023248__Region_2__20250108__200832", fov = "S31542", assay = "Xenium.5k")

# Identify 500 and 5k gene panel samples from metadata (if annotated), else list manually
xenium_5k  <- metadata$XenID_seurat[grepl("5k",  metadata$Xenium_Assay, ignore.case = TRUE)]

# If not in metadata, use your manual groupings:
# xenium_500 <- c("HT1_xen", "HT2a_xen", "HP13b_xen", "NC4_xen") 
# xenium_5k  <- c("HP1", "LP1_xen", "LP2", "HP2", "HT2b", "LT1_xen", "LT2", "HT3_xen")

# Apply gene panel-specific filters
for (xid in names(xenium_objs)) {
    if (xid %in% xenium_5k) {
        xenium_objs[[xid]] <- subset(xenium_objs[[xid]], subset = nCount_Xenium.5k > 0 & nFeature_Xenium.5k > 0)
        xenium_objs[[xid]]$Assay_GenePanel <- "5k"
    }
}

# Optional: Add Assay_GenePanel to metadata for downstream compatibility
metadata$Assay_GenePanel <- NA
metadata$Assay_GenePanel[metadata$XenID_seurat %in% xenium_5k]  <- "5k"



# --- Add percent.mt and cell cycle scoring to Xenium objects (500 and 5k) ---

mt.genes <- c("MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2", "MT-ATP8", "MT-ATP6", "MT-CO3", "MT-ND3", "MT-ND4L", "MT-ND4", "MT-ND5", "MT-ND6", "MT-CYB")
s.genes  <- c("MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8")
g2m.genes <- c("HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA")

for (xid in names(xenium_objs)) {
    if (!(xid %in% xenium_5k)) next  # Only process 5k objects

    obj <- xenium_objs[[xid]]
    DefaultAssay(obj) <- "Xenium.5k"

    gene_names <- rownames(obj)
    avail.mt <- intersect(gene_names, mt.genes)
    if (length(avail.mt) > 0) {
        obj <- PercentageFeatureSet(obj, features = avail.mt, col.name = "percent.mt")
    } else {
        obj$percent.mt <- 0
    }

    obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)

    avail.s   <- intersect(gene_names, s.genes)
    avail.g2m <- intersect(gene_names, g2m.genes)
    obj <- CellCycleScoring(obj, s.features = avail.s, g2m.features = avail.g2m, set.ident = FALSE)
    obj$CC.Difference <- obj$S.Score - obj$G2M.Score

    xenium_objs[[xid]] <- obj
}


# --- VISIUMHD ---
visiumhd_objs <- list()
visiumhd_map <- setNames(metadata$VisHD_seurat, metadata$VisiumHD_seq)
visiumhd_samples <- metadata$VisiumHD_seq[!is.na(metadata$VisiumHD_seq) & metadata$VisiumHD_seq != ""]
for (i in seq_along(visiumhd_samples)) {
  sample_path <- visiumhd_samples[i]
  sample_id   <- visiumhd_map[[sample_path]]
  obj <- Load10X_Spatial(
    data.dir     = file.path("/path/to/local_storage/raid2/project/visiumhd_mnp/customref", sample_path, "outs"),
    slice        = sample_id,
    bin.size     = c(8, 16),
    filter.matrix= TRUE
  )
  visiumhd_objs[[sample_id]] <- obj
}
# ---- Filter VisiumHD objects on Spatial.008um bin ----
for (sample_id in names(visiumhd_objs)) {
  obj <- visiumhd_objs[[sample_id]]
  # Only keep spots with nCount > 0 and nFeature > 0 in the 8um bin
  if ("Spatial.008um" %in% Assays(obj)) {
    obj <- subset(obj, subset = nCount_Spatial.008um > 0 & nFeature_Spatial.008um > 0)
    DefaultAssay(obj) <- "Spatial.008um"
    visiumhd_objs[[sample_id]] <- obj
  } else {
    warning(paste("Sample", sample_id, "does not contain Spatial.008um assay."))
  }
}
# ---- QC chunk: add percent.mt and cell cycle scoring to VisiumHD objects ----

mt.genes   <- c("MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2", "MT-ATP8", "MT-ATP6", "MT-CO3", "MT-ND3", "MT-ND4L", "MT-ND4", "MT-ND5", "MT-ND6", "MT-CYB")
s.genes    <- c("MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8")
g2m.genes  <- c("HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA")

for (sample_id in names(visiumhd_objs)) {
  obj <- visiumhd_objs[[sample_id]]
  # Set default assay to Spatial.008um for QC
  DefaultAssay(obj) <- "Spatial.008um"
  gene_names <- rownames(obj)
  # Add percent.mt
  avail.mt   <- intersect(gene_names, mt.genes)
  if (length(avail.mt) > 0) {
    obj <- PercentageFeatureSet(obj, features = avail.mt, col.name = "percent.mt")
  } else {
    obj$percent.mt <- 0
  }
  # Normalize for cell cycle
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  # Add cell cycle scores
  avail.s    <- intersect(gene_names, s.genes)
  avail.g2m  <- intersect(gene_names, g2m.genes)
  obj <- CellCycleScoring(obj, s.features = avail.s, g2m.features = avail.g2m, set.ident = FALSE)
  obj$CC.Difference <- obj$S.Score - obj$G2M.Score
  visiumhd_objs[[sample_id]] <- obj
}



# --- ANNOTATE METADATA TO EACH OBJECT ---
add_metadata <- function(obj, meta_row) {
  for (col in colnames(meta_row)) {
    obj[[col]] <- meta_row[[col]]
  }
  obj
}

for (id in names(xenium_objs)) {
  idx <- which(metadata$XenID_seurat == id)
  if (length(idx) > 0) xenium_objs[[id]] <- add_metadata(xenium_objs[[id]], metadata[idx, ])
}
for (id in names(visiumhd_objs)) {
  idx <- which(metadata$VisHD_seurat == id)
  if (length(idx) > 0) visiumhd_objs[[id]] <- add_metadata(visiumhd_objs[[id]], metadata[idx, ])
}


# --- FINAL OUTPUT ---
seurat_all <- c(xenium_objs, visiumhd_objs)
# All objects have harmonized metadata and are ready for merging/analysis

# rm(list = ls()[!ls() %in% c("seurat_all")])
gc()


# (Optional) save as RDS
saveRDS(seurat_all, file="all_spatial_seurat_objects_HP2.rds")


library(Seurat)

# This merges all the objects in the list
seurat_all_merged <- merge(
  x = seurat_all[[1]],
  y = seurat_all[-1],
  add.cell.ids = names(seurat_all), # Optional: use names as prefixes for cell barcodes
  project = "spatial_mnp"
)

saveRDS(seurat_all_merged, file="all_spatial_seurat_obj_HP2.rds")






















library(Seurat)
library(spacexr) # for SpatialRNA
library(purrr)   # for list operations

# Output lists
query_list <- list()
coords_list <- list()

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
visHD_samples  <- c("HP2_hd")
xenium_samples <- c("HP2_xen")
results_dir <- "/path/to/local_home/visiumhd_mnp/results/seurat_all_v1"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# Confirm groups
cat("Xenium: ", paste(xenium_samples, collapse = ", "), "\n")
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
run_rctd_batch(xenium_samples,     doublet_mode = "full",    label = "Xenium")
run_rctd_batch(visHD_samples,      doublet_mode = "full",    label = "Visium HD")
#run_rctd_batch(visium_v1_samples,  doublet_mode = "doublet", label = "Visium v1/v2")



library(Seurat)
library(spacexr)

rctd_result_files <- list.files(results_dir, pattern = "\\.rds$", full.names = FALSE)
samples <- gsub("^rctd_|\\.rds$", "", rctd_result_files[grep("^rctd_", rctd_result_files)])

filter_log <- data.frame(Sample = character(), Total = integer(), Annotated = integer(),
                         Unannotated = integer(), Filtered = integer(), FilterType = character(),
                         stringsAsFactors = FALSE)

# Helper for universal barcode harmonization
clean_barcode <- function(x) gsub(".*_", "", gsub("-[12]$", "", x))

for (sample in samples) {
    message("Annotating and filtering: ", sample)
    rctd_obj <- readRDS(file.path(results_dir, paste0("rctd_", sample, ".rds")))
    seurat_obj <- seurat_all[[sample]]
    weights.df <- rctd_obj@results$weights

    # ----- CHECK: skip sample if weights is invalid or empty -----
    if (!is.data.frame(weights.df) || nrow(weights.df) == 0) {
        warning(paste("No weights data for sample", sample, "- skipping"))
        filter_log <- rbind(filter_log, data.frame(
            Sample = sample, Total = ncol(seurat_obj), Annotated = 0, Unannotated = ncol(seurat_obj),
            Filtered = 0, FilterType = "no annotation", stringsAsFactors = FALSE
        ))
        next
    }
    # ----- END CHECK -----

    # --- BARCODE HARMONIZATION ---
    weights_barcodes_clean <- clean_barcode(rownames(weights.df))
    seurat_barcodes_clean  <- clean_barcode(colnames(seurat_obj))
    rownames(weights.df)   <- weights_barcodes_clean
    names(seurat_barcodes_clean) <- colnames(seurat_obj) # mapping to original barcodes

    # --- Assign main cell type and score ---
    main_cell_type_idx <- apply(weights.df, 1, which.max)
    main_cell_type     <- colnames(weights.df)[main_cell_type_idx]
    max_score          <- apply(weights.df, 1, max)

    # Merge into one dataframe for AddMetaData
    weights_annot_df <- data.frame(
        barcode = rownames(weights.df),
        predicted.celltype = main_cell_type,
        predicted.celltype.score = max_score,
        weights.df, # All cell type proportions (optional, for downstream)
        stringsAsFactors = FALSE
    )
    rownames(weights_annot_df) <- weights_annot_df$barcode

    # Build aligned annotation dataframe in order of seurat_obj columns
    ann_aligned <- weights_annot_df[match(seurat_barcodes_clean, weights_barcodes_clean), , drop=FALSE]
    rownames(ann_aligned) <- colnames(seurat_obj)

    # Add all RCTD metadata (NAs for unmatched barcodes)
    seurat_obj <- AddMetaData(seurat_obj, metadata = ann_aligned)

    # Log counts
    total_n <- ncol(seurat_obj)
    annotated_n <- sum(!is.na(seurat_obj$predicted.celltype))
    unannotated_n <- sum(is.na(seurat_obj$predicted.celltype))

    # Save full annotated Seurat object
    saveRDS(seurat_obj, file.path(results_dir, paste0("seurat_", sample, "_annotated.rds")))

    # Filtering by celltype confidence score if present
    if ("predicted.celltype.score" %in% colnames(seurat_obj@meta.data)) {
        obj_filtered <- subset(seurat_obj, subset = predicted.celltype.score > 0.5)
        filter_type <- "score > 0.5"
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



HP2_hd <- readRDS("rctd_HP2_hd.rds") ## this is  11.5 gb
HP2_xen <- readRDS("rctd_HP2_xen.rds") ## this is  11.5 gb
HP2_hd_pre <- readRDS("rctd_HP2_hd_pre.rds") ## this is  11.5 gb
HP2_xen_pre <- readRDS("rctd_HP2_xen_pre.rds") ## this is  11.5 gb




cd /path/to/local_storage/raid2/project/visiumhd_mnp/
mkdir hp2_named
cp /path/to/local_home/visiumhd_mnp/results/seurat_all_v1/ .
