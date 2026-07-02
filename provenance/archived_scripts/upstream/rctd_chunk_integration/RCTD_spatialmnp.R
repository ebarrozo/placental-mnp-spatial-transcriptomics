## RCTD_spatialmnp.R
# 7.29.25
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




seurat_all <- readRDS("/path/to/local_storage/raid2/project/visiumhd_mnp/all_spatial_seurat_objects.rds") ## this is  11.5 gb
# obj <- readRDS("/path/to/local_storage/raid2/project/visiumhd_mnp/all_spatial_seurat_obj.rds") ## this is  11.5 gb

library(Seurat)
library(spacexr) # for SpatialRNA
library(purrr)   # for list operations

# Output lists
query_list <- list()
coords_list <- list()

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
    query <- SpatialRNA(coords, counts, colSums(counts))
    query_list[[sample]] <- query
    coords_list[[sample]] <- coords
}

reference.ref <- readRDS("/path/to/local_storage/raid2/project/spatial_MNP/placental_ref_Jul_2025.rds") ## this is  11.5 gb
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



# 1. Define your sample groups
#visHD_samples  <- c("HP2", "HT3", "HT2_hd", "LP1_hd")
#xenium_samples <- c("HP1", "HP2", "HT1_xen", "HT2a_xen", "HT2b", "HT3_xen", "LP1_xen", "LP2", "LT1_xen", "LT2", "HP13b_xen", "NC4_xen")
#visium_v1_samples <- c("HT1_v1", "HT2_v1", "HT3_v1", "LP1_v1", "LT1_v1", "HP13b_v1", "NC1a", "NC1b", "NC1c", "ND5", "SP11", "SP8", "SP9", "SP10", "ND6", "ND7", "SP12", "HP13a", "NC2", "NC3", "HP14", "NC4_v1")
#visHD_samples  <- c()
visium_v1_samples <- c("HP13b_v1", "NC1a", "NC1b", 
	"NC1c", "ND5", "SP11", "SP8", "SP9", 
	"SP10", "ND6", "ND7", "SP12", "HP13a", 
	"NC2", "NC3", "NC4_v1")
#xenium_samples <- c("HT2b", "HT3_xen", "LP2", "LT1_xen", "LT2")

# [1] "LP1_v1" "LT1_v1" "HP14" missing in seurat_all? yes

# Confirm groups
#cat("Xenium: ", paste(xenium_samples, collapse = ", "), "\n")
#cat("Visium HD: ", paste(visHD_samples, collapse = ", "), "\n")
cat("Visium v1/v2: ", paste(visium_v1_samples, collapse = ", "), "\n")

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
#run_rctd_batch(visHD_samples,      doublet_mode = "full",    label = "Visium HD")
run_rctd_batch(visium_v1_samples,  doublet_mode = "doublet", label = "Visium v1/v2")


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
        filter_log <- rbind(filter_log, data.frame(
            Sample = sample, Total = ncol(seurat_obj), Annotated = 0, Unannotated = ncol(seurat_obj),
            Filtered = 0, FilterType = "no annotation", stringsAsFactors = FALSE
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

}

# Save logs
write.csv(filter_log, file.path(results_dir, "rctd_filtering_log.csv"), row.names = FALSE)



cat(sprintf("%s: %d of %d cells annotated (%.1f%%)\n", sample, annotated_n, total_n, 100 * annotated_n / total_n))

pdf(file.path(results_dir, paste0("rctd_", sample, "_umap_predicted_celltype.pdf")))
print(DimPlot(seurat_obj, group.by = "predicted.celltype", label = TRUE))
dev.off()




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






















seurat_HP1_annotated.rds

HP1.obj <- readRDS("seurat_HP1_annotated.rds") ## this is  11.5 gb
Idents(HP1.obj) <- "predicted.celltype"
levels(HP1.obj)
table(Idents(HP1.obj))
ncol(HP1.obj) # 43376













# --- UNIVERSAL BARCODE HARMONIZATION ---
clean_barcode <- function(x) gsub(".*_", "", gsub("-[12]$", "", x))
annotation_barcodes_clean <- clean_barcode(rownames(annotations.df))
seurat_barcodes_clean <- clean_barcode(colnames(seurat_obj))
annot_map <- annotations.df$first_type
names(annot_map) <- annotation_barcodes_clean
preds_for_seurat <- annot_map[seurat_barcodes_clean]
num_na <- sum(is.na(preds_for_seurat))
if (num_na > 0) {
    warning(paste("Sample", sample, ":", num_na, "barcodes not found in annotation after harmonization."))
    # preds_for_seurat[is.na(preds_for_seurat)] <- "Unknown"
}
seurat_obj$predicted.celltype <- preds_for_seurat



preds_for_seurat <- rep(NA_character_, ncol(seurat_obj))
names(preds_for_seurat) <- colnames(seurat_obj)
preds_for_seurat[names(annotations)] <- as.character(annotations)
seurat_obj$predicted.celltype <- preds_for_seurat
Idents(seurat_obj) <- "predicted.celltype.score"
levels(seurat_obj)
seurat_obj$predicted.celltype <- factor(seurat_obj$predicted.celltype)
seurat_obj$predicted.celltype <- droplevels(seurat_obj$predicted.celltype)
table(Idents(seurat_obj))







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
reference.ref <- Reference(counts, cluster, nUMI) ## 4.1 GB with whole ref
#> reference <- Reference(counts, cluster, nUMI) ## 4.1 GB
#Warning message:
#In Reference(counts, cluster, nUMI) :
#  Reference: number of cells per cell type is 122659, larger than maximum allowable of 10000. Downsampling number of cells to: 10000
rm(list = ls()[!ls() %in% c("obj", "seurat_all", "reference", "query")])
gc()

# run RCTD
RCTD <- create.RCTD(query, reference.ref, max_cores = 1)
saveRDS(RCTD, "/path/to/local_storage/raid2/project/visiumhd_mnp/RCTD.rds") # 4.1GB

rm(list = ls()[!ls() %in% c("obj", "RCTD")])
gc()

?run.RCTD()
	## doublet mode? or subcellular for xenium/visiumhd/visiumv1?
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
obj <- AddMetaData(obj, metadata = RCTD@results$results_df)
annotations.df <- RCTD@results$results_df
annotations <- annotations.df$first_type
names(annotations) <- rownames(annotations.df)
obj$predicted.celltype <- annotations

saveRDS(obj, "/path/to/local_home/spatial_mnp/results/seurat_v1/xen_obj_annot.rds")## this is  gb

# keep.cells <- Cells(obj)[!is.na(obj$predicted.celltype)]
ncol(obj) # 1833667
obj_subset <- subset(obj, subset = predicted.cell_type_expanded.score > 0.5)
ncol(obj_subset) # 1422064 (77.6% over 50% confidence)
table(obj.integrated$orig.ident, obj.integrated$predicted.cell_type_expanded)
t1<- table(obj.integrated$orig.ident, obj.integrated$predicted.cell_type_expanded)
write.csv(t1, "obj-orig.ident_predicted.cell_type_expanded.score_over0.5.csv", row.names = T)





############################### ############################### ###############################
###############################      UMAP with predicted.celltype         ###############################
############################### ############################### ###############################
# obj <- readRDS("/path/to/local_home/spatial_mnp/results/seurat_v1/xen_obj_annot.rds") ## this is  gb

DefaultAssay(obj) <- "Xenium.500"
obj <- NormalizeData(obj)
obj <- ScaleData(obj)
obj <- FindVariableFeatures(obj, verbose=TRUE)
#top.genes <- obj@assays[["Xenium.500"]]@meta.data[["var.features"]]
# View(top.genes)
#top.genes <- top.genes[!is.na(top.genes)]
#obj <- SCTransform(obj, assay = "Xenium.500") ## this will take a minute
obj <- RunPCA(obj, npcs = 30, features = rownames(obj)) ## uses all available features, since some objects have 500 and some have 5000?
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
obj <- FindClusters(obj, resolution = 0.6)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, reduction.name = "pca.umap")


rm(list = ls()[!ls() %in% c("obj", "mypal3")])
gc()
#obj <- IntegrateLayers(
#  object = obj, method = RPCAIntegration, assay="Xenium.500",
#  orig.reduction = "pca", new.reduction = "integrated.rpca",
#  verbose = TRUE
#obj <- FindNeighbors(obj, assay = "sketch", reduction = "integrated.rpca", dims = 1:50)
#obj <- FindClusters(obj, cluster.name = "seurat_cluster.sketched", resolution = 0.6)
#obj <- RunUMAP(obj, reduction = "integrated.rpca", dims = 1:30, reduction.name = "integrated.rpca")
      ################# There are no layers to integrate #################
mypal3 <-get_palette("igv",40)
p2<-DimPlot(obj, reduction = "pca.umap", group.by = "predicted.celltype", alpha = 0.1, label = FALSE,raster=T, cols = mypal3) + ggtitle("Predicted Cell Types (1,833,667 Cells)") + NoAxes()
p2
ggsave("xenium.obj_annot_pca.umap.png", plot = p2, device = "png", width = 8, height = 5, units = "in", dpi = 400)
p2<-DimPlot(obj, reduction = "pca.umap", group.by = "MNPBurden", alpha = 0.1, label = FALSE,raster=T, cols = mypal3) + ggtitle("MNP Burden (Median threshold)") + NoAxes()
p2
ggsave("xenium.obj_MNPBurden_pca.umap.png", plot = p2, device = "png", width = 8, height = 5, units = "in", dpi = 400)
p2<-DimPlot(obj, reduction = "pca.umap", group.by = "GestationalAge", alpha = 0.1, label = FALSE,raster=T, cols = mypal3) + ggtitle("Gestational Age at Delivery") + NoAxes()
p2
ggsave("xenium.obj_GA_pca.umap.png", plot = p2, device = "png", width = 8, height = 5, units = "in", dpi = 400)
p2<-DimPlot(obj, reduction = "pca.umap", group.by = "MNP_GA", alpha = 0.1, label = FALSE,raster=T, cols = mypal3) + ggtitle("MNP and GA") + NoAxes()
p2
ggsave("xenium.obj_MNP_GA_pca.umap.png", plot = p2, device = "png", width = 8, height = 5, units = "in", dpi = 400)
p2<-DimPlot(obj, reduction = "pca.umap", group.by = "Pregnancy_Infections", alpha = 0.1, label = FALSE,raster=T, cols = mypal3) + ggtitle("Pregnancy Infections") + NoAxes()
p2


ggsave("xenium.obj_Pregnancy_Infections_pca.umap.png", plot = p2, device = "png", width = 8, height = 5, units = "in", dpi = 400)
p2<-DimPlot(obj, reduction = "pca.umap", group.by = "Fetal_sex", alpha = 0.1, label = FALSE,raster=T, cols = mypal3) + ggtitle("Fetal_sex") + NoAxes()
p2
ggsave("xenium.obj_Fetal_sex_pca.umap.png", plot = p2, device = "png", width = 8, height = 5, units = "in", dpi = 400)
p2<-DimPlot(obj, reduction = "pca.umap", group.by = "Ethnicity_Race", alpha = 0.1, label = FALSE,raster=T, cols = mypal3) + ggtitle("Ethnicity_Race") + NoAxes()
p2
ggsave("xenium.obj_Ethnicity_Race_pca.umap.png", plot = p2, device = "png", width = 8, height = 5, units = "in", dpi = 400)
p2<-DimPlot(obj, reduction = "pca.umap", group.by = "orig.ident", alpha = 0.1, label = FALSE,raster=T, cols = mypal3) + ggtitle("SampleID") + NoAxes()
p2
ggsave("xenium.obj_SampleID_pca.umap.png", plot = p2, device = "png", width = 8, height = 5, units = "in", dpi = 400)


Idents(obj) <- "orig.ident"  
cell.types.counts <- table(Idents(obj))
write.table(cell.types.counts, "orig.ident_counts_obj.txt", sep="\t")
################ This is the final UMAP $$$$$$$$$$$$$$$$$$
obj <- DietSeurat(obj, layers="counts", dimreducs="pca.umap")
obj <- JoinLayers(obj) ## join layers if wanting to run differential expression downstream
# saveRDS(obj, "/path/to/local_home/spatial_mnp/results/seurat_v1/xenium.obj_annot_UMAP.rds")


############################### ############################### ###############################
###############################      Run Banksy to annotate banksy niches         ###############################
############################### ############################### ###############################
library(SeuratWrappers)
library(Banksy)

obj <- RunBanksy(obj,
  lambda = 0.8, # Influence of the neighborhood. Larger values yield more spatially coherent domains
  verbose = TRUE,
  assay = "Xenium.500", slot = "data", features = "variable",
  k_geom = 50   # Local neighborhood size. Larger values will yield larger domains
)

DefaultAssay(obj) <- "BANKSY"
obj <- RunPCA(obj, assay = "BANKSY", reduction.name = "pca.banksy", features = rownames(obj), npcs = 30)
obj <- FindNeighbors(obj, reduction = "pca.banksy", dims = 1:20)
obj <- FindClusters(obj, cluster.name = "banksy_cluster", resolution = 0.25)
      # 83 clusters resolution 0.5, ; dims 1:20 res=0.25 32 clusters

Idents(obj) <- "banksy_cluster"
p1 <- DimPlot(obj, reduction = "pca.umap", label = T, raster=T, cols = mypal3)  + ggtitle("banksy_clusters") + NoAxes()+ NoLegend()
p1
ggsave("banksy_clusters_obj.png", plot = p1, device = "png", width = 8, height = 5, units = "in", dpi = 400)


plot_cell_types <- function(data, label) {
  p <- ggplot(data, aes(x = get(label), y = n, fill = predicted.celltype)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_text(aes(label = ifelse(n >= min_count_to_show_label, predicted.celltype, "")), position = position_stack(vjust = 0.5), size = 2) +
    xlab(label) +
    ylab("# of Spots") +
    ggtitle(paste0("Distribution of Cell Types across ", label)) +
    theme_minimal()
}

cell_type_banksy_counts <- obj[[]] %>%
  dplyr::filter(predicted.celltype %in% predicted.celltype) %>%
  dplyr::count(predicted.celltype, banksy_cluster)

min_count_to_show_label <- 20

p <- plot_cell_types(cell_type_banksy_counts, "banksy_cluster")
p
ggsave("obj_cell_type_banksy_counts.png", plot = p, device = "png", width = 10, height = 8, units = "in", dpi = 400)

############################### ############################### ###############################
###############################      Annotate Niches         ###############################
############################### ############################### ###############################
#fov.list <- c("S10963", "S11231", "S12129", "S13164a", "S13164b", "S23b", "S26", "S30578", "S3065", "S31175", "S31542", "S9380")
fov.list <- obj@images
# obj <- BuildNicheAssay(object = obj, fov = "S13164b", group.by = "predicted.celltype")
obj <- BuildNicheAssay(object = obj, fov = fov.list, group.by = "predicted.celltype")



table(obj$predicted.celltype, obj$niches)
t1<- table(obj$predicted.celltype, obj$niches)
write.csv(t1, "obj-predicted.celltype_niches_counts.csv", row.names = T)
table(obj$predicted.celltype, obj$niches)
t1<- table(obj$predicted.celltype, obj$banksy_cluster)
write.csv(t1, "obj-predicted.celltype_banksy_cluster_counts.csv", row.names = T)
t1<- table(obj$predicted.celltype, obj$MNPBurden)
write.csv(t1, "obj-predicted.celltype_MNPBurden_counts.csv", row.names = T)
t1<- table(obj$predicted.celltype, obj$MNP_GA)
write.csv(t1, "obj-predicted.celltype_MNP_GA_counts.csv", row.names = T)
t1<- table(obj$predicted.celltype, obj$Pregnancy_Infections
write.csv(t1, "obj-predicted.celltype_Pregnancy_Infections_counts.csv", row.names = T)

############################### ############################### ###############################
###############################      Save Annotations for Xenium Explorer         ###############################
############################### ############################### ###############################

# save the annotations for uploading into Xenium Explorer
Idents(obj) <- "predicted.celltype"

predicted.celltype = as.data.frame(Idents(obj))
predicted.celltype$cell_id = (rownames(predicted.celltype))
predicted.celltype$group = (predicted.celltype[,1])

predicted.celltype = predicted.celltype[,2:3]
write.csv(predicted.celltype, "obj_predicted.celltype_forXeniumExplorer.csv", quote=FALSE, row.names = F)


orig.ident = as.data.frame(Idents(obj))
orig.ident$cell_id = (rownames(orig.ident))
orig.ident$group = (orig.ident[,1])

orig.ident = orig.ident[,2:3]
write.csv(orig.ident, "obj_orig.ident_forXeniumExplorer.csv", quote=FALSE, row.names = F)

# save the annotations for uploading into Xenium Explorer
Idents(obj) <- "niches"

niches = as.data.frame(Idents(obj))
niches$cell_id = (rownames(niches))
niches$group = (niches[,1])

niches = niches[,2:3]
write.csv(niches, "obj_niches_forXeniumExplorer.csv", quote=FALSE, row.names = F)

# save the annotations for uploading into Xenium Explorer
Idents(obj) <- "seurat_cluster.projected"

seurat_cluster.projected = as.data.frame(Idents(obj))
seurat_cluster.projected$cell_id = (rownames(seurat_cluster.projected))
seurat_cluster.projected$group = (seurat_cluster.projected[,1])

seurat_cluster.projected = seurat_cluster.projected[,2:3]
write.csv(seurat_cluster.projected, "obj_seurat_cluster.projected_forXeniumExplorer.csv", quote=FALSE, row.names = F)

# save the annotations for uploading into Xenium Explorer
Idents(obj) <- "banksy_cluster"

banksy_cluster = as.data.frame(Idents(obj))
banksy_cluster$cell_id = (rownames(banksy_cluster))
banksy_cluster$group = (banksy_cluster[,1])
banksy_cluster = banksy_cluster[,2:3]
write.csv(banksy_cluster, "obj_banksy_cluster_forXeniumExplorer.csv", quote=FALSE, row.names = F)

# saveRDS(obj, "/path/to/local_home/spatial_mnp/results/xenium5k/obj_annot.rds")
rm(list = ls()[!ls() %in% c("obj", "mypal3")])




