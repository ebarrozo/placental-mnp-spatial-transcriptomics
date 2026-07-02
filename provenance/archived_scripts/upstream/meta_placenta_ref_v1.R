## meta_placenta_ref_v1.R
## Enrico Barrozo, Ph.D.
## Assistant Professor, Department of OBGYN, Baylor College of Medicine and Texas Children’s Hospital

## download the files from czi in terminal
## in terminal cd /path/to/local_home/spatial_mnp/scripts
## bash wget_czi_placenta.sh

## Load session information
sessionInfo()

## Set working directory and seed
setwd("/path/to/local_home/iut/results")
set.seed(1)
options(future.globals.maxSize = 300e9) ## 300 gb ram


## Load required libraries
library(dplyr)
library(Seurat)  # Seurat_5.2.1
library(sctransform)
library(ggplot2)
library(cowplot)
library(patchwork)
library(glmGamPoi)
library(BPCells)
library(presto)
library(ggpubr)
library(units)
library(hdf5r)
library(ggrepel)
mypal3 <-get_palette("igv",35)
library(SoupX)
library(HDF5Array) # BiocManager::install(c("HDF5Array"))
library(SingleCellExperiment) #BiocManager::install(c("SingleCellExperiment")
library(zellkonverter) # BiocManager::install(c("zellkonverter"))
library(EnsDb.Hsapiens.v86)    # For gene annotation; convert Ensembl IDs to gene symbols
library(AnnotationDbi)  

setwd("/path/to/local_home/spatial_mnp/results")
dir.create("placental_ref")
setwd("/path/to/local_home/spatial_mnp/results/placental_ref")


########################################################
# Load the .h5ad file into an AnnData object
adata <- readH5AD("/path/to/local_home/spatial_mnp/data/czi/50bbe1a2-5f27-47f5-a809-046459a4ae5e.h5ad")
	## this is a SingleCellExperiments Object

## inspect 
str(adata, max.level = 2)
metadata(adata)
assayNames(adata)
str(assay(adata, "X"))	## 255,901 cells
## counts matrix found in X

library(Matrix)  # Load sparse matrix
counts_matrix <- assay(adata, "X")

# Check a small subset
counts_matrix[1:5, 1:5]

## create seurat object
seurat_obj <- CreateSeuratObject(counts = counts_matrix)

# Check the Seurat object
seurat_obj

## Add metadata
seurat_obj@meta.data <- as.data.frame(colData(adata))

## save seurat object
saveRDS(seurat_obj, file = "Hoo.2024.CellSystems.obj.rds") ## 6.7 GB

########################################################
# Load the .h5ad file into an AnnData object
adata <- readH5AD("/path/to/local_home/spatial_mnp/data/czi/b4e813a4-e9e0-4735-9c4b-fac499ac1696.h5ad")
	## this is a SingleCellExperiments Object

## inspect 
str(adata, max.level = 2)
metadata(adata)
assayNames(adata)
str(assay(adata, "X"))	## 255,901 cells
## counts matrix found in X

library(Matrix)  # Load sparse matrix
counts_matrix <- assay(adata, "X")

# Check a small subset
counts_matrix[1:5, 1:5]

## create seurat object
seurat_obj <- CreateSeuratObject(counts = counts_matrix)

# Check the Seurat object
seurat_obj

## Add metadata
seurat_obj@meta.data <- as.data.frame(colData(adata))

## save seurat object
saveRDS(seurat_obj, file = "Arutyunyan.2023.Nature.1.obj.rds") ## 6.7 GB

########################################################
# Load the .h5ad file into an AnnData object
adata <- readH5AD("/path/to/local_home/spatial_mnp/data/czi/b36f6787-e0fe-429c-a024-5a8d6915e781.h5ad")
	## this is a SingleCellExperiments Object

## inspect 
str(adata, max.level = 2)
metadata(adata)
assayNames(adata)
str(assay(adata, "X"))	## 255,901 cells
## counts matrix found in X

library(Matrix)  # Load sparse matrix
counts_matrix <- assay(adata, "X")

# Check a small subset
counts_matrix[1:5, 1:5]

## create seurat object
seurat_obj <- CreateSeuratObject(counts = counts_matrix)

# Check the Seurat object
seurat_obj

## Add metadata
seurat_obj@meta.data <- as.data.frame(colData(adata))

## save seurat object
saveRDS(seurat_obj, file = "Arutyunyan.2023.Nature.2.obj.rds") ## 6.7 GB


########################################################
# Load the .h5ad file into an AnnData object
adata <- readH5AD("/path/to/local_home/spatial_mnp/data/czi/1d0c5d7d-5450-4387-8063-c665e4bed731.h5ad")
	## this is a SingleCellExperiments Object

## inspect 
str(adata, max.level = 2)
metadata(adata)
assayNames(adata)
str(assay(adata, "X"))	## 255,901 cells
## counts matrix found in X

library(Matrix)  # Load sparse matrix
counts_matrix <- assay(adata, "X")

# Check a small subset
counts_matrix[1:5, 1:5]

## create seurat object
seurat_obj <- CreateSeuratObject(counts = counts_matrix)

# Check the Seurat object
seurat_obj

## Add metadata
seurat_obj@meta.data <- as.data.frame(colData(adata))

## save seurat object
saveRDS(seurat_obj, file = "Arutyunyan.2023.Nature.3.obj.rds") ## 6.7 GB


########################################################
# Load the .h5ad file into an AnnData object
rm(adata)
adata <- readH5AD("/path/to/local_home/spatial_mnp/data/czi/5b8941a9-fef8-4348-876c-c333a87c0c3f.h5ad")
	## this is a SingleCellExperiments Object

## inspect 
str(adata, max.level = 2)
metadata(adata)
assayNames(adata)
str(assay(adata, "X"))	## 255,901 cells
## counts matrix found in X

library(Matrix)  # Load sparse matrix
counts_matrix <- assay(adata, "X")

# Check a small subset
counts_matrix[1:5, 1:5]

## create seurat object
seurat_obj <- CreateSeuratObject(counts = counts_matrix)

# Check the Seurat object
seurat_obj

## Add metadata
seurat_obj@meta.data <- as.data.frame(colData(adata))

## save seurat object
saveRDS(seurat_obj, file = "Cao.2020.Science.1.obj.rds") ## 6.7 GB

########################################################
# Load the .h5ad file into an AnnData object
rm(adata)
adata <- readH5AD("/path/to/local_home/spatial_mnp/data/czi/35476cf7-35f1-4daf-8507-734403d09c27.h5ad")
	## this is a SingleCellExperiments Object

## inspect 
str(adata, max.level = 2)
metadata(adata)
assayNames(adata)
str(assay(adata, "X"))	## 255,901 cells
## counts matrix found in X

library(Matrix)  # Load sparse matrix
counts_matrix <- assay(adata, "X")

# Check a small subset
counts_matrix[1:5, 1:5]

## create seurat object
seurat_obj <- CreateSeuratObject(counts = counts_matrix)

# Check the Seurat object
seurat_obj

## Add metadata
seurat_obj@meta.data <- as.data.frame(colData(adata))

## save seurat object
saveRDS(seurat_obj, file = "Cao.2020.Science.2.obj.rds") ## 6.7 GB

########################################################
# Load the .h5ad file into an AnnData object
rm(adata)
adata <- readH5AD("/path/to/local_home/spatial_mnp/data/czi/ae0c62a1-a30f-4033-97d2-0edb2e146c53.h5ad")
	## this is a SingleCellExperiments Object

## inspect 
str(adata, max.level = 2)
metadata(adata)
assayNames(adata)
str(assay(adata, "X"))	## 255,901 cells
## counts matrix found in X

library(Matrix)  # Load sparse matrix
counts_matrix <- assay(adata, "X")

# Check a small subset
counts_matrix[1:5, 1:5]

## create seurat object
seurat_obj <- CreateSeuratObject(counts = counts_matrix)

# Check the Seurat object
seurat_obj

## Add metadata
seurat_obj@meta.data <- as.data.frame(colData(adata))

## save seurat object
saveRDS(seurat_obj, file = "Han.2020.Nature.obj.rds") ## 6.7 GB


########################################################
# Load the .h5ad file into an AnnData object
adata <- readH5AD("/path/to/local_home/spatial_mnp/data/czi/01b1a201-b5cf-4ffd-8e3b-284aeb5f11c5.h5ad")
	## this is a SingleCellExperiments Object

## inspect 
str(adata, max.level = 2)
metadata(adata)
assayNames(adata)
str(assay(adata, "X"))	## 255,901 cells
## counts matrix found in X

library(Matrix)  # Load sparse matrix
counts_matrix <- assay(adata, "X")

# Check a small subset
counts_matrix[1:5, 1:5]

## create seurat object
seurat_obj <- CreateSeuratObject(counts = counts_matrix)

# Check the Seurat object
seurat_obj

## Add metadata
seurat_obj@meta.data <- as.data.frame(colData(adata))

## save seurat object
saveRDS(seurat_obj, file = "VentoTormo.2018.Nature.obj.rds") ## 6.7 GB



########################################################

## Merge Seurat objects and convert Ensembl gene IDs to gene symbols
# Print session information
sessionInfo()

# Set working directory (adjust if needed)
setwd("/path/to/local_home/spatial_mnp/results/placental_ref")
set.seed(1)

# -------------------------------------------------------------------
# Install necessary libraries if not already installed:
#
# For CRAN packages:
# install.packages("Seurat")
#
# For Bioconductor packages:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("EnsDb.Hsapiens.v86")
# BiocManager::install("AnnotationDbi")
# -------------------------------------------------------------------

library(Seurat)               # For Seurat objects and merging
library(EnsDb.Hsapiens.v86)    # For gene annotation; convert Ensembl IDs to gene symbols
library(AnnotationDbi)         # For mapping IDs


# List of saved Seurat object filenames
#seurat_files <- c("Hoo.2024.CellSystems.obj.rds",
 #                 "Arutyunyan.2023.Nature.1.obj.rds",
 #                 "Arutyunyan.2023.Nature.2.obj.rds",
 #                 "Arutyunyan.2023.Nature.3.obj.rds",
 #                 "Cao.2020.Science.1.obj.rds",
 #                 "Cao.2020.Science.2.obj.rds",
 #                 "Han.2020.Nature.obj.rds",
 #                 "VentoTormo.2018.obj.rds")


# Load the Seurat objects into a list
#seurat_list <- lapply(seurat_files, readRDS) ## 63 gb

## clean environment, removing everything but the merged object
# rm(list = ls()[!ls() %in% c("seurat_list")])

## make sure each object has layers joined before merging
Hoo.obj <- readRDS("/path/to/local_home/spatial_mnp/results/placental_ref/Hoo.2024.CellSystems.obj.rds") ## this is 8.7 gb
Arutyunyan.1.obj <- readRDS("/path/to/local_home/spatial_mnp/results/placental_ref/Arutyunyan.2023.Nature.1.obj.rds") ## this is 2.7 gb
#Arutyunyan.2.obj <- readRDS("/path/to/local_home/spatial_mnp/results/placental_ref/Arutyunyan.2023.Nature.2.obj.rds") ## this is 9.1 gb
#Arutyunyan.3.obj <- readRDS("/path/to/local_home/spatial_mnp/results/placental_ref/Arutyunyan.2023.Nature.3.obj.rds") ## this is 1.2 gb
Cao.obj <- readRDS("/path/to/local_home/spatial_mnp/results/placental_ref/Cao.2020.Science.2.obj.rds") ## this is 7.8 gb
Han.obj <- readRDS("/path/to/local_home/spatial_mnp/results/placental_ref/Han.2020.Nature.obj.rds") ## this is 3.8 gb
#Vento.obj <- readRDS("/path/to/local_home/spatial_mnp/results/placental_ref/VentoTormo.2018.Nature.obj.rds") ## this is 2.1 gb

## running at 202 gb

Idents(Hoo.obj) <- "cell_type"
Hoo.obj <- JoinLayers(Hoo.obj)
gc()

Arutyunyan.3.obj
Idents(Arutyunyan.3.obj) <- "cell_type"
Arutyunyan.3.obj <- JoinLayers(Arutyunyan.3.obj)
gc()
Idents(Arutyunyan.2.obj) <- "cell_type"
Arutyunyan.2.obj <- JoinLayers(Arutyunyan.2.obj)
rm(Arutyunyan.2.obj)
rm(Arutyunyan.3.obj)

## use Arutyunyan.1
Idents(Arutyunyan.1.obj) <- "cell_type"
Idents(Arutyunyan.1.obj) <- "dataset"
levels(Arutyunyan.1.obj)
Idents(Arutyunyan.1.obj) <- "development_stage"
levels(Arutyunyan.1.obj)

Arutyunyan.1.obj <- JoinLayers(Arutyunyan.1.obj)
gc()
## includes vento-tormo 2018
rm(Vento.obj)

ncol(Arutyunyan.1.obj) # 325665

Idents(Cao.obj) <- "disease"
Cao.obj <- subset(x = Cao.obj, idents = c("normal"))
Idents(Cao.obj) <- "tissue"
levels(Cao.obj)
#[1] "intestine"     "stomach"       "heart"         "eye"          
# [5] "pancreas"      "muscle organ"  "telencephalon" "placenta"     
# [9] "cerebellum"    "lung"          "spleen"        "liver"        
#[13] "kidney"        "adrenal gland" "thymus"
Cao.obj <- subset(x = Cao.obj, idents = c("placenta"))
Idents(Cao.obj) <- "cell_type"
Cao.obj <- JoinLayers(Cao.obj)
gc()
ncol(Cao.obj) # 7374


Idents(Han.obj) <- "tissue_type"
Han.obj <- subset(x = Han.obj, idents = c("tissue"))
Idents(Han.obj) <- "development_stage"
levels(Han.obj)
Han.obj <- subset(x = Han.obj, idents = c("10th week post-fertilization stage", "11th week post-fertilization stage", "12th week post-fertilization stage", "13th week post-fertilization stage", "14th week post-fertilization stage", "12th week post-fertilization stage"))
Idents(Han.obj) <- "cell_type"
Han.obj <- JoinLayers(Han.obj)
gc()
ncol(Han.obj) # 188673




# Merge the objects into one Seurat object
merged_obj <- merge(x = Hoo.obj, y = c(Arutyunyan.1.obj, Cao.obj, Han.obj), merge.data = TRUE, project = "czi_ref")
ncol(merged_obj) # 680690 19.4 GB , 4 counts matrices, 52 cell_types
merged_obj <- JoinLayers(merged_obj) ## 1 count matrix- however R doubled in size while doing it so be careful



## clean environment, removing everything but the merged object
rm(list = ls()[!ls() %in% c("merged_obj")])
saveRDS(merged_obj, file = "placental_ref.obj.rds") ## 62.5 gb GB

# merged_obj <- readRDS("/path/to/local_home/iut/results/seurat_analysis_v1/placental_ref/placental_ref.obj.rds") ## this is 62.5 gb
Idents(merged_obj) <- "cell_type"  
cell.types.counts <- table(Idents(merged_obj))
write.table(cell.types.counts, "cell.types_counts_merged_obj.txt", sep="\t")
Idents(merged_obj) <- "development_stage"  
cell.types.counts <- table(Idents(merged_obj))
write.table(cell.types.counts, "development_stage_counts_merged_obj.txt", sep="\t")
Idents(merged_obj) <- "donor_id"  
cell.types.counts <- table(Idents(merged_obj))
write.table(cell.types.counts, "donor_id_counts_merged_obj.txt", sep="\t")
Idents(merged_obj) <- "tissue_type"  
cell.types.counts <- table(Idents(merged_obj))
write.table(cell.types.counts, "tissue_type_counts_merged_obj.txt", sep="\t")
Idents(merged_obj) <- "disease"  
cell.types.counts <- table(Idents(merged_obj))
write.table(cell.types.counts, "disease_counts_merged_obj.txt", sep="\t")
Idents(merged_obj) <- "cell_type"  
rm(list = ls()[!ls() %in% c("merged_obj")])


# merged_obj@assays[["RNA"]]@layers[["counts.1"]]
# Convert Ensembl gene IDs to gene symbols
gene_ids <- rownames(merged_obj)  #38,168 genes

# Map Ensembl IDs to gene symbols using EnsDb.Hsapiens.v86
gene_symbols <- mapIds(EnsDb.Hsapiens.v86,
                       keys = gene_ids,
                       column = "SYMBOL",
                       keytype = "GENEID",
                       multiVals = "first")
      #Warning message:
      #Unable to map 2661 of 54471 requested IDs. 

# For genes that fail conversion (NA), retain the original Ensembl ID
new_gene_names <- ifelse(is.na(gene_symbols), gene_ids, gene_symbols)

## @RNA@layers [["counts.1"]] through counts.8 for each object. need to be joined for tru merge


# merged_obj <- JoinLayers(merged_obj)  ## 290 gb before crashing R

# Update gene names in the default RNA assay
DefaultAssay(merged_obj) <- "RNA"
rownames(merged_obj[["RNA"]]) <- new_gene_names
if (!is.null(merged_obj[["RNA"]])) {
  rownames(merged_obj[["RNA"]]) <- new_gene_names
}

# Also update the Seurat object rownames (they are typically synced with the default assay)
rownames(merged_obj) <- new_gene_names

# Check the updated gene names (first few entries)
head(rownames(merged_obj))


#for (layer in layer_names) {
    # Get the number of genes in the current layer
    num_genes_layer <- nrow(merged_obj@assays[["RNA"]]@layers[[layer]])
    
    # Ensure `new_gene_names` has the correct length for this specific layer
    if (length(new_gene_names) >= num_genes_layer) {
        # Subset new_gene_names to match the number of genes in this layer
        layer_gene_names <- new_gene_names[1:num_genes_layer] 
        
        # Assign updated gene names
        rownames(merged_obj@assays[["RNA"]]@layers[[layer]]) <- layer_gene_names
        
        message(paste("Updated gene names for layer:", layer, "with", num_genes_layer, "genes."))
    } else {
        warning(paste("Skipping layer:", layer, "- Gene name mismatch (Expected:", num_genes_layer, ", Got:", length(new_gene_names), ")"))
    }
}


# lapply(layer_names, function(layer) nrow(merged_obj@assays[["RNA"]]@layers[[layer]]))



## save seurat object
saveRDS(merged_obj, file = "placental_ref_named.obj.rds") ## 19.4 GB




# merged_obj <- readRDS("/path/to/local_home/spatial_mnp/results/seurat_analysis_v1/placental_ref/placental_ref.obj.rds") ## this is 62.5 gb
seurat.object <- readRDS("/path/to/local_home/spatial_mnp/results/xenium5k/reference_placenta.rds") ## this is 27 gb
ncol(seurat.object)	#  229186

DefaultAssay(seurat.object) <- "RNA"
ncol(seurat.object)	#273944


Idents(seurat.object) <- "Type"s
seurat.object$cell_type <- Idents(seurat.object)
Idents(seurat.object) <- "cell_type"

seurat.object@images <- list() # from 27gb to 12.9 gb
seurat.object@dimreducs <- list() # from 27gb to 12.9 gb



library(Matrix)

# Get count matrices from each Seurat object
counts_merged <- GetAssayData(merged_obj, layer = "counts")
counts_seurat <- GetAssayData(seurat.object, layer = "counts")
# Create a complete gene list (meta.genes) from both datasets
meta.genes <- union(rownames(counts_merged), rownames(counts_seurat))
missing_genes_merged <- setdiff(meta.genes, rownames(counts_merged))

# Create a zero-filled sparse matrix for missing genes
zero_matrix_merged <- Matrix(0, nrow = length(missing_genes_merged), ncol = ncol(counts_merged), sparse = TRUE)
rownames(zero_matrix_merged) <- missing_genes_merged
colnames(zero_matrix_merged) <- colnames(counts_merged)

# Append zero-filled genes
counts_merged_expanded <- rbind(counts_merged, zero_matrix_merged)
missing_genes_seurat <- setdiff(meta.genes, rownames(counts_seurat))

# Create a zero-filled sparse matrix for missing genes
zero_matrix_seurat <- Matrix(0, nrow = length(missing_genes_seurat), ncol = ncol(counts_seurat), sparse = TRUE)
rownames(zero_matrix_seurat) <- missing_genes_seurat
colnames(zero_matrix_seurat) <- colnames(counts_seurat)

# Append zero-filled genes
counts_seurat_expanded <- rbind(counts_seurat, zero_matrix_seurat)
# Ensure row order matches before merging
counts_merged_expanded <- counts_merged_expanded[meta.genes, ]
counts_seurat_expanded <- counts_seurat_expanded[meta.genes, ]

# Merge the count matrices by columns (cells)
merged_counts <- cbind(counts_merged_expanded, counts_seurat_expanded)


rownames(merged_counts) <- gsub("_", "-", rownames(merged_counts))
rownames(merged_counts) <- make.unique(rownames(merged_counts))

merged_ref <- CreateSeuratObject(counts = merged_counts, project = "czi_barrozo_ref")

meta_merged <- merged_obj@meta.data
meta_seurat <- seurat.object@meta.data

# Columns missing in meta_merged
missing_in_merged <- setdiff(colnames(meta_seurat), colnames(meta_merged))

# Columns missing in meta_seurat
missing_in_seurat <- setdiff(colnames(meta_merged), colnames(meta_seurat))

print(missing_in_merged)  # Columns present in meta_seurat but not in meta_merged
print(missing_in_seurat)  # Columns present in meta_merged but not in meta_seurat
# Add missing columns to meta_merged
for (col in missing_in_merged) {
  meta_merged[[col]] <- NA
}

# Add missing columns to meta_seurat
for (col in missing_in_seurat) {
  meta_seurat[[col]] <- NA
}
all(colnames(meta_merged) == colnames(meta_seurat))  # Should return TRUE
meta_merged <- meta_merged[, order(colnames(meta_merged))]
meta_seurat <- meta_seurat[, order(colnames(meta_seurat))]
identical(colnames(meta_merged), colnames(meta_seurat))  # Should return TRUE
meta_combined <- rbind(meta_merged, meta_seurat)
rownames(meta_combined) <- colnames(merged_ref)  # Ensure row names match cell barcodes
merged_ref@meta.data <- meta_combined
dim(merged_ref@meta.data)  # Should match the number of cells in merged_ref
head(merged_ref@meta.data)  # Confirm metadata looks correct



Idents(merged_ref) <- "cell_type"
levels(merged_ref)


# Define a named vector for grouping similar cell types
cell_type_map <- c(
  # Fibroblasts
  "fibroblast" = "Fibroblast",
  "stromal cell" = "Fibroblast",
  "embryonic fibroblast" = "Fibroblast",
  "mesenchymal stem cell" = "Fibroblast",
  "Stromal" = "Fibroblast",

  # Trophoblasts
  "placental villous trophoblast" = "Trophoblast",
  "extravillous trophoblast" = "Trophoblast",
  "chorionic trophoblast cell" = "Trophoblast",
  "syncytiotrophoblast cell" = "Trophoblast",
  "trophoblast giant cell" = "Trophoblast",
  "VCT" = "Trophoblast",
  "EVT" = "Trophoblast",
  "SYT" = "Trophoblast",

  # Endothelial Cells
  "endothelial cell" = "Endothelial",
  "endothelial cell of placenta" = "Endothelial",
  "endothelial cell of uterus" = "Endothelial",
  "endothelial cell of lymphatic vessel" = "Endothelial",
  "endothelial cell of vascular tree" = "Endothelial",
  "VEC" = "Endothelial",

  # Macrophages
  "Hofbauer cell" = "Macrophage",
  "macrophage" = "Macrophage",
  "alternatively activated macrophage" = "Macrophage",
  "myeloid cell" = "Macrophage",
  "monocyte" = "Macrophage",
  "Macrophage" = "Macrophage",

  # Dendritic Cells
  "dendritic cell" = "Dendritic Cell",
  "DC" = "Dendritic Cell",
  "DC2" = "Dendritic Cell",

  # Natural Killer (NK) Cells
  "natural killer cell" = "NK Cell",
  "decidual natural killer cell, human" = "NK Cell",
  "NK" = "NK Cell",

  # T Cells
  "T cell" = "T Cell",
  "CD8T" = "T Cell",
  "T" = "T Cell",
  "regulatory T cell" = "T Cell",
  "group 3 innate lymphoid cell, human" = "T Cell",
  "innate lymphoid cell" = "T Cell",

  # B Cells
  "B cell" = "B Cell",
  "B" = "B Cell",
  "plasma cell" = "B Cell",

  # Smooth Muscle
  "smooth muscle cell" = "Smooth Muscle",
  "uterine smooth muscle cell" = "Smooth Muscle",
  "SmoothMuscle" = "Smooth Muscle",
  "PVC" = "Smooth Muscle",

  # Neutrophils & Granulocytes
  "granulocyte" = "Neutrophil",
  "neutrophil" = "Neutrophil",

  # Other Common Categories
  "erythroid lineage cell" = "Erythroid",
  "erythroid progenitor cell" = "Erythroid",
  "megakaryocyte" = "Erythroid",

  "unknown" = "Unknown"
)

# Set cell type as identity
Idents(merged_ref) <- "cell_type"

# Replace cell types based on the mapping
merged_ref$cell_type_merged <- recode(merged_ref$cell_type, !!!cell_type_map)

# Update Idents to reflect the new grouping
Idents(merged_ref) <- merged_ref$cell_type

merged_ref$cell_type_l1 <- Idents(merged_ref)

# Define Level 2 Cell Type Annotation with Updates
cell_type_l2_map <- c(
  # Trophoblasts
  "EVT" = "Extravillous Trophoblast",
  "VCT" = "Villous Cytotrophoblast",
  "SYT" = "Syncytiotrophoblast",
  "placental villous trophoblast" = "Villous Cytotrophoblast",
  "extravillous trophoblast" = "Extravillous Trophoblast",
  "syncytiotrophoblast cell" = "Syncytiotrophoblast",
  "chorionic trophoblast cell" = "Chorionic Trophoblast",
  "trophoblast giant cell" = "Trophoblast Giant Cell",

  # Macrophages
  "Hofbauer cell" = "Hofbauer Cell",
  "macrophage" = "Macrophage",
  "alternatively activated macrophage" = "M2 Macrophage",
  "myeloid cell" = "Macrophage",
  "monocyte" = "Macrophage",
  "Macrophage" = "Macrophage",

  # NK Cells
  "decidual natural killer cell, human" = "dNK",
  "natural killer cell" = "NK",
  "NK" = "NK",

  # T Cells
  "T cell" = "T Cell",
  "CD8T" = "T Cell",
  "T" = "T Cell",
  "regulatory T cell" = "Treg",
  "group 3 innate lymphoid cell, human" = "ILC3",
  "innate lymphoid cell" = "Innate Lymphoid Cell",

  # B Cells
  "B cell" = "B Cell",
  "B" = "B Cell",
  "plasma cell" = "Plasma Cell",

  # Dendritic Cells
  "dendritic cell" = "DC",
  "DC" = "DC",
  "DC2" = "DC2",

  # Endothelial Cells
  "endothelial cell" = "Endothelial Cell",
  "endothelial cell of placenta" = "Placental Endothelial",
  "endothelial cell of uterus" = "Uterine Endothelial",
  "endothelial cell of lymphatic vessel" = "Lymphatic Endothelial",
  "endothelial cell of vascular tree" = "Vascular Endothelial",
  "VEC" = "Endothelial Cell",

  # Fibroblasts
  "fibroblast" = "Fibroblast",
  "stromal cell" = "Fibroblast",
  "embryonic fibroblast" = "Fibroblast",
  "mesenchymal stem cell" = "Fibroblast",
  "Stromal" = "Fibroblast",

  # Smooth Muscle
  "smooth muscle cell" = "Smooth Muscle",
  "uterine smooth muscle cell" = "Uterine Smooth Muscle",
  "SmoothMuscle" = "Smooth Muscle",
  "PVC" = "Smooth Muscle",

  # Neutrophils & Granulocytes
  "granulocyte" = "Neutrophil",
  "neutrophil" = "Neutrophil",

  # Other Specific Cells
  "erythroid lineage cell" = "Erythroid",
  "erythroid progenitor cell" = "Erythroid",
  "megakaryocyte" = "Erythroid",

  # Unknown
  "unknown" = "Unknown"
)
# Set Level 1 cell type as identity first
Idents(merged_ref) <- "cell_type_l1"

# Assign Level 2 cell types
merged_ref$cell_type_l2 <- recode(merged_ref$cell_type_l1, !!!cell_type_l2_map)

# Update Idents to reflect the new Level 2 annotation
Idents(merged_ref) <- merged_ref$cell_type_l2
levels(merged_ref)  # Should reflect refined cell types
Idents(merged_ref) <- "cell_type_l2"

table(Idents(merged_ref))  # Cell count per refined type


Idents(merged_ref) <- "cell_type_l2"
final_placenta_decidua_cells <- c(
  # **Trophoblast Lineages**
  "Villous Cytotrophoblast",
  "Extravillous Trophoblast",
  "Syncytiotrophoblast",
  "Trophoblast Giant Cell",
  "Chorionic Trophoblast",

  # **Endothelial & Vascular Cells**
  "Placental Endothelial",
  "Vascular Endothelial",
  "perivascular cell",  # Lowercase fix

  # **Fibroblasts & Stromal Cells**
  "Fibroblast",

  # **Macrophages**
  "Hofbauer Cell",
  "Macrophage",

  # **Natural Killer Cells**
  "dNK",
  "NK",

  # **Dendritic Cells**
  "DC",
  "DC2",

  # **T Cells**
  "T Cell",
  "Treg",
  "ILC3",

  # **B Cells**
  "B Cell",

  # **Plasma Cells**
  "Plasma Cell"
)


merged_ref_final <- subset(merged_ref, idents = final_placenta_decidua_cells)

levels(merged_ref_final)  # Should list only placenta + decidua cell types
t1 <- table(Idents(merged_ref_final))  # Confirm cell counts
write.table(t1, "final_cell_types_l2.txt", sep="\t")

ncol(merged_ref_final) # 743247
## 20 cell types
merged_ref_final <- JoinLayers(merged_ref_final)


Idents(merged_ref_final) <- "cell_type_l2"

all.genes <- rownames(merged_ref_final)
	## ############# After running these I had counts and data in layers
merged_ref_final <- NormalizeData(merged_ref_final)
DefaultAssay(merged_ref_final) <- "RNA"
merged_ref_final <- FindVariableFeatures(merged_ref_final)

rm(list = ls()[!ls() %in% c("merged_ref_final")])
gc()
#merged_ref_final <- JoinLayers(merged_ref_final)
merged_ref_final2 <- DietSeurat(merged_ref_final, layers="counts")
merged_ref_final2@commands <- list() # from 27gb to 12.9 gb

#saveRDS(merged_ref_final2, file = "merged_ref_final_annot.rds") ##  23.9GB GB
rm(merged_ref_final2)

DefaultAssay(merged_ref_final) <- "RNA"



merged_ref_final <- readRDS("/path/to/local_home/spatial_mnp/results/placental_ref/merged_ref_final_annot.rds") ## this is 8.7 gb
DefaultAssay(merged_ref_final) <- "RNA"
merged_ref_final <- NormalizeData(merged_ref_final)
merged_ref_final <- FindVariableFeatures(merged_ref_final)
merged_ref_final <- ScaleData(merged_ref_final)
# we select 50,0000 cells and create a new 'sketch' assay
merged_ref_final <- SketchData(
  object = merged_ref_final,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)
# switch analysis to sketched cells
DefaultAssay(merged_ref_final) <- "sketch"

# perform clustering workflow
merged_ref_final <- ScaleData(merged_ref_final)
merged_ref_final <- RunPCA(merged_ref_final, assay = "sketch", reduction.name = "pca.sketch")
merged_ref_final <- FindNeighbors(merged_ref_final, assay = "sketch", reduction = "pca.sketch", dims = 1:50)
merged_ref_final <- FindClusters(merged_ref_final, cluster.name = "seurat_cluster.sketched", resolution = 0.6)
merged_ref_final <- RunUMAP(merged_ref_final, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50)
# Projected UMAP embeddings will be stored in merged_ref_final[["umap.sketch"]]
merged_ref_final <- ProjectData(
  object = merged_ref_final,
  assay = "RNA",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)
# We can visualize the clustering results for the sketched cells, as well as the projected clustering results for the full dataset:
DefaultAssay(merged_ref_final) <- "sketch"
Idents(merged_ref_final) <- "seurat_cluster.sketched"
mypal3 <-get_palette("igv",40)
p1 <- DimPlot(merged_ref_final, reduction = "umap.sketch", label = T, raster=T, cols = mypal3)  + ggtitle("Sketched clustering (50,000 cells)") + NoAxes()+ NoLegend()
p1
# switch to full dataset
DefaultAssay(merged_ref_final) <- "RNA"
Idents(merged_ref_final) <- "seurat_cluster.projected"
p2 <- DimPlot(merged_ref_final, reduction = "full.umap.sketch", label = T,raster=T, cols = mypal3) + ggtitle("Projected clustering (full dataset)") + NoAxes() + NoLegend()# + theme(legend.position = "bottom")
p2
p1 | p2
p3 <- p1 | p2
ggsave("seurat_cluster_projections_merged_ref_final.png", plot = p3, device = "png", width = 12, height = 5, units = "in", dpi = 400)

p2 <- DimPlot(merged_ref_final, reduction = "full.umap.sketch", label = T,raster=T, cols = mypal3) + ggtitle("Projected clustering (full dataset)") + NoAxes() + NoLegend()# + theme(legend.position = "bottom")
p1 <- DimPlot(merged_ref_final, reduction = "full.umap.sketch", group.by="cell_type_l2", label = F, raster=T, cols = mypal3)  + ggtitle("Cell Types") + NoAxes()

p1 | p2
p3 <-  p2|p1 
ggsave("seurat_cluster_projections_cell_type_l2_merged_ref_final.png", plot = p3, device = "png", width = 12, height = 5, units = "in", dpi = 400)








merged_ref_final <- CellCycleScoring(merged_ref_final, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = FALSE)
merged_ref_final$CC.Difference <- merged_ref_final$S.Score - merged_ref_final$G2M.Score
merged_ref_final <- PercentageFeatureSet(merged_ref_final, pattern = "^MT-", col.name = "percent.mt")
merged_ref_final <- PercentageFeatureSet(merged_ref_final, pattern = "^RP[SL]", col.name = "percent.ribo")
merged_ref_final <- ScaleData(merged_ref_final, vars.to.regress = c("CC.Difference", "percent.mt", "percent.ribo"))

DefaultAssay(merged_ref_final) <- "RNA"
merged_ref_final <- FindVariableFeatures(merged_ref_final, nfeatures=3000, verbose=TRUE)

top.3k <- merged_ref_final@assays[["RNA"]]@meta.data[["var.features"]]
# View(top.3k)
top.3k <- top.3k[!is.na(top.3k)]

merged_ref_final<- RunPCA(merged_ref_final, features = top.3k)
merged_ref_final <- FindNeighbors(merged_ref_final, features = "top.3k", dims = 1:30)
merged_ref_final <- FindClusters(merged_ref_final, resolution = 0.6)
merged_ref_final <- RunUMAP(merged_ref_final, reduction = "pca", dims = 1:30, reduction.name = "pca.umap.var.reg")
merged_ref_final <- JoinLayers(merged_ref_final)

rm(list = ls()[!ls() %in% c("merged_ref_final", "top.3k", "mypal3")])
gc()
#merged_ref_final <- IntegrateLayers(
#  object = merged_ref_final, method = RPCAIntegration, assay="RNA",
#  orig.reduction = "pca", new.reduction = "integrated.rpca",
#  verbose = TRUE
#)
#merged_ref_final <- FindNeighbors(merged_ref_final, assay = "sketch", reduction = "integrated.rpca", dims = 1:50)
#merged_ref_final <- FindClusters(merged_ref_final, cluster.name = "seurat_cluster.sketched", resolution = 0.6)
#merged_ref_final <- RunUMAP(merged_ref_final, reduction = "integrated.rpca", dims = 1:30, reduction.name = "integrated.rpca")

	################# There are no layers to integrate #################

p2<-DimPlot(merged_ref_final, reduction = "pca.umap.var.reg", group.by = "cell_type_l2", alpha = 0.1, label = FALSE,raster=T, cols = mypal3) + ggtitle("Reference (743,247 Cells)")  + NoAxes()
p2
ggsave("ref_annot_pca.umap.var.reg.png", plot = p2, device = "png", width = 8, height = 5, units = "in", dpi = 400)



Idents(merged_ref_final) <- "orig.ident"  
cell.types.counts <- table(Idents(merged_ref_final))
write.table(cell.types.counts, "orig.ident_counts_merged_ref_final.txt", sep="\t")


################ This is the final UMAP $$$$$$$$$$$$$$$$$$





#p4<-DimPlot(merged_ref_final, reduction = "integrated.rpca", group.by = "cell_type_l2", alpha = 0.1, label = FALSE,raster=T, cols = mypal3) + ggtitle("UMAP (RCPA Integrated") + NoAxes()
#p4
#p3 <-  p2|p1 
#ggsave("ref_annot_umapvsrcpa.png", plot = p5, device = "png", width = 12, height = 5,  units = "in", dpi = 400)

merged_ref_final <- DietSeurat(merged_ref_final, layers="counts")

# saveRDS(merged_ref_final, file = "merged_ref_final_UMAP.rds") ## 65.4 GB
merged_ref_final <- readRDS("/path/to/local_home/spatial_mnp/results/placental_ref/merged_ref_final_UMAP.rds") ## this is 8.7 gb


# Create downsampled merged_ref_final to make visualization easier
DefaultAssay(merged_ref_final) <- "RNA"
Idents(merged_ref_final) <- "cell_type_l2"
merged_ref_final_subset <- subset(merged_ref_final, cells = Cells(merged_ref_final[["RNA"]]), downsample = 50000)
rm(merged_ref_final)
gc()

# Order clusters by similarity
DefaultAssay(merged_ref_final_subset) <- "RNA"
Idents(merged_ref_final_subset) <- "cell_type_l2"
merged_ref_final_subset <- BuildClusterTree(merged_ref_final_subset, assay = "RNA", reduction = "pca", reorder = T)
p3 <- PlotClusterTree(merged_ref_final_subset, direction = "downwards")
ggsave("tree_merged_ref_final.png", plot = p3, device = "png", width = 8, height = 11, units = "in", dpi = 400)

Idents(merged_ref_final_subset) <- "cell_type_l2"
# We can also find and visualize the top gene expression markers for each cluster:
de_markers <- FindAllMarkers(merged_ref_final_subset, assay = "RNA", only.pos = TRUE, min.pct = 0.20, logfc.threshold = 1) # Detects FC >2 in 20% of cells in cluster; for marker genes
write.table(de_markers, "merged_ref_final_de_cell_type_l2_subset_2FC_20perccells.txt", sep="\t")
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene

# Define pattern for lncRNAs, uncharacterized genes, and Ensembl IDs
lncRNA_pattern <- "^RP[0-9]+-|^AC[0-9]+.|^AL[0-9]+.|^LINC|^MIR|^SNHG|^HCG|^CT[CDE]-"  # Updated to remove CTC, CTD genes

# Remove genes matching lncRNA pattern and Ensembl IDs
filtered_top2 <- top2[!grepl(lncRNA_pattern, top2) & !grepl("^ENSG", top2)]

# Ensure unique, clean gene list
unique.top2_filtered <- unique(filtered_top2)


# Generate the DotPlot with enhanced aesthetics
feature.plot <- DotPlot(merged_ref_final_subset, features = unique.top2_filtered) +
  theme_minimal() +  # Clean, minimal theme
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12, face = "bold", color = "black"),  # Black, readable x-axis labels
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),  # Black, bold y-axis labels
    legend.text = element_text(size = 10, color = "black"),  # Black legend text
    legend.title = element_text(size = 12, face = "bold", color = "black")  # Black, bold legend title
  ) +
  scale_x_discrete(name = "Top Marker Genes") +  # Rename x-axis
  scale_y_discrete(name = "Cell Type") +  # Rename y-axis
  scale_color_gradient(low = "blue", high = "red") +  # Adjust color scheme for better contrast
  scale_size(range = c(1, 8)) +  # Increase dot size range for better differentiation
  guides(color = guide_colorbar(title = "Expression Level"), size = guide_legend(title = "Pct Expressed"))

# Save high-resolution plot with black labels
png(file = paste0("de_cell_type_l2_subset_2FC_50perccells-DotPlot.png"), res = 400, width = 3000, height = 2000)
print(feature.plot)  # Print to PNG file
dev.off()

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene

# Define pattern for lncRNAs, uncharacterized genes, and Ensembl IDs
lncRNA_pattern <- "^RP[0-9]+-|^AC[0-9]+.|^AL[0-9]+.|^LINC|^MIR|^SNHG|^HCG|^CT[CDE]-"  # Updated to remove CTC, CTD genes

# Remove genes matching lncRNA pattern and Ensembl IDs
filtered_top2 <- top2[!grepl(lncRNA_pattern, top2) & !grepl("^ENSG", top2)]

# Ensure unique, clean gene list
unique.top2_filtered <- unique(filtered_top2)


# Generate the DotPlot with enhanced aesthetics
library(ggplot2)

feature.plot <- DotPlot(merged_ref_final_subset, features = unique.top2_filtered) +
  theme_minimal() +  # Clean, minimal theme
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = rel(1.2), face = "bold", color = "black"),  # Larger x-axis labels
    axis.text.y = element_text(size = rel(1.2), face = "bold", color = "black"),  # Larger y-axis labels
    legend.text = element_text(size = rel(1), color = "black"),  # Larger legend text
    legend.title = element_text(size = rel(1.2), face = "bold", color = "black"),  # Bold legend title
    legend.spacing = unit(0.3, "cm")  # Adjust legend spacing for better fit
  ) +
  scale_x_discrete(name = "Top Marker Genes") +  # Rename x-axis
  scale_y_discrete(name = "Cell Type") +  # Rename y-axis
  scale_color_gradient(low = "blue", high = "red") +  # Adjust color scheme for better contrast
  scale_size(range = c(2, 10)) +  # Increase dot size range for better visibility
  guides(color = guide_colorbar(title = "Expression Level"), size = guide_legend(title = "Pct Expressed"))

# Save high-resolution plot with improved dimensions
png(file = paste0("de_cell_type_l2_subset_2FC_50perccells-DotPlot_top5.png"), res = 400, width = 5000, height = 2500)
print(feature.plot)  # Print to PNG file
dev.off()




feature.plot <- DotPlot(merged_ref_final_subset, features = unique.top2_filtered)
png(file=paste0("de_cell_type_l2_subset_2FC_50perccells-DotPlot.png"),res=300, width=3000, height=2000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Cell Type")
dev.off()

Idents(merged_ref_final_subset) <- "seurat_cluster.projected"
# We can also find and visualize the top gene expression markers for each cluster:
de_markers <- FindAllMarkers(merged_ref_final_subset, assay = "RNA", only.pos = TRUE, min.pct = 0.20, logfc.threshold = 1) # Detects FC >2 in 20% of cells in cluster; for marker genes
write.table(de_markers, "merged_ref_final_de_seurat_cluster.projected_subset_2FC_20perccells.txt", sep="\t")



# Define pattern for lncRNAs, uncharacterized genes, and Ensembl IDs
lncRNA_pattern <- "^RP[0-9]+-|^AC[0-9]+.|^AL[0-9]+.|^LINC|^MIR|^SNHG|^HCG|^CT[CDE]-"  # Updated to remove CTC, CTD genes

# Remove genes matching lncRNA pattern and Ensembl IDs
filtered_top2 <- top2[!grepl(lncRNA_pattern, top2) & !grepl("^ENSG", top2)]

# Ensure unique, clean gene list
unique.top2_filtered <- unique(filtered_top2)


# Generate the DotPlot with enhanced aesthetics
feature.plot <- DotPlot(merged_ref_final_subset, features = unique.top2_filtered) +
  theme_minimal() +  # Clean, minimal theme
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12, face = "bold", color = "black"),  # Black, readable x-axis labels
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),  # Black, bold y-axis labels
    legend.text = element_text(size = 10, color = "black"),  # Black legend text
    legend.title = element_text(size = 12, face = "bold", color = "black")  # Black, bold legend title
  ) +
  scale_x_discrete(name = "Top Marker Genes") +  # Rename x-axis
  scale_y_discrete(name = "Cluster") +  # Rename y-axis
  scale_color_gradient(low = "blue", high = "red") +  # Adjust color scheme for better contrast
  scale_size(range = c(1, 8)) +  # Increase dot size range for better differentiation
  guides(color = guide_colorbar(title = "Expression Level"), size = guide_legend(title = "Pct Expressed"))

# Save high-resolution plot with black labels
png(file = paste0("de_seurat_cluster.projected_subset_2FC_50perccells-DotPlot.png"), res = 400, width = 3000, height = 2500)
print(feature.plot)  # Print to PNG file
dev.off()





## This will require over 350 gb ram
	## revised with downsampling and densify 
DefaultAssay(merged_ref_final) <- "RNA"
gc()
de_markers <- FindAllMarkers(
  merged_ref_final, 
  slot="data", # "data" : difference in the log of the average exponentiated data, with pseudocount. This adjusts for differences in sequencing depth between cells, and assumes that "data" has been log-normalized.
  only.pos = TRUE,         # Keep only up regulated genes
  min.pct = 0.50,           # threshold for marker genes
  logfc.threshold = 1,     # Detects FC >2
  densify = TRUE ,            # Convert the sparse matrix to a dense form before running the DE test. This can provide speedups but might require higher memory
  max.cells.per.ident = 50000 , ## Down sample each identity class to a max number.
  verbose = TRUE
)
write.table(de_markers, "de_cell-types_2FC_50perccells.txt", sep="\t")

top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(merged_ref_final, features = unique.top2)
png(file=paste0("de_cell_type_l2_2FC_50perccells-DotPlot.png"),res=300, width=3000, height=2000)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Cell Type")
dev.off()


library(future)
# plan("multicore", workers = parallel::detectCores() - 2)  # Reserve 2 cores for system stability
plan("multisession", workers = parallel::detectCores() - 2)
plan("sequential")


