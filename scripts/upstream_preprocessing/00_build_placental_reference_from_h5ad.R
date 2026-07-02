#!/usr/bin/env Rscript
# Build a placental reference Seurat object from public h5ad files.
# This public-facing script is a cleaned scaffold derived from the archived server script.
# It assumes h5ad files have already been downloaded into data/external/placental_reference/.

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(zellkonverter)
  library(Matrix)
  library(dplyr)
})

set.seed(1)
options(future.globals.maxSize = 300e9)

input_manifest <- Sys.getenv("PLACENTAL_REF_MANIFEST", "config/placental_reference_h5ad_manifest.csv")
out_dir <- Sys.getenv("PLACENTAL_REF_OUTDIR", "outputs/upstream/placental_reference")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_manifest)) {
  stop("Missing manifest: ", input_manifest,
       "\nCreate it with columns: dataset_id, dataset_name, local_h5ad, keep_filter_note")
}

manifest <- read.csv(input_manifest, stringsAsFactors = FALSE)
required <- c("dataset_id", "dataset_name", "local_h5ad")
missing <- setdiff(required, colnames(manifest))
if (length(missing)) stop("Manifest is missing columns: ", paste(missing, collapse = ", "))

load_one_h5ad <- function(path, dataset_name) {
  message("Reading ", dataset_name, ": ", path)
  sce <- zellkonverter::readH5AD(path)
  counts <- SummarizedExperiment::assay(sce, "X")
  obj <- CreateSeuratObject(counts = counts, project = dataset_name)
  obj@meta.data <- as.data.frame(SummarizedExperiment::colData(sce))
  obj$reference_dataset <- dataset_name
  obj
}

objs <- lapply(seq_len(nrow(manifest)), function(i) {
  fp <- manifest$local_h5ad[i]
  if (!file.exists(fp)) stop("h5ad not found: ", fp)
  load_one_h5ad(fp, manifest$dataset_name[i])
})

# Harmonize cell-type column name when possible.
objs <- lapply(objs, function(o) {
  ct_candidates <- intersect(c("cell_type", "celltype", "Type", "author_cell_type"), colnames(o@meta.data))
  if (length(ct_candidates) && !"cell_type" %in% colnames(o@meta.data)) {
    o$cell_type <- o@meta.data[[ct_candidates[1]]]
  }
  o
})

ref <- Reduce(function(x, y) merge(x, y = y, merge.data = TRUE, project = "placental_reference"), objs)
ref <- JoinLayers(ref)

saveRDS(ref, file.path(out_dir, "placental_reference_merged_seurat.rds"))

write.csv(as.data.frame(table(ref$reference_dataset)),
          file.path(out_dir, "placental_reference_dataset_counts.csv"), row.names = FALSE)
if ("cell_type" %in% colnames(ref@meta.data)) {
  write.csv(as.data.frame(table(ref$cell_type)),
            file.path(out_dir, "placental_reference_celltype_counts.csv"), row.names = FALSE)
}

sessionInfo()
