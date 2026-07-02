## harmonized_vis_v1.R

set.seed(1) # set seed for reproducibility
options(future.globals.maxSize = 400e9) ## 300 gb ram (out of 440 on server 3; adjust if on different server)

setwd("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

## --- config ---
OUTDIR_MAPPED <- "SCT_chunks_mapped"  # your mapped dir
# Optional: cap cells per chunk for plotting (set to Inf to keep all)
MAX_CELLS_PER_CHUNK <- 30000L   # e.g., keep all if each chunk is 30k; set 5e3/1e4 if you want lighter viz
SEED <- 1

## --- helpers ---
parse_sample_chunk <- function(bn) {
  # e.g., sctrans_HP2_hd_chunk_015_SCTonly_mapped.rds
  x <- sub("^sctrans_(.+?)_chunk_(\\d+)_SCTonly_mapped\\.rds$", "\\1__\\2", bn)
  parts <- strsplit(x, "__", fixed = TRUE)[[1]]
  list(sample = parts[1], chunk = sprintf("%03d", as.integer(parts[2])))
}

# load one mapped file -> returns list(umap=..., pca=..., meta=..., cells=chr)
load_one_mapped <- function(fp, keep_n = Inf) {
  o <- readRDS(fp)
  stopifnot("umap" %in% Reductions(o), "ref.pca" %in% Reductions(o))
  um <- Embeddings(o, "umap")
  rp <- Embeddings(o, "ref.pca")

  # align cells across reductions
  common <- intersect(rownames(um), rownames(rp))
  um <- um[common, , drop = FALSE]
  rp <- rp[common, , drop = FALSE]

  # optional downsample per chunk
  if (is.finite(keep_n) && nrow(um) > keep_n) {
    set.seed(SEED)
    idx <- sample.int(nrow(um), keep_n)
    um <- um[idx, , drop = FALSE]
    rp <- rp[idx, , drop = FALSE]
    common <- rownames(um)
  }

  md <- o@meta.data[common, , drop = FALSE]
  # keep only columns you care about; add more if needed
  keep_meta <- intersect(c("predicted.celltype","predicted.celltype.score"), colnames(md))
  md <- md[, keep_meta, drop = FALSE]

  # add sample / chunk from filename
  id <- parse_sample_chunk(basename(fp))
  md$sample <- id$sample
  md$chunk  <- id$chunk

  list(umap = um, pca = rp, meta = md, cells = common)
}

# Build a minimal Seurat object from big embeddings + meta
build_seurat_skeleton <- function(umap_all, pca_all, meta_all) {
  stopifnot(nrow(umap_all) == nrow(pca_all),
            nrow(umap_all) == nrow(meta_all),
            all(rownames(umap_all) == rownames(pca_all)),
            all(rownames(umap_all) == rownames(meta_all)))

  N <- nrow(umap_all)
  # tiny dummy assay (1 feature x N cells)
  dummy_counts <- Matrix(0, nrow = 1, ncol = N, sparse = TRUE)
  rownames(dummy_counts) <- "dummy"
  colnames(dummy_counts) <- rownames(umap_all)

  obj <- CreateSeuratObject(counts = dummy_counts, assay = "RNA", project = "viz_skeleton")
  obj <- AddMetaData(obj, meta = meta_all)
  DefaultAssay(obj) <- "RNA"

  # attach reductions (coerce to base matrix)
  obj[["ref.pca"]] <- CreateDimReducObject(
    embeddings = as.matrix(pca_all),
    key        = "refPC_",
    assay      = DefaultAssay(obj)
  )
  obj[["umap"]] <- CreateDimReducObject(
    embeddings = as.matrix(umap_all),
    key        = "UMAP_",
    assay      = DefaultAssay(obj)
  )
  obj
}

## --- load mapped files (skip empties/corrupts automatically) ---
mapped_files <- list.files(OUTDIR_MAPPED, pattern = "_mapped\\.rds$", full.names = TRUE)
stopifnot(length(mapped_files) > 0)

# stream through files to limit peak RAM
UMAP_list <- vector("list", length(mapped_files))
PCA_list  <- vector("list", length(mapped_files))
META_list <- vector("list", length(mapped_files))
n_kept    <- integer(length(mapped_files))

cat("Collecting embeddings/meta from ", length(mapped_files), " mapped chunks ...\n", sep = "")
for (i in seq_along(mapped_files)) {
  fp <- mapped_files[[i]]
  id <- parse_sample_chunk(basename(fp))
  cat(sprintf("[%4d/%4d] %s (sample=%s chunk=%s)\n",
              i, length(mapped_files), basename(fp), id$sample, id$chunk))
  x <- try(load_one_mapped(fp, keep_n = MAX_CELLS_PER_CHUNK), silent = TRUE)
  if (inherits(x, "try-error")) {
    message("  -> skipped (read/structure error): ", conditionMessage(attr(x, "condition")))
    next
  }
  UMAP_list[[i]] <- x$umap
  PCA_list[[i]]  <- x$pca
  META_list[[i]] <- x$meta
  n_kept[i]      <- nrow(x$umap)
  rm(x); gc()
}
kept <- which(n_kept > 0L)
stopifnot(length(kept) > 0)

UMAP_all <- do.call(rbind, UMAP_list[kept])
PCA_all  <- do.call(rbind, PCA_list[kept])
META_all <- do.call(rbind, META_list[kept])

cat("Total cells in viz object: ", nrow(UMAP_all), "\n", sep = "")

## --- build the lightweight Seurat viz object ---
viz <- build_seurat_skeleton(UMAP_all, PCA_all, META_all)

# set identity (optional)
if ("predicted.celltype" %in% colnames(viz@meta.data)) {
  Idents(viz) <- "predicted.celltype"
}

## --- example plots (Seurat-native) ---
# 1) UMAP by sample
DimPlot(viz, reduction = "umap", group.by = "sample", raster = TRUE, pt.size = 0.2)

# 2) UMAP by predicted.celltype with labels
if ("predicted.celltype" %in% colnames(viz@meta.data)) {
  DimPlot(viz, reduction = "umap", group.by = "predicted.celltype", label = TRUE, repel = TRUE,
          raster = TRUE, pt.size = 0.2)
}

# 3) FeaturePlot on a metadata column (scores)
if ("predicted.celltype.score" %in% colnames(viz@meta.data)) {
  FeaturePlot(viz, reduction = "umap", features = "predicted.celltype.score", pt.size = 0.2)
}

# 4) VlnPlot of a metadata column by cell type
if (all(c("predicted.celltype","predicted.celltype.score") %in% colnames(viz@meta.data))) {
  VlnPlot(viz, features = "predicted.celltype.score", group.by = "predicted.celltype", pt.size = 0)
}

## --- optional: save the lightweight viz object for reuse ---
saveRDS(viz, file = file.path(OUTDIR_MAPPED, "viz_skeleton_umap_pca_meta.rds"))





load_one_mapped <- function(fp, keep_n = Inf) {
  o <- readRDS(fp)
  stopifnot("umap" %in% Reductions(o), "ref.pca" %in% Reductions(o))
  um <- Embeddings(o, "umap")
  rp <- Embeddings(o, "ref.pca")

  # align
  common <- intersect(rownames(um), rownames(rp))
  um <- um[common, , drop = FALSE]
  rp <- rp[common, , drop = FALSE]

  # (optional) downsample per chunk
  if (is.finite(keep_n) && nrow(um) > keep_n) {
    set.seed(SEED)
    idx <- sample.int(nrow(um), keep_n)
    um <- um[idx, , drop = FALSE]
    rp <- rp[idx, , drop = FALSE]
    common <- rownames(um)
  }

  md <- o@meta.data[common, , drop = FALSE]

  # --- ensure a consistent schema across chunks ---
  if (!"predicted.celltype"       %in% colnames(md)) md$predicted.celltype       <- NA_character_
  if (!"predicted.celltype.score" %in% colnames(md)) md$predicted.celltype.score <- NA_real_

  # coerce types (avoid factors/list-cols)
  md$predicted.celltype       <- as.character(md$predicted.celltype)
  md$predicted.celltype.score <- as.numeric(md$predicted.celltype.score)

  # add sample/chunk from filename
  id <- parse_sample_chunk(basename(fp))
  md$sample <- id$sample
  md$chunk  <- id$chunk

  # enforce column order
  md <- md[, c("predicted.celltype","predicted.celltype.score","sample","chunk"), drop = FALSE]

  list(umap = um, pca = rp, meta = md, cells = common)
}
UMAP_all <- do.call(rbind, UMAP_list[kept])
PCA_all  <- do.call(rbind, PCA_list[kept])
META_all <- do.call(rbind, META_list[kept])

stopifnot(
  identical(rownames(UMAP_all), rownames(PCA_all)),
  identical(rownames(UMAP_all), rownames(META_all))
)

viz <- build_seurat_skeleton(UMAP_all, PCA_all, META_all)







suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

# expect these already exist from your previous step:
# - mapped_files
# - UMAP_list, PCA_list, META_list, kept
# - parse_sample_chunk(), build_seurat_skeleton(), MAX_CELLS_PER_CHUNK, SEED

# 1) Define a canonical meta schema
target_cols <- c("predicted.celltype","predicted.celltype.score","sample","chunk")

# 2) Harmonize each meta df to the schema + row order = UMAP rows
META_list_h <- mapply(function(md, um, fp) {
  rn <- rownames(um)

  # ensure we have a data.frame (not tibble/NULL)
  if (is.null(md) || !is.data.frame(md)) {
    md <- data.frame(row.names = rn)
  } else {
    # keep only rows we have embeddings for
    if (!identical(rownames(md), rn)) md <- md[rn, , drop = FALSE]
  }

  # add missing columns with proper types
  if (!"predicted.celltype"       %in% names(md)) md$predicted.celltype       <- NA_character_
  if (!"predicted.celltype.score" %in% names(md)) md$predicted.celltype.score <- NA_real_
  # enforce base types (avoid factors/list-cols)
  md$predicted.celltype       <- as.character(md$predicted.celltype)
  md$predicted.celltype.score <- as.numeric(md$predicted.celltype.score)

  # sample/chunk from filename (override if missing)
  id <- parse_sample_chunk(basename(fp))
  if (!"sample" %in% names(md)) md$sample <- id$sample
  if (!"chunk"  %in% names(md)) md$chunk  <- id$chunk
  md$sample <- as.character(md$sample)
  md$chunk  <- as.character(md$chunk)

  # enforce column order exactly
  md <- md[, target_cols, drop = FALSE]

  # add a 'cell' column so we can safely rbind via data.table
  md$cell <- rn
  md
}, META_list[kept], UMAP_list[kept], mapped_files[kept], SIMPLIFY = FALSE)

# 3) Combine fast & safely (fills if anything is still missing)
META_all <- data.table::rbindlist(META_list_h, fill = TRUE)
# restore rownames, drop helper 'cell'
rn_all <- META_all$cell; META_all$cell <- NULL
META_all <- as.data.frame(META_all)
rownames(META_all) <- rn_all

# 4) Combine embeddings
UMAP_all <- do.call(rbind, UMAP_list[kept])
PCA_all  <- do.call(rbind, PCA_list[kept])

# 5) Final row alignment check (and enforce)
stopifnot(identical(rownames(UMAP_all), rownames(PCA_all)))
if (!identical(rownames(UMAP_all), rownames(META_all))) {
  META_all <- META_all[rownames(UMAP_all), , drop = FALSE]
}

# 6) Build the lightweight Seurat viz object
viz <- build_seurat_skeleton(UMAP_all, PCA_all, META_all)

# Optional: set identities
if ("predicted.celltype" %in% colnames(viz@meta.data)) {
  Idents(viz) <- "predicted.celltype"
}

# 7) Example plots
DimPlot(viz, reduction = "umap", group.by = "sample", raster = TRUE, pt.size = 0.2)
if ("predicted.celltype" %in% colnames(viz@meta.data)) {
  DimPlot(viz, reduction = "umap", group.by = "predicted.celltype", label = TRUE, repel = TRUE,
          raster = TRUE, pt.size = 0.2)
}
if ("predicted.celltype.score" %in% colnames(viz@meta.data)) {
  FeaturePlot(viz, reduction = "umap", features = "predicted.celltype.score", pt.size = 0.2)
}

# Save for reuse
saveRDS(viz, file = file.path(OUTDIR_MAPPED, "viz_skeleton_umap_pca_meta.rds"))

















build_seurat_skeleton <- function(umap_all, pca_all, meta_all, assay_name = "RNA") {
  stopifnot(
    nrow(umap_all) == nrow(pca_all),
    nrow(umap_all) == nrow(meta_all),
    identical(rownames(umap_all), rownames(pca_all)),
    identical(rownames(umap_all), rownames(meta_all))
  )

  N <- nrow(umap_all)

  # --- key change: use ≥2 dummy features to avoid drop-to-vector bugs
  dummy_counts <- Matrix::Matrix(0, nrow = 2, ncol = N, sparse = TRUE)
  rownames(dummy_counts) <- c("dummy1", "dummy2")
  colnames(dummy_counts) <- rownames(umap_all)

  obj <- Seurat::CreateSeuratObject(counts = dummy_counts, assay = assay_name, project = "viz_skeleton")
  Seurat::DefaultAssay(obj) <- assay_name

  # For Seurat v5 Assay5, ensure a "data" layer exists (some code expects it)
  if (inherits(obj@assays[[assay_name]], "Assay5")) {
    ass <- obj@assays[[assay_name]]
    if (is.null(ass@layers[["data"]])) ass@layers[["data"]] <- ass@layers[["counts"]]
    obj@assays[[assay_name]] <- ass
  }

  # add meta
  obj <- Seurat::AddMetaData(obj, meta = meta_all)

  # attach reductions (coerce to base matrices)
  obj[["ref.pca"]] <- Seurat::CreateDimReducObject(
    embeddings = as.matrix(pca_all),
    key        = "refPC_",
    assay      = assay_name
  )
  obj[["umap"]] <- Seurat::CreateDimReducObject(
    embeddings = as.matrix(umap_all),
    key        = "UMAP_",
    assay      = assay_name
  )

  obj
}
# Combine (as you did after harmonizing META_list_*):
UMAP_all <- do.call(rbind, UMAP_list[kept])
PCA_all  <- do.call(rbind, PCA_list[kept])
# META_all was built via rbindlist(fill=TRUE) and aligned to rownames(UMAP_all)

viz <- build_seurat_skeleton(UMAP_all, PCA_all, META_all)

# sanity
Seurat::Reductions(viz)             # "ref.pca", "umap"
head(Seurat::Embeddings(viz, "umap"))
Seurat::DimPlot(viz, reduction = "umap", group.by = "sample", raster = TRUE, pt.size = 0.2)

if ("predicted.celltype" %in% colnames(viz@meta.data)) {
  Seurat::DimPlot(viz, reduction = "umap", group.by = "predicted.celltype",
                  label = TRUE, repel = TRUE, raster = TRUE, pt.size = 0.2)
}

if ("predicted.celltype.score" %in% colnames(viz@meta.data)) {
  Seurat::FeaturePlot(viz, reduction = "umap", features = "predicted.celltype.score", pt.size = 0.2)
}
options(Seurat.object.validate = FALSE)


ncol(viz) # 3969532 cells
levels(viz) # 28 orig ident

library(ggplot2)
library(cowplot)
library(patchwork)
library(ggrepel)
library(ggpubr)


mypal3 <-get_palette("igv",40)

Idents(viz) <- "predicted.celltype"
p2<-DimPlot(viz, reduction = "umap", group.by = "predicted.celltype", alpha = 0.1, label = FALSE,raster=T, cols = mypal3) + ggtitle("Predicted Cell Types (3,969,532 Cells)") + NoAxes()
p2
ggsave("viz_predicted.predicted.celltype.png", plot = p2, device = "png", width = 8, height = 5, units = "in", dpi = 400)
p2<-DimPlot(viz, reduction = "umap", group.by = "predicted.celltype", alpha = 0.1, label = TRUE, repel=TRUE,raster=T, cols = mypal3) + ggtitle("Predicted Cell Types (3,969,532 Cells)") + NoAxes()
p2
ggsave("viz_predicted.predicted.celltype_labeled.png", plot = p2, device = "png", width = 16, height = 10, units = "in", dpi = 400)


table(viz$predicted.celltype, viz$orig.ident)
t1<- table(viz$predicted.celltype, viz$orig.ident)
write.csv(t1, "viz-predicted.celltype_orig.ident.csv", row.names = T)





# Save for reuse
saveRDS(viz, file = file.path(OUTDIR_MAPPED, "viz_skeleton_umap_pca_meta.rds"))


viz_subset <- subset(x = viz, idents = c("Syncytiotrophoblast"     , "Villous Cytotrophoblast" , "Extravillous Trophoblast" ,"Macrophage"         ,      "B Cell"     ,              "NK Cell"                 
,"T Cell"             ,      "Perivascular Cell"   ,     "Hofbauer Cell"  ,          "Placental Endothelial"  ,  "Treg"    ,                 "Decidual NK Cell"        
, "Fibroblast"         ,      "ILC3"           ,          "Plasma Cell"))

ncol(viz_subset) # 




## reference not predicted.celltypes 

suppressPackageStartupMessages({library(data.table); library(ggplot2); library(scales)})

# --- 1) load reference with its native labels ---
ref_path <- file.path(OUTDIR_MAPPED, "integrated_reference_SCT.rds")
ref <- readRDS(ref_path)

# pick a label column or fall back to Idents
ref_md <- as.data.table(ref@meta.data, keep.rownames = "ref_cell")
label_candidates <- c("cell_type_expanded")
labcol <- intersect(label_candidates, colnames(ref_md))
if (length(labcol)) {
  ref_md[, ref_label := as.character(get(labcol[1]))]
} else {
  ref_md[, ref_label := as.character(Idents(ref))]
}

# --- 2) extract original barcode from viz cell names & merge labels into viz meta (ref rows only) ---
md <- as.data.table(viz@meta.data, keep.rownames = "cell")
md[, sample := as.character(sample)]
is_ref <- md$sample == "placental_ref_Jul_2025"

# original barcode in your naming: sample_cNNN_<orig>
md[, orig := sub("^[^_]+_c\\d+_", "", cell)]

# fast exact join first
setkey(ref_md, ref_cell); setkey(md, orig)
md[is_ref, ref_label := ref_md[md[is_ref]]$ref_label]

# fallback: if many NAs, try matching on a version of orig with last library suffix removed
miss <- is_ref & is.na(md$ref_label)
if (any(miss)) {
  tmp <- copy(ref_md)
  tmp[, ref_cell_trim := sub("_[^_]+$", "", ref_cell)]
  setkey(tmp, ref_cell_trim)
  md[miss, ref_label := tmp[.(sub("_[^_]+$", "", orig[miss]))]$ref_label]
}

# summary of matches
n_ref <- sum(is_ref)
n_match <- sum(is_ref & !is.na(md$ref_label) & md$ref_label != "")
cat(sprintf("Reference labeling: matched %s / %s (%.1f%%)\n",
            comma(n_match), comma(n_ref), 100*n_match/n_ref))

# --- 3) build display columns (use ref labels for ref rows) ---
md[, predicted.celltype := ifelse(is.na(predicted.celltype) | predicted.celltype=="", "NA", predicted.celltype)]
md[, celltype_display := predicted.celltype]
md[is_ref & !is.na(ref_label) & ref_label != "", celltype_display := ref_label]
md[, score_display := predicted.celltype.score]
md[is_ref, score_display := NA_real_]  # reference labels aren't RCTD scores

# push back to viz
rownames_md <- md$cell; md$cell <- NULL
viz@meta.data <- as.data.frame(md); rownames(viz@meta.data) <- rownames_md

# --- 4) recompute composition using celltype_display (stacked % bars incl. NA) ---
tot <- nrow(viz@meta.data)
comp <- as.data.table(viz@meta.data)[, .N, by = .(sample, celltype_display)]
comp[, pct := 100*N/sum(N), by = sample]

# order legend by global abundance; samples by NA%
ct_order <- comp[, .(N=sum(N)), by = celltype_display][order(-N), celltype_display]
na_pct   <- comp[celltype_display=="NA", .(na_pct = pct), by = sample]
s_order  <- na_pct[order(-na_pct), sample]

comp[, celltype_display := factor(celltype_display, levels = ct_order)]
comp[, sample := factor(sample, levels = s_order)]

p <- ggplot(comp, aes(x = sample, y = pct, fill = celltype_display)) +
  geom_col(width = 0.9) +
  coord_flip() +
  labs(x = "Sample", y = "% of cells",
       fill = "Cell type (ref or RCTD)",
       title = sprintf("Cell-type composition per sample (n=%s)", comma(tot))) +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())

h <- max(6, 0.25 * length(unique(comp$sample)))
ggsave("composition_per_sample_stacked_pct_refAware.png", p, width = 14, height = h, dpi = 300)

# --- 5) updated NA quantification (after substituting ref labels) ---
md2 <- as.data.table(viz@meta.data)
md2[, is_na := celltype_display == "NA"]
overall_na <- md2[, mean(is_na)]
cat(sprintf("After ref substitution, NA overall: %.2f%% (%s / %s)\n",
            100*overall_na, comma(sum(md2$is_na)), comma(nrow(md2))))

# also per-sample NA %
na_by_sample <- md2[, .(na_pct = 100*mean(is_na), na = sum(is_na), total = .N), by = sample][order(-na_pct)]
data.table::fwrite(na_by_sample, "NA_fraction_by_sample_after_ref.csv")


suppressPackageStartupMessages({library(data.table); library(Seurat)})

## paths/ids you already use
REF_PATH <- "/path/to/local_storage/raid2/project/spatial_MNP/placental_ref_Jul_2025.rds"
REF_ID   <- "placental_ref_Jul_2025"   # sample name used in viz

# --- 1) load reference & pick its label column (fallback: Idents) ---
ref <- readRDS(REF_PATH)
ref_md <- as.data.table(ref@meta.data, keep.rownames = "ref_cell")

label_candidates <- c("celltype","cell_type","annotation","annot","label",
                      "CellType","broad_celltype","fine_celltype","celltype_final")
labcol <- intersect(label_candidates, colnames(ref_md))
if (length(labcol)) {
  ref_md[, ref_label := as.character(get(labcol[1]))]
} else {
  ref_md[, ref_label := as.character(Idents(ref))]
}
# helpful trims for alternate barcode formats
ref_md[, ref_cell_trim1 := sub("-1$", "", ref_cell)]
ref_md[, ref_cell_trim2 := sub("_[^_]+$", "", ref_cell)]   # drop last “_lib” suffix if present

# --- 2) extract reference rows from viz and reconstruct the original barcode ---
md <- as.data.table(viz@meta.data, keep.rownames = "cell")
md[, sample := as.character(sample)]
is_ref <- md$sample == REF_ID

# cell names look like: "<REF_ID>_cXXX_<orig>"
md[is_ref, orig := sub(paste0("^", REF_ID, "_c\\d+_"), "", cell)]
md[is_ref, orig_trim1 := sub("-1$", "", orig)]
md[is_ref, orig_trim2 := sub("_[^_]+$", "", orig)]

# --- 3) join cascade: exact → trim “-1” → trim last “_lib” ---
setkey(ref_md, ref_cell); setkey(md, orig)
md[is_ref, ref_label := ref_md[md[is_ref]]$ref_label]

miss <- is_ref & is.na(md$ref_label)
if (any(miss)) {
  setkey(ref_md, ref_cell_trim1); setkey(md, orig_trim1)
  md[miss, ref_label := ref_md[md[miss]]$ref_label]
}
miss <- is_ref & is.na(md$ref_label)
if (any(miss)) {
  setkey(ref_md, ref_cell_trim2); setkey(md, orig_trim2)
  md[miss, ref_label := ref_md[md[miss]]$ref_label]
}

# --- 4) write back to viz meta: use ref label for reference rows; keep RCTD for others ---
n_ref    <- sum(is_ref)
n_match  <- sum(is_ref & !is.na(md$ref_label) & md$ref_label != "")
cat(sprintf("Reference labeling matched %d / %d (%.1f%%)\n", n_match, n_ref, 100*n_match/n_ref))

# clean RCTD column
md[, predicted.celltype := ifelse(is.na(predicted.celltype) | predicted.celltype=="", "NA", predicted.celltype)]

# final display + source
md[, celltype_display := predicted.celltype]
md[is_ref & !is.na(ref_label) & ref_label != "", celltype_display := ref_label]
md[, label_source := fifelse(is_ref & !is.na(ref_label) & ref_label != "", "ref",
                             fifelse(predicted.celltype!="NA","rctd","NA"))]

# persist back into viz
rn <- md$cell; md$cell <- NULL
viz@meta.data <- as.data.frame(md); rownames(viz@meta.data) <- rn

# quick sanity
table(viz$sample == REF_ID, is.na(viz$ref_label))  # few TRUE/FALSE in second col is good
table(viz$label_source)                            # breakdown of ref/rctd/NA










mypal3 <-get_palette("igv",40)

Idents(viz) <- "predicted.celltype"
p2<-DimPlot(viz, reduction = "umap", group.by = "predicted.celltype", alpha = 0.1, label = FALSE,raster=T, cols = mypal3) + ggtitle("Predicted Cell Types (3,969,532 Cells)") + NoAxes()
p2
ggsave("viz_predicted.predicted.celltype.png", plot = p2, device = "png", width = 8, height = 5, units = "in", dpi = 1200)


c("HP1"  ,     "HP13a"  ,   "HP13b"  ,   "HP2",       "HT1"    ,   "HT2"   ,    "HT2a"   ,   "HT2b"  ,    "HT3"  ,     "LP1"  ,     "LP2"    ,   "LT1"   ,    "LT2"      ,"NC1a" ,     "NC1b"   ,   "NC1c"   ,   "NC2"  ,     "NC3"   ,    "NC4"  ,     "ND5"    ,   "ND6"   ,    "ND7"    ,   "Reference", "SP10"   ,   "SP11"     , "SP12"    , "SP8"  ,    "SP9")


Idents(viz) <- "orig.ident"
levels(viz)         # 
new.metadata <- c("HP1"  ,     "HP13a"  ,   "HP13b"  ,   "HP2",       "HT1"    ,   "HT2"   ,    "HT2a"   ,   "HT2b"  ,    "HT3"  ,     "LP1"  ,     "LP2"    ,   "LT1"   ,    "LT2"      ,"NC1a" ,     "NC1b"   ,   "NC1c"   ,   "NC2"  ,     "NC3"   ,    "NC4"  ,     "ND5"    ,   "ND6"   ,    "ND7"    ,   "Reference", "SP10"   ,   "SP11"     , "SP12"    , "SP8"  ,    "SP9")
names(new.metadata) <- levels(viz)
viz <- RenameIdents(viz, new.metadata)
viz$orig.ident <- Idents(viz)
Idents(viz) <- "orig.ident"

Idents(viz) <- "orig.ident"
p2<-DimPlot(viz, reduction = "umap", group.by = "orig.ident", alpha = 0.1, label = FALSE,raster=T, cols = mypal3) + ggtitle("Original Identity (3,969,532 Cells)") + NoAxes()
p2
ggsave("viz_predicted.orig.ident.png", plot = p2, device = "png", width = 8, height = 5, units = "in", dpi = 1200)


# Save for reuse
saveRDS(viz, file = file.path(OUTDIR_MAPPED, "viz_skeleton_umap_pca_meta.rds"))





viz_subset <- subset(x = viz, idents = c("Syncytiotrophoblast"     , "Villous Cytotrophoblast" , "Extravillous Trophoblast" ,"Macrophage"         ,      "B Cell"     ,              "NK Cell"                 
,"T Cell"             ,      "Perivascular Cell"   ,     "Hofbauer Cell"  ,          "Placental Endothelial"  ,  "Treg"    ,                 "Decidual NK Cell"        
, "Fibroblast"         ,      "ILC3"           ,          "Plasma Cell"))

ncol(viz_subset) # 
