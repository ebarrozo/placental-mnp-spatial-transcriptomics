
# rebuild_viz_embed_only.R
# Reconstruct a lightweight Seurat object ("viz") from already-mapped chunks,
# using the manifest of chunks that actually landed in viz.
# - No counts or metadata enrichment here (keep it minimal to avoid memory/file issues)
# - Produces an object with a tiny dummy SCT assay + ref.pca and (if present) UMAP reductions
# - Safe on server 3: all paths are local/customizable below

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(stringr)
  library(tibble)
})

## =========================
## Config (EDIT ME)
## =========================
# Where your "chunks_in_viz_detailed.csv" lives (from tonight's run)
CIV_CSV     <- "/path/to/local_home/visiumhd_mnp/docs/chunks_in_viz_detailed.csv"   # e.g., "/path/to/local_home/visiumhd_mnp/docs/chunks_in_viz_detailed.csv"

# Location of mapped chunks (produced by your previous pipeline)
MAPPED_DIR  <- "/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/SCT_chunks_mapped"            # e.g., "/path/to/local_storage/raid3/.../SCT_chunks_mapped"

# Output folder (local to server 3). Will be created if missing.
SAVE_DIR    <- file.path(getwd(), "viz_rebuild")  # change to "/path/to/local_home/visiumhd_mnp/viz"
dir.create(SAVE_DIR, showWarnings = FALSE, recursive = TRUE)

# Output filename
OUT_RDS     <- file.path(SAVE_DIR, sprintf("viz_embed_only_%s.rds", format(Sys.time(), "%Y%m%d_%H%M%S")))

# If TRUE, stop early if any mapped files are missing.
STRICT_MAPPED <- TRUE

## =========================
## Helpers
## =========================
# compute expected mapped filename for a (sample, chunk)
mapped_path <- function(sample, chunk, mapped_dir = MAPPED_DIR) {
  file.path(mapped_dir, sprintf("sctrans_%s_chunk_%03d_SCTonly_mapped.rds", sample, as.integer(chunk)))
}

# Read a mapped chunk safely & extract reductions + minimal meta
read_mapped_chunk <- function(path) {
  o <- readRDS(path)
  # sanity: assay
  if (!"SCT" %in% names(o@assays)) stop("No SCT assay in: ", basename(path))

  # embeddings
  has_pca  <- "ref.pca" %in% Reductions(o)
  has_umap <- "umap" %in% Reductions(o)

  if (!has_pca && !has_umap) stop("No ref.pca or umap reductions found in: ", basename(path))

  E_pca  <- if (has_pca)  Embeddings(o, "ref.pca") else NULL  # cells x dims
  E_umap <- if (has_umap) Embeddings(o, "umap")    else NULL

  # minimal meta: keep only stable columns; add sample + chunk from colnames
  md <- o@meta.data
  # coerce factors/lists to character
  for (i in seq_len(ncol(md))) {
    xi <- md[[i]]
    if (is.factor(xi)) md[[i]] <- as.character(xi)
    else if (is.list(xi)) md[[i]] <- vapply(xi, function(z) if (length(z)) as.character(z)[1] else NA_character_, "")
  }
  # always compute sample & chunk from cell names: <sample>_cNNN_<barcode>
  cn <- rownames(md)
  sample <- sub("_c\\d{3}_.+$", "", cn)
  chunk  <- suppressWarnings(as.integer(sub("^.*_c(\\d{3})_.+$", "\\1", cn)))

  keep_cols <- intersect(colnames(md), c("orig.ident","predicted.celltype"))
  md_out <- md[, keep_cols, drop = FALSE]
  md_out$sample <- sample
  md_out$chunk  <- chunk
  rownames(md_out) <- cn

  # free original object
  rm(o); gc()

  list(pca = E_pca, umap = E_umap, meta = md_out)
}

## =========================
## 0) Load the manifest of chunks that made it into viz
## =========================
if (!file.exists(CIV_CSV)) stop("Manifest CSV not found: ", CIV_CSV)
civ <- read.csv(CIV_CSV, stringsAsFactors = FALSE)

# Expect at least columns: sample, chunk
stopifnot(all(c("sample","chunk") %in% names(civ)))

# Build expected mapped paths; deduplicate sample-chunk pairs
civ_sc <- civ %>% distinct(sample, chunk)
civ_sc$mapped_file <- mapped_path(civ_sc$sample, civ_sc$chunk)

exists_flag <- file.exists(civ_sc$mapped_file)
missing_tbl <- civ_sc[!exists_flag, c("sample","chunk","mapped_file")]
if (nrow(missing_tbl)) {
  message("Missing mapped files: ", nrow(missing_tbl))
  print(utils::head(missing_tbl, 20))
  write.csv(missing_tbl, file.path(SAVE_DIR, "missing_mapped_files.csv"), row.names = FALSE)
  if (STRICT_MAPPED) stop("Abort: some mapped files are missing. See missing_mapped_files.csv")
}

present_tbl <- civ_sc[exists_flag, , drop = FALSE]
message("Found mapped files: ", nrow(present_tbl), " / ", nrow(civ_sc))

## =========================
## 1) Stream mapped chunks → collect embeddings and minimal meta
## =========================
pca_list   <- list()
umap_list  <- list()
meta_list  <- list()

min_pca_dims <- Inf
have_any_umap <- FALSE

pb_every <- max(1L, floor(nrow(present_tbl)/50))  # progress every ~2%
for (i in seq_len(nrow(present_tbl))) {
  path <- present_tbl$mapped_file[i]
  if (!file.exists(path)) next

  if (i %% pb_every == 0) message(sprintf("[%d/%d] %s", i, nrow(present_tbl), basename(path)))

  parts <- read_mapped_chunk(path)

  if (!is.null(parts$pca)) {
    min_pca_dims <- min(min_pca_dims, ncol(parts$pca))
    pca_list[[length(pca_list)+1]] <- parts$pca
  }
  if (!is.null(parts$umap)) {
    have_any_umap <- TRUE
    umap_list[[length(umap_list)+1]] <- parts$umap
  }
  meta_list[[length(meta_list)+1]] <- parts$meta
}

if (!length(meta_list)) stop("No mapped chunks could be read.")

# rbind meta (safe because we converted to simple types)
meta_all <- dplyr::bind_rows(meta_list)

# Align PCA dims across chunks (truncate to common min)
if (is.finite(min_pca_dims)) {
  pca_all <- do.call(rbind, lapply(pca_list, function(M) M[, seq_len(min_pca_dims), drop = FALSE]))
  colnames(pca_all) <- paste0("refPC_", seq_len(ncol(pca_all)))
  # ensure rownames match the meta cells
  stopifnot(nrow(pca_all) == nrow(meta_all))
} else {
  pca_all <- NULL
}

# UMAP (only if present everywhere; otherwise we can still store what we have)
umap_all <- NULL
if (have_any_umap && length(umap_list)) {
  # ensure all UMAPs have same ncol (typically 2)
  umap_dims <- unique(vapply(umap_list, ncol, integer(1)))
  if (length(umap_dims) > 1L) {
    # truncate to min
    md <- min(umap_dims)
    umap_all <- do.call(rbind, lapply(umap_list, function(M) M[, seq_len(md), drop = FALSE]))
  } else {
    umap_all <- do.call(rbind, umap_list)
  }
  colnames(umap_all) <- paste0("UMAP_", seq_len(ncol(umap_all)))
  stopifnot(nrow(umap_all) == nrow(meta_all))
}

## =========================
## 2) Build a minimal Seurat object with a dummy assay
## =========================
cells <- rownames(meta_all)
stopifnot(!anyDuplicated(cells))

# 1×N sparse dummy matrix (feature "DUMMY1"), just to satisfy Seurat
dummy <- Matrix(0, nrow = 1, ncol = length(cells), sparse = TRUE,
                dimnames = list(c("DUMMY1"), cells))

assay <- CreateAssayObject(counts = dummy)
viz   <- CreateSeuratObject(assays = list(SCT = assay), meta.data = meta_all)
DefaultAssay(viz) <- "SCT"

# attach reductions
if (!is.null(pca_all)) {
  viz[["ref.pca"]] <- CreateDimReducObject(embeddings = as.matrix(pca_all),
                                           key = "refPC_", assay = "SCT")
}
if (!is.null(umap_all)) {
  viz[["umap"]]    <- CreateDimReducObject(embeddings = as.matrix(umap_all),
                                           key = "UMAP_", assay = "SCT")
}

# identities default to sample (parsed)
Idents(viz) <- factor(viz$sample)

## =========================
## 3) Save + summary
## =========================
saveRDS(viz, OUT_RDS, compress = "gzip")
cat("\nSaved viz ->", OUT_RDS, "\n")

cat("\nSummary:\n")
cat("  Cells:", ncol(viz), "\n")
cat("  Meta columns:", ncol(viz@meta.data), "\n")
cat("  Reductions:", paste(Reductions(viz), collapse = ", "), "\n")
if ("ref.pca" %in% Reductions(viz)) cat("  ref.pca dims:", ncol(Embeddings(viz, 'ref.pca')), "\n")
if ("umap"   %in% Reductions(viz)) cat("  umap dims:",   ncol(Embeddings(viz, 'umap')),   "\n")

# Done.
