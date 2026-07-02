## obj_integrated_final.R — reproducible SCT merge (Seurat-native, v5-safe)
## Enrico R. Barrozo, Ph.D. — 2025-09-25

setwd("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations")

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject)
  library(Matrix); library(dplyr); library(stringr); library(tibble); library(tools)
})

options(warn = 1)
options(error = function() { traceback(2); dump.frames(to.file = TRUE) })
options(future.globals.maxSize = 300e9)   # adjust if on a different server
options(Seurat.object.validate = FALSE)
set.seed(1)

## ---------- Paths ----------
CHUNK_MAPPED_DIR   <- "SCT_chunks_mapped"  # contains *_mapped.rds
REF_INTEGRATED_RDS <- file.path(CHUNK_MAPPED_DIR, "integrated_reference_SCT.rds")
VIZ_RDS            <- file.path(CHUNK_MAPPED_DIR, "viz_re-metadata.rds")

OUT_MERGED_XEN     <- file.path(CHUNK_MAPPED_DIR, "merged_SCT_Xenium.rds")
OUT_MERGED_HD      <- file.path(CHUNK_MAPPED_DIR, "merged_SCT_VisiumHD.rds")
OUT_MERGED_V1      <- file.path(CHUNK_MAPPED_DIR, "merged_SCT_VisiumV1.rds")
OUT_MERGED_REF     <- file.path(CHUNK_MAPPED_DIR, "merged_SCT_RNAref.rds")               # separate
OUT_MERGED_ALL     <- file.path(CHUNK_MAPPED_DIR, "merged_SCT_allplatforms_noDR.rds")
OUT_MERGED_ALL_DR  <- file.path(CHUNK_MAPPED_DIR, "merged_SCT_allplatforms_withDR.rds")

## ---------- Settings ----------
PLAT_LEVELS         <- c("Xenium","VisiumHD","VisiumV1","RNAref")
REF_ID              <- "placental_ref_Jul_2025"
INCLUDE_REF_IN_FINAL<- FALSE     # keep RNAref out of final cross-platform merge
CKPT_EVERY          <- 2         # save checkpoint every N chunks
RESUME_XEN          <- TRUE
RESUME_HD           <- TRUE
RESUME_V1           <- TRUE
RESUME_REF          <- TRUE
DIMS                <- 1:30

## ---------- Helpers ----------
parse_platform <- function(s) dplyr::case_when(
  stringr::str_ends(s, stringr::regex("_xen$", ignore_case = TRUE)) ~ "Xenium",
  stringr::str_ends(s, stringr::regex("_hd$",  ignore_case = TRUE)) ~ "VisiumHD",
  stringr::str_ends(s, stringr::regex("_v1$",  ignore_case = TRUE)) ~ "VisiumV1",
  TRUE ~ NA_character_
)
manual_map <- c("HP1"="Xenium","HP13a"="VisiumV1","HT2b"="Xenium",
                "HT3"="Xenium","LP2"="Xenium","placental_ref_Jul_2025"="RNAref")

parse_manifest_row <- function(fp) {
  bn <- basename(fp)
  m  <- regexec("^sctrans_(.+)_chunk_(\\d{3})_SCTonly(?:_mapped)?\\.rds$", bn)
  mm <- regmatches(bn, m)[[1]]
  if (length(mm) == 0) stop("Filename not recognized: ", bn)
  sample <- mm[2]; chunk <- as.integer(mm[3])
  plat_guess <- parse_platform(sample)
  platform <- ifelse(!is.na(plat_guess), plat_guess,
                     ifelse(sample %in% names(manual_map), unname(manual_map[sample]), "VisiumV1"))
  tibble(sample = sample, chunk = chunk, file = fp, platform = platform)
}

# keep SCT assay as-is (Assay5-safe), drop heavy slots, sanitize meta
diet_SCT_only <- function(o) {
  DefaultAssay(o) <- "SCT"
  # drop non-SCT assays explicitly
  o@assays <- o@assays["SCT"]
  # drop heavy slots
  o@reductions <- list(); o@graphs <- list(); o@images <- list(); o@tools <- list()
  # sanitize meta
  md <- o@meta.data
  md <- md[match(colnames(o), rownames(md)), , drop = FALSE]
  rownames(md) <- colnames(o)
  keep <- rep(TRUE, ncol(md))
  for (i in seq_len(ncol(md))) {
    xi <- md[[i]]
    if (is.factor(xi)) { md[[i]] <- as.character(xi); next }
    if (is.list(xi) || (is.object(xi) && !is.atomic(xi))) {
      md[[i]] <- vapply(
        xi,
        function(z) if (is.null(z) || !length(z)) NA_character_ else paste(as.character(z), collapse = ";"),
        FUN.VALUE = character(1)
      ); next
    }
    if (!(is.character(xi) || is.numeric(xi) || is.integer(xi) || is.logical(xi))) keep[i] <- FALSE
  }
  o@meta.data <- md[, keep, drop = FALSE]
  if (any(is.na(Idents(o)))) Idents(o) <- factor(rep("unknown", ncol(o)))
  DefaultAssay(o) <- "SCT"
  o
}

# logger (console + file)
mklogger <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  function(...) {
    msg <- paste0("[", format(Sys.time(), "%F %T"), "] ", paste0(..., collapse = ""))
    cat(msg, "\n")
    cat(msg, "\n", file = path, append = TRUE)
  }
}

# streaming Seurat merge with checkpoints + resume
merge_platform <- function(df_subset, out_path, ckpt_every = 5, resume = TRUE) {
  if (!nrow(df_subset)) return(NULL)
  log_path <- sub("\\.rds$", ".log", out_path)
  logit    <- mklogger(log_path)
  hb_path  <- sub("\\.rds$", ".heartbeat", out_path)

  # resume?
  start_i <- 1L; acc <- NULL
  if (resume) {
    base  <- file_path_sans_ext(basename(out_path))
    ckpts <- list.files(dirname(out_path), paste0("^", base, "_ckpt_\\d{3}\\.rds$"), full.names = TRUE)
    if (length(ckpts)) {
      k <- max(as.integer(sub(".*_ckpt_(\\d{3})\\.rds$", "\\1", ckpts)))
      ck <- file.path(dirname(out_path), sprintf("%s_ckpt_%03d.rds", base, k))
      logit(sprintf("Auto-resume from checkpoint %03d (%s)", k, ck))
      acc <- readRDS(ck); DefaultAssay(acc) <- "SCT"
      start_i <- k + 1L
    }
  }

  plat <- as.character(unique(df_subset$platform))
  logit(sprintf("Merging platform: %s  (%d chunks)", plat, nrow(df_subset)))
  logit(sprintf("Writing logs -> %s", log_path))
  writeLines(format(Sys.time(), "%F %T"), hb_path)

  for (i in seq(from = start_i, to = nrow(df_subset))) {
    fp <- df_subset$file[i]; nm <- basename(fp)
    logit(sprintf("[%3d/%3d] load + diet %s", i, nrow(df_subset), nm))
    o <- readRDS(fp); DefaultAssay(o) <- "SCT"
    o <- diet_SCT_only(o)

    if (is.null(acc)) {
      acc <- o
    } else {
      acc <- merge(acc, y = o, merge.data = TRUE)
    }
    rm(o); gc()

    if (i %% ckpt_every == 0 || i == nrow(df_subset)) {
      writeLines(format(Sys.time(), "%F %T"), hb_path)
      ck_path <- sub("\\.rds$", sprintf("_ckpt_%03d.rds", i), out_path)
      saveRDS(acc, ck_path)
      logit(sprintf("…ckpt saved (%03d): %s  | cells=%d | genes=%d", i, ck_path, ncol(acc), nrow(acc)))
    }
  }

  DefaultAssay(acc) <- "SCT"
  saveRDS(acc, out_path)
  logit(sprintf("Saved -> %s | cells=%d | genes=%d", out_path, ncol(acc), nrow(acc)))
  acc
}

# DR projection helpers
project_ref_pca <- function(o, ref_integrated, dims_use = DIMS) {
  stopifnot("pca" %in% Reductions(ref_integrated))
  L <- as.matrix(Loadings(ref_integrated, reduction = "pca"))
  dims_use <- dims_use[dims_use <= ncol(L)]
  feats <- intersect(rownames(L), rownames(o))
  if (length(feats) < 50) stop("Too few overlapping genes with ref PCA (n=", length(feats), ").")
  X <- as.matrix(GetAssayData(o, "SCT", "data")[feats, , drop = FALSE])
  PC <- t(X) %*% L[feats, dims_use, drop = FALSE]
  colnames(PC) <- paste0("refPC_", seq_along(dims_use))
  rownames(PC) <- colnames(o)
  o[["ref.pca"]] <- CreateDimReducObject(embeddings = as.matrix(PC), key = "refPC_", assay = DefaultAssay(o))
  o
}
ensure_umap_model <- function(ref_integrated, dims_use = DIMS) {
  um_model <- tryCatch(ref_integrated[["umap"]]@misc$model, error = function(e) NULL)
  if (is.list(um_model) && !is.null(um_model$nn_index)) return(um_model)
  stopifnot(requireNamespace("uwot", quietly = TRUE))
  Xref <- Embeddings(ref_integrated, "pca")[, dims_use[dims_use <= ncol(Embeddings(ref_integrated,"pca"))], drop = FALSE]
  fit  <- uwot::umap(Xref, n_neighbors = 30, min_dist = 0.3, metric = "euclidean", ret_model = TRUE, verbose = TRUE)
  emb  <- as.matrix(fit$embedding); rownames(emb) <- rownames(Xref); colnames(emb) <- paste0("UMAP_", seq_len(ncol(emb)))
  ref_integrated[["umap"]] <- CreateDimReducObject(embeddings = emb, key = "UMAP_", assay = DefaultAssay(ref_integrated))
  ref_integrated[["umap"]]@misc$model <- fit
  fit
}
project_umap <- function(o, um_model, dims_use = DIMS) {
  stopifnot("ref.pca" %in% Reductions(o))
  coords <- try(Seurat::ProjectUMAP(umap.model = um_model,
                                    query = as.matrix(Embeddings(o, "ref.pca")[, dims_use, drop = FALSE])),
                silent = TRUE)
  UM <- if (!inherits(coords, "try-error")) as.matrix(coords) else {
    stopifnot(requireNamespace("uwot", quietly = TRUE))
    as.matrix(uwot::umap_transform(as.matrix(Embeddings(o,"ref.pca")[, dims_use, drop = FALSE]), model = um_model))
  }
  if (is.null(colnames(UM))) colnames(UM) <- paste0("UMAP_", seq_len(ncol(UM)))
  rownames(UM) <- colnames(o)
  o[["umap"]] <- CreateDimReducObject(embeddings = UM, key = "UMAP_", assay = DefaultAssay(o))
  o
}
prefix_meta_cols <- function(df, prefix) { names(df) <- paste0(prefix, names(df)); df }

## ---------- Build manifest from mapped files ----------
mapped_files <- list.files(CHUNK_MAPPED_DIR, pattern = "_mapped\\.rds$", full.names = TRUE)
stopifnot(length(mapped_files) > 0)

man <- dplyr::bind_rows(lapply(mapped_files, parse_manifest_row)) %>%
  dplyr::mutate(platform = factor(platform, levels = PLAT_LEVELS))
print(table(man$platform))

## ---------- Ensure reference is present and put FIRST in a saved manifest (for audit) ----------
has_ref <- any(man$sample == REF_ID)
cat("Reference present in manifest? ", has_ref, "\n")
man_ref_first <- if (has_ref) {
  bind_rows(
    filter(man, sample == REF_ID) %>% arrange(chunk),
    filter(man, sample != REF_ID)
  )
} else man

write.csv(man_ref_first, file.path(CHUNK_MAPPED_DIR, "chunk_manifest_refFirst.csv"), row.names = FALSE)

## ---------- Split by platform ----------
df_xen <- filter(man, platform == "Xenium")
df_hd  <- filter(man, platform == "VisiumHD")
df_v1  <- filter(man, platform == "VisiumV1")
df_ref <- filter(man, platform == "RNAref")

## ---------- Merge per platform (independent, with resume+ckpt) ----------
if (nrow(df_xen)) { obj_xen <- merge_platform(df_xen, OUT_MERGED_XEN, ckpt_every = CKPT_EVERY, resume = RESUME_XEN); rm(obj_xen); gc() }
if (nrow(df_hd))  { obj_hd  <- merge_platform(df_hd,  OUT_MERGED_HD,  ckpt_every = CKPT_EVERY, resume = RESUME_HD);  rm(obj_hd);  gc() }
if (nrow(df_v1))  { obj_v1  <- merge_platform(df_v1,  OUT_MERGED_V1,  ckpt_every = CKPT_EVERY, resume = RESUME_V1);  rm(obj_v1);  gc() }
if (nrow(df_ref)) { obj_ref <- merge_platform(df_ref, OUT_MERGED_REF, ckpt_every = CKPT_EVERY, resume = RESUME_REF); rm(obj_ref); gc() }

## ---------- Cross-platform merge (exclude RNAref unless toggled) ----------
to_merge <- c(OUT_MERGED_XEN, OUT_MERGED_HD, OUT_MERGED_V1)
if (INCLUDE_REF_IN_FINAL && file.exists(OUT_MERGED_REF)) to_merge <- c(to_merge, OUT_MERGED_REF)
to_merge <- to_merge[file.exists(to_merge)]
stopifnot(length(to_merge) > 0)

# streaming cross-platform merge
log_all <- mklogger(sub("\\.rds$", ".log", OUT_MERGED_ALL))
log_all(sprintf("Cross-platform merging %d platform objects -> %s", length(to_merge), OUT_MERGED_ALL))
acc_all <- NULL
for (p in to_merge) {
  o <- readRDS(p); DefaultAssay(o) <- "SCT"
  o <- diet_SCT_only(o)   # defensive: keep only SCT, drop reductions, sanitize meta
  if (is.null(acc_all)) acc_all <- o else acc_all <- merge(acc_all, y = o, merge.data = TRUE)
  rm(o); gc()
}
DefaultAssay(acc_all) <- "SCT"
saveRDS(acc_all, OUT_MERGED_ALL)
log_all(sprintf("Saved combined (no DR) -> %s | cells=%d | genes=%d", OUT_MERGED_ALL, ncol(acc_all), nrow(acc_all)))

## ---------- DR: project to integrated reference ----------
ref_integrated <- readRDS(REF_INTEGRATED_RDS)
acc_all <- project_ref_pca(acc_all, ref_integrated, dims_use = DIMS)
um_model <- ensure_umap_model(ref_integrated, dims_use = DIMS)
acc_all <- project_umap(acc_all, um_model, dims_use = DIMS)
saveRDS(acc_all, OUT_MERGED_ALL_DR)
log_all(sprintf("Saved combined with DR -> %s", OUT_MERGED_ALL_DR))

## ---------- Metadata graft: keep BOTH (orig_* from merged; viz_* from viz) ----------
orig_meta <- prefix_meta_cols(acc_all@meta.data, "orig_")
acc_all@meta.data <- orig_meta

viz <- readRDS(VIZ_RDS)
viz_meta <- viz@meta.data
common_cells <- intersect(rownames(acc_all@meta.data), rownames(viz_meta))
cat("Cells overlapping merged/viz:", length(common_cells), "of", ncol(acc_all), "\n")
stopifnot(length(common_cells) > 0)

viz_meta <- prefix_meta_cols(viz_meta, "viz_")
# attach columns (only common cells; NA elsewhere)
for (nm in colnames(viz_meta)) {
  acc_all[[nm]] <- NA
  acc_all[[nm]][common_cells] <- viz_meta[common_cells, nm]
}

if ("viz_sample_id_std" %in% colnames(acc_all@meta.data)) {
  Idents(acc_all) <- "viz_sample_id_std"
} else if ("viz_predicted.celltype" %in% colnames(acc_all@meta.data)) {
  Idents(acc_all) <- "viz_predicted.celltype"
}

saveRDS(acc_all, OUT_MERGED_ALL_DR)
cat("Final object saved with DR + prefixed metadata ->", OUT_MERGED_ALL_DR, "\n")

## ---------- Repro: write session info ----------
sink(file.path(CHUNK_MAPPED_DIR, "merge_sessionInfo.txt"))
print(sessionInfo())
sink()

cat("\nDone.\nLog files:\n",
    "  Xenium:   ", sub("\\.rds$", ".log", OUT_MERGED_XEN), "\n",
    "  VisiumHD: ", sub("\\.rds$", ".log", OUT_MERGED_HD), "\n",
    "  VisiumV1: ", sub("\\.rds$", ".log", OUT_MERGED_V1), "\n",
    "  RNAref:   ", sub("\\.rds$", ".log", OUT_MERGED_REF), "\n",
    "  Final:    ", sub("\\.rds$", ".log", OUT_MERGED_ALL), "\n", sep = "")
