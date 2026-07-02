# weee_harmony.R


set.seed(1) # set seed for reproducibility
options(future.globals.maxSize = 400e9) ## 300 gb ram (out of 440 on server 3; adjust if on different server)

library(Seurat)

seurat_all <- readRDS("/path/to/local_storage/raid2/project/visiumhd_mnp/all_spatial_seurat_objects.rds") ## this is  11.5 gb
# Check for duplicate names
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
table(names(seurat_all))

setwd("/path/to/local_storage/raid3/project/seurat_all_ref_merged")



suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(stringr)
})

## --- Config ---
annot_dir <- "/path/to/local_storage/raid3/project/seurat_all_ref_merged"
cols_want <- c("predicted.celltype", "predicted.celltype.score")
save_updated_objects <- TRUE
outdir <- file.path(annot_dir, "with_annotations")
if (save_updated_objects && !dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

stopifnot(exists("seurat_all"), is.list(seurat_all), length(seurat_all) > 0)

## --- Discover annotated files, FIXED parser ---
annot_files <- list.files(annot_dir, pattern = "^seurat_.*_annotated\\.rds$", full.names = TRUE)
stopifnot(length(annot_files) > 0)

get_sample_id_from_file <- function(fname) {
  bn <- basename(fname)
  bn <- sub("^seurat_", "", bn)                  # drop prefix
  bn <- sub("_annotated\\.rds$", "", bn)         # drop suffix (run separately!)
  bn
}
annot_sample_ids <- vapply(annot_files, get_sample_id_from_file, character(1))
annot_map <- setNames(annot_files, annot_sample_ids)

cat("Samples in seurat_all: ", length(seurat_all), "\n", sep = "")
cat("Annotated files parsed: ", length(annot_map), "\n", sep = "")

cat("Unmatched in seurat_all (if any):\n")
print(setdiff(names(seurat_all), names(annot_map)))
cat("Annotated without seurat_all (if any):\n")
print(setdiff(names(annot_map), names(seurat_all)))

## --- Barcode harmonization helpers ---
id_identity <- function(x) x
id_strip_lane <- function(x) gsub("-[0-9]+$", "", x)                       # drop -1/-2
id_drop_prefix <- function(x) sub("^.*_", "", x)                           # drop sample_ prefix
id_drop_prefix_and_lane <- function(x) gsub("-[0-9]+$", "", sub("^.*_", "", x))

choose_best_barcode_mapper <- function(src, ann) {
  strategies <- list(
    identity = id_identity,
    strip_lane = id_strip_lane,
    drop_prefix = id_drop_prefix,
    drop_prefix_and_lane = id_drop_prefix_and_lane
  )
  scores <- sapply(strategies, function(f) length(intersect(f(src), f(ann))))
  best_name <- names(which.max(scores))
  list(fun = strategies[[best_name]], name = best_name, overlap = max(scores))
}

## --- Safe column resolver (exact first, regex fallback) ---
find_col_safe <- function(df, target_name, fallback_regex = NULL) {
  if (target_name %in% colnames(df)) return(target_name)
  if (!is.null(fallback_regex)) {
    hits <- grep(fallback_regex, colnames(df), ignore.case = TRUE, value = TRUE)
    if (length(hits) > 0) return(hits[1])
  }
  NA_character_
}

## --- Transfer only the two columns, conflict-safe ---
transfer_two_cols <- function(target_obj, annotated_obj,
                              cols = c("predicted.celltype","predicted.celltype.score"),
                              prefix_conflicts = "RCTD_") {
  src_cells <- colnames(target_obj)
  ann_cells <- colnames(annotated_obj)

  strat <- choose_best_barcode_mapper(src_cells, ann_cells)
  mapfun <- strat$fun
  src_key <- mapfun(src_cells)
  ann_key <- mapfun(ann_cells)

  m <- match(src_key, ann_key)

  ann_md <- annotated_obj@meta.data
  ann_md$key <- ann_key

  col_ct    <- find_col_safe(ann_md, cols[1], fallback_regex = "predicted.*cell.?type")
  col_score <- find_col_safe(ann_md, cols[2], fallback_regex = "predicted.*(score|prob)")

  out_df <- data.frame(row.names = src_cells)

  if (!is.na(col_ct)) {
    v <- rep(NA, length(src_cells)); hit <- !is.na(m)
    if (any(hit)) v[hit] <- ann_md[[col_ct]][m[hit]]
    nm <- cols[1]; if (nm %in% colnames(target_obj@meta.data)) nm <- paste0(prefix_conflicts, nm)
    out_df[[nm]] <- v
  } else {
    warning("Missing predicted.celltype in annotated meta.")
  }

  if (!is.na(col_score)) {
    v <- rep(NA, length(src_cells)); hit <- !is.na(m)
    if (any(hit)) v[hit] <- ann_md[[col_score]][m[hit]]
    nm <- cols[2]; if (nm %in% colnames(target_obj@meta.data)) nm <- paste0(prefix_conflicts, nm)
    out_df[[nm]] <- v
  } else {
    warning("Missing predicted.celltype.score in annotated meta.")
  }

  target_obj <- AddMetaData(target_obj, out_df)

  list(
    object = target_obj,
    stats = list(
      strategy = strat$name,
      overlap  = strat$overlap,
      n_target = length(src_cells),
      n_ann    = length(ann_cells),
      matched  = sum(!is.na(m)),
      added    = colnames(out_df)
    )
  )
}

## --- Apply to the intersection only ---
transfer_log <- list()
common <- intersect(names(seurat_all), names(annot_map))
cat("Transferring for", length(common), "matched samples.\n")

for (sid in common) {
  cat("\n== Sample:", sid, "==\n")
  ann_obj <- readRDS(annot_map[[sid]])
  res <- transfer_two_cols(seurat_all[[sid]], ann_obj, cols = cols_want, prefix_conflicts = "RCTD_")
  seurat_all[[sid]] <- res$object
  st <- res$stats
  cat(sprintf("  Strategy: %s | overlap=%d | matched=%d/%d | added: %s\n",
              st$strategy, st$overlap, st$matched, st$n_target, paste(st$added, collapse=", ")))
  transfer_log[[sid]] <- data.frame(sample = sid, matched = st$matched, n_target = st$n_target,
                                    strategy = st$strategy, added = paste(st$added, collapse=";"))
  if (save_updated_objects) {
    saveRDS(seurat_all[[sid]], file.path(outdir, paste0("seurat_", sid, "_withPredCols.rds")))
  }
}

## --- Summary ---
if (length(transfer_log) > 0) {
  transfer_summary <- bind_rows(transfer_log)
  write.csv(transfer_summary, file.path(outdir, "predicted_celltype_transfer_summary.csv"), row.names = FALSE)
  cat("\nSummary saved ->", file.path(outdir, "predicted_celltype_transfer_summary.csv"), "\n")
}
cat("\nDone.\n")





## 


