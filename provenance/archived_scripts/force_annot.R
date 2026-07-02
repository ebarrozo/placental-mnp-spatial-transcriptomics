## force_annot.R


setwd("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/SCT_chunks_mapped/fully_annotated")
obj <-readRDS("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/SCT_chunks_mapped/final_sct_atlas.rds") # 7.5 gb object with SCT data, all 3,969,532 cells, and umap


##1.5.26

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(RcppAnnoy)
})

rescue_lowconf_labels_v2 <- function(
  obj,
  group_col = "sample_id_std",
  celltype_col = "predicted.celltype",
  lowconf_tags = c("RCTD_lowconf","RCTD-lowconf"),
  bad_label_regex = "^UMAP_",      # <- key fix
  treat_na_as_lowconf = TRUE,
  exclude_groups = NULL,           # e.g., c("placental_ref")
  reduction = "umap",
  dims = NULL,
  K = 75,
  n_trees = 50,
  support_min = 0.40,
  assign_all_remaining = FALSE,
  chunk = 20000,
  min_hi = 200,
  max_iter = 2,
  seed_support_min = 0.60
) {
  stopifnot(group_col %in% colnames(obj@meta.data))
  stopifnot(celltype_col %in% colnames(obj@meta.data))
  stopifnot(reduction %in% Reductions(obj))

  Eall <- Embeddings(obj, reduction)
  if (is.null(dims)) dims <- seq_len(ncol(Eall))
  dims <- dims[dims <= ncol(Eall)]
  Eall <- Eall[, dims, drop = FALSE]

  ct0 <- obj@meta.data[[celltype_col]]
  grp <- obj@meta.data[[group_col]]

  # preserve original only once (do not overwrite on subsequent passes)
  orig_col <- paste0(celltype_col, "_orig")
  if (!orig_col %in% colnames(obj@meta.data)) {
    obj@meta.data[[orig_col]] <- ct0
  }

  new_lab <- ct0
  support <- rep(NA_real_, length(ct0))

  is_bad_label <- function(x) {
    x <- as.character(x)
    bad <- rep(FALSE, length(x))
    bad <- bad | is.na(x) | x == ""
    bad <- bad | (x %in% lowconf_tags)
    bad <- bad | (!is.null(bad_label_regex) & grepl(bad_label_regex, x))
    if (!treat_na_as_lowconf) bad <- bad & !(is.na(x) | x == "")
    bad
  }

  groups <- unique(as.character(grp))
  groups <- groups[!is.na(groups) & groups != ""]
  if (!is.null(exclude_groups)) groups <- setdiff(groups, exclude_groups)

  message("Groups: ", length(groups), " | reduction=", reduction, " dims=", paste(range(dims), collapse="-"))

  for (g in groups) {
    g_idx <- which(as.character(grp) == g)
    if (!length(g_idx)) next

    Eg <- Eall[g_idx, , drop = FALSE]

    iter <- 1
    changed_any <- TRUE

    while (iter <= max_iter && changed_any) {
      changed_any <- FALSE

      # define hi/lo *based on current* labels (so iterative expansion is real)
      lab_g <- new_lab[g_idx]
      bad_g <- is_bad_label(lab_g)

      hi_local_use <- which(!bad_g)
      lo_local     <- which(bad_g)  # includes lowconf + UMAP_* + NA/blank (if enabled)

      if (length(lo_local) == 0) break
      if (length(hi_local_use) < min_hi) {
        message("Group ", g, " iter ", iter, " skipped: hi_local_use=", length(hi_local_use),
                " < min_hi=", min_hi, " | lo=", length(lo_local))
        break
      }

      hi_labels_use <- as.character(lab_g[hi_local_use])
      # sanity: ensure no bad labels slipped in
      if (any(is_bad_label(hi_labels_use))) {
        stop("Bad labels present in hi seeds for group ", g, " (check bad_label_regex / lowconf_tags).")
      }

      # build annoy index over hi seeds (0-based ids)
      ann <- new(AnnoyEuclidean, ncol(Eg))
      for (i in seq_along(hi_local_use)) {
        ann$addItem(i - 1L, as.numeric(Eg[hi_local_use[i], ]))
      }
      ann$build(n_trees)

      # rescue lows in chunks
      seq_starts <- seq(1, length(lo_local), by = chunk)
      for (s in seq_starts) {
        sel <- lo_local[s:min(s + chunk - 1, length(lo_local))]
        for (jj in sel) {
          v <- as.numeric(Eg[jj, ])

          nn <- tryCatch(
            ann$getNNsByVector(v, K, -1, TRUE),
            error = function(e) ann$getNNsByVectorList(v, K, -1, TRUE)
          )

          nn_i <- nn$item      # 0-based indices into hi_local_use/hi_labels_use
          nn_d <- nn$distance

          w <- 1 / (nn_d + 1e-6)
          labs <- hi_labels_use[nn_i + 1L]

          tab <- tapply(w, labs, sum)
          best <- names(tab)[which.max(tab)]
          sup  <- max(tab) / sum(tab)

          glob_j <- g_idx[jj]
          support[glob_j] <- sup

          if (!is.na(sup) && (sup >= support_min || assign_all_remaining)) {
            if (!identical(new_lab[glob_j], best)) changed_any <- TRUE
            new_lab[glob_j] <- best
          }
        }
      }

      message("Group ", g, " iter ", iter,
              " done | lo=", length(lo_local),
              " | hi_use=", length(hi_local_use),
              " | labels=", length(unique(hi_labels_use)))
      iter <- iter + 1
    }
  }

  # write back
  obj@meta.data[[celltype_col]] <- new_lab
  obj@meta.data[[paste0(celltype_col, "_support")]] <- support

  # make support usable downstream: set support=1 for cells we never had to rescue
  sup_col <- paste0(celltype_col, "_support")
  bad_final <- is_bad_label(obj@meta.data[[celltype_col]])
  obj@meta.data[[sup_col]][is.na(obj@meta.data[[sup_col]]) & !bad_final] <- 1

  obj
}


BAD <- c("UMAP_1","UMAP_2")
LOWCONF_TAGS <- c("RCTD_lowconf","RCTD-lowconf")

obj$predicted.celltype[obj$predicted.celltype %in% BAD] <- LOWCONF_TAGS[1]

# re-run rescue (2-pass strategy like you did)
obj <- rescue_lowconf_labels_v2(
  obj,
  group_col = "sample_id_std",
  reduction = "umap",
  K = 50,
  support_min = 0.60,
  max_iter = 1
)

obj <- rescue_lowconf_labels_v2(
  obj,
  group_col = "sample_id_std",
  reduction = "umap",
  K = 75,
  support_min = 0.40,
  max_iter = 2,
  seed_support_min = 0.60
)

# confirm UMAP labels are gone
table(obj$predicted.celltype %in% BAD, useNA="ifany")
table(obj$predicted.celltype %in% LOWCONF_TAGS, useNA="ifany")

obj <- rescue_lowconf_labels_v2(
  obj,
  group_col = "sample_id_std",
  exclude_groups = "placental_ref",
  reduction = "umap",
  K = 75,
  support_min = 0.40,
  max_iter = 2
)


BAD <- c("UMAP_1","UMAP_2")
LOWCONF_TAGS <- c("RCTD_lowconf","RCTD-lowconf")

obj$predicted.celltype[obj$predicted.celltype %in% BAD] <- LOWCONF_TAGS[1]

table(obj$predicted.celltype %in% c("UMAP_1","UMAP_2"), useNA="ifany")
obj@meta.data |>
  dplyr::count(sample_id_std, is_lowconf = predicted.celltype %in% c("RCTD_lowconf","RCTD-lowconf")) |>
  tidyr::pivot_wider(names_from = is_lowconf, values_from = n, values_fill = 0) |>
  dplyr::mutate(frac_lowconf = `TRUE` / (`TRUE` + `FALSE`)) |>
  dplyr::arrange(dplyr::desc(frac_lowconf)) |>
  print(n = Inf)


table(obj$predicted.celltype %in% LOWCONF_TAGS, useNA="ifany")
summary(obj$predicted.celltype_support)

# After Pass 2: lowconf = 62,660 out of 3,789,532 total → 1.65% of all transcriptomes still lowconf.
# Relative to the original lowconf pool (2,492,192): remaining 62,660 / 2,492,192 = 2.5%. So you rescued ~97.5% of previously lowconf cells.
setwd("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/SCT_chunks_mapped/fully_annotated")
saveRDS(obj, file = "/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/integrated_SCT_final_rescue-RCTD-lowconf_pass2.rds")


# Pass 3 (fill remaining, keep support for downstream filtering)
obj <- rescue_lowconf_labels_v2(
    obj,
    group_col = "sample_id_std",
    exclude_groups = "placental_ref",
    reduction = "umap",
    K = 75,
    support_min = 0.00,
    assign_all_remaining = TRUE,
    max_iter = 1
)
table(obj$predicted.celltype %in% c("UMAP_1","UMAP_2"), useNA="ifany")
obj@meta.data |>
  dplyr::count(sample_id_std, is_lowconf = predicted.celltype %in% c("RCTD_lowconf","RCTD-lowconf")) |>
  tidyr::pivot_wider(names_from = is_lowconf, values_from = n, values_fill = 0) |>
  dplyr::mutate(frac_lowconf = `TRUE` / (`TRUE` + `FALSE`)) |>
  dplyr::arrange(dplyr::desc(frac_lowconf)) |>
  print(n = Inf)
  
table(obj$predicted.celltype %in% LOWCONF_TAGS, useNA="ifany")
summary(obj$predicted.celltype_support)
saveRDS(obj, file = "/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/integrated_SCT_final_rescue-RCTD-lowconf_all.rds")


cp /path/to/collaborator_home/MNPs/Data/mnp_obj_w_nn_predicted_celltype.rds /path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/

mypal3 <-get_palette("igv",40)

## 4) Sanity + quick plot (palette later)
Idents(obj) <- "predicted.celltype"
levels(obj)

t1<- table(obj$predicted.celltype, obj$sample_id_std)
write.csv(t1, "obj-predicted.celltype_sample_id_std.csv", row.names = T)
t1


stopifnot(identical(rownames(Embeddings(obj, "umap")), colnames(obj)))
DimPlot(obj, reduction = "umap", group.by = "predicted.celltype",
        alpha = 0.1, raster = TRUE, label = FALSE) + NoAxes()


p2 <- DimPlot(
  obj, reduction = "umap",
  group.by = "predicted.celltype", alpha = 0.1, raster = TRUE, label = FALSE,
  cols = mypal3
) + ggtitle(sprintf("Predicted Cell Types (%s Cells)", format(ncol(obj), big.mark=",", cols = mypal3))) + ggtitle("Predicted Cell Types (3,969,532 Cells)") + NoAxes()
p2



Idents(obj) <- "predicted.celltype"
p2<-DimPlot(obj, reduction = "umap", group.by = "predicted.celltype", alpha = 0.1, label = FALSE,raster=T, cols = mypal3) + ggtitle("Predicted Cell Types (3,969,532 Cells)") + NoAxes()
p2
ggsave("obj_predicted.celltype.png", plot = p2, device = "png", width = 8, height = 5, units = "in", dpi = 1200)
## Identify the reference cells (prefer barcodes from `ref`; fallback to the known tag)
ref_cells <- if (exists("ref")) {
  intersect(colnames(obj), colnames(ref))          # barcodes from ref present in obj
} else {
  rownames(obj@meta.data)[obj$sample_id_std %in% c("placental_ref_Jul_2025")]
}

stopifnot(length(ref_cells) > 0)

## Ensure the columns exist and are character (avoid factor overhead on 4M rows)
if (!"platform" %in% colnames(obj@meta.data))       obj@meta.data$platform       <- NA_character_
if (!"sample_id_std" %in% colnames(obj@meta.data))  obj@meta.data$sample_id_std  <- NA_character_

obj@meta.data$platform      <- as.character(obj@meta.data$platform)
obj@meta.data$sample_id_std <- as.character(obj@meta.data$sample_id_std)

## Set both fields for reference cells only
obj@meta.data[ref_cells, "platform"]      <- "reference_dataset"
obj@meta.data[ref_cells, "sample_id_std"] <- "reference_dataset"

## Quick sanity
c(length_ref = length(ref_cells),
  n_platform_ref = sum(obj@meta.data[ref_cells, "platform"]      == "reference_dataset"),
  n_sampleid_ref = sum(obj@meta.data[ref_cells, "sample_id_std"] == "reference_dataset"))


p2<-DimPlot(obj, reduction = "umap", group.by = "platform", alpha = 0.1, label = FALSE,raster=T, cols = mypal3) + ggtitle("Platform (3,969,532 cells)") + NoAxes()
ggsave("obj_platform.png", plot = p2, device = "png", width = 8, height = 5, units = "in", dpi = 1200)
p2<-DimPlot(obj, reduction = "umap", split.by = "platform", alpha = 0.1, label = FALSE,raster=T, cols = mypal3) + ggtitle("Platform (3,969,532 cells)") + NoAxes()
ggsave("obj_split.by_platform.png", plot = p2, device = "png", width = 16, height = 5, units = "in", dpi = 1200)

p2<-DimPlot(obj, reduction = "umap", group.by = "sample_id_std", alpha = 0.1, label = FALSE,raster=T, cols = mypal3) + ggtitle("Sample ID_std (3,969,532 cells)") + NoAxes()
ggsave("obj_study_id_std.png", plot = p2, device = "png", width = 8, height = 5, units = "in", dpi = 1200)
p2<-DimPlot(obj, reduction = "umap", group.by = "study_id", alpha = 0.1, label = FALSE,raster=T, cols = mypal3) + ggtitle("Sample ID (3,969,532 cells)") + NoAxes()
ggsave("obj_study_id.png", plot = p2, device = "png", width = 8, height = 5, units = "in", dpi = 1200)

Idents(obj) <- "predicted.celltype"

n_panels <- length(unique(obj$predicted.celltype))
ncol_pl  <- 4
nrow_pl  <- ceiling(n_panels / ncol_pl)

p2 <- DimPlot(
  obj, reduction = "umap",
  #group.by = "predicted.celltype",
  split.by = "predicted.celltype",
  ncol = ncol_pl,           # <-- number of columns
  alpha = 0.1, label = FALSE, raster = TRUE,
  cols = mypal3
) + ggtitle(sprintf("Predicted Cell Types (3,969,532 cells)", format(ncol(obj), big.mark=","))) +
  NoAxes()

# make the export size scale with the grid so panels don’t get squished
ggsave(
  "obj_split.by_predicted.celltype.png", plot = p2, device = "png",
  width  = 3.0 * ncol_pl,     # tweak panel size multiplier as you like
  height = 3.0 * nrow_pl,
  units = "in", dpi = 1200, limitsize = FALSE
)
p2<-DimPlot(obj, reduction = "umap", group.by = "predicted.celltype", alpha = 0.1, label = FALSE,raster=T, cols = mypal3) + ggtitle("Predicted Cell Types (3,969,532 Cells)") + NoAxes()
p2
ggsave("obj_predicted.celltype.png", plot = p2, device = "png", width = 8, height = 5, units = "in", dpi = 1200)



## save
BASE <- "/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/SCT_chunks_mapped"
ANN_DIR <- file.path(BASE, "fully_annotated")
dir.create(ANN_DIR, showWarnings = FALSE, recursive = TRUE)

saveRDS(obj, file.path(ANN_DIR, "final_sct_atlas_fully_annotated.rds"))



setwd("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/SCT_chunks_mapped/fully_annotated")
saveRDS(obj, file = "/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/integrated_SCT_final_forced.annot.rds")




t1<- table(obj$predicted.celltype, obj$sample_id_std)
write.csv(t1, "obj-predicted.celltype_sample_id_std.csv", row.names = T)


###############################################################################################
##1.2.26    This is the one that ran, but named cells as UMAP_1 and UMAP_2 # updated above 
###############################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(RcppAnnoy)
})

rescue_lowconf_labels <- function(
  obj,
  group_col = "sample_id_std",
  celltype_col = "predicted.celltype",
  lowconf_tags = c("RCTD_lowconf","RCTD-lowconf"),
  reduction = "umap",
  dims = NULL,                 # NULL => all dims in reduction
  K = 75,
  n_trees = 50,
  support_min = 0.40,
  assign_all_remaining = FALSE,
  chunk = 20000,
  min_hi = 200,                # skip groups with too few highconf seeds
  max_iter = 2,                # iteratively expand seeds using rescued labels
  seed_support_min = 0.60      # only promote rescued labels as seeds if >= this
) {
  stopifnot(group_col %in% colnames(obj@meta.data))
  stopifnot(celltype_col %in% colnames(obj@meta.data))
  stopifnot(reduction %in% Reductions(obj))

  Eall <- Embeddings(obj, reduction)
  if (is.null(dims)) dims <- seq_len(ncol(Eall))
  dims <- dims[dims <= ncol(Eall)]
  Eall <- Eall[, dims, drop = FALSE]

  ct0 <- obj@meta.data[[celltype_col]]
  grp <- obj@meta.data[[group_col]]

  is_low0 <- !is.na(ct0) & ct0 %in% lowconf_tags
  new_lab <- ct0
  support <- rep(NA_real_, length(ct0))

  # keep provenance
  obj[[paste0(celltype_col, "_orig")]] <- ct0

  groups <- unique(as.character(grp))
  groups <- groups[!is.na(groups) & groups != ""]

  message("Groups: ", length(groups), " | reduction=", reduction, " dims=", paste(range(dims), collapse="-"))

  for (g in groups) {
    g_idx <- which(as.character(grp) == g)
    if (!length(g_idx)) next

    # iterative loop: expand seed set with high-support rescues
    hi_seed <- which(!is_low0 & !is.na(ct0) & ct0 != "" & as.character(grp) == g)

    lo_idx <- which(is_low0 & as.character(grp) == g)

    if (length(lo_idx) == 0) next
    if (length(hi_seed) < min_hi) {
      message("Group ", g, " skipped: hi_seed=", length(hi_seed), " < min_hi=", min_hi,
              " | lo=", length(lo_idx))
      next
    }

    # work in local indices for speed
    Eg <- Eall[g_idx, , drop = FALSE]

    # map global -> local row positions
    g_pos <- setNames(seq_along(g_idx), g_idx)
    hi_local <- unname(g_pos[as.character(hi_seed)])
    lo_local <- unname(g_pos[as.character(lo_idx)])

    # if any NA (should not), drop
    hi_local <- hi_local[!is.na(hi_local)]
    lo_local <- lo_local[!is.na(lo_local)]

    # iterative seed expansion
    iter <- 1
    changed_any <- TRUE
    while (iter <= max_iter && changed_any) {
      changed_any <- FALSE

      hi_labels <- new_lab[g_idx[hi_local]]
      keep_hi <- !is.na(hi_labels) & hi_labels != "" & !(hi_labels %in% lowconf_tags)
      hi_local_use <- hi_local[keep_hi]
      hi_labels_use <- hi_labels[keep_hi]

      if (length(hi_local_use) < min_hi) {
        message("Group ", g, " iter ", iter, " stopped: hi_local_use=", length(hi_local_use),
                " < min_hi=", min_hi)
        break
      }

      # build annoy over hi seeds (0-based item ids)
      ann <- new(AnnoyEuclidean, ncol(Eg))
      for (i in seq_along(hi_local_use)) {
        ann$addItem(i - 1L, as.numeric(Eg[hi_local_use[i], ]))
      }
      ann$build(n_trees)

      # rescue lows in chunks
      seq_starts <- seq(1, length(lo_local), by = chunk)
      for (s in seq_starts) {
        sel <- lo_local[s:min(s + chunk - 1, length(lo_local))]
        for (jj in sel) {
          v <- as.numeric(Eg[jj, ])
          nn <- ann$getNNsByVectorList(v, K, -1, TRUE)  # <-- correct API :contentReference[oaicite:2]{index=2}

          nn_i <- nn$item      # 0-based indices into hi_local_use/hi_labels_use
          nn_d <- nn$distance

          # weighted vote
          w <- 1 / (nn_d + 1e-6)
          labs <- hi_labels_use[nn_i + 1L]

          tab <- tapply(w, labs, sum)
          best <- names(tab)[which.max(tab)]
          sup  <- max(tab) / sum(tab)

          glob_j <- g_idx[jj]
          support[glob_j] <- sup

          if (!is.na(sup) && (sup >= support_min || assign_all_remaining)) {
            if (!identical(new_lab[glob_j], best)) changed_any <- TRUE
            new_lab[glob_j] <- best
          }
        }
      }

      # promote strong rescues into seed set for next iter
      if (iter < max_iter) {
        promotable <- g_idx[lo_local][!is.na(support[g_idx[lo_local]]) & support[g_idx[lo_local]] >= seed_support_min &
                                        !is.na(new_lab[g_idx[lo_local]]) & !(new_lab[g_idx[lo_local]] %in% lowconf_tags)]
        if (length(promotable)) {
          # mark them as no longer lowconf for seed purposes
          is_low0[promotable] <- FALSE
        }
      }

      message("Group ", g, " iter ", iter,
              " done | lo=", length(lo_local),
              " | hi_use=", length(hi_local_use),
              " | labels=", length(unique(hi_labels_use)))
      iter <- iter + 1
    }
  }

  # write back “final” labels into canonical column
  obj@meta.data[[celltype_col]] <- new_lab
  obj@meta.data[[paste0(celltype_col, "_support")]] <- support

  obj
}




LOWCONF_TAGS <- c("RCTD_lowconf","RCTD-lowconf")

# Pass 1
obj <- rescue_lowconf_labels(
  obj,
  group_col = "sample_id_std",
  reduction = "umap",
  K = 50,
  support_min = 0.60,
  assign_all_remaining = FALSE,
  max_iter = 1
)

# Pass 2 (expand)
obj <- rescue_lowconf_labels(
  obj,
  group_col = "sample_id_std",
  reduction = "umap",
  K = 75,
  support_min = 0.40,
  assign_all_remaining = FALSE,
  max_iter = 2,
  seed_support_min = 0.60
)


## only 
LOWCONF_TAGS <- c("RCTD_lowconf","RCTD-lowconf")
print(obj@meta.data |>
  dplyr::count(sample_id_std, is_lowconf = predicted.celltype %in% LOWCONF_TAGS) |>
  tidyr::pivot_wider(names_from = is_lowconf, values_from = n, values_fill = 0),
  n = Inf
)

sample_id_std `FALSE` `TRUE`
   <fct>           <int>  <int>
 1 HP1_xen         42848    356
 2 HP13a_v1          369      0
 3 HP13b_v1         1275     48
 4 HP13b_xen      196580     68
 5 HP2_hd         441572   1941
 6 HP2_xen        130447    126
 7 HT1_v1           1826      0
 8 HT1_xen        196637     11
 9 HT2_hd         487598  14146
10 HT2_v1           1336     11
11 HT2a_xen       295939      3
12 HT2b_xen        99192      9
13 HT3_hd         367294  42969
14 HT3_xen        108269    714
15 LP1_hd         578449     89
16 LP1_xen        186082   1331
17 LP2_xen        163048    451
18 LT1_xen        181648    247
19 LT2_xen         94576    102
20 NC1a_v1           548      0
21 NC1b_v1           454     16
22 NC1c_v1           967      0
23 NC2_v1           2027      0
24 NC3_v1           2601      4
25 NC4_v1           1894      0
26 NC4_xen        134965     18
27 ND5_v1           1505      0
28 ND6_v1            515      0
29 ND7_v1            341      0
30 SP10_v1           918      0
31 SP11_v1          1547      0
32 SP12_v1           669      0
33 SP8_v1           1375      0
34 SP9_v1           1561      0

# After Pass 2: lowconf = 62,660 out of 3,789,532 total → 1.65% of all transcriptomes still lowconf.
# Relative to the original lowconf pool (2,492,192): remaining 62,660 / 2,492,192 = 2.5%. So you rescued ~97.5% of previously lowconf cells.


# Pass 3 (fill remaining, keep support for downstream filtering)
#obj <- rescue_lowconf_labels(
#  obj,
#  group_col = "sample_id_std",
#  reduction = "umap",
#  K = 100,
#  support_min = 0.00,
#  assign_all_remaining = TRUE,
#  max_iter = 1
#)

table(obj$predicted.celltype %in% LOWCONF_TAGS, useNA="ifany")
summary(obj$predicted.celltype_support)





# ===================== CONFIG =====================
BASE <- "/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/SCT_chunks_mapped"

ANN_DIR <- file.path(BASE, "fully_annotated")
dir.create(ANN_DIR, showWarnings = FALSE, recursive = TRUE)

# save the updated annotations once (run this right after rescue passes)
OBJ_RDS <- file.path(ANN_DIR, "final_sct_atlas_fully_annotated.rds")
# saveRDS(obj, OBJ_RDS)

OUTDIR <- file.path(ANN_DIR, "out_mecp2_prelim")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

PLATFORMS_KEEP <- c("Xenium", "VisiumHD", "VisiumV1")

CELLTYPE_COL <- "predicted.celltype"
STUDY_COL    <- "study_id"
PE_COL       <- "Preeclampsia_bin"





obj <- rescue_lowconf_labels(
  obj,
  group_col = "sample_id_std",
  reduction = "umap",
  K = 150,
  support_min = 0.00,
  assign_all_remaining = TRUE,     # hard-fill everything
  max_iter = 1
)


obj$predicted.celltype_final <- obj$predicted.celltype
obj$predicted.celltype_final[!is.na(obj$predicted.celltype_support) &
                             obj$predicted.celltype_support < 0.40] <- NA

table(obj$predicted.celltype %in% LOWCONF_TAGS, useNA="ifany")


## save
BASE <- "/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/SCT_chunks_mapped"
ANN_DIR <- file.path(BASE, "fully_annotated")
dir.create(ANN_DIR, showWarnings = FALSE, recursive = TRUE)

saveRDS(obj, file.path(ANN_DIR, "final_sct_atlas_fully_annotated.rds"))



setwd("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/SCT_chunks_mapped/fully_annotated")
saveRDS(obj, file = "/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/integrated_SCT_final_forced.annot.rds")




t1<- table(obj$predicted.celltype, obj$sample_id_std)
write.csv(t1, "obj-predicted.celltype_sample_id_std.csv", row.names = T)

Idents(obj) <- "predicted.celltype"
levels(obj)
stopifnot(identical(rownames(Embeddings(obj, "umap")), colnames(obj)))
DimPlot(obj, reduction = "umap", group.by = "predicted.celltype",
  alpha = 0.1, raster = TRUE, label = FALSE) + NoAxes()
mypal3 <-get_palette("igv",40)




