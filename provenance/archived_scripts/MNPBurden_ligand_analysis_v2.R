## MNPBurden_ligand_analysis_v2.R
## server3
## Enrico Barrozo, Ph.D.
## Assistant Professor, Department of OB/GYN, Baylor College of Medicine and Texas Children’s Hospital
# 2.23.26

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(RcppAnnoy)
})
suppressPackageStartupMessages({
    library(dplyr); library(stringr); library(tibble); library(Seurat)
})
library(ggplot2)
library(cowplot)
library(patchwork)
library(ggrepel)
library(ggpubr)
library(ggplot2)
setwd("/path/to/local_home/visiumhd_mnp/results/audra")
obj <-readRDS("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/SCT_chunks_mapped/audra/mnp_obj_w_nn_predicted_celltype.rds") # 7.1 gb object with SCT data [no counts], all 3,969,532 cells, and umap


mypal3 <-get_palette("igv",40)

## 4) Sanity + quick plot (palette later)
Idents(obj) <- "nn.predicted.celltype"
levels(obj)

t1<- table(obj$nn.predicted.celltype, obj$sample_id_std)
write.csv(t1, "obj-nn.predicted.celltype_sample_id_std.csv", row.names = T)
t1
Idents(obj) <- "predicted.celltype"
levels(obj)
t1<- table(obj$nn.predicted.celltype, obj$sample_id_std)
write.csv(t1, "obj-predicted.celltype_sample_id_std.csv", row.names = T)
#Idents(obj) <- "cell.type"
#levels(obj)
t1<- table(obj$cell_type, obj$sample_id_std)
write.csv(t1, "obj-cell_type_by_sample_id_std.csv", row.names = T)
t1
Idents(obj) <- "nn.predicted.celltype"
levels(obj)

stopifnot(identical(rownames(Embeddings(obj, "umap")), colnames(obj)))
DimPlot(obj, reduction = "umap", group.by = "nn.predicted.celltype",
        alpha = 0.1, raster = TRUE, label = FALSE) + NoAxes()


p2 <- DimPlot(
  obj, reduction = "umap",
  group.by = "nn.predicted.celltype", alpha = 0.1, raster = TRUE, label = FALSE,
  cols = mypal3
) + ggtitle(sprintf("Predicted Cell Types (%s Cells)", format(ncol(obj), big.mark=",", cols = mypal3))) + ggtitle("Predicted Cell Types (3,969,532 Cells)") + NoAxes()
p2



Idents(obj) <- "nn.predicted.celltype"
p2<-DimPlot(obj, reduction = "umap", group.by = "nn.predicted.celltype", alpha = 0.1, label = FALSE,raster=T, cols = mypal3) + ggtitle("Predicted Cell Types (3,969,532 Cells)") + NoAxes()
p2
ggsave("obj_nn.predicted.celltype.png", plot = p2, device = "png", width = 8, height = 5, units = "in", dpi = 1200)



# ------------------------------------------------------------
# Add placental_ref (sample_id_std) cell_type labels into nn.predicted.celltype
# ------------------------------------------------------------

library(Seurat)
library(dplyr)

md <- obj[[]]  # metadata as data.frame

# sanity checks
stopifnot(all(c("cell_type", "nn.predicted.celltype", "sample_id_std") %in% colnames(md)))

# define which samples are "placental_ref"
# OPTION 1 (recommended): sample_id_std is literally "placental_ref" for those cells
is_ref <- md$sample_id_std %in% "placental_ref"

# If instead sample_id_std contains multiple IDs and you have a vector of ref sample IDs, use:
# ref_ids <- c("placental_ref_01", "placental_ref_02")  # <-- your IDs here
# is_ref <- md$sample_id_std %in% ref_ids

# -------------------------
# A) Overwrite nn.predicted.celltype for placental_ref cells (even if currently non-NA)
# -------------------------
obj$nn.predicted.celltype <- as.character(obj$nn.predicted.celltype)
obj$cell_type             <- as.character(obj$cell_type)

obj$nn.predicted.celltype[is_ref] <- obj$cell_type[is_ref]

# Optional: make it a factor after editing (keeps levels clean)
obj$nn.predicted.celltype <- factor(obj$nn.predicted.celltype)

# -------------------------
# Quick QA
# -------------------------
cat("placental_ref cells:", sum(is_ref), "\n")
cat("placental_ref nn.predicted.celltype NA after:", sum(is.na(obj$nn.predicted.celltype[is_ref])), "\n")
print(table(
  cell_type = obj$cell_type[is_ref],
  nn_pred   = obj$nn.predicted.celltype[is_ref],
  useNA = "ifany"
))



Idents(obj) <- "nn.predicted.celltype"
p2<-DimPlot(obj, reduction = "umap", group.by = "nn.predicted.celltype", alpha = 0.1, label = FALSE,raster=T, cols = mypal3) + ggtitle("Predicted Cell Types (3,969,532 Cells)") + NoAxes()
p2
ggsave("obj_nn.predicted.celltype.png", plot = p2, device = "png", width = 8, height = 5, units = "in", dpi = 1200)




# ------------------------------------------------------------
# Reorder nn.predicted.celltype levels (group similar types) + keep RCTD_lowconf last
# ------------------------------------------------------------
library(Seurat)
library(ggplot2)

Idents(obj) <- "nn.predicted.celltype"

# Current levels (what exists in the object)
lvl <- levels(factor(obj$nn.predicted.celltype))

# Desired order (edit this list to taste)
desired <- c(
  # Immune (lymphoid)
  "T Cell", "Treg", "NK Cell", "Decidual NK Cell", "ILC3",
  "B Cell", "Plasma Cell",

  # Immune (myeloid)
  "Macrophage", "Hofbauer Cell",

  # Stromal / vascular
  "Fibroblast", "Perivascular Cell", "Placental Endothelial",

  # Trophoblast
  "Extravillous Trophoblast", "Villous Cytotrophoblast", "Syncytiotrophoblast",

  # Low-confidence last
  "RCTD_lowconf"
)

# Keep only levels that actually exist; append any missing/unexpected types at the end,
# but still force RCTD_lowconf to be last if present.
desired_present <- desired[desired %in% lvl]
extras <- setdiff(lvl, desired_present)

# If RCTD_lowconf exists but got swept into extras, pull it out and re-append last
extras_no_lowconf <- setdiff(extras, "RCTD_lowconf")
new_levels <- c(
  setdiff(desired_present, "RCTD_lowconf"),
  extras_no_lowconf,
  intersect("RCTD_lowconf", lvl)  # adds it only if present
)

obj$nn.predicted.celltype <- factor(obj$nn.predicted.celltype, levels = new_levels)
Idents(obj) <- "nn.predicted.celltype"  # refresh identities to match

# ------------------------------------------------------------
# Colors: use a NAMED vector so reordering does not scramble colors
# ------------------------------------------------------------
mypal3 <- get_palette("igv", 40)
cell_levels <- levels(obj$nn.predicted.celltype)

# assign colors in level order (truncate palette if fewer than 40)
col_map <- setNames(mypal3[seq_along(cell_levels)], cell_levels)

# ------------------------------------------------------------
# Plot + save
# ------------------------------------------------------------
p2 <- DimPlot(
  obj,
  reduction = "umap",
  group.by  = "nn.predicted.celltype",
  alpha     = 0.1,
  label     = FALSE,
  raster    = TRUE,
  cols      = col_map
) +
  ggtitle("Predicted Cell Types (3,969,532 Cells)") +
  NoAxes()

p2

ggsave(
  "obj_nn.predicted.celltype.png",
  plot   = p2,
  device = "png",
  width  = 8,
  height = 5,
  units  = "in",
  dpi    = 1200
)


library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# ---- choose your subject column name here ----
subject_col <- "sample_id_std"   

md <- obj[[]]

stopifnot(subject_col %in% colnames(md))
stopifnot("nn.predicted.celltype" %in% colnames(md))

# Ensure order is what you want (uses your factor levels)
md$nn.predicted.celltype <- factor(
  md$nn.predicted.celltype,
  levels = levels(obj$nn.predicted.celltype)
)

# Recreate the same named color map if it is not already in your environment
# (If you already have col_map from your UMAP code, you can skip this block.)
if (!exists("col_map")) {
  mypal3 <- get_palette("igv", 40)
  cell_levels <- levels(md$nn.predicted.celltype)
  col_map <- setNames(mypal3[seq_along(cell_levels)], cell_levels)
}

df_prop <- md %>%
  transmute(
    subject = .data[[subject_col]],
    celltype = nn.predicted.celltype
  ) %>%
  filter(!is.na(subject), !is.na(celltype)) %>%
  count(subject, celltype, name = "n") %>%
  group_by(subject) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Optional: order subjects by a biologically relevant cell type proportion (e.g., EVT)
# target_ct <- "Extravillous Trophoblast"
# subj_order <- df_prop %>% filter(celltype == target_ct) %>% arrange(desc(prop)) %>% pull(subject)
# df_prop$subject <- factor(df_prop$subject, levels = subj_order)

p_bar <- ggplot(df_prop, aes(x = subject, y = prop, fill = celltype)) +
  geom_col(width = 0.95, color = NA) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = col_map, drop = FALSE) +
  labs(x = "Subject", y = "Proportion of cells", fill = "Cell type") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )

p_bar

ggsave(
  "obj_nn_pred_celltype_prop_by_subject.png",
  plot = p_bar,
  device = "png",
  width = 10,
  height = 6,
  units = "in",
  dpi = 600
)
















## Identify the reference cells (prefer barcodes from `ref`; fallback to the known tag)
ref_cells <- if (exists("ref")) {
  intersect(colnames(obj), colnames(ref))          # barcodes from ref present in obj
} else {
  rownames(obj@meta.data)[obj$sample_id_std %in% c("placental_ref")]
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

Idents(obj) <- "nn.predicted.celltype"

n_panels <- length(unique(obj$nn.predicted.celltype))
ncol_pl  <- 4
nrow_pl  <- ceiling(n_panels / ncol_pl)

p2 <- DimPlot(
  obj, reduction = "umap",
  #group.by = "nn.predicted.celltype",
  split.by = "nn.predicted.celltype",
  ncol = ncol_pl,           # <-- number of columns
  alpha = 0.1, label = FALSE, raster = TRUE,
  cols = mypal3
) + ggtitle(sprintf("Predicted Cell Types (3,969,532 cells)", format(ncol(obj), big.mark=","))) +
  NoAxes()

# make the export size scale with the grid so panels don’t get squished
ggsave(
  "obj_split.by_nn.predicted.celltype.png", plot = p2, device = "png",
  width  = 3.0 * ncol_pl,     # tweak panel size multiplier as you like
  height = 3.0 * nrow_pl,
  units = "in", dpi = 1200, limitsize = FALSE
)
p2<-DimPlot(obj, reduction = "umap", group.by = "nn.predicted.celltype", alpha = 0.1, label = FALSE,raster=T, cols = mypal3) + ggtitle("Predicted Cell Types (3,969,532 Cells)") + NoAxes()
p2
ggsave("obj_nn.predicted.celltype.png", plot = p2, device = "png", width = 8, height = 5, units = "in", dpi = 1200)









# ============================================================
# CellChat drop-in: High vs Low within Term/Preterm
# Robust to CellChat API differences (subsetData keep.signaling)
# Writes: RDS + CSV tables + PDF plots (+ optional per-platform runs)
# ============================================================

set.seed(1)
options(future.globals.maxSize = 300e9)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr); library(tidyr); library(purrr); library(stringr); library(tibble); library(readr)
  library(Matrix)
  library(CellChat)
})

# -------------------- CONFIG --------------------
out_root <- "MNP_interactomes_CellChat_v3"

run_per_platform <- TRUE   # you said "sure" -> set TRUE; set FALSE to skip
platform_levels <- c("Xenium","VisiumHD","VisiumV1")

min_cells_ct   <- 50   # minimum cells per celltype within a condition to keep that celltype
min_cells_comm <- 10   # minimum cells per group for communication filtering

# celltype order + colors (yours)
ord_celltypes <- c(
  "Syncytiotrophoblast","Villous Cytotrophoblast","Extravillous Trophoblast",
  "Fibroblast","Perivascular Cell","Placental Endothelial",
  "Hofbauer Cell","Macrophage","NK Cell","Decidual NK Cell","ILC3",
  "Plasma Cell","B Cell","T Cell","Treg","RCTD_lowconf"
)

ct_colors <- c(
  "Syncytiotrophoblast"     = "#5050FF",
  "Villous Cytotrophoblast" = "#CE3D32",
  "Extravillous Trophoblast"= "#749B58",
  "Fibroblast"              = "#F0E685",
  "Perivascular Cell"       = "#466983",
  "Placental Endothelial"   = "#BA6338",
  "Hofbauer Cell"           = "#5DB1DD",
  "Macrophage"              = "#802268",
  "NK Cell"                 = "#6BD76B",
  "Decidual NK Cell"        = "#D595A7",
  "ILC3"                    = "#924822",
  "Plasma Cell"             = "#837B8D",
  "B Cell"                  = "#C75127",
  "T Cell"                  = "#D58F5C",
  "Treg"                    = "#7A65A5",
  "RCTD_lowconf"            = "#BFBFBF"
)

safe_dir <- function(path) { dir.create(path, showWarnings = FALSE, recursive = TRUE); path }
`%||%` <- function(a,b) if (!is.null(a)) a else b

# -------------------- HELPERS --------------------
use_SCT <- function(seu) {
  if ("SCT" %in% names(seu@assays)) DefaultAssay(seu) <- "SCT"
  seu
}

factorize_ct <- function(seu) {
  md <- seu[[]]
  stopifnot("nn.predicted.celltype" %in% colnames(md))
  md$nn.predicted.celltype <- factor(as.character(md$nn.predicted.celltype), levels = ord_celltypes)
  seu@meta.data <- md
  Idents(seu) <- "nn.predicted.celltype"
  seu
}

get_data_mat <- function(seu, assay = DefaultAssay(seu)) {
  lay <- tryCatch(SeuratObject::Layers(seu[[assay]]), error = function(e) NULL)
  if (!is.null(lay)) GetAssayData(seu, assay = assay, layer = "data") else GetAssayData(seu, assay = assay, slot = "data")
}

parse_subject_platform <- function(seu) {
  md <- seu[[]]
  stopifnot("sample_id_std" %in% colnames(md))
  md <- md %>%
    mutate(
      subject = if_else(str_detect(sample_id_std, "_"),
                        str_replace(sample_id_std, "_.*$", ""),
                        sample_id_std),
      platform_parsed = if_else(str_detect(sample_id_std, "_"),
                                str_replace(sample_id_std, "^.*_", ""),
                                "ref"),
      platform_parsed = recode(platform_parsed,
                               "xen"="Xenium","hd"="VisiumHD","v1"="VisiumV1",
                               .default = platform_parsed)
    )
  seu@meta.data <- md
  seu
}

prep_obj_for_cellchat <- function(seu) {
  seu <- use_SCT(seu)
  seu <- parse_subject_platform(seu)

  md <- seu[[]]
  need <- c("platform","MNPBurden","Delivery_class","nn.predicted.celltype","sample_id_std")
  miss <- setdiff(need, colnames(md))
  if (length(miss)) stop("Missing metadata columns: ", paste(miss, collapse=", "))

  seu <- subset(
    seu,
    subset =
      platform %in% platform_levels &
      sample_id_std != "placental_ref" &
      !is.na(MNPBurden) & MNPBurden %in% c("High","Low") &
      !is.na(Delivery_class) & Delivery_class %in% c("Term","Preterm")
  )

  seu <- factorize_ct(seu)
  seu
}

ct_cols_for <- function(cellchat_obj) {
  grps <- levels(cellchat_obj@idents)
  cols <- ct_colors[grps]; cols[is.na(cols)] <- "#808080"
  as.character(cols)
}

# Robust wrapper for subsetData across CellChat versions
subsetData_robust <- function(cellchat_obj) {
  # Newer/older versions differ; safest is call without extra args.
  # If your version supports keep.signaling or other args, you can add logic here.
  subsetData(cellchat_obj)
}

# Robust wrapper for extracting comm table across CellChat versions
extract_comm_robust <- function(cellchat_obj) {
  # Most versions: subsetCommunication(cellchat)
  out <- tryCatch({
    subsetCommunication(cellchat_obj)
  }, error = function(e) {
    # fallback: try slot directly (not ideal but better than nothing)
    tryCatch(as.data.frame(cellchat_obj@net$weight), error = function(e2) NULL)
  })
  out
}

# -------------------- RUN ONE CONDITION --------------------
run_cellchat_one <- function(seu, delivery, burden, outdir,
                             min_cells_ct = 50, min_cells_comm = 10) {

  sub <- subset(seu, subset = Delivery_class == delivery & MNPBurden == burden)
  if (ncol(sub) == 0) stop("No cells for ", delivery, "_", burden)

  # drop tiny cell types
  tab <- table(sub$nn.predicted.celltype)
  keep_ct <- names(tab)[tab >= min_cells_ct]
  sub <- subset(sub, subset = nn.predicted.celltype %in% keep_ct)

  if (length(keep_ct) < 2) stop("Too few cell types >= ", min_cells_ct, " cells in ", delivery, "_", burden)

  data.input <- get_data_mat(sub, assay = DefaultAssay(sub))
  meta.input <- sub[[]]

  cellchat <- createCellChat(object = data.input, meta = meta.input, group.by = "nn.predicted.celltype")
  cellchat@DB <- CellChatDB.human

  # API-safe: do not pass keep.signaling
  cellchat <- subsetData_robust(cellchat)

  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, raw.use = FALSE)
  cellchat <- filterCommunication(cellchat, min.cells = min_cells_comm)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat)

  od <- safe_dir(outdir)
  saveRDS(cellchat, file = file.path(od, "cellchat.rds"))

  # ----- tables -----
  write.csv(as.data.frame(cellchat@net$weight), file.path(od, "net_weight_matrix.csv"), row.names = TRUE)
  write.csv(as.data.frame(cellchat@net$count),  file.path(od, "net_count_matrix.csv"),  row.names = TRUE)

  comm <- extract_comm_robust(cellchat)
  if (!is.null(comm)) {
    # subsetCommunication returns LR-level table in most versions
    tryCatch(
      readr::write_csv(as_tibble(comm), file.path(od, "communications_table.csv")),
      error = function(e) {
        # if it's a matrix-like object
        write.csv(as.data.frame(comm), file.path(od, "communications_table.csv"), row.names = TRUE)
      }
    )
  }

  # ----- plots -----
  cols <- ct_cols_for(cellchat)
  pdf(file.path(od, "net_circle_weight.pdf"), width = 8, height = 8)
  netVisual_circle(cellchat@net$weight, color.use = cols,
                   vertex.weight = as.numeric(table(cellchat@idents)),
                   label.edge = FALSE)
  dev.off()

  pdf(file.path(od, "net_circle_count.pdf"), width = 8, height = 8)
  netVisual_circle(cellchat@net$count, color.use = cols,
                   vertex.weight = as.numeric(table(cellchat@idents)),
                   label.edge = FALSE)
  dev.off()

  pdf(file.path(od, "net_heatmap.pdf"), width = 9, height = 7)
  netVisual_heatmap(cellchat, color.use = cols)
  dev.off()

  invisible(cellchat)
}

# -------------------- RUN ALL 4 CONDITIONS (POOLED OR PLATFORM-SPLIT) --------------------
run_cellchat_all_conditions <- function(seu, out_root, prefix = "pooled") {
  cc_root <- safe_dir(file.path(out_root, "CellChat", prefix))

  conds <- tidyr::crossing(
    Delivery_class = c("Term","Preterm"),
    MNPBurden = c("Low","High")
  ) %>% mutate(cond = paste(Delivery_class, MNPBurden, sep = "_"))

  cc_list <- list()
  for (i in seq_len(nrow(conds))) {
    delivery <- conds$Delivery_class[i]
    burden   <- conds$MNPBurden[i]
    cond     <- conds$cond[i]
    message("[CellChat] Running: ", prefix, " :: ", cond)

    cc_list[[cond]] <- run_cellchat_one(
      seu, delivery = delivery, burden = burden,
      outdir = file.path(cc_root, cond),
      min_cells_ct = min_cells_ct,
      min_cells_comm = min_cells_comm
    )
  }

  saveRDS(cc_list, file.path(cc_root, "cellchat_all_conditions.list.rds"))
  cc_list
}

# -------------------- COMPARISON (HIGH VS LOW) WITH STATS + OUTPUTS --------------------
compare_high_low <- function(cc_low, cc_high, outdir, label = "High_vs_Low") {
  od <- safe_dir(outdir)

  merged <- mergeCellChat(list(Low = cc_low, High = cc_high), add.names = c("Low","High"))
  saveRDS(merged, file.path(od, "cellchat_merged.rds"))

  # ---- basic stats tables (supported in most versions) ----
  # compareInteractions returns summary; capture if it prints
  stats_count <- tryCatch(compareInteractions(merged, group = c(1,2), measure = "count"), error = function(e) NULL)
  stats_weight <- tryCatch(compareInteractions(merged, group = c(1,2), measure = "weight"), error = function(e) NULL)

  # write summary (even if NULL, keep a note)
  writeLines(c(
    paste0("label: ", label),
    paste0("compareInteractions(count): ", if (is.null(stats_count)) "FAILED/UNSUPPORTED" else "OK"),
    paste0("compareInteractions(weight): ", if (is.null(stats_weight)) "FAILED/UNSUPPORTED" else "OK")
  ), file.path(od, "compareInteractions_status.txt"))

  # ---- differential network plots ----
  pdf(file.path(od, "diff_circle_count.pdf"), width = 8, height = 8)
  netVisual_diffInteraction(merged, weight.scale = TRUE, measure = "count")
  dev.off()

  pdf(file.path(od, "diff_circle_weight.pdf"), width = 8, height = 8)
  netVisual_diffInteraction(merged, weight.scale = TRUE, measure = "weight")
  dev.off()

  # ---- ranking plot; try stats if supported ----
  # do.stat may not exist / may be heavy; try TRUE then fallback FALSE
  pdf(file.path(od, "rankNet_comparison.pdf"), width = 10, height = 6)
  ok <- tryCatch({
    rankNet(merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)
    TRUE
  }, error = function(e) FALSE)
  if (!ok) rankNet(merged, mode = "comparison", stacked = TRUE, do.stat = FALSE)
  dev.off()

  # ---- export merged comm tables per condition (LR-level) ----
  comm_low  <- tryCatch(subsetCommunication(merged, sources.use = NULL, targets.use = NULL, slot.name = "net", thresh = 0), error = function(e) NULL)
  # Some versions differ; so we also export the individual comm tables we already saved in each condition dir.
  # If comm_low exists and is a data.frame-like object, write it:
  if (!is.null(comm_low)) {
    tryCatch(
      readr::write_csv(as_tibble(comm_low), file.path(od, "merged_subsetCommunication.csv")),
      error = function(e) write.csv(as.data.frame(comm_low), file.path(od, "merged_subsetCommunication.csv"), row.names = TRUE)
    )
  }

  invisible(merged)
}

# ============================================================
# RUN
# ============================================================
obj_lr <- prep_obj_for_cellchat(obj)
safe_dir(out_root)

# pooled (all platforms together)
cc_list_pooled <- run_cellchat_all_conditions(obj_lr, out_root = out_root, prefix = "pooled")

# comparisons pooled
comp_root <- safe_dir(file.path(out_root, "CellChat_compare", "pooled"))
merged_term_pooled <- compare_high_low(cc_list_pooled[["Term_Low"]],    cc_list_pooled[["Term_High"]],
                                       outdir = file.path(comp_root, "Term_High_vs_Low"),
                                       label = "pooled_Term_High_vs_Low")
merged_pre_pooled  <- compare_high_low(cc_list_pooled[["Preterm_Low"]], cc_list_pooled[["Preterm_High"]],
                                       outdir = file.path(comp_root, "Preterm_High_vs_Low"),
                                       label = "pooled_Preterm_High_vs_Low")

# optional: per-platform runs
if (run_per_platform) {
  for (plat in platform_levels) {
    message("[CellChat] PLATFORM split: ", plat)
    seu_p <- subset(obj_lr, subset = platform == plat)
    if (ncol(seu_p) == 0) next

    prefix <- paste0("platform_", plat)
    cc_list_plat <- run_cellchat_all_conditions(seu_p, out_root = out_root, prefix = prefix)

    comp_root_p <- safe_dir(file.path(out_root, "CellChat_compare", prefix))
    compare_high_low(cc_list_plat[["Term_Low"]],    cc_list_plat[["Term_High"]],
                     outdir = file.path(comp_root_p, "Term_High_vs_Low"),
                     label = paste0(prefix, "_Term_High_vs_Low"))
    compare_high_low(cc_list_plat[["Preterm_Low"]], cc_list_plat[["Preterm_High"]],
                     outdir = file.path(comp_root_p, "Preterm_High_vs_Low"),
                     label = paste0(prefix, "_Preterm_High_vs_Low"))

    saveRDS(cc_list_plat, file.path(out_root, "CellChat", prefix, "cellchat_all_conditions.list.rds"))
  }
}

message("DONE. Outputs under: ", normalizePath(out_root))
























# ============================================================
# FULL DROP-IN: CellChat (and wrappers for High vs Low within Term/Preterm)
# Robust to CellChat API differences + droplevels fixes + outputs + visuals + stats tables
# Assumes Seurat object in memory: obj
# Uses metadata: nn.predicted.celltype, MNPBurden, Delivery_class, platform, sample_id_std
# ============================================================

set.seed(1)
options(future.globals.maxSize = 300e9)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr); library(tidyr); library(purrr); library(stringr); library(tibble); library(readr)
  library(Matrix)
  library(CellChat)
})

# -------------------- CONFIG --------------------
out_root <- "MNP_interactomes_CellChat_v4"

run_per_platform <- TRUE
platform_levels <- c("Xenium","VisiumHD","VisiumV1")  # keep only these

min_cells_ct   <- 50   # min cells per celltype within a condition
min_cells_comm <- 10   # min cells per group for communication filtering

# optional: cap cells per celltype to reduce dense allocations (set NULL to disable)
cap_cells_per_celltype <- NULL  # e.g., 50000

# celltype order + colors (your palette; RCTD_lowconf last)
ord_celltypes <- c(
  "Syncytiotrophoblast","Villous Cytotrophoblast","Extravillous Trophoblast",
  "Fibroblast","Perivascular Cell","Placental Endothelial",
  "Hofbauer Cell","Macrophage","NK Cell","Decidual NK Cell","ILC3",
  "Plasma Cell","B Cell","T Cell","Treg","RCTD_lowconf"
)

ct_colors <- c(
  "Syncytiotrophoblast"     = "#5050FF",
  "Villous Cytotrophoblast" = "#CE3D32",
  "Extravillous Trophoblast"= "#749B58",
  "Fibroblast"              = "#F0E685",
  "Perivascular Cell"       = "#466983",
  "Placental Endothelial"   = "#BA6338",
  "Hofbauer Cell"           = "#5DB1DD",
  "Macrophage"              = "#802268",
  "NK Cell"                 = "#6BD76B",
  "Decidual NK Cell"        = "#D595A7",
  "ILC3"                    = "#924822",
  "Plasma Cell"             = "#837B8D",
  "B Cell"                  = "#C75127",
  "T Cell"                  = "#D58F5C",
  "Treg"                    = "#7A65A5",
  "RCTD_lowconf"            = "#BFBFBF"
)

# -------------------- UTIL --------------------
safe_dir <- function(path) { dir.create(path, showWarnings = FALSE, recursive = TRUE); path }
`%||%` <- function(a,b) if (!is.null(a)) a else b

use_SCT <- function(seu) {
  if ("SCT" %in% names(seu@assays)) DefaultAssay(seu) <- "SCT"
  seu
}

get_data_mat <- function(seu, assay = DefaultAssay(seu)) {
  lay <- tryCatch(SeuratObject::Layers(seu[[assay]]), error = function(e) NULL)
  m <- if (!is.null(lay)) GetAssayData(seu, assay = assay, layer = "data") else GetAssayData(seu, assay = assay, slot = "data")
  if (!inherits(m, "dgCMatrix")) m <- as(m, "dgCMatrix")
  m
}

parse_subject_platform <- function(seu) {
  md <- seu[[]]
  stopifnot("sample_id_std" %in% colnames(md))
  md <- md %>%
    mutate(
      subject = if_else(str_detect(sample_id_std, "_"),
                        str_replace(sample_id_std, "_.*$", ""),
                        sample_id_std),
      platform_parsed = if_else(str_detect(sample_id_std, "_"),
                                str_replace(sample_id_std, "^.*_", ""),
                                "ref"),
      platform_parsed = recode(platform_parsed,
                               "xen"="Xenium","hd"="VisiumHD","v1"="VisiumV1",
                               .default = platform_parsed)
    )
  seu@meta.data <- md
  seu
}

refactor_celltypes <- function(seu, levels_master) {
  md <- seu[[]]
  md$nn.predicted.celltype <- as.character(md$nn.predicted.celltype)
  observed <- unique(md$nn.predicted.celltype)
  observed <- observed[!is.na(observed)]
  lev <- c(intersect(levels_master, observed), setdiff(observed, levels_master))
  md$nn.predicted.celltype <- factor(md$nn.predicted.celltype, levels = lev)
  md$nn.predicted.celltype <- droplevels(md$nn.predicted.celltype)
  seu@meta.data <- md
  Idents(seu) <- "nn.predicted.celltype"
  seu
}

ct_cols_for <- function(cellchat_obj) {
  grps <- levels(cellchat_obj@idents)
  cols <- ct_colors[grps]; cols[is.na(cols)] <- "#808080"
  as.character(cols)
}

# CellChat API robustness
subsetData_robust <- function(cellchat_obj) {
  # do NOT pass keep.signaling to support your installed CellChat
  subsetData(cellchat_obj)
}

extract_comm_robust <- function(cellchat_obj) {
  # most versions support subsetCommunication
  out <- tryCatch(subsetCommunication(cellchat_obj), error = function(e) NULL)
  out
}

# optional downsample per celltype
downsample_per_celltype <- function(seu, cap) {
  if (is.null(cap)) return(seu)
  md <- seu[[]] %>% mutate(cell = rownames(.), ct = as.character(nn.predicted.celltype))
  keep_cells <- md %>%
    filter(!is.na(ct)) %>%
    group_by(ct) %>%
    slice_sample(n = min(n(), cap)) %>%
    pull(cell)
  subset(seu, cells = keep_cells)
}

# -------------------- PREP OBJECT --------------------
prep_obj_for_cellchat <- function(seu) {
  seu <- use_SCT(seu)
  seu <- parse_subject_platform(seu)

  md <- seu[[]]
  need <- c("platform","MNPBurden","Delivery_class","nn.predicted.celltype","sample_id_std")
  miss <- setdiff(need, colnames(md))
  if (length(miss)) stop("Missing metadata columns: ", paste(miss, collapse=", "))

  # keep only spatial platforms; drop placental_ref; keep High/Low and Term/Preterm
  seu <- subset(
    seu,
    subset =
      platform %in% platform_levels &
      sample_id_std != "placental_ref" &
      !is.na(MNPBurden) & MNPBurden %in% c("High","Low") &
      !is.na(Delivery_class) & Delivery_class %in% c("Term","Preterm")
  )

  # enforce celltype order + droplevels
  seu <- refactor_celltypes(seu, ord_celltypes)
  seu
}

# -------------------- RUN ONE CONDITION --------------------
run_cellchat_one <- function(seu, delivery, burden, outdir,
                             min_cells_ct = 50, min_cells_comm = 10,
                             cap_cells_per_celltype = NULL) {

  od <- safe_dir(outdir)

  sub <- subset(seu, subset = Delivery_class == delivery & MNPBurden == burden)
  if (ncol(sub) == 0) stop("No cells for ", delivery, "_", burden)

  # optional downsample (reduces dense allocations)
  sub <- downsample_per_celltype(sub, cap_cells_per_celltype)

  # drop tiny celltypes using character (pre-factor)
  sub$nn.predicted.celltype <- as.character(sub$nn.predicted.celltype)
  tab <- table(sub$nn.predicted.celltype, useNA = "no")
  keep_ct <- names(tab)[tab >= min_cells_ct]
  sub <- subset(sub, subset = nn.predicted.celltype %in% keep_ct)
  if (length(keep_ct) < 2) stop("Too few cell types >= ", min_cells_ct, " cells in ", delivery, "_", burden)

  # CRITICAL: refactor + droplevels AFTER filtering
  sub <- refactor_celltypes(sub, ord_celltypes)

  # QA outputs
  qa_levels <- tibble(level = levels(sub$nn.predicted.celltype),
                      present = level %in% unique(as.character(sub$nn.predicted.celltype)))
  write_csv(qa_levels, file.path(od, "idents_levels_present.csv"))
  write_csv(
    as_tibble(table(sub$nn.predicted.celltype), .name_repair = "unique") %>%
      rename(celltype = 1, n_cells = 2),
    file.path(od, "idents_cell_counts.csv")
  )

  data.input <- get_data_mat(sub, assay = DefaultAssay(sub))
  meta.input <- sub[[]] %>%
    mutate(nn.predicted.celltype = droplevels(nn.predicted.celltype))

  cellchat <- createCellChat(object = data.input, meta = meta.input, group.by = "nn.predicted.celltype")
  cellchat@DB <- CellChatDB.human

  # ensure idents are clean
  cellchat@idents <- droplevels(cellchat@idents)

  cellchat <- subsetData_robust(cellchat)

  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  # drop unused again right before computeCommunProb
  cellchat@idents <- droplevels(cellchat@idents)

  cellchat <- computeCommunProb(cellchat, raw.use = FALSE)
  cellchat <- filterCommunication(cellchat, min.cells = min_cells_comm)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat)

  saveRDS(cellchat, file.path(od, "cellchat.rds"))

  # Tables
  write.csv(as.data.frame(cellchat@net$weight), file.path(od, "net_weight_matrix.csv"), row.names = TRUE)
  write.csv(as.data.frame(cellchat@net$count),  file.path(od, "net_count_matrix.csv"),  row.names = TRUE)

  comm <- extract_comm_robust(cellchat)
  if (!is.null(comm)) write_csv(as_tibble(comm), file.path(od, "communications_table.csv"))

  # Plots
  cols <- ct_cols_for(cellchat)

  pdf(file.path(od, "net_circle_weight.pdf"), width = 8, height = 8)
  netVisual_circle(cellchat@net$weight, color.use = cols,
                   vertex.weight = as.numeric(table(cellchat@idents)),
                   label.edge = FALSE)
  dev.off()

  pdf(file.path(od, "net_circle_count.pdf"), width = 8, height = 8)
  netVisual_circle(cellchat@net$count, color.use = cols,
                   vertex.weight = as.numeric(table(cellchat@idents)),
                   label.edge = FALSE)
  dev.off()

  pdf(file.path(od, "net_heatmap.pdf"), width = 9, height = 7)
  netVisual_heatmap(cellchat, color.use = cols)
  dev.off()

  invisible(cellchat)
}

# -------------------- RUN ALL 4 CONDITIONS --------------------
run_cellchat_all_conditions <- function(seu, out_root, prefix = "pooled") {
  cc_root <- safe_dir(file.path(out_root, "CellChat", prefix))

  conds <- tidyr::crossing(
    Delivery_class = c("Term","Preterm"),
    MNPBurden = c("Low","High")
  ) %>% mutate(cond = paste(Delivery_class, MNPBurden, sep = "_"))

  cc_list <- list()
  for (i in seq_len(nrow(conds))) {
    delivery <- conds$Delivery_class[i]
    burden   <- conds$MNPBurden[i]
    cond     <- conds$cond[i]

    message("[CellChat] Running: ", prefix, " :: ", cond)

    cc_list[[cond]] <- run_cellchat_one(
      seu, delivery = delivery, burden = burden,
      outdir = file.path(cc_root, cond),
      min_cells_ct = min_cells_ct,
      min_cells_comm = min_cells_comm,
      cap_cells_per_celltype = cap_cells_per_celltype
    )
  }

  saveRDS(cc_list, file.path(cc_root, "cellchat_all_conditions.list.rds"))
  cc_list
}

# -------------------- COMPARE HIGH VS LOW (TERM/PRETERM) + OUTPUTS --------------------
compare_high_low <- function(cc_low, cc_high, outdir, label) {
  od <- safe_dir(outdir)

  merged <- mergeCellChat(list(Low = cc_low, High = cc_high), add.names = c("Low","High"))
  saveRDS(merged, file.path(od, "cellchat_merged.rds"))

  # diff plots
  pdf(file.path(od, "diff_circle_count.pdf"), width = 8, height = 8)
  netVisual_diffInteraction(merged, weight.scale = TRUE, measure = "count")
  dev.off()

  pdf(file.path(od, "diff_circle_weight.pdf"), width = 8, height = 8)
  netVisual_diffInteraction(merged, weight.scale = TRUE, measure = "weight")
  dev.off()

  # rankNet (try stats; fall back)
  pdf(file.path(od, "rankNet_comparison.pdf"), width = 10, height = 6)
  ok <- tryCatch({
    rankNet(merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)
    TRUE
  }, error = function(e) FALSE)
  if (!ok) rankNet(merged, mode = "comparison", stacked = TRUE, do.stat = FALSE)
  dev.off()

  # Summary stats tables (version dependent)
  # Try compareInteractions and capture outputs in a text file
  sink(file.path(od, "compareInteractions_summary.txt"))
  cat("Label:", label, "\n\n")
  cat("compareInteractions (count):\n")
  print(tryCatch(compareInteractions(merged, group = c(1,2), measure = "count"), error = function(e) e))
  cat("\ncompareInteractions (weight):\n")
  print(tryCatch(compareInteractions(merged, group = c(1,2), measure = "weight"), error = function(e) e))
  sink()

  # Export a merged communications table if supported
  comm_merged <- tryCatch(subsetCommunication(merged), error = function(e) NULL)
  if (!is.null(comm_merged)) write_csv(as_tibble(comm_merged), file.path(od, "merged_subsetCommunication.csv"))

  invisible(merged)
}

# ============================================================
# RUN EVERYTHING
# ============================================================
safe_dir(out_root)
obj_lr <- prep_obj_for_cellchat(obj)

# pooled run
cc_list_pooled <- run_cellchat_all_conditions(obj_lr, out_root = out_root, prefix = "pooled")

# pooled comparisons
comp_root <- safe_dir(file.path(out_root, "CellChat_compare", "pooled"))
merged_term_pooled <- compare_high_low(
  cc_list_pooled[["Term_Low"]], cc_list_pooled[["Term_High"]],
  outdir = file.path(comp_root, "Term_High_vs_Low"),
  label = "pooled_Term_High_vs_Low"
)
merged_preterm_pooled <- compare_high_low(
  cc_list_pooled[["Preterm_Low"]], cc_list_pooled[["Preterm_High"]],
  outdir = file.path(comp_root, "Preterm_High_vs_Low"),
  label = "pooled_Preterm_High_vs_Low"
)

# optional per-platform runs + comparisons
if (run_per_platform) {
  for (plat in platform_levels) {
    message("[CellChat] PLATFORM split: ", plat)
    seu_p <- subset(obj_lr, subset = platform == plat)
    if (ncol(seu_p) == 0) next

    prefix <- paste0("platform_", plat)
    cc_list_plat <- run_cellchat_all_conditions(seu_p, out_root = out_root, prefix = prefix)

    comp_root_p <- safe_dir(file.path(out_root, "CellChat_compare", prefix))
    compare_high_low(
      cc_list_plat[["Term_Low"]], cc_list_plat[["Term_High"]],
      outdir = file.path(comp_root_p, "Term_High_vs_Low"),
      label = paste0(prefix, "_Term_High_vs_Low")
    )
    compare_high_low(
      cc_list_plat[["Preterm_Low"]], cc_list_plat[["Preterm_High"]],
      outdir = file.path(comp_root_p, "Preterm_High_vs_Low"),
      label = paste0(prefix, "_Preterm_High_vs_Low")
    )
  }
}

message("DONE. Outputs under: ", normalizePath(out_root))
























# ============================================================
# FULLY REVISED DROP-IN: run_cellchat_smoke_forceLR_v8_SCTdata
# Uses: sub@assays[["SCT"]]@data  (your request)
# Notes:
#  - No CellChat_input assay created.
#  - Uses SCT *data* (log-normalized) as the expression matrix for CellChat.
#  - Standardizes feature names (ENSG->SYMBOL, uppercase, collapses duplicates).
#  - Does NOT pre-intersect to DB genes (lets CellChat handle it).
#  - Forces LRsig from DB restricted to genes present.
#  - Prints diagnostics to console and writes key files.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr); library(tibble); library(readr)
  library(Matrix)
  library(CellChat)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

safe_dir <- function(path) { dir.create(path, showWarnings = FALSE, recursive = TRUE); path }

# Seurat v5-safe data accessor; falls back to slot for older objects
get_sct_data <- function(seu) {
  if (!"SCT" %in% names(seu@assays)) stop("No SCT assay found.")
  tryCatch(
    GetAssayData(seu, assay = "SCT", layer = "data"),
    error = function(e) GetAssayData(seu, assay = "SCT", slot = "data")
  )
}

downsample_per_celltype <- function(seu, cap, ct_col = "nn.predicted.celltype") {
  if (is.null(cap)) return(seu)
  md <- seu[[]]
  ct <- as.character(md[[ct_col]])
  names(ct) <- rownames(md)
  ct <- ct[!is.na(ct)]
  sp <- split(names(ct), ct)
  keep <- unlist(lapply(sp, function(cells) if (length(cells) <= cap) cells else sample(cells, cap)),
                 use.names = FALSE)
  subset(seu, cells = keep)
}

map_mixed_ids_to_symbol <- function(mat) {
  if (!inherits(mat, "dgCMatrix")) mat <- as(mat, "dgCMatrix")
  rn <- rownames(mat)
  is_ensg <- grepl("^ENSG", rn)

  rn_std <- rn
  rn_std[!is_ensg] <- toupper(rn_std[!is_ensg])

  if (any(is_ensg)) {
    ensg <- rn[is_ensg]
    ensg_clean <- sub("\\..*$", "", ensg)
    ensg_clean <- sub("\\|.*$", "", ensg)

    sym <- suppressMessages(
      AnnotationDbi::mapIds(
        org.Hs.eg.db,
        keys = unique(ensg_clean),
        column = "SYMBOL",
        keytype = "ENSEMBL",
        multiVals = "first"
      )
    )
    mapped <- toupper(as.character(sym[ensg_clean]))
    rn_std[is_ensg] <- ifelse(!is.na(mapped) & mapped != "", mapped, toupper(ensg))
  }

  rownames(mat) <- rn_std

  # collapse duplicates (sparse-safe)
  if (any(duplicated(rownames(mat)))) {
    grp <- factor(rownames(mat), levels = unique(rownames(mat)))
    G <- Matrix::sparse.model.matrix(~ 0 + grp)
    mat <- Matrix::t(G) %*% mat
    rownames(mat) <- levels(grp)
    mat <- as(mat, "dgCMatrix")
  }
  mat
}

subsetData_perm <- function(cc) {
  fm <- names(formals(CellChat::subsetData))
  args <- list(object = cc)
  if ("thresh" %in% fm) args$thresh <- 0
  if ("thresh.pc" %in% fm) args$thresh.pc <- 0
  do.call(CellChat::subsetData, args)
}

force_LRsig_from_DB <- function(cc, genes_present, DB = CellChatDB.human) {
  genes_present <- unique(toupper(genes_present))
  lr0 <- DB$interaction %>%
    transmute(
      interaction_name = interaction_name,
      pathway_name     = pathway_name,
      ligand           = toupper(ligand),
      receptor         = toupper(receptor),
      interaction_name_2 = if ("interaction_name_2" %in% colnames(DB$interaction)) interaction_name_2 else interaction_name
    ) %>%
    filter(ligand %in% genes_present, receptor %in% genes_present)

  cc@LR$LRsig <- as.data.frame(lr0)
  cc
}

# quick overlap diagnostics vs DB
db_overlap_diag <- function(gene_names, DB = CellChatDB.human) {
  mg <- unique(toupper(gene_names))
  dbL <- unique(toupper(DB$interaction$ligand))
  dbR <- unique(toupper(DB$interaction$receptor))
  tibble(
    n_genes = length(mg),
    overlap_ligands = sum(dbL %in% mg),
    overlap_receptors = sum(dbR %in% mg),
    n_db_ligands = length(dbL),
    n_db_receptors = length(dbR)
  )
}

run_cellchat_smoke_forceLR_v8_SCTdata <- function(
  seu,
  delivery = "Preterm", burden = "High",
  outdir,
  cap_per_ct = 10000,
  min_cells_ct = 10,
  min_cells_comm = 1,
  restrict_ct = c("Syncytiotrophoblast","Villous Cytotrophoblast","Extravillous Trophoblast",
                  "Macrophage","NK Cell","T Cell","B Cell"),
  type_use = c("median","thresholdedMean","truncatedMean","triMean"),
  raw_use = FALSE
) {
  type_use <- match.arg(type_use)
  od <- safe_dir(outdir)
  message("[SMOKE forceLR v8 | SCT@data] ", delivery, "_", burden, " -> ", od)

  message("[v8] Subsetting...")
  sub <- subset(seu, subset = Delivery_class == delivery & MNPBurden == burden)
  sub <- subset(sub, subset = nn.predicted.celltype != "RCTD_lowconf")
  sub <- subset(sub, subset = nn.predicted.celltype %in% restrict_ct)

  sub <- downsample_per_celltype(sub, cap_per_ct)

  tab <- table(as.character(sub$nn.predicted.celltype))
  keep <- names(tab)[tab >= min_cells_ct]
  sub <- subset(sub, subset = nn.predicted.celltype %in% keep)

  write_csv(
    tibble(celltype = names(table(sub$nn.predicted.celltype)),
           n_cells  = as.integer(table(sub$nn.predicted.celltype))),
    file.path(od, "idents_cell_counts.csv")
  )

  message("[v8] Pulling SCT@data...")
  data.input <- get_sct_data(sub)
  if (!inherits(data.input, "dgCMatrix")) data.input <- as(data.input, "dgCMatrix")

  message("[v8] Mapping ENSG->SYMBOL and collapsing duplicates...")
  data.mapped <- map_mixed_ids_to_symbol(data.input)
  data.mapped <- data.mapped[!grepl("^ENSG", rownames(data.mapped)), , drop = FALSE]

  diag0 <- tibble(n_cells = ncol(data.mapped), n_genes_mapped = nrow(data.mapped))
  diag1 <- db_overlap_diag(rownames(data.mapped), CellChatDB.human)
  write_csv(diag0, file.path(od, "diag_matrix_dims.csv"))
  write_csv(diag1, file.path(od, "diag_db_overlap.csv"))

  message("[v8 DIAG] n_cells=", diag0$n_cells, " | n_genes_mapped=", diag0$n_genes_mapped)
  message("[v8 DIAG] DB overlap: ligands=", diag1$overlap_ligands, "/", diag1$n_db_ligands,
          " | receptors=", diag1$overlap_receptors, "/", diag1$n_db_receptors)

  meta.input <- sub[[]] %>%
    mutate(nn.predicted.celltype = droplevels(factor(as.character(nn.predicted.celltype))))

  message("[v8] Creating CellChat...")
  cc <- createCellChat(object = data.mapped, meta = meta.input, group.by = "nn.predicted.celltype")
  cc@DB <- CellChatDB.human
  cc@idents <- droplevels(cc@idents)

  message("[v8] subsetData (permissive)...")
  cc <- subsetData_perm(cc)

  genes_present <- rownames(cc@data)
  message("[v8 DIAG] genes_present_in_cc@data=", length(genes_present))

  message("[v8] Forcing LRsig from DB...")
  cc <- force_LRsig_from_DB(cc, genes_present = genes_present, DB = CellChatDB.human)
  message("[v8 DIAG] Forced LRsig rows: ", nrow(cc@LR$LRsig))
  writeLines(paste0("Forced LRsig rows: ", nrow(cc@LR$LRsig)), file.path(od, "forced_LRsig_n.txt"))
  write_csv(as_tibble(head(cc@LR$LRsig, 50)), file.path(od, "forced_LRsig_head.csv"))

  if (nrow(cc@LR$LRsig) == 0) {
    message("[v8 FAIL] Forced LRsig is 0. This is a gene-symbol / DB overlap problem. See diag_db_overlap.csv.")
    writeLines("Forced LRsig is 0. Check diag_db_overlap.csv.", file.path(od, "FAIL_forced_LRsig_zero.txt"))
    saveRDS(cc, file.path(od, "cellchat_FAIL_forced_LRsig_zero.rds"))
    return(invisible(NULL))
  }

  message("[v8] identifyOverExpressedGenes/Interactions...")
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)

  message("[v8] Running computeCommunProb... raw.use=", raw_use)
  cc <- tryCatch({
    computeCommunProb(
      object = cc,
      type = type_use,
      trim = 0,
      LR.use = cc@LR$LRsig,      # use all forced LR pairs
      raw.use = raw_use,         # you can toggle this if needed
      population.size = FALSE
    )
  }, error = function(e) {
    message("[v8 FAIL] computeCommunProb error: ", conditionMessage(e))
    saveRDS(cc, file.path(od, "cellchat_FAIL_computeCommunProb.rds"))
    writeLines(conditionMessage(e), file.path(od, "FAIL_computeCommunProb.txt"))
    stop(e)
  })

  message("[v8] Downstream aggregation...")
  cc <- filterCommunication(cc, min.cells = min_cells_comm)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)

  comm <- tryCatch(subsetCommunication(cc), error = function(e) NULL)
  if (!is.null(comm) && is.data.frame(comm)) write_csv(as_tibble(comm), file.path(od, "communications_table.csv"))

  saveRDS(cc, file.path(od, "cellchat_SMOKE_forceLR_OK.rds"))
  message("[v8 OK] Saved: cellchat_SMOKE_forceLR_OK.rds")
  invisible(cc)
}



cc_smoke_force <- run_cellchat_smoke_forceLR_v8_SCTdata(
  obj_lr, "Preterm", "High",
  outdir = "MNP_interactomes_CellChat_v12/SMOKETEST/Preterm_High_forceLR_v8_SCTdata_rawTRUE",
  type_use = "truncatedMean",
  raw_use = TRUE
)




















## 3 days later


# ============================================================
# FULL DROP-IN: CellChat full dataset (SCT@data + raw.use=TRUE)
# Robust + diagnostics + pooled + per-platform + comparisons + visuals
# Metadata required:
#   nn.predicted.celltype, MNPBurden, Delivery_class, platform, sample_id_std
# ============================================================

set.seed(1)
options(future.globals.maxSize = 300e9)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr); library(tidyr); library(stringr); library(tibble); library(readr)
  library(Matrix)
  library(CellChat)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# -------------------- CONFIG --------------------
out_root <- "MNP_interactomes_CellChat_v8_full"

run_per_platform <- TRUE
platform_levels <- c("Xenium","VisiumHD","VisiumV1")  # keep only these

min_cells_ct   <- 50   # min cells per celltype within a condition
min_cells_comm <- 10   # min cells per group for communication filtering

cap_cells_per_celltype <- NULL  # e.g., 50000 (set NULL to disable)

ord_celltypes <- c(
  "Syncytiotrophoblast","Villous Cytotrophoblast","Extravillous Trophoblast",
  "Fibroblast","Perivascular Cell","Placental Endothelial",
  "Hofbauer Cell","Macrophage","NK Cell","Decidual NK Cell","ILC3",
  "Plasma Cell","B Cell","T Cell","Treg","RCTD_lowconf"
)

ct_colors <- c(
  "Syncytiotrophoblast"     = "#5050FF",
  "Villous Cytotrophoblast" = "#CE3D32",
  "Extravillous Trophoblast"= "#749B58",
  "Fibroblast"              = "#F0E685",
  "Perivascular Cell"       = "#466983",
  "Placental Endothelial"   = "#BA6338",
  "Hofbauer Cell"           = "#5DB1DD",
  "Macrophage"              = "#802268",
  "NK Cell"                 = "#6BD76B",
  "Decidual NK Cell"        = "#D595A7",
  "ILC3"                    = "#924822",
  "Plasma Cell"             = "#837B8D",
  "B Cell"                  = "#C75127",
  "T Cell"                  = "#D58F5C",
  "Treg"                    = "#7A65A5",
  "RCTD_lowconf"            = "#BFBFBF"
)

# -------------------- UTIL --------------------
safe_dir <- function(path) { dir.create(path, showWarnings = FALSE, recursive = TRUE); path }
`%||%` <- function(a,b) if (!is.null(a)) a else b

# Seurat v5-safe accessor for SCT@data
get_sct_data <- function(seu) {
  if (!"SCT" %in% names(seu@assays)) stop("No SCT assay found.")
  tryCatch(
    GetAssayData(seu, assay = "SCT", layer = "data"),
    error = function(e) GetAssayData(seu, assay = "SCT", slot = "data")
  )
}

# downsample per celltype (optional)
downsample_per_celltype <- function(seu, cap, ct_col = "nn.predicted.celltype") {
  if (is.null(cap)) return(seu)
  md <- seu[[]]
  ct <- as.character(md[[ct_col]])
  names(ct) <- rownames(md)
  ct <- ct[!is.na(ct)]
  sp <- split(names(ct), ct)
  keep <- unlist(lapply(sp, function(cells) if (length(cells) <= cap) cells else sample(cells, cap)),
                 use.names = FALSE)
  subset(seu, cells = keep)
}

# parse subject/platform from sample_id_std
parse_subject_platform <- function(seu) {
  md <- seu[[]]
  stopifnot("sample_id_std" %in% colnames(md))
  md <- md %>%
    mutate(
      subject = if_else(str_detect(sample_id_std, "_"),
                        str_replace(sample_id_std, "_.*$", ""),
                        sample_id_std),
      platform_parsed = if_else(str_detect(sample_id_std, "_"),
                                str_replace(sample_id_std, "^.*_", ""),
                                "ref"),
      platform_parsed = recode(platform_parsed,
                               "xen"="Xenium","hd"="VisiumHD","v1"="VisiumV1",
                               .default = platform_parsed)
    )
  seu@meta.data <- md
  seu
}

# enforce celltype order + droplevels
refactor_celltypes <- function(seu, levels_master) {
  md <- seu[[]]
  md$nn.predicted.celltype <- as.character(md$nn.predicted.celltype)
  observed <- unique(md$nn.predicted.celltype)
  observed <- observed[!is.na(observed)]
  lev <- c(intersect(levels_master, observed), setdiff(observed, levels_master))
  md$nn.predicted.celltype <- factor(md$nn.predicted.celltype, levels = lev)
  md$nn.predicted.celltype <- droplevels(md$nn.predicted.celltype)
  seu@meta.data <- md
  Idents(seu) <- "nn.predicted.celltype"
  seu
}

ct_cols_for <- function(cellchat_obj) {
  grps <- levels(cellchat_obj@idents)
  cols <- ct_colors[grps]; cols[is.na(cols)] <- "#808080"
  as.character(cols)
}

# Mixed IDs -> SYMBOL (ENSEMBL->SYMBOL), uppercase, collapse duplicates (sparse-safe)
map_mixed_ids_to_symbol <- function(mat) {
  if (!inherits(mat, "dgCMatrix")) mat <- as(mat, "dgCMatrix")
  rn <- rownames(mat)
  is_ensg <- grepl("^ENSG", rn)

  rn_std <- rn
  rn_std[!is_ensg] <- toupper(rn_std[!is_ensg])

  if (any(is_ensg)) {
    ensg <- rn[is_ensg]
    ensg_clean <- sub("\\..*$", "", ensg)
    ensg_clean <- sub("\\|.*$", "", ensg)

    sym <- suppressMessages(
      AnnotationDbi::mapIds(
        org.Hs.eg.db,
        keys = unique(ensg_clean),
        column = "SYMBOL",
        keytype = "ENSEMBL",
        multiVals = "first"
      )
    )
    mapped <- toupper(as.character(sym[ensg_clean]))
    rn_std[is_ensg] <- ifelse(!is.na(mapped) & mapped != "", mapped, toupper(ensg))
  }

  rownames(mat) <- rn_std

  if (any(duplicated(rownames(mat)))) {
    grp <- factor(rownames(mat), levels = unique(rownames(mat)))
    G <- Matrix::sparse.model.matrix(~ 0 + grp)
    mat <- Matrix::t(G) %*% mat
    rownames(mat) <- levels(grp)
    mat <- as(mat, "dgCMatrix")
  }
  mat
}

# subsetData permissive (API-safe)
subsetData_perm <- function(cc) {
  fm <- names(formals(CellChat::subsetData))
  args <- list(object = cc)
  if ("thresh" %in% fm) args$thresh <- 0
  if ("thresh.pc" %in% fm) args$thresh.pc <- 0
  do.call(CellChat::subsetData, args)
}

# Force LRsig from DB restricted to genes present
force_LRsig_from_DB <- function(cc, genes_present, DB = CellChatDB.human) {
  genes_present <- unique(toupper(genes_present))
  lr0 <- DB$interaction %>%
    transmute(
      interaction_name   = interaction_name,
      pathway_name       = pathway_name,
      ligand             = toupper(ligand),
      receptor           = toupper(receptor),
      interaction_name_2 = if ("interaction_name_2" %in% colnames(DB$interaction)) interaction_name_2 else interaction_name
    ) %>%
    filter(ligand %in% genes_present, receptor %in% genes_present)
  cc@LR$LRsig <- as.data.frame(lr0)
  cc
}

# DB overlap diag (console + file)
db_overlap_diag <- function(gene_names, DB = CellChatDB.human) {
  mg  <- unique(toupper(gene_names))
  dbL <- unique(toupper(DB$interaction$ligand))
  dbR <- unique(toupper(DB$interaction$receptor))
  tibble(
    n_genes = length(mg),
    overlap_ligands = sum(dbL %in% mg),
    overlap_receptors = sum(dbR %in% mg),
    n_db_ligands = length(dbL),
    n_db_receptors = length(dbR)
  )
}

# -------------------- PREP OBJECT --------------------
prep_obj_for_cellchat <- function(seu) {
  if (!"SCT" %in% names(seu@assays)) stop("Object must contain SCT assay.")
  DefaultAssay(seu) <- "SCT"
  seu <- parse_subject_platform(seu)

  md <- seu[[]]
  need <- c("platform","MNPBurden","Delivery_class","nn.predicted.celltype","sample_id_std")
  miss <- setdiff(need, colnames(md))
  if (length(miss)) stop("Missing metadata columns: ", paste(miss, collapse=", "))

  seu <- subset(
    seu,
    subset =
      platform %in% platform_levels &
      sample_id_std != "placental_ref" &
      !is.na(MNPBurden) & MNPBurden %in% c("High","Low") &
      !is.na(Delivery_class) & Delivery_class %in% c("Term","Preterm")
  )

  seu <- refactor_celltypes(seu, ord_celltypes)
  seu
}

# -------------------- RUN ONE CONDITION --------------------
run_cellchat_one_v8 <- function(seu, delivery, burden, outdir,
                                min_cells_ct = 50, min_cells_comm = 10,
                                cap_cells_per_celltype = NULL,
                                raw_use = TRUE,
                                type_use = c("median","thresholdedMean","truncatedMean","triMean")) {
  type_use <- match.arg(type_use)
  od <- safe_dir(outdir)

  message("[CellChat] ", delivery, "_", burden, " -> ", od)

  sub <- subset(seu, subset = Delivery_class == delivery & MNPBurden == burden)
  if (ncol(sub) == 0) stop("No cells for ", delivery, "_", burden)

  sub <- downsample_per_celltype(sub, cap_cells_per_celltype)

  # drop tiny celltypes
  sub$nn.predicted.celltype <- as.character(sub$nn.predicted.celltype)
  tab <- table(sub$nn.predicted.celltype, useNA = "no")
  keep_ct <- names(tab)[tab >= min_cells_ct]
  sub <- subset(sub, subset = nn.predicted.celltype %in% keep_ct)
  if (length(keep_ct) < 2) stop("Too few cell types >= ", min_cells_ct, " cells in ", delivery, "_", burden)

  sub <- refactor_celltypes(sub, ord_celltypes)

  # QA
  write_csv(
    tibble(celltype = names(table(sub$nn.predicted.celltype)),
           n_cells  = as.integer(table(sub$nn.predicted.celltype))),
    file.path(od, "idents_cell_counts.csv")
  )

  # Pull SCT@data and harmonize gene IDs
  data.input <- get_sct_data(sub)
  if (!inherits(data.input, "dgCMatrix")) data.input <- as(data.input, "dgCMatrix")
  data.mapped <- map_mixed_ids_to_symbol(data.input)
  data.mapped <- data.mapped[!grepl("^ENSG", rownames(data.mapped)), , drop = FALSE]

  # Diagnostics
  diag_dims <- tibble(n_cells = ncol(data.mapped), n_genes_mapped = nrow(data.mapped))
  diag_db   <- db_overlap_diag(rownames(data.mapped), CellChatDB.human)
  write_csv(diag_dims, file.path(od, "diag_matrix_dims.csv"))
  write_csv(diag_db,   file.path(od, "diag_db_overlap.csv"))
  message("[DIAG] n_cells=", diag_dims$n_cells,
          " | n_genes_mapped=", diag_dims$n_genes_mapped,
          " | DB overlap lig=", diag_db$overlap_ligands, "/", diag_db$n_db_ligands,
          " rec=", diag_db$overlap_receptors, "/", diag_db$n_db_receptors)

  meta.input <- sub[[]] %>% mutate(nn.predicted.celltype = droplevels(factor(as.character(nn.predicted.celltype))))

  cc <- createCellChat(object = data.mapped, meta = meta.input, group.by = "nn.predicted.celltype")
  cc@DB <- CellChatDB.human
  cc@idents <- droplevels(cc@idents)

  cc <- subsetData_perm(cc)

  genes_present <- rownames(cc@data)
  cc <- force_LRsig_from_DB(cc, genes_present = genes_present, DB = CellChatDB.human)

  message("[DIAG] Forced LRsig rows: ", nrow(cc@LR$LRsig))
  writeLines(paste0("Forced LRsig rows: ", nrow(cc@LR$LRsig)), file.path(od, "forced_LRsig_n.txt"))
  write_csv(as_tibble(head(cc@LR$LRsig, 50)), file.path(od, "forced_LRsig_head.csv"))

  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)

  # CRITICAL: raw.use=TRUE (your smoke test)
  cc <- tryCatch({
    computeCommunProb(cc, raw.use = raw_use, type = type_use, trim = 0, LR.use = cc@LR$LRsig, population.size = FALSE)
  }, error = function(e) {
    saveRDS(cc, file.path(od, "cellchat_FAIL_computeCommunProb.rds"))
    writeLines(conditionMessage(e), file.path(od, "FAIL_computeCommunProb.txt"))
    stop(e)
  })

  cc <- filterCommunication(cc, min.cells = min_cells_comm)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc)

  saveRDS(cc, file.path(od, "cellchat.rds"))

  # Tables
  write.csv(as.data.frame(cc@net$weight), file.path(od, "net_weight_matrix.csv"), row.names = TRUE)
  write.csv(as.data.frame(cc@net$count),  file.path(od, "net_count_matrix.csv"),  row.names = TRUE)

  comm <- tryCatch(subsetCommunication(cc), error = function(e) NULL)
  if (!is.null(comm) && is.data.frame(comm)) write_csv(as_tibble(comm), file.path(od, "communications_table.csv"))

  # Visuals
  cols <- ct_cols_for(cc)

  pdf(file.path(od, "net_circle_weight.pdf"), width = 8, height = 8)
  netVisual_circle(cc@net$weight, color.use = cols,
                   vertex.weight = as.numeric(table(cc@idents)),
                   label.edge = FALSE)
  dev.off()

  pdf(file.path(od, "net_circle_count.pdf"), width = 8, height = 8)
  netVisual_circle(cc@net$count, color.use = cols,
                   vertex.weight = as.numeric(table(cc@idents)),
                   label.edge = FALSE)
  dev.off()

  pdf(file.path(od, "net_heatmap_weight.pdf"), width = 9, height = 7)
  netVisual_heatmap(cc, measure = "weight", color.use = cols)
  dev.off()

  pdf(file.path(od, "net_heatmap_count.pdf"), width = 9, height = 7)
  netVisual_heatmap(cc, measure = "count", color.use = cols)
  dev.off()

  # Optional: signaling role / centrality plots (version-dependent but usually present)
  pdf(file.path(od, "signalingRole_network.pdf"), width = 10, height = 6)
  tryCatch(netAnalysis_signalingRole_network(cc, color.use = cols), error = function(e) NULL)
  dev.off()

  pdf(file.path(od, "signalingRole_heatmap.pdf"), width = 10, height = 8)
  tryCatch(netAnalysis_signalingRole_heatmap(cc, pattern = "all"), error = function(e) NULL)
  dev.off()

  invisible(cc)
}

# -------------------- RUN ALL 4 CONDITIONS --------------------
run_cellchat_all_conditions <- function(seu, out_root, prefix = "pooled",
                                       raw_use = TRUE, type_use = "truncatedMean") {
  cc_root <- safe_dir(file.path(out_root, "CellChat", prefix))

  conds <- tidyr::crossing(
    Delivery_class = c("Term","Preterm"),
    MNPBurden = c("Low","High")
  ) %>% mutate(cond = paste(Delivery_class, MNPBurden, sep = "_"))

  cc_list <- list()
  for (i in seq_len(nrow(conds))) {
    delivery <- conds$Delivery_class[i]
    burden   <- conds$MNPBurden[i]
    cond     <- conds$cond[i]

    cc_list[[cond]] <- run_cellchat_one_v8(
      seu, delivery = delivery, burden = burden,
      outdir = file.path(cc_root, cond),
      min_cells_ct = min_cells_ct,
      min_cells_comm = min_cells_comm,
      cap_cells_per_celltype = cap_cells_per_celltype,
      raw_use = raw_use,
      type_use = type_use
    )
  }

  saveRDS(cc_list, file.path(cc_root, "cellchat_all_conditions.list.rds"))
  cc_list
}

# -------------------- COMPARE HIGH VS LOW --------------------
compare_high_low <- function(cc_low, cc_high, outdir, label) {
  od <- safe_dir(outdir)

  merged <- mergeCellChat(list(Low = cc_low, High = cc_high), add.names = c("Low","High"))
  saveRDS(merged, file.path(od, "cellchat_merged.rds"))

  pdf(file.path(od, "diff_circle_count.pdf"), width = 8, height = 8)
  netVisual_diffInteraction(merged, weight.scale = TRUE, measure = "count")
  dev.off()

  pdf(file.path(od, "diff_circle_weight.pdf"), width = 8, height = 8)
  netVisual_diffInteraction(merged, weight.scale = TRUE, measure = "weight")
  dev.off()

  pdf(file.path(od, "rankNet_comparison.pdf"), width = 10, height = 6)
  ok <- tryCatch({
    rankNet(merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)
    TRUE
  }, error = function(e) FALSE)
  if (!ok) rankNet(merged, mode = "comparison", stacked = TRUE, do.stat = FALSE)
  dev.off()

  sink(file.path(od, "compareInteractions_summary.txt"))
  cat("Label:", label, "\n\n")
  cat("compareInteractions(count):\n")
  print(tryCatch(compareInteractions(merged, group = c(1,2), measure = "count"), error = function(e) e))
  cat("\ncompareInteractions(weight):\n")
  print(tryCatch(compareInteractions(merged, group = c(1,2), measure = "weight"), error = function(e) e))
  sink()

  comm_merged <- tryCatch(subsetCommunication(merged), error = function(e) NULL)
  if (!is.null(comm_merged) && is.data.frame(comm_merged)) {
    write_csv(as_tibble(comm_merged), file.path(od, "merged_subsetCommunication.csv"))
  }

  invisible(merged)
}

# ============================================================
# RUN EVERYTHING
# ============================================================
safe_dir(out_root)
obj_lr <- prep_obj_for_cellchat(obj)

# pooled run (all platforms together)
cc_list_pooled <- run_cellchat_all_conditions(obj_lr, out_root = out_root, prefix = "pooled",
                                              raw_use = TRUE, type_use = "truncatedMean")

# pooled comparisons
comp_root <- safe_dir(file.path(out_root, "CellChat_compare", "pooled"))
merged_term_pooled <- compare_high_low(
  cc_list_pooled[["Term_Low"]], cc_list_pooled[["Term_High"]],
  outdir = file.path(comp_root, "Term_High_vs_Low"),
  label = "pooled_Term_High_vs_Low"
)
merged_preterm_pooled <- compare_high_low(
  cc_list_pooled[["Preterm_Low"]], cc_list_pooled[["Preterm_High"]],
  outdir = file.path(comp_root, "Preterm_High_vs_Low"),
  label = "pooled_Preterm_High_vs_Low"
)

# optional per-platform runs + comparisons
if (run_per_platform) {
  for (plat in platform_levels) {
    message("[CellChat] PLATFORM split: ", plat)
    seu_p <- subset(obj_lr, subset = platform == plat)
    if (ncol(seu_p) == 0) next

    prefix <- paste0("platform_", plat)
    cc_list_plat <- run_cellchat_all_conditions(seu_p, out_root = out_root, prefix = prefix,
                                                raw_use = TRUE, type_use = "truncatedMean")

    comp_root_p <- safe_dir(file.path(out_root, "CellChat_compare", prefix))
    compare_high_low(
      cc_list_plat[["Term_Low"]], cc_list_plat[["Term_High"]],
      outdir = file.path(comp_root_p, "Term_High_vs_Low"),
      label = paste0(prefix, "_Term_High_vs_Low")
    )
    compare_high_low(
      cc_list_plat[["Preterm_Low"]], cc_list_plat[["Preterm_High"]],
      outdir = file.path(comp_root_p, "Preterm_High_vs_Low"),
      label = paste0(prefix, "_Preterm_High_vs_Low")
    )
  }
}

message("DONE. Outputs under: ", normalizePath(out_root))

## crashed

# ============================================================
# FULL DROP-IN v9: CellChat full dataset (SCT@data + raw.use=TRUE)
# Scales to big n by capping per cell type + robust resume + per-condition tryCatch
# ============================================================

set.seed(1)
options(future.globals.maxSize = 300e9)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr); library(tidyr); library(stringr); library(tibble); library(readr)
  library(Matrix)
  library(CellChat)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# -------------------- CONFIG --------------------
out_root <- "MNP_interactomes_CellChat_v9_full"

run_per_platform <- TRUE
platform_levels <- c("Xenium","VisiumHD","VisiumV1")

min_cells_ct   <- 50
min_cells_comm <- 10

# CRITICAL FOR SCALE:
cap_cells_per_celltype <- 25000   # <- start here; if stable, try 50000

drop_lowconf <- TRUE  # drop RCTD_lowconf everywhere
resume_mode  <- TRUE  # skip completed conditions

ord_celltypes <- c(
  "Syncytiotrophoblast","Villous Cytotrophoblast","Extravillous Trophoblast",
  "Fibroblast","Perivascular Cell","Placental Endothelial",
  "Hofbauer Cell","Macrophage","NK Cell","Decidual NK Cell","ILC3",
  "Plasma Cell","B Cell","T Cell","Treg","RCTD_lowconf"
)

ct_colors <- c(
  "Syncytiotrophoblast"     = "#5050FF",
  "Villous Cytotrophoblast" = "#CE3D32",
  "Extravillous Trophoblast"= "#749B58",
  "Fibroblast"              = "#F0E685",
  "Perivascular Cell"       = "#466983",
  "Placental Endothelial"   = "#BA6338",
  "Hofbauer Cell"           = "#5DB1DD",
  "Macrophage"              = "#802268",
  "NK Cell"                 = "#6BD76B",
  "Decidual NK Cell"        = "#D595A7",
  "ILC3"                    = "#924822",
  "Plasma Cell"             = "#837B8D",
  "B Cell"                  = "#C75127",
  "T Cell"                  = "#D58F5C",
  "Treg"                    = "#7A65A5",
  "RCTD_lowconf"            = "#BFBFBF"
)

# -------------------- UTIL --------------------
safe_dir <- function(path) { dir.create(path, showWarnings = FALSE, recursive = TRUE); path }

get_sct_data <- function(seu) {
  if (!"SCT" %in% names(seu@assays)) stop("No SCT assay found.")
  tryCatch(
    GetAssayData(seu, assay = "SCT", layer = "data"),
    error = function(e) GetAssayData(seu, assay = "SCT", slot = "data")
  )
}

downsample_per_celltype <- function(seu, cap, ct_col = "nn.predicted.celltype") {
  if (is.null(cap)) return(seu)
  md <- seu[[]]
  ct <- as.character(md[[ct_col]])
  names(ct) <- rownames(md)
  ct <- ct[!is.na(ct)]
  sp <- split(names(ct), ct)
  keep <- unlist(lapply(sp, function(cells) if (length(cells) <= cap) cells else sample(cells, cap)),
                 use.names = FALSE)
  subset(seu, cells = keep)
}

parse_subject_platform <- function(seu) {
  md <- seu[[]]
  stopifnot("sample_id_std" %in% colnames(md))
  md <- md %>%
    mutate(
      subject = if_else(str_detect(sample_id_std, "_"),
                        str_replace(sample_id_std, "_.*$", ""),
                        sample_id_std),
      platform_parsed = if_else(str_detect(sample_id_std, "_"),
                                str_replace(sample_id_std, "^.*_", ""),
                                "ref"),
      platform_parsed = recode(platform_parsed,
                               "xen"="Xenium","hd"="VisiumHD","v1"="VisiumV1",
                               .default = platform_parsed)
    )
  seu@meta.data <- md
  seu
}

refactor_celltypes <- function(seu, levels_master) {
  md <- seu[[]]
  md$nn.predicted.celltype <- as.character(md$nn.predicted.celltype)
  observed <- unique(md$nn.predicted.celltype)
  observed <- observed[!is.na(observed)]
  lev <- c(intersect(levels_master, observed), setdiff(observed, levels_master))
  md$nn.predicted.celltype <- factor(md$nn.predicted.celltype, levels = lev)
  md$nn.predicted.celltype <- droplevels(md$nn.predicted.celltype)
  seu@meta.data <- md
  Idents(seu) <- "nn.predicted.celltype"
  seu
}

ct_cols_for <- function(cellchat_obj) {
  grps <- levels(cellchat_obj@idents)
  cols <- ct_colors[grps]; cols[is.na(cols)] <- "#808080"
  as.character(cols)
}

map_mixed_ids_to_symbol <- function(mat) {
  if (!inherits(mat, "dgCMatrix")) mat <- as(mat, "dgCMatrix")
  rn <- rownames(mat)
  is_ensg <- grepl("^ENSG", rn)

  rn_std <- rn
  rn_std[!is_ensg] <- toupper(rn_std[!is_ensg])

  if (any(is_ensg)) {
    ensg <- rn[is_ensg]
    ensg_clean <- sub("\\..*$", "", ensg)
    ensg_clean <- sub("\\|.*$", "", ensg)

    sym <- suppressMessages(
      AnnotationDbi::mapIds(
        org.Hs.eg.db,
        keys = unique(ensg_clean),
        column = "SYMBOL",
        keytype = "ENSEMBL",
        multiVals = "first"
      )
    )
    mapped <- toupper(as.character(sym[ensg_clean]))
    rn_std[is_ensg] <- ifelse(!is.na(mapped) & mapped != "", mapped, toupper(ensg))
  }

  rownames(mat) <- rn_std

  if (any(duplicated(rownames(mat)))) {
    grp <- factor(rownames(mat), levels = unique(rownames(mat)))
    G <- Matrix::sparse.model.matrix(~ 0 + grp)
    mat <- Matrix::t(G) %*% mat
    rownames(mat) <- levels(grp)
    mat <- as(mat, "dgCMatrix")
  }
  mat
}

subsetData_perm <- function(cc) {
  fm <- names(formals(CellChat::subsetData))
  args <- list(object = cc)
  if ("thresh" %in% fm) args$thresh <- 0
  if ("thresh.pc" %in% fm) args$thresh.pc <- 0
  do.call(CellChat::subsetData, args)
}

force_LRsig_from_DB <- function(cc, genes_present, DB = CellChatDB.human) {
  genes_present <- unique(toupper(genes_present))
  lr0 <- DB$interaction %>%
    transmute(
      interaction_name   = interaction_name,
      pathway_name       = pathway_name,
      ligand             = toupper(ligand),
      receptor           = toupper(receptor),
      interaction_name_2 = if ("interaction_name_2" %in% colnames(DB$interaction)) interaction_name_2 else interaction_name
    ) %>%
    filter(ligand %in% genes_present, receptor %in% genes_present)
  cc@LR$LRsig <- as.data.frame(lr0)
  cc
}

db_overlap_diag <- function(gene_names, DB = CellChatDB.human) {
  mg  <- unique(toupper(gene_names))
  dbL <- unique(toupper(DB$interaction$ligand))
  dbR <- unique(toupper(DB$interaction$receptor))
  tibble(
    n_genes = length(mg),
    overlap_ligands = sum(dbL %in% mg),
    overlap_receptors = sum(dbR %in% mg),
    n_db_ligands = length(dbL),
    n_db_receptors = length(dbR)
  )
}

# Safe writer for small matrices (group x group)
write_small_mat_csv <- function(x, path) {
  x2 <- tryCatch(as.matrix(x), error = function(e) x)
  write.csv(as.data.frame(x2), path, row.names = TRUE)
}

# -------------------- PREP OBJECT --------------------
prep_obj_for_cellchat <- function(seu) {
  if (!"SCT" %in% names(seu@assays)) stop("Object must contain SCT assay.")
  DefaultAssay(seu) <- "SCT"
  seu <- parse_subject_platform(seu)

  md <- seu[[]]
  need <- c("platform","MNPBurden","Delivery_class","nn.predicted.celltype","sample_id_std")
  miss <- setdiff(need, colnames(md))
  if (length(miss)) stop("Missing metadata columns: ", paste(miss, collapse=", "))

  seu <- subset(
    seu,
    subset =
      platform %in% platform_levels &
      sample_id_std != "placental_ref" &
      !is.na(MNPBurden) & MNPBurden %in% c("High","Low") &
      !is.na(Delivery_class) & Delivery_class %in% c("Term","Preterm")
  )

  if (drop_lowconf) {
    seu <- subset(seu, subset = nn.predicted.celltype != "RCTD_lowconf")
  }

  seu <- refactor_celltypes(seu, ord_celltypes)
  seu
}

# -------------------- RUN ONE CONDITION --------------------
run_cellchat_one_v9 <- function(seu, delivery, burden, outdir,
                                min_cells_ct = 50, min_cells_comm = 10,
                                cap_cells_per_celltype = 25000,
                                raw_use = TRUE,
                                type_use = c("median","thresholdedMean","truncatedMean","triMean")) {
  type_use <- match.arg(type_use)
  od <- safe_dir(outdir)
  done_file <- file.path(od, "cellchat.rds")

  if (resume_mode && file.exists(done_file)) {
    message("[CellChat] SKIP (exists): ", delivery, "_", burden)
    return(readRDS(done_file))
  }

  message("[CellChat] ", delivery, "_", burden, " -> ", od)

  sub <- subset(seu, subset = Delivery_class == delivery & MNPBurden == burden)
  if (ncol(sub) == 0) stop("No cells for ", delivery, "_", burden)

  # cap (CRITICAL for stability)
  sub <- downsample_per_celltype(sub, cap_cells_per_celltype)

  # filter tiny groups
  sub$nn.predicted.celltype <- as.character(sub$nn.predicted.celltype)
  tab <- table(sub$nn.predicted.celltype, useNA = "no")
  keep_ct <- names(tab)[tab >= min_cells_ct]
  sub <- subset(sub, subset = nn.predicted.celltype %in% keep_ct)
  if (length(keep_ct) < 2) stop("Too few cell types >= ", min_cells_ct, " cells in ", delivery, "_", burden)

  sub <- refactor_celltypes(sub, ord_celltypes)

  write_csv(
    tibble(celltype = names(table(sub$nn.predicted.celltype)),
           n_cells  = as.integer(table(sub$nn.predicted.celltype))),
    file.path(od, "idents_cell_counts.csv")
  )

  # SCT@data -> mapped symbols
  data.input <- get_sct_data(sub)
  if (!inherits(data.input, "dgCMatrix")) data.input <- as(data.input, "dgCMatrix")
  data.mapped <- map_mixed_ids_to_symbol(data.input)
  data.mapped <- data.mapped[!grepl("^ENSG", rownames(data.mapped)), , drop = FALSE]

  diag_dims <- tibble(n_cells = ncol(data.mapped), n_genes_mapped = nrow(data.mapped))
  diag_db   <- db_overlap_diag(rownames(data.mapped), CellChatDB.human)
  write_csv(diag_dims, file.path(od, "diag_matrix_dims.csv"))
  write_csv(diag_db,   file.path(od, "diag_db_overlap.csv"))

  message("[DIAG] n_cells=", diag_dims$n_cells,
          " | n_genes_mapped=", diag_dims$n_genes_mapped,
          " | DB overlap lig=", diag_db$overlap_ligands, "/", diag_db$n_db_ligands,
          " rec=", diag_db$overlap_receptors, "/", diag_db$n_db_receptors)

  meta.input <- sub[[]] %>% mutate(nn.predicted.celltype = droplevels(factor(as.character(nn.predicted.celltype))))

  cc <- createCellChat(object = data.mapped, meta = meta.input, group.by = "nn.predicted.celltype")
  cc@DB <- CellChatDB.human
  cc@idents <- droplevels(cc@idents)

  cc <- subsetData_perm(cc)

  cc <- force_LRsig_from_DB(cc, genes_present = rownames(cc@data), DB = CellChatDB.human)
  message("[DIAG] Forced LRsig rows: ", nrow(cc@LR$LRsig))
  writeLines(paste0("Forced LRsig rows: ", nrow(cc@LR$LRsig)), file.path(od, "forced_LRsig_n.txt"))

  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)

  # CRITICAL
  cc <- computeCommunProb(cc, raw.use = raw_use, type = type_use, trim = 0, LR.use = cc@LR$LRsig, population.size = FALSE)

  cc <- filterCommunication(cc, min.cells = min_cells_comm)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc)

  saveRDS(cc, done_file)

  # Tables (small group x group)
  write_small_mat_csv(cc@net$weight, file.path(od, "net_weight_matrix.csv"))
  write_small_mat_csv(cc@net$count,  file.path(od, "net_count_matrix.csv"))

  comm <- tryCatch(subsetCommunication(cc), error = function(e) NULL)
  if (!is.null(comm) && is.data.frame(comm)) write_csv(as_tibble(comm), file.path(od, "communications_table.csv"))

  # Visuals
  cols <- ct_cols_for(cc)

  pdf(file.path(od, "net_circle_weight.pdf"), width = 8, height = 8)
  netVisual_circle(cc@net$weight, color.use = cols,
                   vertex.weight = as.numeric(table(cc@idents)),
                   label.edge = FALSE)
  dev.off()

  pdf(file.path(od, "net_circle_count.pdf"), width = 8, height = 8)
  netVisual_circle(cc@net$count, color.use = cols,
                   vertex.weight = as.numeric(table(cc@idents)),
                   label.edge = FALSE)
  dev.off()

  pdf(file.path(od, "net_heatmap_weight.pdf"), width = 9, height = 7)
  netVisual_heatmap(cc, measure = "weight", color.use = cols)
  dev.off()

  pdf(file.path(od, "net_heatmap_count.pdf"), width = 9, height = 7)
  netVisual_heatmap(cc, measure = "count", color.use = cols)
  dev.off()

  pdf(file.path(od, "signalingRole_network.pdf"), width = 10, height = 6)
  tryCatch(netAnalysis_signalingRole_network(cc, color.use = cols), error = function(e) NULL)
  dev.off()

  pdf(file.path(od, "signalingRole_heatmap.pdf"), width = 10, height = 8)
  tryCatch(netAnalysis_signalingRole_heatmap(cc, pattern = "all"), error = function(e) NULL)
  dev.off()

  invisible(cc)
}

# -------------------- RUN ALL 4 CONDITIONS --------------------
run_cellchat_all_conditions_v9 <- function(seu, out_root, prefix = "pooled",
                                          raw_use = TRUE, type_use = "truncatedMean") {
  cc_root <- safe_dir(file.path(out_root, "CellChat", prefix))

  conds <- tidyr::crossing(
    Delivery_class = c("Term","Preterm"),
    MNPBurden = c("Low","High")
  ) %>% mutate(cond = paste(Delivery_class, MNPBurden, sep = "_"))

  cc_list <- list()
  for (i in seq_len(nrow(conds))) {
    delivery <- conds$Delivery_class[i]
    burden   <- conds$MNPBurden[i]
    cond     <- conds$cond[i]

    cc_list[[cond]] <- tryCatch(
      run_cellchat_one_v9(
        seu, delivery = delivery, burden = burden,
        outdir = file.path(cc_root, cond),
        min_cells_ct = min_cells_ct,
        min_cells_comm = min_cells_comm,
        cap_cells_per_celltype = cap_cells_per_celltype,
        raw_use = raw_use,
        type_use = type_use
      ),
      error = function(e) {
        od <- safe_dir(file.path(cc_root, cond))
        writeLines(conditionMessage(e), file.path(od, "FAIL_condition.txt"))
        message("[CellChat FAIL] ", cond, " :: ", conditionMessage(e))
        NULL
      }
    )
  }

  saveRDS(cc_list, file.path(cc_root, "cellchat_all_conditions.list.rds"))
  cc_list
}

# -------------------- COMPARE HIGH VS LOW --------------------
compare_high_low <- function(cc_low, cc_high, outdir, label) {
  od <- safe_dir(outdir)

  if (is.null(cc_low) || is.null(cc_high)) {
    writeLines("One or both CellChat objects are NULL (failed earlier).", file.path(od, "FAIL_compare_missing.txt"))
    return(invisible(NULL))
  }

  merged <- mergeCellChat(list(Low = cc_low, High = cc_high), add.names = c("Low","High"))
  saveRDS(merged, file.path(od, "cellchat_merged.rds"))

  pdf(file.path(od, "diff_circle_count.pdf"), width = 8, height = 8)
  netVisual_diffInteraction(merged, weight.scale = TRUE, measure = "count")
  dev.off()

  pdf(file.path(od, "diff_circle_weight.pdf"), width = 8, height = 8)
  netVisual_diffInteraction(merged, weight.scale = TRUE, measure = "weight")
  dev.off()

  pdf(file.path(od, "rankNet_comparison.pdf"), width = 10, height = 6)
  ok <- tryCatch({
    rankNet(merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)
    TRUE
  }, error = function(e) FALSE)
  if (!ok) rankNet(merged, mode = "comparison", stacked = TRUE, do.stat = FALSE)
  dev.off()

  sink(file.path(od, "compareInteractions_summary.txt"))
  cat("Label:", label, "\n\n")
  cat("compareInteractions(count):\n")
  print(tryCatch(compareInteractions(merged, group = c(1,2), measure = "count"), error = function(e) e))
  cat("\ncompareInteractions(weight):\n")
  print(tryCatch(compareInteractions(merged, group = c(1,2), measure = "weight"), error = function(e) e))
  sink()

  comm_merged <- tryCatch(subsetCommunication(merged), error = function(e) NULL)
  if (!is.null(comm_merged) && is.data.frame(comm_merged)) {
    write_csv(as_tibble(comm_merged), file.path(od, "merged_subsetCommunication.csv"))
  }

  invisible(merged)
}

# ============================================================
# RUN
# ============================================================
safe_dir(out_root)
obj_lr <- prep_obj_for_cellchat(obj)

# pooled run
cc_list_pooled <- run_cellchat_all_conditions_v9(obj_lr, out_root = out_root, prefix = "pooled",
                                                 raw_use = TRUE, type_use = "truncatedMean")

# pooled comparisons
comp_root <- safe_dir(file.path(out_root, "CellChat_compare", "pooled"))
merged_term_pooled <- compare_high_low(
  cc_list_pooled[["Term_Low"]], cc_list_pooled[["Term_High"]],
  outdir = file.path(comp_root, "Term_High_vs_Low"),
  label = "pooled_Term_High_vs_Low"
)
merged_preterm_pooled <- compare_high_low(
  cc_list_pooled[["Preterm_Low"]], cc_list_pooled[["Preterm_High"]],
  outdir = file.path(comp_root, "Preterm_High_vs_Low"),
  label = "pooled_Preterm_High_vs_Low"
)

# optional per-platform runs
if (run_per_platform) {
  for (plat in platform_levels) {
    message("[CellChat] PLATFORM split: ", plat)
    seu_p <- subset(obj_lr, subset = platform == plat)
    if (ncol(seu_p) == 0) next

    prefix <- paste0("platform_", plat)
    cc_list_plat <- run_cellchat_all_conditions_v9(seu_p, out_root = out_root, prefix = prefix,
                                                   raw_use = TRUE, type_use = "truncatedMean")

    comp_root_p <- safe_dir(file.path(out_root, "CellChat_compare", prefix))
    compare_high_low(
      cc_list_plat[["Term_Low"]], cc_list_plat[["Term_High"]],
      outdir = file.path(comp_root_p, "Term_High_vs_Low"),
      label = paste0(prefix, "_Term_High_vs_Low")
    )
    compare_high_low(
      cc_list_plat[["Preterm_Low"]], cc_list_plat[["Preterm_High"]],
      outdir = file.path(comp_root_p, "Preterm_High_vs_Low"),
      label = paste0(prefix, "_Preterm_High_vs_Low")
    )
  }
}

message("DONE. Outputs under: ", normalizePath(out_root))



run_cellchat_one_v9(
  obj_lr, "Preterm", "High",
  outdir = file.path(out_root, "SMOKETEST_Preterm_High_cap25k"),
  cap_cells_per_celltype = 25000,
  raw_use = TRUE,
  type_use = "truncatedMean",
  min_cells_ct = 50,
  min_cells_comm = 10
)






align_cellchat_groups <- function(cc1, cc2) {
  g1 <- levels(cc1@idents)
  g2 <- levels(cc2@idents)
  g_keep <- intersect(g1, g2)

  if (length(g_keep) < 2) stop("Too few overlapping cell groups to compare.")

  # subset data + nets to common groups
  keep_net <- function(cc, g_keep) {
    # idents
    cc@idents <- droplevels(factor(cc@idents, levels = g_keep))

    # net slots (weight/count are matrices)
    if (!is.null(cc@net$weight)) cc@net$weight <- cc@net$weight[g_keep, g_keep, drop = FALSE]
    if (!is.null(cc@net$count))  cc@net$count  <- cc@net$count[g_keep, g_keep, drop = FALSE]

    # pathway-level netP is more complex; safest: recompute after restricting idents
    cc
  }

  cc1a <- keep_net(cc1, g_keep)
  cc2a <- keep_net(cc2, g_keep)

  list(cc1 = cc1a, cc2 = cc2a, groups = g_keep)
}

compare_high_low_safe <- function(cc_low, cc_high, outdir, label) {
  od <- safe_dir(outdir)

  ali <- align_cellchat_groups(cc_low, cc_high)
  cc_low2  <- ali$cc1
  cc_high2 <- ali$cc2
  writeLines(paste("Common groups:", paste(ali$groups, collapse = ", ")),
             file.path(od, "common_groups.txt"))

  merged <- mergeCellChat(list(Low = cc_low2, High = cc_high2), add.names = c("Low","High"))
  saveRDS(merged, file.path(od, "cellchat_merged.rds"))

  pdf(file.path(od, "diff_circle_count.pdf"), width = 8, height = 8)
  netVisual_diffInteraction(merged, weight.scale = TRUE, measure = "count")
  dev.off()

  pdf(file.path(od, "diff_circle_weight.pdf"), width = 8, height = 8)
  netVisual_diffInteraction(merged, weight.scale = TRUE, measure = "weight")
  dev.off()

  pdf(file.path(od, "rankNet_comparison.pdf"), width = 10, height = 6)
  ok <- tryCatch({
    rankNet(merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)
    TRUE
  }, error = function(e) FALSE)
  if (!ok) rankNet(merged, mode = "comparison", stacked = TRUE, do.stat = FALSE)
  dev.off()

  sink(file.path(od, "compareInteractions_summary.txt"))
  cat("Label:", label, "\n\n")
  cat("compareInteractions(count):\n")
  print(tryCatch(compareInteractions(merged, group = c(1,2), measure = "count"), error = function(e) e))
  cat("\ncompareInteractions(weight):\n")
  print(tryCatch(compareInteractions(merged, group = c(1,2), measure = "weight"), error = function(e) e))
  sink()

  comm_merged <- tryCatch(subsetCommunication(merged), error = function(e) NULL)
  if (!is.null(comm_merged) && is.data.frame(comm_merged)) {
    write_csv(as_tibble(comm_merged), file.path(od, "merged_subsetCommunication.csv"))
  }

  invisible(merged)
}

merged_term_pooled <- compare_high_low_safe(
  cc_list_pooled[["Term_Low"]], cc_list_pooled[["Term_High"]],
  outdir = file.path(comp_root, "Term_High_vs_Low"),
  label = "pooled_Term_High_vs_Low"
)






























# ============================================================
# FULL DROP-IN v10: CellChat rerun with min_cells_ct=500 + PNG(300dpi) + UMAP colors
# - Uses obj@assays[["SCT"]]@data
# - raw.use=TRUE (required for your data)
# - min_cells_ct=500 (drops rare types like B cells everywhere)
# - Saves plots as PNG 300 dpi (white background) using ct_colors
# - Robust comparisons: aligns common groups before merge to avoid non-conformable arrays
# ============================================================
options(mc.cores = 1)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
set.seed(1)
options(future.globals.maxSize = 300e9)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr); library(tidyr); library(stringr); library(tibble); library(readr)
  library(Matrix)
  library(CellChat)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# -------------------- CONFIG --------------------
out_root <- "MNP_interactomes_CellChat_v10_full_min500_png"

run_per_platform <- TRUE
platform_levels <- c("Xenium","VisiumHD","VisiumV1")

min_cells_ct   <- 500   # <- your preference; forces B cells out if rare
min_cells_comm <- 10

cap_cells_per_celltype <- 10000   # stability for big conditions; tune up/down

resume_mode <- TRUE   # skip condition if cellchat.rds exists

# celltype order + colors (your UMAP palette)
ord_celltypes <- c(
  "Syncytiotrophoblast","Villous Cytotrophoblast","Extravillous Trophoblast",
  "Fibroblast","Perivascular Cell","Placental Endothelial",
  "Hofbauer Cell","Macrophage","NK Cell","Decidual NK Cell","ILC3",
  "Plasma Cell","B Cell","T Cell","Treg","RCTD_lowconf"
)

ct_colors <- c(
  "Syncytiotrophoblast"     = "#5050FF",
  "Villous Cytotrophoblast" = "#CE3D32",
  "Extravillous Trophoblast"= "#749B58",
  "Fibroblast"              = "#F0E685",
  "Perivascular Cell"       = "#466983",
  "Placental Endothelial"   = "#BA6338",
  "Hofbauer Cell"           = "#5DB1DD",
  "Macrophage"              = "#802268",
  "NK Cell"                 = "#6BD76B",
  "Decidual NK Cell"        = "#D595A7",
  "ILC3"                    = "#924822",
  "Plasma Cell"             = "#837B8D",
  "B Cell"                  = "#C75127",
  "T Cell"                  = "#D58F5C",
  "Treg"                    = "#7A65A5",
  "RCTD_lowconf"            = "#BFBFBF"
)

# -------------------- UTIL --------------------
safe_dir <- function(path) { dir.create(path, showWarnings = FALSE, recursive = TRUE); path }

# Seurat v5-safe accessor for SCT@data
get_sct_data <- function(seu) {
  if (!"SCT" %in% names(seu@assays)) stop("No SCT assay found.")
  tryCatch(
    GetAssayData(seu, assay = "SCT", layer = "data"),
    error = function(e) GetAssayData(seu, assay = "SCT", slot = "data")
  )
}

downsample_per_celltype <- function(seu, cap, ct_col = "nn.predicted.celltype") {
  if (is.null(cap)) return(seu)
  md <- seu[[]]
  ct <- as.character(md[[ct_col]])
  names(ct) <- rownames(md)
  ct <- ct[!is.na(ct)]
  sp <- split(names(ct), ct)
  keep <- unlist(lapply(sp, function(cells) if (length(cells) <= cap) cells else sample(cells, cap)),
                 use.names = FALSE)
  subset(seu, cells = keep)
}

parse_subject_platform <- function(seu) {
  md <- seu[[]]
  stopifnot("sample_id_std" %in% colnames(md))
  md <- md %>%
    mutate(
      subject = if_else(str_detect(sample_id_std, "_"),
                        str_replace(sample_id_std, "_.*$", ""),
                        sample_id_std),
      platform_parsed = if_else(str_detect(sample_id_std, "_"),
                                str_replace(sample_id_std, "^.*_", ""),
                                "ref"),
      platform_parsed = recode(platform_parsed,
                               "xen"="Xenium","hd"="VisiumHD","v1"="VisiumV1",
                               .default = platform_parsed)
    )
  seu@meta.data <- md
  seu
}

refactor_celltypes <- function(seu, levels_master) {
  md <- seu[[]]
  md$nn.predicted.celltype <- as.character(md$nn.predicted.celltype)
  observed <- unique(md$nn.predicted.celltype)
  observed <- observed[!is.na(observed)]
  lev <- c(intersect(levels_master, observed), setdiff(observed, levels_master))
  md$nn.predicted.celltype <- factor(md$nn.predicted.celltype, levels = lev)
  md$nn.predicted.celltype <- droplevels(md$nn.predicted.celltype)
  seu@meta.data <- md
  Idents(seu) <- "nn.predicted.celltype"
  seu
}

ct_cols_for <- function(cellchat_obj) {
  grps <- levels(cellchat_obj@idents)
  cols <- ct_colors[grps]; cols[is.na(cols)] <- "#808080"
  as.character(cols)
}

# Mixed IDs -> SYMBOL, uppercase, collapse duplicates (sparse-safe)
map_mixed_ids_to_symbol <- function(mat) {
  if (!inherits(mat, "dgCMatrix")) mat <- as(mat, "dgCMatrix")
  rn <- rownames(mat)
  is_ensg <- grepl("^ENSG", rn)

  rn_std <- rn
  rn_std[!is_ensg] <- toupper(rn_std[!is_ensg])

  if (any(is_ensg)) {
    ensg <- rn[is_ensg]
    ensg_clean <- sub("\\..*$", "", ensg)
    ensg_clean <- sub("\\|.*$", "", ensg)

    sym <- suppressMessages(
      AnnotationDbi::mapIds(
        org.Hs.eg.db,
        keys = unique(ensg_clean),
        column = "SYMBOL",
        keytype = "ENSEMBL",
        multiVals = "first"
      )
    )
    mapped <- toupper(as.character(sym[ensg_clean]))
    rn_std[is_ensg] <- ifelse(!is.na(mapped) & mapped != "", mapped, toupper(ensg))
  }

  rownames(mat) <- rn_std

  if (any(duplicated(rownames(mat)))) {
    grp <- factor(rownames(mat), levels = unique(rownames(mat)))
    G <- Matrix::sparse.model.matrix(~ 0 + grp)
    mat <- Matrix::t(G) %*% mat
    rownames(mat) <- levels(grp)
    mat <- as(mat, "dgCMatrix")
  }
  mat
}

subsetData_perm <- function(cc) {
  fm <- names(formals(CellChat::subsetData))
  args <- list(object = cc)
  if ("thresh" %in% fm) args$thresh <- 0
  if ("thresh.pc" %in% fm) args$thresh.pc <- 0
  do.call(CellChat::subsetData, args)
}

force_LRsig_from_DB <- function(cc, genes_present, DB = CellChatDB.human) {
  genes_present <- unique(toupper(genes_present))
  lr0 <- DB$interaction %>%
    transmute(
      interaction_name   = interaction_name,
      pathway_name       = pathway_name,
      ligand             = toupper(ligand),
      receptor           = toupper(receptor),
      interaction_name_2 = if ("interaction_name_2" %in% colnames(DB$interaction)) interaction_name_2 else interaction_name
    ) %>%
    filter(ligand %in% genes_present, receptor %in% genes_present)
  cc@LR$LRsig <- as.data.frame(lr0)
  cc
}

db_overlap_diag <- function(gene_names, DB = CellChatDB.human) {
  mg  <- unique(toupper(gene_names))
  dbL <- unique(toupper(DB$interaction$ligand))
  dbR <- unique(toupper(DB$interaction$receptor))
  tibble(
    n_genes = length(mg),
    overlap_ligands = sum(dbL %in% mg),
    overlap_receptors = sum(dbR %in% mg),
    n_db_ligands = length(dbL),
    n_db_receptors = length(dbR)
  )
}

# ---- PNG wrappers (300 dpi) ----
png_open <- function(path, w_in = 8, h_in = 8, res = 300) {
  png(filename = path, width = w_in, height = h_in, units = "in", res = res, bg = "white")
}
png_close <- function() dev.off()

# Robust plotting: protect against "heatmap based on a single object" chatter/errors
plot_safe <- function(expr, on_fail = NULL) {
  tryCatch(expr, error = function(e) { message("[PLOT FAIL] ", conditionMessage(e)); on_fail })
}

# -------------------- PREP OBJECT --------------------
prep_obj_for_cellchat <- function(seu) {
  if (!"SCT" %in% names(seu@assays)) stop("Object must contain SCT assay.")
  DefaultAssay(seu) <- "SCT"
  seu <- parse_subject_platform(seu)

  md <- seu[[]]
  need <- c("platform","MNPBurden","Delivery_class","nn.predicted.celltype","sample_id_std")
  miss <- setdiff(need, colnames(md))
  if (length(miss)) stop("Missing metadata columns: ", paste(miss, collapse=", "))

  seu <- subset(
    seu,
    subset =
      platform %in% platform_levels &
      sample_id_std != "placental_ref" &
      !is.na(MNPBurden) & MNPBurden %in% c("High","Low") &
      !is.na(Delivery_class) & Delivery_class %in% c("Term","Preterm")
  )

  seu <- refactor_celltypes(seu, ord_celltypes)
  seu
}

# -------------------- RUN ONE CONDITION --------------------
run_cellchat_one_v10 <- function(seu, delivery, burden, outdir,
                                 min_cells_ct = 500, min_cells_comm = 10,
                                 cap_default = 25000,
                                 cap_term_high = 10000,
                                 raw_use = TRUE,
                                 type_use = c("median","thresholdedMean","truncatedMean","triMean")) {

  type_use <- match.arg(type_use)
  od <- safe_dir(outdir)
  done_file <- file.path(od, "cellchat.rds")

  if (resume_mode && file.exists(done_file)) {
    message("[CellChat] SKIP (exists): ", delivery, "_", burden)
    return(readRDS(done_file))
  }

  # conditional cap
  cap_use <- if (delivery == "Term" && burden == "High") cap_term_high else cap_default
  message("[CellChat] ", delivery, "_", burden, " cap=", cap_use, " -> ", od)

  sub <- subset(seu, subset = Delivery_class == delivery & MNPBurden == burden)
  if (ncol(sub) == 0) stop("No cells for ", delivery, "_", burden)

  # drop lowconf always
  sub <- subset(sub, subset = nn.predicted.celltype != "RCTD_lowconf")

  # cap for stability
  sub <- downsample_per_celltype(sub, cap_use)

  # filter celltypes by min_cells_ct
  sub$nn.predicted.celltype <- as.character(sub$nn.predicted.celltype)
  tab <- table(sub$nn.predicted.celltype, useNA = "no")
  keep_ct <- names(tab)[tab >= min_cells_ct]
  sub <- subset(sub, subset = nn.predicted.celltype %in% keep_ct)
  if (length(keep_ct) < 2) stop("Too few cell types >= ", min_cells_ct, " cells in ", delivery, "_", burden)

  sub <- refactor_celltypes(sub, ord_celltypes)

  write_csv(
    tibble(celltype = names(table(sub$nn.predicted.celltype)),
           n_cells  = as.integer(table(sub$nn.predicted.celltype))),
    file.path(od, "idents_cell_counts.csv")
  )

  data.input <- get_sct_data(sub)
  if (!inherits(data.input, "dgCMatrix")) data.input <- as(data.input, "dgCMatrix")
  data.mapped <- map_mixed_ids_to_symbol(data.input)
  data.mapped <- data.mapped[!grepl("^ENSG", rownames(data.mapped)), , drop = FALSE]

  diag_dims <- tibble(n_cells = ncol(data.mapped), n_genes_mapped = nrow(data.mapped))
  diag_db   <- db_overlap_diag(rownames(data.mapped), CellChatDB.human)
  write_csv(diag_dims, file.path(od, "diag_matrix_dims.csv"))
  write_csv(diag_db,   file.path(od, "diag_db_overlap.csv"))

  message("[DIAG] n_cells=", diag_dims$n_cells,
          " | n_genes_mapped=", diag_dims$n_genes_mapped,
          " | DB overlap lig=", diag_db$overlap_ligands, "/", diag_db$n_db_ligands,
          " rec=", diag_db$overlap_receptors, "/", diag_db$n_db_receptors)

  meta.input <- sub[[]] %>% mutate(nn.predicted.celltype = droplevels(factor(as.character(nn.predicted.celltype))))

  cc <- createCellChat(object = data.mapped, meta = meta.input, group.by = "nn.predicted.celltype")
  cc@DB <- CellChatDB.human
  cc@idents <- droplevels(cc@idents)

  cc <- subsetData_perm(cc)
  cc <- force_LRsig_from_DB(cc, genes_present = rownames(cc@data), DB = CellChatDB.human)
  message("[DIAG] Forced LRsig rows: ", nrow(cc@LR$LRsig))
  writeLines(paste0("Forced LRsig rows: ", nrow(cc@LR$LRsig)), file.path(od, "forced_LRsig_n.txt"))

  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)

  run_compute <- function(cc_obj) {
    computeCommunProb(cc_obj, raw.use = raw_use, type = type_use, trim = 0,
                      LR.use = cc_obj@LR$LRsig, population.size = FALSE)
  }

  cc <- tryCatch(run_compute(cc), error = function(e) {
    # one retry with smaller cap if we hit the hard crash path
    message("[WARN] computeCommunProb failed once: ", conditionMessage(e))
    message("[RETRY] Rebuilding with cap=5000")
    sub2 <- downsample_per_celltype(sub, 5000)
    data2 <- map_mixed_ids_to_symbol(get_sct_data(sub2))
    data2 <- data2[!grepl("^ENSG", rownames(data2)), , drop=FALSE]
    meta2 <- sub2[[]] %>% mutate(nn.predicted.celltype = droplevels(factor(as.character(nn.predicted.celltype))))
    cc2 <- createCellChat(object = data2, meta = meta2, group.by = "nn.predicted.celltype")
    cc2@DB <- CellChatDB.human
    cc2@idents <- droplevels(cc2@idents)
    cc2 <- subsetData_perm(cc2)
    cc2 <- force_LRsig_from_DB(cc2, genes_present = rownames(cc2@data), DB = CellChatDB.human)
    cc2 <- identifyOverExpressedGenes(cc2)
    cc2 <- identifyOverExpressedInteractions(cc2)
    run_compute(cc2)
  })

  cc <- filterCommunication(cc, min.cells = min_cells_comm)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc)

  saveRDS(cc, done_file)

  # PNGs as in v10
  cols <- ct_cols_for(cc)

  png_open <- function(path, w=8, h=8) png(path, width=w, height=h, units="in", res=300, bg="white")

  png_open(file.path(od, "net_circle_weight.png"), 8, 8)
  plot_safe(netVisual_circle(cc@net$weight, color.use = cols,
                             vertex.weight = as.numeric(table(cc@idents)),
                             label.edge = FALSE))
  dev.off()

  png_open(file.path(od, "net_circle_count.png"), 8, 8)
  plot_safe(netVisual_circle(cc@net$count, color.use = cols,
                             vertex.weight = as.numeric(table(cc@idents)),
                             label.edge = FALSE))
  dev.off()

  png_open(file.path(od, "net_heatmap_weight.png"), 9, 7)
  plot_safe(netVisual_heatmap(cc, measure = "weight", color.use = cols))
  dev.off()

  png_open(file.path(od, "net_heatmap_count.png"), 9, 7)
  plot_safe(netVisual_heatmap(cc, measure = "count", color.use = cols))
  dev.off()

  invisible(cc)
}

# -------------------- RUN ALL 4 CONDITIONS --------------------
run_cellchat_all_conditions_v10 <- function(seu, out_root, prefix = "pooled",
                                           raw_use = TRUE, type_use = "truncatedMean") {
  cc_root <- safe_dir(file.path(out_root, "CellChat", prefix))

  conds <- tidyr::crossing(
    Delivery_class = c("Term","Preterm"),
    MNPBurden = c("Low","High")
  ) %>% mutate(cond = paste(Delivery_class, MNPBurden, sep = "_"))

  cc_list <- list()
  for (i in seq_len(nrow(conds))) {
    delivery <- conds$Delivery_class[i]
    burden   <- conds$MNPBurden[i]
    cond     <- conds$cond[i]

    cc_list[[cond]] <- tryCatch(
      run_cellchat_one_v10(
        seu, delivery = delivery, burden = burden,
        outdir = file.path(cc_root, cond),
        min_cells_ct = min_cells_ct,
        min_cells_comm = min_cells_comm,
        cap_cells_per_celltype = cap_cells_per_celltype,
        raw_use = raw_use,
        type_use = type_use
      ),
      error = function(e) {
        od <- safe_dir(file.path(cc_root, cond))
        writeLines(conditionMessage(e), file.path(od, "FAIL_condition.txt"))
        message("[CellChat FAIL] ", cond, " :: ", conditionMessage(e))
        NULL
      }
    )
  }

  saveRDS(cc_list, file.path(cc_root, "cellchat_all_conditions.list.rds"))
  cc_list
}

# -------------------- COMPARE HIGH VS LOW (robust) --------------------
align_cellchat_groups <- function(cc1, cc2) {
  g1 <- levels(cc1@idents)
  g2 <- levels(cc2@idents)
  g_keep <- intersect(g1, g2)
  if (length(g_keep) < 2) stop("Too few overlapping cell groups to compare.")
  cc1@idents <- droplevels(factor(cc1@idents, levels = g_keep))
  cc2@idents <- droplevels(factor(cc2@idents, levels = g_keep))
  if (!is.null(cc1@net$weight)) cc1@net$weight <- cc1@net$weight[g_keep, g_keep, drop = FALSE]
  if (!is.null(cc2@net$weight)) cc2@net$weight <- cc2@net$weight[g_keep, g_keep, drop = FALSE]
  if (!is.null(cc1@net$count))  cc1@net$count  <- cc1@net$count[g_keep, g_keep, drop = FALSE]
  if (!is.null(cc2@net$count))  cc2@net$count  <- cc2@net$count[g_keep, g_keep, drop = FALSE]
  list(cc1 = cc1, cc2 = cc2, groups = g_keep)
}

# helper: export top delta edges
# -------------------- HOTFIX: export_delta_edges + compare_high_low --------------------
export_delta_edges <- function(cc_low, cc_high, od, top_n=50) {

  edge_long <- function(M) {
    df <- as.data.frame(as.table(as.matrix(M)), stringsAsFactors = FALSE)
    # Normalize the value column name across R versions
    # as.table -> Freq
    names(df)[1:2] <- c("source","target")
    if ("Freq" %in% names(df)) names(df)[names(df)=="Freq"] <- "value"
    if (!("value" %in% names(df))) stop("edge_long: missing 'value' (expected 'Freq' or 'value').")
    df
  }

  Wl <- edge_long(cc_low@net$weight)  %>% dplyr::rename(weight_low=value)
  Wh <- edge_long(cc_high@net$weight) %>% dplyr::rename(weight_high=value)
  Cl <- edge_long(cc_low@net$count)   %>% dplyr::rename(count_low=value)
  Ch <- edge_long(cc_high@net$count)  %>% dplyr::rename(count_high=value)

  out <- Wh %>%
    dplyr::left_join(Wl, by=c("source","target")) %>%
    dplyr::mutate(delta_weight = weight_high - weight_low) %>%
    dplyr::left_join(Ch, by=c("source","target")) %>%
    dplyr::left_join(Cl, by=c("source","target")) %>%
    dplyr::mutate(delta_count = count_high - count_low)

  readr::write_csv(out, file.path(od, "edge_delta_full.csv"))
  readr::write_csv(out %>% dplyr::arrange(dplyr::desc(delta_weight)) %>% dplyr::slice_head(n=top_n),
                   file.path(od, paste0("top_", top_n, "_edges_by_delta_weight.csv")))
  readr::write_csv(out %>% dplyr::arrange(dplyr::desc(delta_count)) %>% dplyr::slice_head(n=top_n),
                   file.path(od, paste0("top_", top_n, "_edges_by_delta_count.csv")))
}

align_cellchat_groups <- function(cc1, cc2) {
  g_keep <- intersect(levels(cc1@idents), levels(cc2@idents))
  if (length(g_keep) < 2) stop("Too few overlapping groups.")

  cc1@idents <- droplevels(factor(cc1@idents, levels=g_keep))
  cc2@idents <- droplevels(factor(cc2@idents, levels=g_keep))

  cc1@net$weight <- cc1@net$weight[g_keep, g_keep, drop=FALSE]
  cc2@net$weight <- cc2@net$weight[g_keep, g_keep, drop=FALSE]
  cc1@net$count  <- cc1@net$count[g_keep, g_keep, drop=FALSE]
  cc2@net$count  <- cc2@net$count[g_keep, g_keep, drop=FALSE]

  list(cc1=cc1, cc2=cc2, groups=g_keep)
}

compare_high_low <- function(cc_low, cc_high, outdir, label) {
  od <- safe_dir(outdir)
  if (is.null(cc_low) || is.null(cc_high)) {
    writeLines("Missing CellChat objects (NULL).", file.path(od, "FAIL_compare_missing.txt"))
    return(invisible(NULL))
  }

  ali <- align_cellchat_groups(cc_low, cc_high)
  cc_low2 <- ali$cc1; cc_high2 <- ali$cc2
  writeLines(paste("Common groups:", paste(ali$groups, collapse=", ")), file.path(od, "common_groups.txt"))

  merged <- mergeCellChat(list(Low=cc_low2, High=cc_high2), add.names=c("Low","High"))
  saveRDS(merged, file.path(od, "cellchat_merged.rds"))

  png(file.path(od, "diff_circle_count.png"), width=8, height=8, units="in", res=300, bg="white")
  try(netVisual_diffInteraction(merged, weight.scale=TRUE, measure="count"), silent=TRUE)
  dev.off()

  png(file.path(od, "diff_circle_weight.png"), width=8, height=8, units="in", res=300, bg="white")
  try(netVisual_diffInteraction(merged, weight.scale=TRUE, measure="weight"), silent=TRUE)
  dev.off()

  png(file.path(od, "rankNet_comparison.png"), width=10, height=6, units="in", res=300, bg="white")
  ok <- tryCatch({ rankNet(merged, mode="comparison", stacked=TRUE, do.stat=TRUE); TRUE },
                 error=function(e) FALSE)
  if (!ok) rankNet(merged, mode="comparison", stacked=TRUE, do.stat=FALSE)
  dev.off()

  sink(file.path(od, "compareInteractions_summary.txt"))
  cat("Label:", label, "\n\n")
  cat("compareInteractions(count):\n")
  print(tryCatch(compareInteractions(merged, group=c(1,2), measure="count"), error=function(e) e))
  cat("\ncompareInteractions(weight):\n")
  print(tryCatch(compareInteractions(merged, group=c(1,2), measure="weight"), error=function(e) e))
  sink()

  export_delta_edges(cc_low2, cc_high2, od, top_n=50)

  invisible(merged)
}

# ============================================================
# RUN
# ============================================================
safe_dir(out_root)
obj_lr <- prep_obj_for_cellchat(obj)

# pooled (all platforms)
cc_list_pooled <- run_cellchat_all_conditions_v10(obj_lr, out_root = out_root, prefix = "pooled",
                                                  raw_use = TRUE, type_use = "truncatedMean")

# pooled comparisons
comp_root <- safe_dir(file.path(out_root, "CellChat_compare", "pooled"))

merged_term_pooled <- compare_high_low(
  cc_list_pooled[["Term_Low"]], cc_list_pooled[["Term_High"]],
  outdir=file.path(comp_root, "Term_High_vs_Low"),
  label="pooled_Term_High_vs_Low"
)

# optional per-platform
if (run_per_platform) {
  for (plat in platform_levels) {
    message("[CellChat] PLATFORM split: ", plat)
    seu_p <- subset(obj_lr, subset = platform == plat)
    if (ncol(seu_p) == 0) next

    prefix <- paste0("platform_", plat)
    cc_list_plat <- run_cellchat_all_conditions_v10(seu_p, out_root = out_root, prefix = prefix,
                                                    raw_use = TRUE, type_use = "truncatedMean")

    comp_root_p <- safe_dir(file.path(out_root, "CellChat_compare", prefix))
    compare_high_low_v10(
      cc_list_plat[["Term_Low"]], cc_list_plat[["Term_High"]],
      outdir = file.path(comp_root_p, "Term_High_vs_Low"),
      label = paste0(prefix, "_Term_High_vs_Low")
    )
    compare_high_low_v10(
      cc_list_plat[["Preterm_Low"]], cc_list_plat[["Preterm_High"]],
      outdir = file.path(comp_root_p, "Preterm_High_vs_Low"),
      label = paste0(prefix, "_Preterm_High_vs_Low")
    )
  }
}

message("DONE. Outputs under: ", normalizePath(out_root))

















# ============================================================
# DROP-IN ADD-ON (Publication): CellRep-style Alluvial A/B
# - Enforces your preferred pathway order + group annotations
# - Uses your ct_colors (UMAP-matching)
# - Saves PNG 300 dpi
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(ggalluvial)
  library(CellChat)
})

# ---- YOUR CELLTYPE ORDER + COLORS ----
ord_celltypes <- c(
  "Syncytiotrophoblast","Villous Cytotrophoblast","Extravillous Trophoblast",
  "Fibroblast","Perivascular Cell","Placental Endothelial",
  "Hofbauer Cell","Macrophage","NK Cell","Decidual NK Cell","ILC3",
  "Plasma Cell","B Cell","T Cell","Treg","RCTD_lowconf"
)

ct_colors <- c(
  "Syncytiotrophoblast"     = "#5050FF",
  "Villous Cytotrophoblast" = "#CE3D32",
  "Extravillous Trophoblast"= "#749B58",
  "Fibroblast"              = "#F0E685",
  "Perivascular Cell"       = "#466983",
  "Placental Endothelial"   = "#BA6338",
  "Hofbauer Cell"           = "#5DB1DD",
  "Macrophage"              = "#802268",
  "NK Cell"                 = "#6BD76B",
  "Decidual NK Cell"        = "#D595A7",
  "ILC3"                    = "#924822",
  "Plasma Cell"             = "#837B8D",
  "B Cell"                  = "#C75127",
  "T Cell"                  = "#D58F5C",
  "Treg"                    = "#7A65A5",
  "RCTD_lowconf"            = "#BFBFBF"
)

# ---- YOUR PREFERRED PATHWAY ORDER (ECM -> adhesion -> chemokines -> cytokines) ----
pathways_keep <- c(
  "COLLAGEN","LAMININ","FN1","THBS","TENASCIN","SPP1",
  "ICAM","VCAM",
  "CCL","CXCL",
  "IL1","TNF"
)

# ---- helper ----
png300 <- function(path, w=10, h=4) png(path, width=w, height=h, units="in", res=300, bg="white")

.get_comm_long <- function(cc, pathways = NULL) {
  comm <- tryCatch(subsetCommunication(cc), error = function(e) NULL)
  if (is.null(comm) || !is.data.frame(comm) || nrow(comm) == 0) {
    stop("subsetCommunication(cc) returned nothing. Ensure you ran computeCommunProbPathway + aggregateNet.")
  }
  if (!("pathway_name" %in% names(comm)) && ("pathway" %in% names(comm))) {
    comm$pathway_name <- comm$pathway
  }
  needed <- c("source","target","pathway_name","prob")
  miss <- setdiff(needed, names(comm))
  if (length(miss)) stop("subsetCommunication missing columns: ", paste(miss, collapse=", "))

  comm <- comm %>%
    mutate(
      source = as.character(source),
      target = as.character(target),
      pathway_name = as.character(pathway_name),
      prob = as.numeric(prob)
    ) %>%
    filter(is.finite(prob), prob > 0)

  if (!is.null(pathways)) comm <- comm %>% filter(pathway_name %in% pathways)
  comm
}

build_alluvial_tables <- function(cc, pathways) {
  comm <- .get_comm_long(cc, pathways)

  # Panel A: sender participation per pathway (sum prob over receivers)
  A <- comm %>%
    group_by(sender = source, pathway = pathway_name) %>%
    summarise(strength = sum(prob, na.rm=TRUE), .groups="drop") %>%
    filter(strength > 0)

  # Panel B: receiver participation per pathway (sum prob over senders)
  B <- comm %>%
    group_by(pathway = pathway_name, receiver = target) %>%
    summarise(strength = sum(prob, na.rm=TRUE), .groups="drop") %>%
    filter(strength > 0)

  A$sender  <- factor(A$sender,  levels = intersect(ord_celltypes, unique(A$sender)))
  B$receiver<- factor(B$receiver,levels = intersect(ord_celltypes, unique(B$receiver)))
  A$pathway <- factor(A$pathway, levels = pathways)
  B$pathway <- factor(B$pathway, levels = pathways)

  list(A=A, B=B)
}

# -------------------- FIXED GROUP ANNOTATION (drop-in) --------------------
.pathway_group_labels <- function(pathways) {
  # indices assume pathways are exactly the 12 in pathways_keep, in that order
  tibble(
    group = c("ECM / Matricellular","Adhesion","Chemokines","Cytokines"),
    start = c(1, 7, 9, 11),
    end   = c(6, 8, 10, 12)
  ) %>%
    mutate(mid = (start + end) / 2)
}

plot_alluvial_A <- function(dfA, pathways, title = NULL) {
  ann <- .pathway_group_labels(pathways)

  ggplot(dfA, aes(axis1 = sender, axis2 = pathway, y = strength)) +
    geom_alluvium(aes(fill = sender), width = 0.15, alpha = 0.85) +
    geom_stratum(width = 0.15, color = "grey30", fill = "grey90") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_fill_manual(values = ct_colors, drop = FALSE) +
    scale_x_discrete(limits = c("Sender cell", "Signaling pathway"), expand = c(.05, .05)) +
    labs(title = title, x = NULL, y = "Interaction strength (sum prob)") +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    # group labels above the pathway axis (x = 2)
    geom_text(
      data = ann,
      aes(x = 2, y = Inf, label = group),
      inherit.aes = FALSE,
      vjust = 1.4,
      size = 3.2
    )
}

make_cellrep_alluvial_AB <- function(
  cc,
  outdir,
  condition_label = "PTL",
  pathways = pathways_keep,
  width_in = 10,
  height_in = 4
) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  tabs <- build_alluvial_tables(cc, pathways)

  write_csv(tabs$A, file.path(outdir, paste0("alluvial_A_sender_by_pathway_", condition_label, ".csv")))
  write_csv(tabs$B, file.path(outdir, paste0("alluvial_B_pathway_by_receiver_", condition_label, ".csv")))
  writeLines(pathways, file.path(outdir, paste0("pathways_used_", condition_label, ".txt")))

  pA <- plot_alluvial_A(tabs$A, pathways, title = paste0("A  ", condition_label, ": Sender roles"))
  png300(file.path(outdir, paste0("Fig_A_alluvial_sender_", condition_label, ".png")), w=width_in, h=height_in)
  print(pA)
  dev.off()

  pB <- plot_alluvial_B(tabs$B, pathways, title = paste0("B  ", condition_label, ": Receiver roles"))
  png300(file.path(outdir, paste0("Fig_B_alluvial_receiver_", condition_label, ".png")), w=width_in, h=height_in)
  print(pB)
  dev.off()

  invisible(list(pathways = pathways, A = tabs$A, B = tabs$B))
}

# ============================================================
# RUN (example)
# ============================================================
cc_control <- readRDS("/path/to/local_home/visiumhd_mnp/results/audra/MNP_interactomes_CellChat_v12_full_min500_png/CellChat/pooled/Preterm_Low/cellchat.rds")
cc_ptl     <- readRDS("/path/to/local_home/visiumhd_mnp/results/audra/MNP_interactomes_CellChat_v12_full_min500_png/CellChat/pooled/Preterm_High/cellchat.rds")

unique(subsetCommunication(cc_ptl)$pathway_name) |> sort()


out_ab <- file.path(out_root, "Figure_addons", "alluvial_AB_core12_ordered")
make_cellrep_alluvial_AB(cc_control, outdir = out_ab, condition_label="Control",
                         pathways = pathways_keep, width_in=8, height_in=8)
make_cellrep_alluvial_AB(cc_ptl, outdir = out_ab, condition_label="PTL",
                         pathways = pathways_keep, width_in=8, height_in=8)



