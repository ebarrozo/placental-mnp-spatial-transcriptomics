## nichenet_v1.R
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
setwd("/home/ebarrozo/visiumhd_mnp/results/audra")
obj <-readRDS("/media/jochum00/Aagaard_Raid3/ebarrozo/seurat_all_ref_merged/with_annotations/SCT_chunks_mapped/audra/mnp_obj_w_nn_predicted_celltype.rds") # 7.1 gb object with SCT data [no counts], all 3,969,532 cells, and umap


mypal3 <-get_palette("igv",40)

## 4) Sanity + quick plot (palette later)
Idents(obj) <- "nn.predicted.celltype"
levels(obj)

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


# ============================================================
# DROP-IN: NicheNet (sender -> receiver) for Control vs PTL
# - Works directly from Seurat object (SCT assay recommended)
# - Produces:
#   (1) per receiver: top ligands (PTL vs Control)
#   (2) ligand->target heatmap table
#   (3) sender ligand expression table
#   (4) optional differential expression summaries
# Notes:
# - This is the "classic" NicheNet workflow (nichenetr package).
# - You MUST install/load: nichenetr, Seurat, dplyr, tibble, readr, Matrix
# - Assumes metadata columns: Delivery_class, MNPBurden, nn.predicted.celltype
# - You define:
#   senders, receivers, condition columns, and which two groups to compare.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr); library(tibble); library(readr)
  library(Matrix)
  library(nichenetr)
})

# -------------------- USER CONFIG --------------------
out_root <- "MNP_interactomes_NicheNet_v1"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

assay_use <- "SCT"         # use SCT for integrated objects
slot_use  <- "data"        # NicheNet expects log-normalized expression

# Define groups
group_col_delivery <- "Delivery_class"   # Term/Preterm
group_col_burden   <- "MNPBurden"        # High/Low
celltype_col       <- "nn.predicted.celltype"

# Choose comparison (example: PTL=Preterm_High vs Control=Preterm_Low)
delivery_use <- "Preterm"
condition_A_name <- "Control"
condition_B_name <- "PTL"
condition_A_filter <- list(Delivery_class = "Preterm", MNPBurden = "Low")
condition_B_filter <- list(Delivery_class = "Preterm", MNPBurden = "High")

# Sender/receiver sets (edit to match your filtered min_cells_ct=500 regime)
senders   <- c("Macrophage","NK Cell","T Cell","Extravillous Trophoblast","Villous Cytotrophoblast","Syncytiotrophoblast")
receivers <- c("Extravillous Trophoblast","Villous Cytotrophoblast","Syncytiotrophoblast","Macrophage","NK Cell")

min_cells_per_group <- 200   # NicheNet is sensitive; raise if noisy
expr_pct_sender <- 0.05      # ligand expressed in >=5% of sender cells
expr_pct_receiver <- 0.05    # receptor expressed in >=5% of receiver cells

# DE settings (keep conservative)
logfc_threshold <- 0.10
p_adj_threshold <- 0.05
top_ligands <- 30
top_targets <- 200

# -------------------- HELPERS --------------------
safe_dir <- function(path) { dir.create(path, showWarnings = FALSE, recursive = TRUE); path }

get_mat <- function(seu, assay = assay_use, slot = slot_use) {
  m <- tryCatch(GetAssayData(seu, assay = assay, layer = slot),
                error = function(e) GetAssayData(seu, assay = assay, slot = slot))
  if (!inherits(m, "dgCMatrix")) m <- as(m, "dgCMatrix")
  m
}

subset_by_filters <- function(seu, filters) {
  # filters: named list of meta filters
  expr <- paste0(names(filters), " == '", unlist(filters), "'", collapse = " & ")
  subset(seu, subset = eval(parse(text = expr)))
}

min_cells_filter <- function(seu, celltype_col, keep_types, min_cells) {
  md <- seu[[]]
  md[[celltype_col]] <- as.character(md[[celltype_col]])
  md <- md[!is.na(md[[celltype_col]]), , drop = FALSE]
  md <- md[md[[celltype_col]] %in% keep_types, , drop = FALSE]
  tab <- table(md[[celltype_col]])
  keep <- names(tab)[tab >= min_cells]
  subset(seu, cells = rownames(md)[md[[celltype_col]] %in% keep])
}

pct_expr_gene <- function(expr_mat, cells_idx) {
  Matrix::rowMeans(expr_mat[, cells_idx, drop = FALSE] > 0)
}

# -------------------- LOAD NICHENET PRIORS --------------------
# built-in datasets
# -------------------- LOAD NICHENET PRIORS (version-robust) --------------------
suppressPackageStartupMessages(library(nichenetr))
lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
library(jsonlite)

rec_id <- "7074291"
rec <- fromJSON(paste0("https://zenodo.org/api/records/", rec_id))

read_zenodo_rds <- function(key, rec){
  stopifnot(key %in% rec$files$key)
  url <- rec$files$links$self[match(key, rec$files$key)]
  tmp <- tempfile(fileext = ".rds")
  download.file(url, tmp, mode = "wb", quiet = TRUE)
  readRDS(tmp)
}

# Human (the ones your stopifnot wants)
lr_network           <- read_zenodo_rds("lr_network_human_21122021.rds", rec)
weighted_networks    <- read_zenodo_rds("weighted_networks_nsga2r_final.rds", rec)
ligand_target_matrix <- read_zenodo_rds("ligand_target_matrix_nsga2r_final.rds", rec)
ligand_tf_matrix     <- read_zenodo_rds("ligand_tf_matrix_nsga2r_final.rds", rec)

# sanity
stopifnot(exists("lr_network"), exists("weighted_networks"),
          exists("ligand_target_matrix"), exists("ligand_tf_matrix"))

# quick structure checks
str(lr_network)
str(weighted_networks, max.level = 1)
dim(ligand_target_matrix)
dim(ligand_tf_matrix)
signaling_network <- read_zenodo_rds("signaling_network_human_21122021.rds", rec)
gr_network        <- read_zenodo_rds("gr_network_human_21122021.rds", rec)

str(signaling_network, max.level = 1)
str(gr_network, max.level = 1)

# sanity (optional)
stopifnot(exists("lr_network"), exists("weighted_networks"),
          exists("ligand_target_matrix"), exists("ligand_tf_matrix"))

# -------------------- PREPARE SEURAT SUBSETS --------------------
# Subset to the two conditions
obj_A <- subset_by_filters(obj_lr, condition_A_filter)
obj_B <- subset_by_filters(obj_lr, condition_B_filter)

# Keep only desired cell types and ensure each has enough cells (per condition)
obj_A <- min_cells_filter(obj_A, celltype_col, unique(c(senders, receivers)), min_cells_per_group)
obj_B <- min_cells_filter(obj_B, celltype_col, unique(c(senders, receivers)), min_cells_per_group)

# -------------------- EXPRESSION MATRICES --------------------
expr_A <- get_mat(obj_A)
expr_B <- get_mat(obj_B)

meta_A <- obj_A[[]] %>% mutate(cell = rownames(.))
meta_B <- obj_B[[]] %>% mutate(cell = rownames(.))

# Use gene symbols that NicheNet expects
genes <- intersect(rownames(expr_A), rownames(expr_B))
expr_A <- expr_A[genes, , drop = FALSE]
expr_B <- expr_B[genes, , drop = FALSE]

# -------------------- MAIN LOOP: per receiver --------------------
out_main <- safe_dir(file.path(out_root, paste0(delivery_use, "_", condition_B_name, "_vs_", condition_A_name)))

receiver_results <- list()

for (receiver in receivers) {
  message("[NicheNet] Receiver: ", receiver)

  # cells for this receiver in each condition
  recv_cells_A <- meta_A$cell[meta_A[[celltype_col]] == receiver]
  recv_cells_B <- meta_B$cell[meta_B[[celltype_col]] == receiver]

  if (length(recv_cells_A) < min_cells_per_group || length(recv_cells_B) < min_cells_per_group) {
    message("  skip (too few receiver cells)")
    next
  }

  # define expressed genes in receiver (for receptor selection and background)
  recv_expr_pct_A <- pct_expr_gene(expr_A, which(colnames(expr_A) %in% recv_cells_A))
  recv_expr_pct_B <- pct_expr_gene(expr_B, which(colnames(expr_B) %in% recv_cells_B))
  expressed_genes_receiver <- names(which(pmax(recv_expr_pct_A, recv_expr_pct_B) >= expr_pct_receiver))

  # DE in receiver: PTL vs Control (targets)
  # Using Seurat DE for robustness; uses SCT assay by default if set
  Idents(obj_A) <- celltype_col
  Idents(obj_B) <- celltype_col

  # Build a combined object for DE convenience
  comb <- merge(obj_A, y = obj_B, add.cell.ids = c(condition_A_name, condition_B_name))
  comb$Condition <- ifelse(grepl(paste0("^", condition_A_name, "_"), colnames(comb)),
                           condition_A_name, condition_B_name)
  Idents(comb) <- "Condition"

  # restrict to receiver cells only
  comb_recv <- subset(comb, cells = colnames(comb)[comb[[celltype_col]] == receiver])

  # FindMarkers: PTL vs Control in receiver
  de <- tryCatch(
    FindMarkers(comb_recv, ident.1 = condition_B_name, ident.2 = condition_A_name,
                assay = assay_use, slot = slot_use, logfc.threshold = 0, min.pct = 0.05),
    error = function(e) NULL
  )

  if (is.null(de) || nrow(de) == 0) {
    message("  skip (no DE results)")
    next
  }

  de <- de %>%
    rownames_to_column("gene") %>%
    mutate(gene = toupper(gene))

  # NicheNet expects these columns (Seurat uses avg_log2FC or avg_logFC depending version)
  fc_col <- dplyr::case_when(
    "avg_log2FC" %in% names(de) ~ "avg_log2FC",
    "avg_logFC"  %in% names(de) ~ "avg_logFC",
    TRUE ~ NA_character_
  )
  if (is.na(fc_col)) stop("Could not find avg_log2FC/avg_logFC in DE table for receiver ", receiver)

  de_sig <- de %>%
    filter(.data[[fc_col]] >= logfc_threshold, p_val_adj <= p_adj_threshold) %>%
    arrange(desc(.data[[fc_col]]))

  if (nrow(de_sig) == 0) {
    message("  skip (no significant up targets in receiver)")
    next
  }

  geneset_oi <- intersect(de_sig$gene, rownames(ligand_target_matrix))
  geneset_oi <- head(geneset_oi, top_targets)

  # ---- define potential ligands from senders ----
  sender_cells_A <- meta_A %>% filter(.data[[celltype_col]] %in% senders) %>% pull(cell)
  sender_cells_B <- meta_B %>% filter(.data[[celltype_col]] %in% senders) %>% pull(cell)

  sender_expr_pct_A <- pct_expr_gene(expr_A, which(colnames(expr_A) %in% sender_cells_A))
  sender_expr_pct_B <- pct_expr_gene(expr_B, which(colnames(expr_B) %in% sender_cells_B))
  expressed_genes_sender <- names(which(pmax(sender_expr_pct_A, sender_expr_pct_B) >= expr_pct_sender))

  # Receiver expressed receptors
  receptors <- lr_network %>% pull(to) %>% unique()
  expressed_receptors <- intersect(receptors, expressed_genes_receiver)

  # Ligands expressed by senders
  ligands <- lr_network %>% pull(from) %>% unique()
  expressed_ligands <- intersect(ligands, expressed_genes_sender)

  # Potential ligands are those with at least one expressed receptor in receiver
  lr_pairs_filtered <- lr_network %>%
    filter(from %in% expressed_ligands, to %in% expressed_receptors)

  potential_ligands <- lr_pairs_filtered %>% pull(from) %>% unique()

  if (length(potential_ligands) == 0) {
    message("  skip (no potential ligands after receptor/ligand expression filtering)")
    next
  }

  # ---- NicheNet ligand activity ----
  ligand_activities <- nichenetr::predict_ligand_activities(
    geneset = geneset_oi,
    background_expressed_genes = intersect(expressed_genes_receiver, rownames(ligand_target_matrix)),
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = potential_ligands
  ) %>% arrange(desc(aupr_corrected))

  top_ligs <- head(ligand_activities$test_ligand, top_ligands)
get_links_safe <- function(top_ligs, geneset_oi, ligand_target_matrix, n = 250) {
  f <- nichenetr::get_weighted_ligand_target_links
  fn <- names(formals(f))

  # candidate argument names used across nichenetr versions
  lig_arg <- intersect(
    c("ligands", "ligands_oi", "ligand_oi", "ligand", "test_ligands", "ligands_of_interest"),
    fn
  )
  if (length(lig_arg) == 0) {
    stop("Could not find a ligand argument in get_weighted_ligand_target_links(). Formals: ",
         paste(fn, collapse = ", "))
  }
  lig_arg <- lig_arg[1]

  # build call using only arguments that exist in this version
  call_list <- list()
  if ("ligand_target_matrix" %in% fn) call_list$ligand_target_matrix <- ligand_target_matrix
  if ("geneset" %in% fn)             call_list$geneset <- geneset_oi
  call_list[[lig_arg]] <- top_ligs
  if ("n" %in% fn)                   call_list$n <- n

  do.call(f, call_list)
}
  # ligand -> target links for top ligands
# top ligands
top_ligs <- head(ligand_activities$test_ligand, top_ligands)

# ligand -> target links (your nichenetr version expects `ligand` and works per-ligand)
ligand_target_links_df <- dplyr::bind_rows(lapply(top_ligs, function(lg) {
  nichenetr::get_weighted_ligand_target_links(
    ligand = lg,
    geneset = geneset_oi,
    ligand_target_matrix = ligand_target_matrix,
    n = 250
  )
})) %>%
  dplyr::filter(!is.na(target), !is.na(weight))

  # sender expression summary for top ligands (PTL vs Control)
  # compute avg expr in each sender celltype for each condition
  sender_expr_summary <- function(expr_mat, meta_df, ligs) {
    out <- lapply(unique(meta_df[[celltype_col]][meta_df[[celltype_col]] %in% senders]), function(ct) {
      cells <- meta_df$cell[meta_df[[celltype_col]] == ct]
      if (!length(cells)) return(NULL)
      idx <- which(colnames(expr_mat) %in% cells)
      mu <- Matrix::rowMeans(expr_mat[intersect(ligs, rownames(expr_mat)), idx, drop=FALSE])
      tibble(celltype = ct, ligand = names(mu), mean_expr = as.numeric(mu))
    })
    bind_rows(out)
  }
  sender_A <- sender_expr_summary(expr_A, meta_A, top_ligs) %>% mutate(Condition = condition_A_name)
  sender_B <- sender_expr_summary(expr_B, meta_B, top_ligs) %>% mutate(Condition = condition_B_name)
  sender_expr_tbl <- bind_rows(sender_A, sender_B) %>%
    tidyr::pivot_wider(names_from = Condition, values_from = mean_expr, values_fill = 0) %>%
    mutate(delta = .data[[condition_B_name]] - .data[[condition_A_name]]) %>%
    arrange(desc(delta))

  # write outputs per receiver
  od <- safe_dir(file.path(out_main, paste0("receiver_", gsub("[^A-Za-z0-9]+","_", receiver))))
  write_csv(ligand_activities, file.path(od, "ligand_activities_full.csv"))
  write_csv(ligand_activities %>% slice_head(n = top_ligands), file.path(od, paste0("top_", top_ligands, "_ligands.csv")))
  write_csv(ligand_target_links_df, file.path(od, "ligand_target_links_topLigands.csv"))
  write_csv(sender_expr_tbl, file.path(od, "sender_ligand_expression_delta.csv"))
  write_csv(de_sig, file.path(od, "receiver_DE_up_targets.csv"))

  receiver_results[[receiver]] <- list(
    ligand_activities = ligand_activities,
    top_ligands = top_ligs,
    ligand_target_links = ligand_target_links_df,
    sender_expr = sender_expr_tbl,
    receiver_de_up = de_sig
  )
}

saveRDS(receiver_results, file.path(out_main, "nichenet_receiver_results.rds"))

message("DONE. NicheNet outputs under: ", normalizePath(out_main))
































suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(readr)
  library(Matrix)
  library(nichenetr)
  library(jsonlite)
  library(tidyr)
})

# -------------------- USER CONFIG --------------------
out_root <- "MNP_interactomes_NicheNet_v2"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

assay_use <- "SCT"
slot_use  <- "data"  # log-normalized

group_col_delivery <- "Delivery_class"
group_col_burden   <- "MNPBurden"
celltype_col       <- "nn.predicted.celltype"

delivery_use <- "Preterm"
condition_A_name <- "Control"
condition_B_name <- "PTL"
condition_A_filter <- list(Delivery_class = "Preterm", MNPBurden = "Low")
condition_B_filter <- list(Delivery_class = "Preterm", MNPBurden = "High")

senders   <- c("Macrophage","NK Cell","T Cell","Extravillous Trophoblast","Villous Cytotrophoblast","Syncytiotrophoblast")
receivers <- c("Extravillous Trophoblast","Villous Cytotrophoblast","Syncytiotrophoblast","Macrophage","NK Cell")

min_cells_per_group <- 200
expr_pct_sender <- 0.05
expr_pct_receiver <- 0.05

logfc_threshold <- 0.10
p_adj_threshold <- 0.05
top_ligands <- 30
top_targets <- 200

# -------------------- HELPERS --------------------
safe_dir <- function(path) { dir.create(path, showWarnings = FALSE, recursive = TRUE); path }

get_mat <- function(seu, assay = assay_use, slot = slot_use) {
  m <- tryCatch(GetAssayData(seu, assay = assay, layer = slot),
                error = function(e) GetAssayData(seu, assay = assay, slot = slot))
  if (!inherits(m, "dgCMatrix")) m <- as(m, "dgCMatrix")
  m
}

subset_by_filters_safe <- function(seu, filters) {
  md <- seu[[]]
  missing_cols <- setdiff(names(filters), colnames(md))
  if (length(missing_cols)) stop("Missing metadata cols: ", paste(missing_cols, collapse=", "))
  keep <- rep(TRUE, nrow(md))
  for (nm in names(filters)) {
    keep <- keep & as.character(md[[nm]]) %in% as.character(filters[[nm]])
  }
  subset(seu, cells = rownames(md)[keep])
}

min_cells_filter <- function(seu, celltype_col, keep_types, min_cells) {
  md <- seu[[]]
  md[[celltype_col]] <- as.character(md[[celltype_col]])
  md <- md[!is.na(md[[celltype_col]]), , drop = FALSE]
  md <- md[md[[celltype_col]] %in% keep_types, , drop = FALSE]
  tab <- table(md[[celltype_col]])
  keep <- names(tab)[tab >= min_cells]
  subset(seu, cells = rownames(md)[md[[celltype_col]] %in% keep])
}

pct_expr_gene <- function(expr_mat, cells) {
  idx <- which(colnames(expr_mat) %in% cells)
  if (!length(idx)) return(setNames(rep(0, nrow(expr_mat)), rownames(expr_mat)))
  Matrix::rowMeans(expr_mat[, idx, drop = FALSE] > 0)
}

get_fc_col <- function(de_tbl) {
  if ("avg_log2FC" %in% names(de_tbl)) return("avg_log2FC")
  if ("avg_logFC"  %in% names(de_tbl)) return("avg_logFC")
  stop("Could not find avg_log2FC/avg_logFC in DE table.")
}

# -------------------- LOAD NICHENET PRIORS (Zenodo 7074291) --------------------
rec_id <- "7074291"
rec <- fromJSON(paste0("https://zenodo.org/api/records/", rec_id))

read_zenodo_rds <- function(key, rec){
  stopifnot(key %in% rec$files$key)
  url <- rec$files$links$self[match(key, rec$files$key)]
  tmp <- tempfile(fileext = ".rds")
  download.file(url, tmp, mode = "wb", quiet = TRUE)
  readRDS(tmp)
}

lr_network           <- read_zenodo_rds("lr_network_human_21122021.rds", rec)
ligand_target_matrix <- read_zenodo_rds("ligand_target_matrix_nsga2r_final.rds", rec)

stopifnot(exists("lr_network"), exists("ligand_target_matrix"))

ligands_all   <- unique(lr_network$from)
receptors_all <- unique(lr_network$to)

# -------------------- PREPARE SEURAT SUBSETS --------------------
# obj_lr must exist in your environment; if your object is named obj, set obj_lr <- obj
stopifnot(exists("obj_lr"))

obj_A <- subset_by_filters_safe(obj_lr, condition_A_filter)
obj_B <- subset_by_filters_safe(obj_lr, condition_B_filter)

obj_A <- min_cells_filter(obj_A, celltype_col, unique(c(senders, receivers)), min_cells_per_group)
obj_B <- min_cells_filter(obj_B, celltype_col, unique(c(senders, receivers)), min_cells_per_group)

expr_A <- get_mat(obj_A)
expr_B <- get_mat(obj_B)

meta_A <- obj_A[[]] %>% mutate(cell = rownames(.), celltype = as.character(.data[[celltype_col]]))
meta_B <- obj_B[[]] %>% mutate(cell = rownames(.), celltype = as.character(.data[[celltype_col]]))

# Align gene universe
genes <- intersect(rownames(expr_A), rownames(expr_B))
expr_A <- expr_A[genes, , drop = FALSE]
expr_B <- expr_B[genes, , drop = FALSE]

# -------------------- BUILD COMBINED OBJECT ONCE (for DE) --------------------
comb <- merge(obj_A, y = obj_B, add.cell.ids = c(condition_A_name, condition_B_name))
comb$Condition <- ifelse(grepl(paste0("^", condition_A_name, "_"), colnames(comb)),
                         condition_A_name, condition_B_name)

# Force celltype and condition to plain character
comb[[celltype_col]] <- as.character(comb[[celltype_col]][,1])
comb$Condition <- as.character(comb$Condition)

Idents(comb) <- "Condition"

out_main <- safe_dir(file.path(out_root, paste0(delivery_use, "_", condition_B_name, "_vs_", condition_A_name)))
receiver_results <- list()

# Precompute sender expressed genes (union across conditions) once
sender_cells_A <- meta_A %>% filter(celltype %in% senders) %>% pull(cell)
sender_cells_B <- meta_B %>% filter(celltype %in% senders) %>% pull(cell)
sender_expr_pct_A <- pct_expr_gene(expr_A, sender_cells_A)
sender_expr_pct_B <- pct_expr_gene(expr_B, sender_cells_B)
expressed_genes_sender <- names(which(pmax(sender_expr_pct_A, sender_expr_pct_B) >= expr_pct_sender))
expressed_ligands <- intersect(ligands_all, expressed_genes_sender)

for (receiver in receivers) {
  message("[NicheNet v2] Receiver: ", receiver)

  receiver_safe <- gsub("[^A-Za-z0-9]+","_", receiver)
  od <- safe_dir(file.path(out_main, paste0("receiver_", receiver_safe)))

  tryCatch({

    recv_cells_A <- meta_A$cell[meta_A$celltype == receiver]
    recv_cells_B <- meta_B$cell[meta_B$celltype == receiver]

    if (length(recv_cells_A) < min_cells_per_group || length(recv_cells_B) < min_cells_per_group) {
      message("  skip (too few receiver cells)")
      next
    }

    # expressed genes in receiver (union across A/B)
    recv_expr_pct_A <- pct_expr_gene(expr_A, recv_cells_A)
    recv_expr_pct_B <- pct_expr_gene(expr_B, recv_cells_B)
    expressed_genes_receiver <- names(which(pmax(recv_expr_pct_A, recv_expr_pct_B) >= expr_pct_receiver))

    expressed_receptors <- intersect(receptors_all, expressed_genes_receiver)

    # ligand-receptor filter
    lr_pairs_filtered <- lr_network %>%
      dplyr::filter(from %in% expressed_ligands, to %in% expressed_receptors)

    potential_ligands <- unique(lr_pairs_filtered$from)
    if (!length(potential_ligands)) {
      message("  skip (no potential ligands after expression filtering)")
      next
    }

    # ---- DE in receiver on combined object (only receiver cells) ----
    comb_recv_cells <- colnames(comb)[comb[[celltype_col]][,1] == receiver]
    comb_recv <- subset(comb, cells = comb_recv_cells)

    de <- FindMarkers(
      comb_recv,
      ident.1 = condition_B_name, ident.2 = condition_A_name,
      assay = assay_use, slot = slot_use,
      logfc.threshold = 0, min.pct = 0.05
    )

    if (is.null(de) || nrow(de) == 0) {
      message("  skip (no DE results)")
      next
    }

    de <- de %>%
      tibble::rownames_to_column("gene") %>%
      dplyr::mutate(gene = toupper(gene))

    fc_col <- get_fc_col(de)

    de_sig <- de %>%
      dplyr::filter(.data[[fc_col]] >= logfc_threshold, p_val_adj <= p_adj_threshold) %>%
      dplyr::arrange(dplyr::desc(.data[[fc_col]]))

    if (nrow(de_sig) == 0) {
      message("  skip (no significant up targets in receiver)")
      next
    }

    geneset_oi <- intersect(de_sig$gene, rownames(ligand_target_matrix))
    geneset_oi <- head(geneset_oi, top_targets)

    # ---- ligand activity ----
    ligand_activities <- nichenetr::predict_ligand_activities(
      geneset = geneset_oi,
      background_expressed_genes = intersect(expressed_genes_receiver, rownames(ligand_target_matrix)),
      ligand_target_matrix = ligand_target_matrix,
      potential_ligands = potential_ligands
    ) %>% dplyr::arrange(dplyr::desc(aupr_corrected))

    top_ligs <- head(ligand_activities$test_ligand, top_ligands)

    # ---- ligand -> target links (your version: ligand=, run per ligand) ----
    ligand_target_links_df <- dplyr::bind_rows(lapply(top_ligs, function(lg) {
      nichenetr::get_weighted_ligand_target_links(
        ligand = lg,
        geneset = geneset_oi,
        ligand_target_matrix = ligand_target_matrix,
        n = 250
      )
    })) %>%
      dplyr::filter(!is.na(target), !is.na(weight))

    # ---- sender expression deltas for top ligands ----
    sender_expr_summary <- function(expr_mat, meta_df, ligs) {
      out <- lapply(unique(meta_df$celltype[meta_df$celltype %in% senders]), function(ct) {
        cells <- meta_df$cell[meta_df$celltype == ct]
        if (!length(cells)) return(NULL)
        idx <- which(colnames(expr_mat) %in% cells)
        g <- intersect(ligs, rownames(expr_mat))
        if (!length(g) || !length(idx)) return(NULL)
        mu <- Matrix::rowMeans(expr_mat[g, idx, drop = FALSE])
        tibble::tibble(celltype = ct, ligand = g, mean_expr = as.numeric(mu))
      })
      dplyr::bind_rows(out)
    }

    sender_A <- sender_expr_summary(expr_A, meta_A, top_ligs) %>% dplyr::mutate(Condition = condition_A_name)
    sender_B <- sender_expr_summary(expr_B, meta_B, top_ligs) %>% dplyr::mutate(Condition = condition_B_name)

    sender_expr_tbl <- dplyr::bind_rows(sender_A, sender_B) %>%
      tidyr::pivot_wider(names_from = Condition, values_from = mean_expr, values_fill = 0) %>%
      dplyr::mutate(delta = .data[[condition_B_name]] - .data[[condition_A_name]]) %>%
      dplyr::arrange(dplyr::desc(delta))

    # ---- write outputs ----
    readr::write_csv(ligand_activities, file.path(od, "ligand_activities_full.csv"))
    readr::write_csv(dplyr::slice_head(ligand_activities, n = top_ligands),
                     file.path(od, paste0("top_", top_ligands, "_ligands.csv")))
    readr::write_csv(ligand_target_links_df, file.path(od, "ligand_target_links_topLigands.csv"))
    readr::write_csv(sender_expr_tbl, file.path(od, "sender_ligand_expression_delta.csv"))
    readr::write_csv(de_sig, file.path(od, "receiver_DE_up_targets.csv"))

    receiver_results[[receiver]] <- list(
      ligand_activities = ligand_activities,
      top_ligands = top_ligs,
      ligand_target_links = ligand_target_links_df,
      sender_expr = sender_expr_tbl,
      receiver_de_up = de_sig
    )

  }, error = function(e) {
    writeLines(conditionMessage(e), file.path(od, "ERROR.txt"))
    message("  ERROR (saved): ", conditionMessage(e))
    NULL
  })
}

saveRDS(receiver_results, file.path(out_main, "nichenet_receiver_results.rds"))
message("DONE. NicheNet outputs under: ", normalizePath(out_main))















suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(Matrix)
})

# -------------------- CONFIG --------------------
# set this to your actual run folder
out_main <- "/home/ebarrozo/visiumhd_mnp/results/audra/MNP_interactomes_NicheNet_v2/Preterm_PTL_vs_Control"

fig_dir <- file.path(out_main, "FIGS")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# plotting knobs
topN_bar <- 20          # top ligands shown in bar plot
top_targets_heat <- 40  # how many targets to show in ligand->target heatmap
min_link_weight <- 0    # filter ligand-target links by weight (try 0.02 if too dense)

# -------------------- HELPERS --------------------
receiver_dirs <- list.dirs(out_main, recursive = FALSE, full.names = TRUE)
receiver_dirs <- receiver_dirs[grepl("receiver_", basename(receiver_dirs))]

receiver_name <- function(path) {
  x <- sub("^receiver_", "", basename(path))
  gsub("_", " ", x)
}

read_if_exists <- function(path) if (file.exists(path)) suppressMessages(read_csv(path)) else NULL

save_png <- function(p, path, w=7, h=5, dpi=300) {
  ggsave(path, p, width = w, height = h, dpi = dpi, bg = "white")
}

# safe percent rescale
rescale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) return(rep(0, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

# -------------------- 1) PER-RECEIVER FIGURES --------------------
per_receiver_summary <- list()

for (rd in receiver_dirs) {
  recv <- receiver_name(rd)
  recv_safe <- basename(rd)

  lig_full <- read_if_exists(file.path(rd, "ligand_activities_full.csv"))
  links    <- read_if_exists(file.path(rd, "ligand_target_links_topLigands.csv"))
  sender   <- read_if_exists(file.path(rd, "sender_ligand_expression_delta.csv"))
  de_up    <- read_if_exists(file.path(rd, "receiver_DE_up_targets.csv"))

  if (is.null(lig_full)) next

  lig_full <- lig_full %>%
    mutate(test_ligand = as.character(test_ligand)) %>%
    arrange(desc(aupr_corrected))

  # ---------- (A) Top ligands bar ----------
  lig_top <- lig_full %>% slice_head(n = topN_bar) %>%
    mutate(test_ligand = factor(test_ligand, levels = rev(test_ligand)))

  p_bar <- ggplot(lig_top, aes(x = test_ligand, y = aupr_corrected)) +
    geom_col() +
    coord_flip() +
    labs(title = paste0(recv, ": Top ligands (AUPR_corrected)"),
         x = NULL, y = "AUPR_corrected") +
    theme_classic(base_size = 12)

  save_png(p_bar, file.path(fig_dir, paste0(recv_safe, "_topLigands_bar.png")), w=7, h=5)

  # ---------- (B) Ligand -> target heatmap (tile) ----------
  # links has columns ligand/target/weight (from your fixed per-ligand bind_rows)
  if (!is.null(links) && nrow(links) > 0) {
    links2 <- links %>%
      mutate(ligand = as.character(ligand),
             target = as.character(target),
             weight = as.numeric(weight)) %>%
      filter(!is.na(target), !is.na(weight), weight >= min_link_weight)

    # choose targets to show: top by max weight among top ligands
    top_ligs <- lig_full %>% slice_head(n = topN_bar) %>% pull(test_ligand)
    links2 <- links2 %>% filter(ligand %in% top_ligs)

    if (nrow(links2) > 0) {
      top_tg <- links2 %>%
        group_by(target) %>%
        summarise(mx = max(weight, na.rm = TRUE), .groups="drop") %>%
        arrange(desc(mx)) %>%
        slice_head(n = top_targets_heat) %>%
        pull(target)

      links2 <- links2 %>% filter(target %in% top_tg)

      # order axes
      lig_order <- lig_full %>% filter(test_ligand %in% unique(links2$ligand)) %>%
        arrange(desc(aupr_corrected)) %>% pull(test_ligand)
      tg_order <- links2 %>%
        group_by(target) %>% summarise(mx=max(weight, na.rm=TRUE), .groups="drop") %>%
        arrange(mx) %>% pull(target)

      links2$ligand <- factor(links2$ligand, levels = lig_order)
      links2$target <- factor(links2$target, levels = tg_order)

      p_lt <- ggplot(links2, aes(x = ligand, y = target, fill = weight)) +
        geom_tile() +
        coord_cartesian(expand = FALSE) +
        labs(title = paste0(recv, ": Ligand → target weights"),
             x = "Ligand", y = "Target", fill = "Weight") +
        theme_classic(base_size = 10) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.text.y = element_text(size = 7))

      save_png(p_lt, file.path(fig_dir, paste0(recv_safe, "_ligandTarget_heat.png")), w=10, h=8)
    }
  }

  # ---------- (C) Sender delta heatmap (senders x ligands) ----------
  # sender table columns: celltype, ligand, Control, PTL, delta
  if (!is.null(sender) && nrow(sender) > 0) {
    sender2 <- sender %>%
      mutate(celltype = as.character(celltype),
             ligand = as.character(ligand),
             delta = as.numeric(delta)) %>%
      filter(celltype %in% unique(celltype))

    # focus on top ligands from activity
    top_ligs2 <- lig_full %>% slice_head(n = topN_bar) %>% pull(test_ligand)
    sender2 <- sender2 %>% filter(ligand %in% top_ligs2)

    # order
    sender2$ligand <- factor(sender2$ligand, levels = rev(top_ligs2))
    sender2$celltype <- factor(sender2$celltype, levels = sort(unique(sender2$celltype)))

    p_sender <- ggplot(sender2, aes(x = ligand, y = celltype, fill = delta)) +
      geom_tile() +
      coord_cartesian(expand = FALSE) +
      labs(title = paste0(recv, ": Sender expression Δ (PTL – Control)"),
           x = "Ligand", y = NULL, fill = "Δ expr") +
      theme_classic(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    save_png(p_sender, file.path(fig_dir, paste0(recv_safe, "_senderDelta_heat.png")), w=9, h=4.5)
  }

  # ---------- (D) Evidence dot plot: activity vs max sender delta ----------
  if (!is.null(sender) && nrow(sender) > 0) {
    max_delta <- sender %>%
      group_by(ligand) %>%
      summarise(max_sender_delta = max(delta, na.rm = TRUE), .groups="drop")

    dot <- lig_full %>%
      inner_join(max_delta, by = c("test_ligand" = "ligand")) %>%
      mutate(rank = row_number()) %>%
      slice_head(n = 50)

    p_dot <- ggplot(dot, aes(x = aupr_corrected, y = max_sender_delta, label = test_ligand)) +
      geom_point() +
      geom_text(check_overlap = TRUE, size = 3, nudge_y = 0.01) +
      labs(title = paste0(recv, ": Ligand evidence (activity vs sender Δ)"),
           x = "AUPR_corrected (ligand activity)", y = "max sender Δ expr (PTL–Control)") +
      theme_classic(base_size = 12)

    save_png(p_dot, file.path(fig_dir, paste0(recv_safe, "_evidence_dot.png")), w=7, h=5)
  }

  # store for cross-receiver
  per_receiver_summary[[recv]] <- lig_full %>%
    select(test_ligand, aupr_corrected, auroc, pearson, aupr) %>%
    mutate(receiver = recv)
}

# -------------------- 2) CROSS-RECEIVER SUMMARIES --------------------
all_lig <- bind_rows(per_receiver_summary)

if (nrow(all_lig) > 0) {

  # (E) Ligand recurrence: how many receivers is each ligand top-N?
  topN_cross <- 30
  lig_rec <- all_lig %>%
    group_by(receiver) %>%
    arrange(desc(aupr_corrected), .by_group = TRUE) %>%
    slice_head(n = topN_cross) %>%
    ungroup() %>%
    count(test_ligand, name = "n_receivers") %>%
    arrange(desc(n_receivers), test_ligand)

  p_rec <- ggplot(lig_rec, aes(x = reorder(test_ligand, n_receivers), y = n_receivers)) +
    geom_col() +
    coord_flip() +
    labs(title = paste0("Ligand recurrence across receivers (top ", topN_cross, " per receiver)"),
         x = NULL, y = "# receivers") +
    theme_classic(base_size = 12)

  save_png(p_rec, file.path(fig_dir, "ALL_ligand_recurrence_bar.png"), w=8, h=10)

  # (F) Receiver x ligand dot heatmap (AUPR_corrected)
  # choose ligands that appear in >=2 receivers, cap to top 40 by recurrence
  lig_keep <- lig_rec %>% filter(n_receivers >= 2) %>% slice_head(n = 40) %>% pull(test_ligand)

  mat_df <- all_lig %>%
    filter(test_ligand %in% lig_keep) %>%
    select(receiver, test_ligand, aupr_corrected)

  # fill missing combos with 0
  mat_df <- mat_df %>%
    tidyr::complete(receiver, test_ligand, fill = list(aupr_corrected = 0))

  # order ligands by recurrence then mean activity
  lig_order <- mat_df %>%
    group_by(test_ligand) %>%
    summarise(mu = mean(aupr_corrected), .groups="drop") %>%
    left_join(lig_rec, by = "test_ligand") %>%
    arrange(desc(n_receivers), desc(mu)) %>%
    pull(test_ligand)

  recv_order <- unique(mat_df$receiver)

  mat_df$test_ligand <- factor(mat_df$test_ligand, levels = lig_order)
  mat_df$receiver <- factor(mat_df$receiver, levels = recv_order)

  p_mat <- ggplot(mat_df, aes(x = test_ligand, y = receiver, size = aupr_corrected, fill = aupr_corrected)) +
    geom_point(shape = 21, color = "black", stroke = 0.2) +
    labs(title = "Receiver × ligand activity (AUPR_corrected)",
         x = "Ligand", y = "Receiver", size = "AUPR", fill = "AUPR") +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_png(p_mat, file.path(fig_dir, "ALL_receiver_ligand_dot.png"), w=12, h=5.5)

  # Save the cross tables too
  write_csv(all_lig, file.path(fig_dir, "ALL_ligand_activities_long.csv"))
  write_csv(lig_rec, file.path(fig_dir, "ALL_ligand_recurrence.csv"))
}

message("DONE. Figures written to: ", normalizePath(fig_dir))













suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(readr)
    library(stringr)
    library(ggplot2)
    library(patchwork)
    library(cowplot)
})

# -------------------- CONFIG --------------------
base_dir <- "/home/ebarrozo/visiumhd_mnp/results/audra/MNP_interactomes_NicheNet_v2/Preterm_PTL_vs_Control"
fig_dir  <- file.path(base_dir, "figs_v1")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

top_ligands_plot <- 25      # bars per receiver plot
top_targets_plot <- 40      # targets to show in ligand->target heatmap
top_links_per_ligand <- 20  # cap targets per ligand for heatmap readability

# -------------------- HELPERS --------------------
list_receiver_dirs <- function(base_dir) {
    dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
    dirs[grepl("/receiver_", dirs)]
}

receiver_name_from_dir <- function(d) {
    sub("^receiver_", "", basename(d))
}

read_if_exists <- function(path) {
    if (!file.exists(path)) return(NULL)
    readr::read_csv(path, show_col_types = FALSE)
}

safe_ggsave <- function(p, filename, w = 10, h = 6, dpi = 300) {
    ggplot2::ggsave(filename = filename, plot = p, width = w, height = h, dpi = dpi, limitsize = FALSE)
}

# -------------------- LOAD ALL RECEIVERS --------------------
recv_dirs <- list_receiver_dirs(base_dir)
stopifnot(length(recv_dirs) > 0)

all_ligacts <- list()
all_sender  <- list()
all_links   <- list()

for (d in recv_dirs) {
    rname <- receiver_name_from_dir(d)
    
    ligacts <- read_if_exists(file.path(d, "ligand_activities_full.csv"))
    sender  <- read_if_exists(file.path(d, "sender_ligand_expression_delta.csv"))
    links   <- read_if_exists(file.path(d, "ligand_target_links_topLigands.csv"))
    de_up   <- read_if_exists(file.path(d, "receiver_DE_up_targets.csv")) # optional for annotation later
    
    if (!is.null(ligacts)) all_ligacts[[rname]] <- ligacts %>% mutate(receiver = rname)
    if (!is.null(sender))  all_sender[[rname]]  <- sender  %>% mutate(receiver = rname)
    if (!is.null(links))   all_links[[rname]]   <- links   %>% mutate(receiver = rname)
}

ligacts_all <- bind_rows(all_ligacts)
sender_all  <- bind_rows(all_sender)
links_all   <- bind_rows(all_links)

# Sanity
stopifnot(nrow(ligacts_all) > 0)
message("Loaded: ligand_activities rows=", nrow(ligacts_all),
        " | sender_expr rows=", nrow(sender_all),
        " | links rows=", nrow(links_all))

# -------------------- (1) Ligand activity barplots per receiver --------------------
make_ligand_activity_bar <- function(df, receiver, n = 25) {
    dd <- df %>%
        filter(receiver == !!receiver) %>%
        arrange(desc(aupr_corrected)) %>%
        slice_head(n = n) %>%
        mutate(test_ligand = factor(test_ligand, levels = rev(test_ligand)))
    
    ggplot(dd, aes(x = test_ligand, y = aupr_corrected)) +
        geom_col() +
        coord_flip() +
        labs(title = paste0("NicheNet ligand activity (", receiver, ")"),
             x = NULL, y = "AUPR_corrected") +
        theme_bw(base_size = 12)
}

for (r in unique(ligacts_all$receiver)) {
    p <- make_ligand_activity_bar(ligacts_all, r, n = top_ligands_plot)
    safe_ggsave(p, file.path(fig_dir, paste0("ligand_activity_bar_", r, ".png")), w = 8, h = 7)
}

# Also a combined faceted view (top 15 each; fast overview)
ligacts_top <- ligacts_all %>%
    group_by(receiver) %>%
    arrange(desc(aupr_corrected)) %>%
    slice_head(n = 15) %>%
    ungroup() %>%
    mutate(test_ligand = reorder_within(test_ligand, aupr_corrected, receiver))

p_facet <- ggplot(ligacts_top, aes(x = test_ligand, y = aupr_corrected)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~ receiver, scales = "free_y") +
    scale_x_reordered() +
    labs(title = "Top ligand activities per receiver", x = NULL, y = "AUPR_corrected") +
    theme_bw(base_size = 11)

safe_ggsave(p_facet, file.path(fig_dir, "ligand_activity_bar_facet.png"), w = 14, h = 10)

# Helpers for reorder_within (avoids extra deps)
reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
    new_x <- paste(x, within, sep = sep)
    stats::reorder(new_x, by, FUN = fun)
}
scale_x_reordered <- function(..., sep = "___") {
    ggplot2::scale_x_discrete(labels = function(x) gsub(paste0(sep, ".*$"), "", x), ...)
}

# -------------------- (2) Sender × Ligand delta heatmap --------------------
# sender_all has columns: celltype, ligand, Control, PTL, delta, receiver (from your pipeline)
make_sender_delta_heat <- function(df, receiver, n_ligs = 25) {
    dd <- df %>%
        filter(receiver == !!receiver) %>%
        group_by(ligand) %>%
        summarize(max_abs_delta = max(abs(delta), na.rm = TRUE), .groups = "drop") %>%
        arrange(desc(max_abs_delta)) %>%
        slice_head(n = n_ligs) %>%
        pull(ligand)
    
    hm <- df %>%
        filter(receiver == !!receiver, ligand %in% dd) %>%
        mutate(ligand = factor(ligand, levels = dd),
               celltype = factor(celltype, levels = sort(unique(celltype)))) %>%
        select(celltype, ligand, delta)
    
    ggplot(hm, aes(x = ligand, y = celltype, fill = delta)) +
        geom_tile() +
        coord_fixed() +
        labs(title = paste0("Sender expression Δ (PTL − Control): ", receiver),
             x = NULL, y = NULL) +
        theme_bw(base_size = 11) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

if (nrow(sender_all) > 0) {
    for (r in unique(sender_all$receiver)) {
        p <- make_sender_delta_heat(sender_all, r, n_ligs = min(top_ligands_plot, length(unique(sender_all$ligand))))
        safe_ggsave(p, file.path(fig_dir, paste0("sender_delta_heat_", r, ".png")), w = 12, h = 6)
    }
}

# -------------------- (3) Ligand → Target weighted heatmap --------------------
# links_all columns: ligand, target, weight, receiver
make_ligand_target_heat <- function(df_links, receiver, top_ligs = 20, top_targets = 40, top_links_per_lig = 20) {
    dd <- df_links %>%
        filter(receiver == !!receiver) %>%
        filter(!is.na(target), !is.na(weight)) %>%
        group_by(ligand) %>%
        summarize(max_w = max(weight, na.rm = TRUE), .groups = "drop") %>%
        arrange(desc(max_w)) %>%
        slice_head(n = top_ligs) %>%
        pull(ligand)
    
    sub <- df_links %>%
        filter(receiver == !!receiver, ligand %in% dd) %>%
        filter(!is.na(target), !is.na(weight)) %>%
        group_by(ligand) %>%
        arrange(desc(weight)) %>%
        slice_head(n = top_links_per_lig) %>%
        ungroup()
    
    # choose top targets globally within this receiver subset
    targ_keep <- sub %>%
        group_by(target) %>%
        summarize(max_w = max(weight, na.rm = TRUE), .groups = "drop") %>%
        arrange(desc(max_w)) %>%
        slice_head(n = top_targets) %>%
        pull(target)
    
    hm <- sub %>%
        filter(target %in% targ_keep) %>%
        mutate(ligand = factor(ligand, levels = dd),
               target = factor(target, levels = rev(targ_keep)))
    
    ggplot(hm, aes(x = ligand, y = target, fill = weight)) +
        geom_tile() +
        labs(title = paste0("Ligand → target links (weights): ", receiver),
             x = NULL, y = NULL) +
        theme_bw(base_size = 11) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

if (nrow(links_all) > 0) {
    for (r in unique(links_all$receiver)) {
        p <- make_ligand_target_heat(links_all, r,
                                     top_ligs = min(20, length(unique(links_all$ligand))),
                                     top_targets = top_targets_plot,
                                     top_links_per_lig = top_links_per_ligand)
        safe_ggsave(p, file.path(fig_dir, paste0("ligand_target_heat_", r, ".png")), w = 12, h = 10)
    }
}

# -------------------- (4) Cross-receiver summary tables + plots --------------------
# Which ligands recur across receivers (top 30 by aupr_corrected per receiver)
lig_summary <- ligacts_all %>%
    group_by(receiver) %>%
    arrange(desc(aupr_corrected)) %>%
    slice_head(n = 30) %>%
    ungroup() %>%
    group_by(test_ligand) %>%
    summarize(
        n_receivers = n_distinct(receiver),
        mean_auprc = mean(aupr_corrected, na.rm = TRUE),
        max_auprc  = max(aupr_corrected, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    arrange(desc(n_receivers), desc(mean_auprc))

write_csv(lig_summary, file.path(fig_dir, "summary_ligand_recurrence.csv"))

p_rec <- lig_summary %>%
    slice_head(n = 40) %>%
    mutate(test_ligand = factor(test_ligand, levels = rev(test_ligand))) %>%
    ggplot(aes(x = test_ligand, y = n_receivers)) +
    geom_col() +
    coord_flip() +
    labs(title = "Ligand recurrence across receivers (top-30 lists)",
         x = NULL, y = "# receivers where ligand is top-30") +
    theme_bw(base_size = 12)

safe_ggsave(p_rec, file.path(fig_dir, "summary_ligand_recurrence_bar.png"), w = 8, h = 10)

# Sender delta summary: max delta per ligand across senders/receivers
if (nrow(sender_all) > 0) {
    sender_summary <- sender_all %>%
        group_by(receiver, ligand) %>%
        summarize(max_delta = max(delta, na.rm = TRUE), .groups = "drop")
    
    write_csv(sender_summary, file.path(fig_dir, "summary_sender_max_delta_by_receiver.csv"))
}

message("DONE. Figures in: ", normalizePath(fig_dir))

























# Viz module inspired by:
# https://raw.githubusercontent.com/saeyslab/nichenetr/master/vignettes/seurat_steps_prioritization.md

suppressPackageStartupMessages({
  library(nichenetr)
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(RColorBrewer)
})

# --- assumes these already exist from your run ---
# obj_lr (Seurat), lr_network, ligand_target_matrix
# condition_A_filter, condition_B_filter, condition_A_name, condition_B_name
# senders, receivers, celltype_col, assay_use, slot_use
# out_main

viz_dir <- file.path(out_main, "FIGS_prioritization")
dir.create(viz_dir, showWarnings = FALSE, recursive = TRUE)

# helper: safe subset by meta filters (no eval/parse)
subset_by_filters_safe <- function(seu, filters) {
  md <- seu[[]]
  missing_cols <- setdiff(names(filters), colnames(md))
  if (length(missing_cols)) stop("Missing metadata cols: ", paste(missing_cols, collapse=", "))
  keep <- rep(TRUE, nrow(md))
  for (nm in names(filters)) keep <- keep & as.character(md[[nm]]) %in% as.character(filters[[nm]])
  subset(seu, cells = rownames(md)[keep])
}

# Prepare A/B objects once
obj_A <- subset_by_filters_safe(obj_lr, condition_A_filter)
obj_B <- subset_by_filters_safe(obj_lr, condition_B_filter)

# Add Condition column + merge once for convenience
comb <- merge(obj_A, y = obj_B, add.cell.ids = c(condition_A_name, condition_B_name))
comb$Condition <- ifelse(grepl(paste0("^", condition_A_name, "_"), colnames(comb)),
                         condition_A_name, condition_B_name)

# Make sure identities are your celltype labels
Idents(comb) <- celltype_col

# enforce uppercase consistency (recommended)
lr_network$from <- toupper(lr_network$from)
lr_network$to   <- toupper(lr_network$to)
colnames(ligand_target_matrix) <- toupper(colnames(ligand_target_matrix))
rownames(ligand_target_matrix) <- toupper(rownames(ligand_target_matrix))

# Color map across all involved cell types (senders + receivers)
all_ct <- sort(unique(c(senders, receivers)))
pal <- suppressWarnings(brewer.pal(max(3, min(12, length(all_ct))), "Set3"))
celltype_colors <- rep(pal, length.out = length(all_ct)) |> setNames(all_ct)

# main loop
for (receiver in receivers) {
  message("[Prioritization viz] Receiver: ", receiver)
  od <- file.path(viz_dir, paste0("receiver_", gsub("[^A-Za-z0-9]+","_", receiver)))
  dir.create(od, showWarnings = FALSE, recursive = TRUE)

  # expressed genes receiver/senders (use nichenetr helpers on Seurat object)
  expressed_genes_receiver <- get_expressed_genes(receiver, comb, pct = 0.05)
  list_expressed_sender <- senders %>% unique() %>% lapply(get_expressed_genes, comb, pct = 0.05)
  expressed_genes_sender <- unique(unlist(list_expressed_sender))

  all_ligands   <- unique(lr_network$from)
  all_receptors <- unique(lr_network$to)

  expressed_ligands   <- intersect(all_ligands, toupper(expressed_genes_sender))
  expressed_receptors <- intersect(all_receptors, toupper(expressed_genes_receiver))

  lr_network_filtered <- lr_network %>%
    dplyr::select(from, to) %>%
    dplyr::filter(from %in% expressed_ligands, to %in% expressed_receptors) %>%
    dplyr::distinct(from, to)

  if (nrow(lr_network_filtered) == 0) {
    writeLines("No LR pairs after expression filtering.", file.path(od, "SKIP.txt"))
    next
  }

  # Receiver DE: PTL vs Control within receiver cells
  seurat_obj_receiver <- subset(comb, idents = receiver)
  de_receiver <- tryCatch(
    FindMarkers(
      object = seurat_obj_receiver,
      ident.1 = condition_B_name, ident.2 = condition_A_name,
      group.by = "Condition",
      min.pct = 0.05, logfc.threshold = 0,
      assay = assay_use, slot = slot_use
    ) %>% rownames_to_column("gene") %>% mutate(gene = toupper(gene)),
    error = function(e) NULL
  )
  if (is.null(de_receiver) || nrow(de_receiver) == 0) {
    writeLines("Receiver DE failed or empty.", file.path(od, "SKIP.txt"))
    next
  }

  # geneset of interest (up in PTL)
  fc_col <- if ("avg_log2FC" %in% names(de_receiver)) "avg_log2FC" else if ("avg_logFC" %in% names(de_receiver)) "avg_logFC" else NA
  if (is.na(fc_col)) stop("No avg_log2FC/avg_logFC column in DE.")

  geneset_oi <- de_receiver %>%
    dplyr::filter(p_val_adj <= 0.05, .data[[fc_col]] >= 0.10) %>%
    dplyr::arrange(dplyr::desc(.data[[fc_col]])) %>%
    dplyr::pull(gene) %>%
    intersect(rownames(ligand_target_matrix))

  if (length(geneset_oi) < 20) {
    writeLines("Too few DE up targets for geneset_oi.", file.path(od, "SKIP.txt"))
    next
  }

  background_expressed_genes <- intersect(toupper(expressed_genes_receiver), rownames(ligand_target_matrix))

  # Ligand activity
  potential_ligands <- lr_network_filtered %>% pull(from) %>% unique()
  ligand_activities <- predict_ligand_activities(
    geneset = geneset_oi,
    background_expressed_genes = background_expressed_genes,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = potential_ligands
  ) %>% arrange(desc(aupr_corrected))

  # Info tables + prioritization table (this enables circos + mushroom)
  info_tables <- generate_info_tables(
    seurat_obj = comb,
    celltype_colname = celltype_col,
    senders_oi = senders,
    receivers_oi = receiver,
    lr_network = lr_network_filtered,
    condition_colname = "Condition",
    condition_oi = condition_B_name,
    condition_reference = condition_A_name,
    scenario = "case_control"
  )

  prior_table <- generate_prioritization_tables(
    sender_receiver_info = info_tables$sender_receiver_info %>% distinct(),
    sender_receiver_de   = info_tables$sender_receiver_de,
    ligand_activities    = ligand_activities,
    lr_condition_de      = info_tables$lr_condition_de %>% distinct(),
    scenario = "case_control"
  )

  write_csv(prior_table, file.path(od, "prior_table.csv"))
  write_csv(ligand_activities, file.path(od, "ligand_activities.csv"))

  # ---- Circos plot (top LR pairs by prioritization_score) ----
  prior_table_oi <- prior_table %>% slice_max(prioritization_score, n = 50)
  circ <- make_circos_lr(
    prior_table_oi,
    colors_sender = celltype_colors,
    colors_receiver = celltype_colors
  )
  ggsave(file.path(od, "circos_top50.png"), circ, width = 8, height = 8, dpi = 300, bg = "white")

  # ---- Mushroom plot (receiver-specific) ----
  mush <- make_mushroom_plot(
    prior_table %>% filter(receiver == !!receiver),
    top_n = 30,
    true_color_range = TRUE,
    show_rankings = TRUE,
    show_all_datapoints = TRUE
  ) + theme(axis.title.x = element_text(hjust = 0.25))

  ggsave(file.path(od, "mushroom_top30.png"), mush, width = 10, height = 5.5, dpi = 300, bg = "white")
}

message("DONE. Prioritization figures under: ", normalizePath(viz_dir))




















suppressPackageStartupMessages({
  library(nichenetr)
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(RColorBrewer)
})

# ---------- CONFIG ----------
viz_dir <- file.path(out_main, "FIGS_prioritization")
dir.create(viz_dir, showWarnings = FALSE, recursive = TRUE)

# downsample per (Condition x celltype) to make nichenetr info tables stable
cap_per_ct_per_cond <- 5000   # start here; can go 2000 if still heavy

# ---------- HELPERS ----------
downsample_cond_celltype <- function(seu, cond_col = "Condition", ct_col = celltype_col, cap = 5000) {
  md <- seu[[]]
  md[[cond_col]] <- as.character(md[[cond_col]])
  md[[ct_col]]   <- as.character(md[[ct_col]])
  md <- md[!is.na(md[[cond_col]]) & !is.na(md[[ct_col]]), , drop = FALSE]

  keep <- unlist(lapply(split(rownames(md), paste(md[[cond_col]], md[[ct_col]], sep="__")), function(cells){
    if (length(cells) <= cap) return(cells)
    sample(cells, cap)
  }), use.names = FALSE)

  subset(seu, cells = keep)
}

# Safe color map
all_ct <- sort(unique(c(senders, receivers)))
pal <- suppressWarnings(brewer.pal(max(3, min(12, length(all_ct))), "Set3"))
celltype_colors <- rep(pal, length.out = length(all_ct)) |> setNames(all_ct)

# ---------- PREP ----------
# You already have comb built; if not, rebuild it once as before.
# CRITICAL: make sure these are plain vectors in meta
comb$Condition <- as.character(comb$Condition)
comb[[celltype_col]] <- as.character(comb[[celltype_col]][,1])

# Downsample ONCE for all receivers (massive stability improvement)
comb_ds <- downsample_cond_celltype(comb, cond_col = "Condition", ct_col = celltype_col, cap = cap_per_ct_per_cond)

# Re-set identities
Idents(comb_ds) <- celltype_col

# enforce uppercase consistency once
lr_network$from <- toupper(lr_network$from)
lr_network$to   <- toupper(lr_network$to)
colnames(ligand_target_matrix) <- toupper(colnames(ligand_target_matrix))
rownames(ligand_target_matrix) <- toupper(rownames(ligand_target_matrix))

seuratObj = comb_ds




args(nichenetr::generate_info_tables)
# ---------- MAIN LOOP ----------
for (receiver in receivers) {
  message("[Prioritization viz DS] Receiver: ", receiver)
  od <- file.path(viz_dir, paste0("receiver_", gsub("[^A-Za-z0-9]+","_", receiver)))
  dir.create(od, showWarnings = FALSE, recursive = TRUE)

  # expressed genes receiver/senders on DOWN-SAMPLED object
  expressed_genes_receiver <- get_expressed_genes(receiver, comb_ds, pct = 0.05)
  list_expressed_sender <- senders %>% unique() %>% lapply(get_expressed_genes, comb_ds, pct = 0.05)
  expressed_genes_sender <- unique(unlist(list_expressed_sender))

  expressed_ligands   <- intersect(unique(lr_network$from), toupper(expressed_genes_sender))
  expressed_receptors <- intersect(unique(lr_network$to),  toupper(expressed_genes_receiver))

  lr_network_filtered <- lr_network %>%
    dplyr::select(from, to) %>%
    dplyr::filter(from %in% expressed_ligands, to %in% expressed_receptors) %>%
    dplyr::distinct(from, to)

  if (nrow(lr_network_filtered) == 0) {
    writeLines("No LR pairs after expression filtering.", file.path(od, "SKIP.txt"))
    next
  }

  # Receiver DE on downsampled receiver cells (Condition-aware)
  seurat_obj_receiver <- subset(comb_ds, idents = receiver)

  de_receiver <- tryCatch(
    FindMarkers(
      object = seurat_obj_receiver,
      group.by = "Condition",
      ident.1 = condition_B_name, ident.2 = condition_A_name,
      min.pct = 0.05, logfc.threshold = 0,
      assay = assay_use, slot = slot_use
    ) %>% rownames_to_column("gene") %>% mutate(gene = toupper(gene)),
    error = function(e) NULL
  )

  if (is.null(de_receiver) || nrow(de_receiver) == 0) {
    writeLines("Receiver DE failed or empty.", file.path(od, "SKIP.txt"))
    next
  }

  fc_col <- if ("avg_log2FC" %in% names(de_receiver)) "avg_log2FC" else if ("avg_logFC" %in% names(de_receiver)) "avg_logFC" else NA
  if (is.na(fc_col)) stop("No avg_log2FC/avg_logFC column in DE.")

  geneset_oi <- de_receiver %>%
    dplyr::filter(p_val_adj <= 0.05, .data[[fc_col]] >= 0.10) %>%
    dplyr::arrange(dplyr::desc(.data[[fc_col]])) %>%
    dplyr::pull(gene) %>%
    intersect(rownames(ligand_target_matrix))

  if (length(geneset_oi) < 20) {
    writeLines("Too few DE up targets for geneset_oi.", file.path(od, "SKIP.txt"))
    next
  }

  background_expressed_genes <- intersect(toupper(expressed_genes_receiver), rownames(ligand_target_matrix))

  potential_ligands <- unique(lr_network_filtered$from)
  ligand_activities <- predict_ligand_activities(
    geneset = geneset_oi,
    background_expressed_genes = background_expressed_genes,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = potential_ligands
  ) %>% arrange(desc(aupr_corrected))

  # ---- This is the step that crashed before; run on DOWN-SAMPLED object ----
  info_tables <- tryCatch(
  generate_info_tables(
    seuratObj = comb_ds,
    celltype_colname = celltype_col,
    senders_oi = senders,
    receivers_oi = receiver,
    lr_network = lr_network_filtered,
    condition_colname = "Condition",
    condition_oi = condition_B_name,
    condition_reference = condition_A_name,
    scenario = "case_control"
  ),
  error = function(e) {
    writeLines(paste0("generate_info_tables failed: ", conditionMessage(e)), file.path(od, "ERROR_info_tables.txt"))
    NULL
  }
)
  if (is.null(info_tables)) next

  prior_table <- tryCatch(
    generate_prioritization_tables(
      sender_receiver_info = info_tables$sender_receiver_info %>% distinct(),
      sender_receiver_de   = info_tables$sender_receiver_de,
      ligand_activities    = ligand_activities,
      lr_condition_de      = info_tables$lr_condition_de %>% distinct(),
      scenario = "case_control"
    ),
    error = function(e) {
      writeLines(paste0("generate_prioritization_tables failed: ", conditionMessage(e)), file.path(od, "ERROR_prior_tables.txt"))
      NULL
    }
  )
  if (is.null(prior_table) || nrow(prior_table) == 0) {
    writeLines("prior_table empty.", file.path(od, "SKIP.txt"))
    next
  }

  write_csv(prior_table, file.path(od, "prior_table.csv"))
  write_csv(ligand_activities, file.path(od, "ligand_activities.csv"))

  # Circos (top 50)
  prior_table_oi <- prior_table %>% slice_max(prioritization_score, n = 50)
  circ <- make_circos_lr(prior_table_oi, colors_sender = celltype_colors, colors_receiver = celltype_colors)
  ggsave(file.path(od, "circos_top50.png"), circ, width = 8, height = 8, dpi = 300, bg = "white")

  # Mushroom
  mush <- make_mushroom_plot(
    prior_table %>% filter(receiver == !!receiver),
    top_n = 30, true_color_range = TRUE, show_rankings = TRUE, show_all_datapoints = TRUE
  )
  ggsave(file.path(od, "mushroom_top30.png"), mush, width = 10, height = 5.5, dpi = 300, bg = "white")
}

message("DONE. Prioritization DS figures under: ", normalizePath(viz_dir))
























##########

debug_receiver <- function(receiver_name, seu_use = comb_ds) {
  message("==== ", receiver_name, " ====")

  expressed_genes_receiver <- get_expressed_genes(receiver_name, seu_use, pct = 0.05)
  list_expressed_sender <- senders %>% unique() %>% lapply(get_expressed_genes, seu_use, pct = 0.05)
  expressed_genes_sender <- unique(unlist(list_expressed_sender))

  expressed_ligands   <- intersect(unique(lr_network$from), toupper(expressed_genes_sender))
  expressed_receptors <- intersect(unique(lr_network$to),  toupper(expressed_genes_receiver))

  lr_network_filtered <- lr_network %>%
    dplyr::select(from, to) %>%
    dplyr::filter(from %in% expressed_ligands, to %in% expressed_receptors) %>%
    dplyr::distinct(from, to)

  cat("receiver expressed genes:", length(expressed_genes_receiver), "\n")
  cat("sender expressed genes:", length(expressed_genes_sender), "\n")
  cat("expressed ligands:", length(expressed_ligands), "\n")
  cat("expressed receptors:", length(expressed_receptors), "\n")
  cat("filtered LR pairs:", nrow(lr_network_filtered), "\n")

  seurat_obj_receiver <- subset(seu_use, idents = receiver_name)

  cat("\nCondition counts in receiver:\n")
  print(table(seurat_obj_receiver$Condition))

  de_receiver <- FindMarkers(
    object = seurat_obj_receiver,
    group.by = "Condition",
    ident.1 = condition_B_name, ident.2 = condition_A_name,
    min.pct = 0.05, logfc.threshold = 0,
    assay = assay_use, slot = slot_use
  ) %>%
    rownames_to_column("gene") %>%
    mutate(gene = toupper(gene))

  fc_col <- if ("avg_log2FC" %in% names(de_receiver)) "avg_log2FC" else "avg_logFC"

  cat("\nDE rows total:", nrow(de_receiver), "\n")
  cat("adj.P <= 0.05:", sum(de_receiver$p_val_adj <= 0.05, na.rm = TRUE), "\n")
  cat("adj.P <= 0.05 & logFC >= 0.10:",
      sum(de_receiver$p_val_adj <= 0.05 & de_receiver[[fc_col]] >= 0.10, na.rm = TRUE), "\n")
  cat("adj.P <= 0.05 & abs(logFC) >= 0.10:",
      sum(de_receiver$p_val_adj <= 0.05 & abs(de_receiver[[fc_col]]) >= 0.10, na.rm = TRUE), "\n")

  geneset_up <- de_receiver %>%
    filter(p_val_adj <= 0.05, .data[[fc_col]] >= 0.10) %>%
    pull(gene) %>%
    intersect(rownames(ligand_target_matrix))

  geneset_abs <- de_receiver %>%
    filter(p_val_adj <= 0.05, abs(.data[[fc_col]]) >= 0.10) %>%
    pull(gene) %>%
    intersect(rownames(ligand_target_matrix))

  cat("geneset_up after LT-matrix intersect:", length(geneset_up), "\n")
  cat("geneset_abs after LT-matrix intersect:", length(geneset_abs), "\n")

  invisible(list(
    de_receiver = de_receiver,
    lr_network_filtered = lr_network_filtered,
    geneset_up = geneset_up,
    geneset_abs = geneset_abs
  ))
}

dbg_nk  <- debug_receiver("NK Cell")
dbg_mac <- debug_receiver("Macrophage")


## now use full dataset for macs
debug_receiver <- function(receiver_name, seu_use = comb) {
  message("==== ", receiver_name, " ====")

  expressed_genes_receiver <- get_expressed_genes(receiver_name, seu_use, pct = 0.05)
  list_expressed_sender <- senders %>% unique() %>% lapply(get_expressed_genes, seu_use, pct = 0.05)
  expressed_genes_sender <- unique(unlist(list_expressed_sender))

  expressed_ligands   <- intersect(unique(lr_network$from), toupper(expressed_genes_sender))
  expressed_receptors <- intersect(unique(lr_network$to),  toupper(expressed_genes_receiver))

  lr_network_filtered <- lr_network %>%
    dplyr::select(from, to) %>%
    dplyr::filter(from %in% expressed_ligands, to %in% expressed_receptors) %>%
    dplyr::distinct(from, to)

  cat("receiver expressed genes:", length(expressed_genes_receiver), "\n")
  cat("sender expressed genes:", length(expressed_genes_sender), "\n")
  cat("expressed ligands:", length(expressed_ligands), "\n")
  cat("expressed receptors:", length(expressed_receptors), "\n")
  cat("filtered LR pairs:", nrow(lr_network_filtered), "\n")

  seurat_obj_receiver <- subset(seu_use, idents = receiver_name)

  cat("\nCondition counts in receiver:\n")
  print(table(seurat_obj_receiver$Condition))

  de_receiver <- FindMarkers(
    object = seurat_obj_receiver,
    group.by = "Condition",
    ident.1 = condition_B_name, ident.2 = condition_A_name,
    min.pct = 0.05, logfc.threshold = 0,
    assay = assay_use, slot = slot_use
  ) %>%
    rownames_to_column("gene") %>%
    mutate(gene = toupper(gene))

  fc_col <- if ("avg_log2FC" %in% names(de_receiver)) "avg_log2FC" else "avg_logFC"

  cat("\nDE rows total:", nrow(de_receiver), "\n")
  cat("adj.P <= 0.05:", sum(de_receiver$p_val_adj <= 0.05, na.rm = TRUE), "\n")
  cat("adj.P <= 0.05 & logFC >= 0.10:",
      sum(de_receiver$p_val_adj <= 0.05 & de_receiver[[fc_col]] >= 0.10, na.rm = TRUE), "\n")
  cat("adj.P <= 0.05 & abs(logFC) >= 0.10:",
      sum(de_receiver$p_val_adj <= 0.05 & abs(de_receiver[[fc_col]]) >= 0.10, na.rm = TRUE), "\n")

  geneset_up <- de_receiver %>%
    filter(p_val_adj <= 0.05, .data[[fc_col]] >= 0.10) %>%
    pull(gene) %>%
    intersect(rownames(ligand_target_matrix))

  geneset_abs <- de_receiver %>%
    filter(p_val_adj <= 0.05, abs(.data[[fc_col]]) >= 0.10) %>%
    pull(gene) %>%
    intersect(rownames(ligand_target_matrix))

  cat("geneset_up after LT-matrix intersect:", length(geneset_up), "\n")
  cat("geneset_abs after LT-matrix intersect:", length(geneset_abs), "\n")

  invisible(list(
    de_receiver = de_receiver,
    lr_network_filtered = lr_network_filtered,
    geneset_up = geneset_up,
    geneset_abs = geneset_abs
  ))
}
dbg_mac <- debug_receiver("Macrophage")






geneset_oi <- de_receiver %>%
  dplyr::filter(p_val_adj <= 0.10, .data[[fc_col]] >= 0.05) %>%
  dplyr::arrange(dplyr::desc(.data[[fc_col]])) %>%
  dplyr::pull(gene) %>%
  intersect(rownames(ligand_target_matrix))

if (length(geneset_oi) < 10) {
  writeLines("Too few DE up targets for geneset_oi.", file.path(od, "SKIP.txt"))
  next
}


de_receiver %>%
  filter(p_val_adj <= 0.05) %>%
  summarise(
    up = sum(.data[[fc_col]] >= 0.10, na.rm = TRUE),
    down = sum(.data[[fc_col]] <= -0.10, na.rm = TRUE)
  )





comb_ds2[["RNA"]] <- CreateAssayObject(data = sct_data)
DefaultAssay(comb_ds2) <- "RNA"

DefaultAssay(comb_ds2) <- "SCT"
sct_data <- GetAssayData(comb_ds2, assay = "SCT", slot = "data")  # dgCMatrix log1p

# sparse-safe pseudo-counts: operate only on non-zero entries
sct_counts <- sct_data
sct_counts@x <- expm1(sct_counts@x)  # base expm1 on numeric vector
sct_counts@x[sct_counts@x < 0] <- 0
sct_counts@x <- as.numeric(round(sct_counts@x))

comb_ds2[["RNA"]] <- CreateAssayObject(counts = sct_counts)
comb_ds2[["RNA"]]@data <- sct_data
DefaultAssay(comb_ds2) <- "RNA"

Assays(comb_ds2)
DefaultAssay(comb_ds2)
dim(GetAssayData(comb_ds2, assay="RNA", slot="data"))
dim(GetAssayData(comb_ds2, assay="RNA", slot="counts"))






receiver_use <- "Macrophage"
od <- file.path(viz_dir, "receiver_Macrophage")
dir.create(od, showWarnings = FALSE, recursive = TRUE)

# 1) expressed genes on DS object
expressed_genes_receiver <- get_expressed_genes(receiver_use, comb_ds2, pct = 0.05)
list_expressed_sender <- senders %>% unique() %>% lapply(get_expressed_genes, comb_ds2, pct = 0.05)
expressed_genes_sender <- unique(unlist(list_expressed_sender))

expressed_ligands   <- intersect(unique(lr_network$from), toupper(expressed_genes_sender))
expressed_receptors <- intersect(unique(lr_network$to),  toupper(expressed_genes_receiver))

lr_network_filtered <- lr_network %>%
  dplyr::select(from, to) %>%
  dplyr::filter(from %in% expressed_ligands, to %in% expressed_receptors) %>%
  dplyr::distinct(from, to)

# 2) DE on FULL object
seurat_obj_receiver <- subset(comb, idents = receiver_use)

de_receiver <- FindMarkers(
  object = seurat_obj_receiver,
  group.by = "Condition",
  ident.1 = condition_B_name, ident.2 = condition_A_name,
  min.pct = 0.05, logfc.threshold = 0,
  assay = assay_use, slot = slot_use
) %>%
  rownames_to_column("gene") %>%
  mutate(gene = toupper(gene))

fc_col <- if ("avg_log2FC" %in% names(de_receiver)) "avg_log2FC" else "avg_logFC"

p_cut <- 0.10
fc_cut <- 0.05
min_targets <- 10

geneset_oi <- de_receiver %>%
  dplyr::filter(p_val_adj <= p_cut, .data[[fc_col]] >= fc_cut) %>%
  dplyr::arrange(dplyr::desc(.data[[fc_col]])) %>%
  dplyr::pull(gene) %>%
  intersect(rownames(ligand_target_matrix))

length(geneset_oi)
stopifnot(length(geneset_oi) >= min_targets)

# 3) ligand activities
background_expressed_genes <- intersect(toupper(expressed_genes_receiver), rownames(ligand_target_matrix))
potential_ligands <- unique(lr_network_filtered$from)

ligand_activities_mac <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
) %>%
  arrange(desc(aupr_corrected))

# 4) info tables
info_tables_mac <- nichenetr::generate_info_tables(
  seuratObj = comb_ds2,
  celltype_colname = celltype_col,
  senders_oi = senders,
  receivers_oi = receiver_use,
  lr_network = lr_network_filtered,
  condition_colname = "Condition",
  condition_oi = condition_B_name,
  condition_reference = condition_A_name,
  scenario = "case_control"
)

# sanity check
unique(info_tables_mac$sender_receiver_info$receiver)

# 5) prioritization table
prior_table_mac <- generate_prioritization_tables(
  sender_receiver_info = info_tables_mac$sender_receiver_info %>% distinct(),
  sender_receiver_de   = info_tables_mac$sender_receiver_de,
  ligand_activities    = ligand_activities_mac,
  lr_condition_de      = info_tables_mac$lr_condition_de %>% distinct(),
  scenario = "case_control"
)

unique(prior_table_mac$receiver)

# 6) plot from macrophage-specific table
prior_table_oi <- prior_table_mac %>%
  dplyr::filter(.data$receiver == receiver_use) %>%
  slice_max(prioritization_score, n = 50)

circ <- make_circos_lr(
  prior_table_oi,
  colors_sender = celltype_colors,
  colors_receiver = celltype_colors
)

ggsave(
  file.path(od, "circos_top50_macrophage.png"),
  circ,
  width = 8, height = 8, dpi = 300, bg = "white"
)

mush <- make_mushroom_plot(
  prior_table_mac %>% dplyr::filter(.data$receiver == receiver_use),
  top_n = 30,
  true_color_range = TRUE,
  show_rankings = TRUE,
  show_all_datapoints = TRUE
)

ggsave(
  file.path(od, "mushroom_top30_macrophage.png"),
  mush,
  width = 10, height = 5.5, dpi = 300, bg = "white"
)




# same palette as UMAP
mypal3 <- get_palette("igv", 40)
cell_levels <- levels(obj$nn.predicted.celltype)
col_map <- setNames(mypal3[seq_along(cell_levels)], cell_levels)

# optional: inspect what names are actually needed here
sort(unique(c(prior_table_mac$sender, prior_table_mac$receiver)))

# keep only colors for levels present in the circos table
circos_col_map <- col_map[names(col_map) %in% unique(c(prior_table_mac$sender, prior_table_mac$receiver))]

prior_table_oi <- prior_table_mac %>%
  dplyr::filter(.data$receiver == receiver_use) %>%
  slice_max(prioritization_score, n = 50)

circ <- make_circos_lr(
  prior_table_oi,
  colors_sender   = circos_col_map,
  colors_receiver = circos_col_map
)
circ2 <- circ +
  theme(
    plot.margin = margin(10, 80, 10, 10),
    legend.margin = margin(0, 10, 0, 0),
    legend.box.margin = margin(0, 10, 0, 0)
  )

ggsave(
  file.path(od, "circos_top50_macrophage_igv.png"),
  circ2,
  width = 11, height = 8,
  dpi = 300, bg = "white"
)






receiver_use <- "NK Cell"
od <- file.path(viz_dir, "receiver_NK Cell")
dir.create(od, showWarnings = FALSE, recursive = TRUE)

# 1) expressed genes on DS object
expressed_genes_receiver <- get_expressed_genes(receiver_use, comb_ds2, pct = 0.05)
list_expressed_sender <- senders %>% unique() %>% lapply(get_expressed_genes, comb_ds2, pct = 0.05)
expressed_genes_sender <- unique(unlist(list_expressed_sender))

expressed_ligands   <- intersect(unique(lr_network$from), toupper(expressed_genes_sender))
expressed_receptors <- intersect(unique(lr_network$to),  toupper(expressed_genes_receiver))

lr_network_filtered <- lr_network %>%
  dplyr::select(from, to) %>%
  dplyr::filter(from %in% expressed_ligands, to %in% expressed_receptors) %>%
  dplyr::distinct(from, to)

# 2) DE on FULL object
seurat_obj_receiver <- subset(comb, idents = receiver_use)

de_receiver <- FindMarkers(
  object = seurat_obj_receiver,
  group.by = "Condition",
  ident.1 = condition_B_name, ident.2 = condition_A_name,
  min.pct = 0.05, logfc.threshold = 0,
  assay = assay_use, slot = slot_use
) %>%
  rownames_to_column("gene") %>%
  mutate(gene = toupper(gene))

fc_col <- if ("avg_log2FC" %in% names(de_receiver)) "avg_log2FC" else "avg_logFC"

p_cut <- 0.10
fc_cut <- 0.05
min_targets <- 10

geneset_oi <- de_receiver %>%
  dplyr::filter(p_val_adj <= p_cut, .data[[fc_col]] >= fc_cut) %>%
  dplyr::arrange(dplyr::desc(.data[[fc_col]])) %>%
  dplyr::pull(gene) %>%
  intersect(rownames(ligand_target_matrix))

length(geneset_oi)
stopifnot(length(geneset_oi) >= min_targets)

# 3) ligand activities
background_expressed_genes <- intersect(toupper(expressed_genes_receiver), rownames(ligand_target_matrix))
potential_ligands <- unique(lr_network_filtered$from)

ligand_activities_mac <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
) %>%
  arrange(desc(aupr_corrected))

# 4) info tables
info_tables_mac <- nichenetr::generate_info_tables(
  seuratObj = comb_ds2,
  celltype_colname = celltype_col,
  senders_oi = senders,
  receivers_oi = receiver_use,
  lr_network = lr_network_filtered,
  condition_colname = "Condition",
  condition_oi = condition_B_name,
  condition_reference = condition_A_name,
  scenario = "case_control"
)

# sanity check
unique(info_tables_mac$sender_receiver_info$receiver)

# 5) prioritization table
prior_table_mac <- generate_prioritization_tables(
  sender_receiver_info = info_tables_mac$sender_receiver_info %>% distinct(),
  sender_receiver_de   = info_tables_mac$sender_receiver_de,
  ligand_activities    = ligand_activities_mac,
  lr_condition_de      = info_tables_mac$lr_condition_de %>% distinct(),
  scenario = "case_control"
)

unique(prior_table_mac$receiver)

# 6) plot from NK Cell-specific table
prior_table_oi <- prior_table_mac %>%
  dplyr::filter(.data$receiver == receiver_use) %>%
  slice_max(prioritization_score, n = 50)

circ <- make_circos_lr(
  prior_table_oi,
  colors_sender = celltype_colors,
  colors_receiver = celltype_colors
)

ggsave(
  file.path(od, "circos_top50_NK Cell.png"),
  circ,
  width = 8, height = 8, dpi = 300, bg = "white"
)

mush <- make_mushroom_plot(
  prior_table_mac %>% dplyr::filter(.data$receiver == receiver_use),
  top_n = 30,
  true_color_range = TRUE,
  show_rankings = TRUE,
  show_all_datapoints = TRUE
)

ggsave(
  file.path(od, "mushroom_top30_NK Cell.png"),
  mush,
  width = 10, height = 5.5, dpi = 300, bg = "white"
)




# same palette as UMAP
mypal3 <- get_palette("igv", 40)
cell_levels <- levels(obj$nn.predicted.celltype)
col_map <- setNames(mypal3[seq_along(cell_levels)], cell_levels)

# optional: inspect what names are actually needed here
sort(unique(c(prior_table_mac$sender, prior_table_mac$receiver)))

# keep only colors for levels present in the circos table
circos_col_map <- col_map[names(col_map) %in% unique(c(prior_table_mac$sender, prior_table_mac$receiver))]

prior_table_oi <- prior_table_mac %>%
  dplyr::filter(.data$receiver == receiver_use) %>%
  slice_max(prioritization_score, n = 50)

circ <- make_circos_lr(
  prior_table_oi,
  colors_sender   = circos_col_map,
  colors_receiver = circos_col_map
)
circ2 <- circ +
  theme(
    plot.margin = margin(10, 80, 10, 10),
    legend.margin = margin(0, 10, 0, 0),
    legend.box.margin = margin(0, 10, 0, 0)
  )

ggsave(
  file.path(od, "circos_top50_NK Cell_igv.png"),
  circ2,
  width = 11, height = 8,
  dpi = 300, bg = "white"
)