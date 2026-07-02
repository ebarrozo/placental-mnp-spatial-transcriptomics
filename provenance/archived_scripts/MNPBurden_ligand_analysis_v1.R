## MNPBurden_ligand_analysis_v1.R
## server3
## Enrico Barrozo, Ph.D.
## Assistant Professor, Department of OB/GYN, Baylor College of Medicine and Texas Children’s Hospital
# 10.7.25



set.seed(1) # set seed for reproducibility
options(future.globals.maxSize = 300e9) ## 300 gb ram (out of 440 on server 3; adjust if on different server)
suppressPackageStartupMessages({
    library(dplyr); library(stringr); library(tibble); library(Seurat)
})
library(ggplot2)
library(cowplot)
library(patchwork)
library(ggrepel)
library(ggpubr)
library(ggplot2)
mypal3 <-get_palette("igv",40)

# #5050FF, #CE3D32, #749B58, #F0E685, #466983, #BA6338, #5DB1DD, #802268, #6BD76B, #D595A7, #924822, #837B8D, #C75127, #D58F5C, #7A65A5 



# NicheNet (GitHub)
#install.packages("devtools")
#devtools::install_github("saeyslab/nichenetr")

# CellChat (GitHub)
#install.packages("devtools")
#BiocManager::install("BiocNeighbors")
#devtools::install_github("sqjin/CellChat")   # see notes about NMF/circlize deps



library(Seurat); library(CellChat); library(nichenetr)

obj <-readRDS("/media/jochum00/Aagaard_Raid3/ebarrozo/seurat_all_ref_merged/with_annotations/SCT_chunks_mapped/final_sct_atlas.rds") # 7.5 gb object with SCT data, all 3,969,532 cells, and umap

## remove RNA reference data 
Idents(obj) <- "platform"
ncol(obj)
levels(obj)
DefaultAssay(obj) <- "SCT"
obj <- subset(x = obj, idents = c("Xenium","VisiumHD","VisiumV1"))
ncol(obj)

## remove RNA reference data 
Idents(obj) <- "MNPBurden"
ncol(obj)
levels(obj)
DefaultAssay(obj) <- "SCT"
obj2 <- subset(x = obj, idents = c("High","Low"))
ncol(obj2) # 3439267

obj <- subset(x = obj, idents = c("High","Low"))
ncol(obj) 


# mnp_cellchat_nichenet_dropin.R
# Drop-in utilities to run CellChat + NicheNet using your existing Seurat object `obj`
# and metadata columns: predicted.celltype, MNPBurden, Delivery_class, study_id.
# Uses SCT assay if present. Saves RDS objects and CSV/PNG summaries.
# mnp_cellchat_nichenet_dropin.R  (robust priors + Seurat v5 + colors)
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr); library(tidyr); library(purrr); library(readr); library(stringr); library(tibble)
  library(Matrix)
  library(CellChat)
  library(nichenetr)
  library(ggplot2)
})

# ===================== CONFIG =====================
Y_GENES <- c("DDX3Y","KDM5D","RPS4Y1","EIF1AY","UTY","USP9Y","ZFY","TXLNGY","TMSB4Y",
             "PRKY","PCDH11Y","NLGN4Y","AMELY")
lnc_pattern <- "^RP[0-9]+-|^AC[0-9]+\\.|^AL[0-9]+\\.|^LINC|^MIR|^SNHG|^HCG|^CT[CDE]-|^ENSG"

ord_celltypes <- c(
  "Syncytiotrophoblast","Villous Cytotrophoblast","Extravillous Trophoblast",
  "Fibroblast","Perivascular Cell","Placental Endothelial",
  "Hofbauer Cell","Macrophage","NK Cell","Decidual NK Cell","ILC3",
  "Plasma Cell","B Cell","T Cell","Treg","RCTD_lowconf"
)

# ===================== CELL TYPE COLORS =====================
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
.ct_colors_for_cellchat <- function(cellchat_obj) {
  grps <- levels(cellchat_obj@idents)
  cols <- ct_colors[grps]; cols[is.na(cols)] <- "#808080"
  as.character(cols)
}

# ===================== HELPERS =====================
`%||%` <- function(a, b) if (!is.null(a)) a else b
.use_SCT <- function(o) { if ("SCT" %in% names(o@assays)) DefaultAssay(o) <- "SCT"; o }
.check_meta <- function(o) {
  need <- c("predicted.celltype","MNPBurden","Delivery_class","study_id")
  miss <- setdiff(need, colnames(o@meta.data))
  if (length(miss)) stop("Missing meta columns: ", paste(miss, collapse=", "))
  invisible(TRUE)
}
factorize_ct <- function(o) {
  md <- o@meta.data
  if ("predicted.celltype" %in% names(md)) {
    md$predicted.celltype <- factor(as.character(md$predicted.celltype), levels = ord_celltypes)
    o@meta.data <- md
  }
  o
}
get_data_mat <- function(seu, assay = DefaultAssay(seu)) {
  lay <- tryCatch(SeuratObject::Layers(seu[[assay]]), error = function(e) NULL)
  if (!is.null(lay)) GetAssayData(seu, assay = assay, layer = "data") else GetAssayData(seu, assay = assay, slot = "data")
}
get_expressed_genes <- function(seu, celltype, pct = 0.10, assay = DefaultAssay(seu)) {
  md <- seu@meta.data
  cells <- rownames(md)[as.character(md$predicted.celltype) == celltype]
  if (length(cells) == 0) return(character(0))
  mat <- get_data_mat(seu, assay = assay)[, cells, drop = FALSE]
  rownames(mat)[Matrix::rowMeans(mat > 0) >= pct]
}
get_expressed_union <- function(seu, celltypes, pct = 0.10, assay = DefaultAssay(seu)) {
  unique(unlist(lapply(celltypes, function(ct) get_expressed_genes(seu, ct, pct = pct, assay = assay))))
}
safe_dir <- function(path) { dir.create(path, showWarnings = FALSE, recursive = TRUE); path }
str_sanitize <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)

# ===================== Robust NicheNet priors =====================

nn_load_priors <- function() {
  lr <- ltm <- wn <- NULL

  # 1) Try to load bundled datasets first
  for (pkg in c("nichenetr","nichenet")) {
    try(utils::data("weighted_networks", package = pkg, envir = environment()), silent = TRUE)
    try(utils::data("ligand_target_matrix", package = pkg, envir = environment()), silent = TRUE)
    try(utils::data("lr_network", package = pkg, envir = environment()), silent = TRUE)
  }
  if (exists("weighted_networks", inherits = FALSE)) wn  <- get("weighted_networks")
  if (exists("ligand_target_matrix", inherits = FALSE)) ltm <- get("ligand_target_matrix")
  if (exists("lr_network", inherits = FALSE)) lr <- get("lr_network")

  # 2) If missing, construct programmatically using exported constructors
  exps <- try(getNamespaceExports("nichenetr"), silent = TRUE); if (inherits(exps, "try-error")) exps <- character(0)
  if (is.null(wn) && "construct_weighted_networks" %in% exps) {
    # Pass species only if supported by this build
    cw <- nichenetr::construct_weighted_networks
    sp_ok <- FALSE
    fm <- try(formals(cw), silent = TRUE)
    if (!inherits(fm, "try-error") && !is.null(fm) && "species" %in% names(fm)) sp_ok <- TRUE
    wn <- if (sp_ok) cw(species = "human") else cw()
  }
  if (is.null(ltm) && "construct_ligand_target_matrix" %in% exps && !is.null(wn)) {
    ltm <- nichenetr::construct_ligand_target_matrix(weighted_networks = wn)
  }
  if (is.null(lr) && "get_weighted_ligand_receptor_links" %in% exps && !is.null(wn)) {
    lr_tbl <- nichenetr::get_weighted_ligand_receptor_links(weighted_networks = wn)
    if (!is.null(lr_tbl) && nrow(lr_tbl)) {
      lr <- lr_tbl %>% dplyr::transmute(from = ligand, to = receptor,
                                        database = if ("database" %in% colnames(lr_tbl)) database else "constructed")
    }
  }

  if (is.null(ltm) || is.null(wn)) stop("NicheNet priors not available: failed to build 'weighted_networks' and/or 'ligand_target_matrix'.")
  if (is.null(lr)) lr <- tibble::tibble(from = character(0), to = character(0), database = character(0))

  list(lr_network = lr, ligand_target_matrix = ltm, weighted_networks = wn)
}



# ===================== NICHE NET runner =====================
read_de_genes <- function(de_dir, ct, q = 0.10, lfc = 0.25, direction = c("UP","DOWN")) {
  direction <- match.arg(direction)
  fp <- file.path(de_dir, sprintf("pb_DE_%s.tsv", str_sanitize(ct)))
  if (!file.exists(fp)) return(character(0))
  de <- suppressMessages(readr::read_tsv(fp, show_col_types = FALSE))
  de <- dplyr::filter(de, !grepl(lnc_pattern, gene, perl = TRUE), !(gene %in% Y_GENES))
  sig <- dplyr::filter(de, adj.P.Val < q, abs(logFC) >= lfc)
  if (direction == "UP") sig <- dplyr::filter(sig, logFC > 0) else sig <- dplyr::filter(sig, logFC < 0)
  unique(sig$gene)
}

run_nichenet_mnp <- function(seu, delivery_class = c("Preterm","Term"),
                             receiver_cts, sender_cts,
                             de_dir, outdir,
                             pct_expr = 0.10, top_ligands = 25, top_targets = 200) {
  delivery_class <- match.arg(delivery_class)
  seu <- .use_SCT(seu); .check_meta(seu); seu <- factorize_ct(seu)

  sub <- subset(seu, subset = Delivery_class == delivery_class)
  if (ncol(sub) == 0) stop("No cells for Delivery_class == ", delivery_class)

  assay_use <- DefaultAssay(sub)
  expressed_receivers <- get_expressed_union(sub, receiver_cts, pct = pct_expr, assay = assay_use)
  expressed_senders   <- get_expressed_union(sub, sender_cts,   pct = pct_expr, assay = assay_use)

  geneset <- unique(unlist(lapply(receiver_cts, function(ct) read_de_genes(de_dir, ct, q = 0.10, lfc = 0.25, direction = "UP"))))
  if (!length(geneset)) message("NicheNet: no receiver DE genes found at q<=0.10; consider lowering threshold or using ranked set.")

  priors <- nn_load_priors()
  lr_network <- priors$lr_network
  ligand_target_matrix <- priors$ligand_target_matrix

  ligands_all <- if (!is.null(rownames(ligand_target_matrix))) rownames(ligand_target_matrix) else unique(ligand_target_matrix$ligand)
  ligands <- intersect(ligands_all, expressed_senders)

  receptors_all <- if ("to" %in% colnames(lr_network)) unique(lr_network$to) else character(0)
  receptors <- intersect(receptors_all, expressed_receivers)

  ligand_activities <- nichenetr::predict_ligand_activities(
    geneset                    = geneset,
    background_expressed_genes = expressed_receivers,
    ligand_target_matrix       = ligand_target_matrix
  ) %>% dplyr::arrange(dplyr::desc(pearson)) %>% dplyr::filter(test_ligand %in% ligands) %>% dplyr::slice_head(n = top_ligands)

  best_ligands <- ligand_activities$test_ligand
  l2t <- nichenetr::get_weighted_ligand_target_links(best_ligands, ligand_target_matrix, n = top_targets)

  if (length(receptors)) {
    lr_pairs <- dplyr::filter(lr_network, from %in% best_ligands, to %in% receptors)
  } else {
    lr_pairs <- tibble::tibble(from = character(0), to = character(0))
  }

  od <- safe_dir(file.path(outdir, paste0("NicheNet_", delivery_class)))
  readr::write_csv(ligand_activities, file.path(od, "ligand_activities.csv"))
  readr::write_csv(l2t,              file.path(od, "ligand_target_links.csv"))
  readr::write_csv(lr_pairs,         file.path(od, "ligand_receptor_pairs_expressed.csv"))
  writeLines(best_ligands,           file.path(od, "top_ligands.txt"))
  writeLines(geneset,                file.path(od, "receiver_geneset_used.txt"))

  if (nrow(ligand_activities)) {
    p <- ggplot(ligand_activities, aes(x = reorder(test_ligand, pearson), y = pearson)) +
      geom_col() + coord_flip() + theme_classic() +
      labs(x = "Ligand", y = "Activity (Pearson)", title = paste("NicheNet ligand activities —", delivery_class, paste(receiver_cts, collapse = "+")))
    ggsave(file.path(od, "ligand_activities.barplot.png"), p, width = 7, height = 6, dpi = 300, bg = "white")
  }

  if (nrow(l2t)) {
    tgt_keep <- l2t %>% dplyr::group_by(target) %>% dplyr::summarise(score = sum(weight), .groups="drop") %>% dplyr::arrange(dplyr::desc(score)) %>% dplyr::slice_head(n = min(50, n())) %>% dplyr::pull(target)
    ht <- dplyr::filter(l2t, target %in% tgt_keep)
    mat <- tidyr::pivot_wider(ht, names_from = ligand, values_from = weight, values_fill = 0)
    mat_df <- as.data.frame(mat); rownames(mat_df) <- mat_df$target; mat_df$target <- NULL
    mat_m <- as.matrix(mat_df)
    mat_m <- apply(mat_m, 2, function(x) { s <- stats::sd(x); if (!is.finite(s) || s == 0) x else (x - mean(x))/s })
    write.csv(mat_m, file.path(od, "ligand_target_heatmap_matrix.csv"))
  }
  invisible(list(ligand_activities = ligand_activities, ligand_target = l2t, lr_pairs = lr_pairs, geneset = geneset))
}

# ===================== CELLCHAT runner =====================
run_cellchat_mnp <- function(seu, delivery_class = c("Preterm","Term"),
                             outdir, min_cells = 50, min_cells_group = 10,
                             keep_signaling = NULL,
                             focus_sources = NULL, focus_targets = NULL) {
  delivery_class <- match.arg(delivery_class)
  seu <- .use_SCT(seu); .check_meta(seu); seu <- factorize_ct(seu)

  sub <- subset(seu, subset = Delivery_class == delivery_class)
  if (ncol(sub) == 0) stop("No cells for Delivery_class == ", delivery_class)

  tab <- table(sub$predicted.celltype)
  keep_ct <- names(tab)[tab >= min_cells]
  sub <- subset(sub, subset = predicted.celltype %in% keep_ct)
  if (length(keep_ct) < 2) stop("Too few cell types with >= ", min_cells, " cells in ", delivery_class)

  data.input <- get_data_mat(sub, assay = DefaultAssay(sub))
  meta.input <- sub@meta.data

  cellchat <- createCellChat(object = data.input, meta = meta.input, group.by = "predicted.celltype")
  cellchat@DB <- CellChatDB.human

  if (!is.null(keep_signaling)) {
    avail <- unique(CellChatDB.human$interaction$pathway_name)
    keep_signaling <- intersect(keep_signaling, avail)
  }

  cellchat <- subsetData(cellchat, keep.signaling = keep_signaling)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, raw.use = FALSE)
  cellchat <- filterCommunication(cellchat, min.cells = min_cells_group)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat)

  od <- safe_dir(file.path(outdir, paste0("CellChat_", delivery_class)))
  saveRDS(cellchat, file = file.path(od, "cellchat.rds"))
  net_df <- as.data.frame(cellchat@net$weight)
  write.csv(net_df, file.path(od, "net_weight_matrix.csv"), row.names = TRUE)

  pdf(file.path(od, "net_circle_all.pdf"), width = 8, height = 8)
  netVisual_circle(cellchat@net$weight, color.use = .ct_colors_for_cellchat(cellchat),
                   vertex.weight = as.numeric(table(cellchat@idents)), label.edge = FALSE)
  dev.off()

  if (!is.null(focus_sources) || !is.null(focus_targets)) {
    sources.use <- intersect(focus_sources %||% keep_ct, keep_ct)
    targets.use <- intersect(focus_targets %||% keep_ct, keep_ct)
    if (length(sources.use) && length(targets.use)) {
      pdf(file.path(od, "net_circle_focus.pdf"), width = 8, height = 8)
      netVisual_circle(cellchat@net$weight, color.use = .ct_colors_for_cellchat(cellchat),
                       sources.use = sources.use, targets.use = targets.use, label.edge = TRUE)
      dev.off()
    }
  }

  pdf(file.path(od, "net_heatmap.pdf"), width = 9, height = 7)
  netVisual_heatmap(cellchat, color.use = .ct_colors_for_cellchat(cellchat))
  dev.off()

  invisible(cellchat)
}

# ===================== WRAPPERS =====================
run_ptb_interactome <- function(seu,
                                de_dir = "pb_Preterm_MNPBurden_limma_SCT",
                                out_root = "MNP_interactomes") {
  receiver_cts <- c("Macrophage","Hofbauer Cell")
  all_cts <- sort(unique(as.character(seu$predicted.celltype)))
  sender_cts <- setdiff(all_cts, c(receiver_cts, "RCTD_lowconf", NA))
  outdir <- safe_dir(out_root)

  message("[PTB] NicheNet ...")
  nn <- run_nichenet_mnp(seu, delivery_class = "Preterm",
                         receiver_cts = receiver_cts, sender_cts = sender_cts,
                         de_dir = de_dir, outdir = file.path(outdir, "PTB"))

  message("[PTB] CellChat ...")
  cc <- run_cellchat_mnp(seu, delivery_class = "Preterm",
                         outdir = file.path(outdir, "PTB"),
                         keep_signaling = NULL,
                         focus_sources = c("Macrophage","Hofbauer Cell"),
                         focus_targets = c("Syncytiotrophoblast","Extravillous Trophoblast","Villous Cytotrophoblast",
                                           "Perivascular Cell","Placental Endothelial","Fibroblast"))
  invisible(list(nichenet = nn, cellchat = cc))
}

run_term_interactome <- function(seu,
                                 de_dir = "pb_Term_MNPBurden_limma_SCT",
                                 out_root = "MNP_interactomes") {
  receiver_cts <- c("NK Cell","Decidual NK Cell")
  all_cts <- sort(unique(as.character(seu$predicted.celltype)))
  sender_cts <- setdiff(all_cts, c(receiver_cts, "RCTD_lowconf", NA))
  outdir <- safe_dir(out_root)

  message("[Term] NicheNet ...")
  nn <- run_nichenet_mnp(seu, delivery_class = "Term",
                         receiver_cts = receiver_cts, sender_cts = sender_cts,
                         de_dir = de_dir, outdir = file.path(outdir, "Term"))

  message("[Term] CellChat ...")
  cc <- run_cellchat_mnp(seu, delivery_class = "Term",
                         outdir = file.path(outdir, "Term"),
                         keep_signaling = NULL,
                         focus_sources = c("NK Cell","Decidual NK Cell"),
                         focus_targets = c("Fibroblast","Perivascular Cell","Placental Endothelial",
                                           "Extravillous Trophoblast","Villous Cytotrophoblast","Syncytiotrophoblast"))
  invisible(list(nichenet = nn, cellchat = cc))
}



# ===================== EXAMPLE USAGE =====================
# Assuming `obj` is in memory:
res_ptb  <- run_ptb_interactome(obj,  de_dir = "/home/ebarrozo/visiumhd_mnp/results/seurat_all_v1/seurat_all_ref_merged/de_MNPBurden_Delivery_class_bulk_sensitive_go/pb_Preterm_MNPBurden_limma_SCT", out_root = "MNP_interactomes")
res_term <- run_term_interactome(obj, de_dir = "/home/ebarrozo/visiumhd_mnp/results/seurat_all_v1/seurat_all_ref_merged/de_MNPBurden_Delivery_class_bulk_sensitive_go/pb_Term_MNPBurden_limma_SCT",   out_root = "MNP_interactomes")
# Look in MNP_interactomes/PTB/* and MNP_interactomes/Term/* for CSV/PNG/RDS outputs.
