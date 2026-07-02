## cellphonedb_v1.R
## server3
## Enrico Barrozo, Ph.D.
## 2.27.26

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(readr)
  library(stringr)
})

# -------------------- USER CONFIG --------------------
setwd("/home/ebarrozo/visiumhd_mnp/results/audra")

obj <- readRDS("/media/jochum00/Aagaard_Raid3/ebarrozo/seurat_all_ref_merged/with_annotations/SCT_chunks_mapped/audra/mnp_obj_w_nn_predicted_celltype.rds")

out_root <- "MNP_interactomes_CellPhoneDB_v1"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

assay_use <- "SCT"
slot_use  <- "data"  # we'll derive pseudo-counts from SCT@data

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

# Downsample per celltype per condition to keep CPDB feasible (huge object)
cap_cells_per_celltype <- 3000

# CPDB params
cpdb_db <- "cellphonedb.zip"   # you will create this once (see notes below)
threads <- 16
iterations <- 1000             # permutation iterations
threshold_pct <- 0.10          # gene expression threshold (CPDB)
result_precision <- 3          # rounding

# -------------------- HELPERS --------------------
safe_dir <- function(path) { dir.create(path, showWarnings = FALSE, recursive = TRUE); path }

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
downsample_by_celltype <- function(seu, celltype_col, cap) {
  md <- seu[[]]
  md[[celltype_col]] <- as.character(md[[celltype_col]])
  keep_cells <- unlist(lapply(split(rownames(md), md[[celltype_col]]), function(cells){
    if (length(cells) <= cap) return(cells)
    sample(cells, cap)
  }), use.names = FALSE)
  subset(seu, cells = keep_cells)
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

downsample_by_celltype <- function(seu, celltype_col, cap) {
  md <- seu[[]]
  md[[celltype_col]] <- as.character(md[[celltype_col]])
  keep_cells <- unlist(lapply(split(rownames(md), md[[celltype_col]]), function(cells){
    if (length(cells) <= cap) return(cells)
    sample(cells, cap)
  }), use.names = FALSE)
  subset(seu, cells = keep_cells)
}

# get matrix from Seurat (SCT@data is log1p normalized)
get_mat <- function(seu, assay = assay_use, slot = slot_use) {
  m <- tryCatch(GetAssayData(seu, assay = assay, layer = slot),
                error = function(e) GetAssayData(seu, assay = assay, slot = slot))
  if (!inherits(m, "dgCMatrix")) m <- as(m, "dgCMatrix")
  m
}

# derive pseudo-counts from log1p data: counts ~= expm1(data)
# clamp small numerical negatives to 0 and round to integers
logdata_to_pseudocounts <- function(logmat) {
  if (!inherits(logmat, "dgCMatrix")) logmat <- as(logmat, "dgCMatrix")

  x <- logmat
  # operate only on non-zero entries (keeps it sparse)
  x@x <- expm1(x@x)          # base expm1 on numeric vector
  x@x[x@x < 0] <- 0
  x@x <- as.numeric(round(x@x))
  x
}

# Write CPDB input files:
# - meta.tsv: columns "Cell" and "cell_type"
# - counts.tsv: genes as rows, cells as columns (tab-delimited), integer counts
write_cpdb_inputs <- function(expr_counts, meta_df, celltype_col, outdir) {
  meta <- meta_df %>%
    transmute(Cell = cell, cell_type = .data[[celltype_col]]) %>%
    mutate(cell_type = gsub("[^A-Za-z0-9_]+", "_", cell_type))
  write_tsv(meta, file.path(outdir, "meta.tsv"))

  counts <- expr_counts[, meta$Cell, drop = FALSE]

  df <- as.data.frame(as.matrix(counts))  # OK only if you downsample enough
  df <- tibble::rownames_to_column(df, "gene")
  write_tsv(df, file.path(outdir, "counts.tsv"))
}

# Run CellPhoneDB statistical analysis (requires CPDB installed in environment)
cpdb_python <- "/home/ebarrozo/miniconda3/envs/r-reticulate/envs/cpdb/bin/python3"
cpdb_runner <- "/home/ebarrozo/visiumhd_mnp/results/audra/cpdb_run_statistical.py"

cpdb_db <- "/full/path/to/cellphonedb.zip"  # keep as full path
run_cpdb <- function(input_dir, outdir, py = "python3", runner_py,
                     db_path, threads, iterations, threshold_pct) {
  safe_dir(outdir)

  cmd <- sprintf(
    "%s %s --meta %s --counts %s --outdir %s --db %s --iterations %d --threshold %s --threads %d",
    shQuote(py),
    shQuote(runner_py),
    shQuote(file.path(input_dir, "meta.tsv")),
    shQuote(file.path(input_dir, "counts.tsv")),
    shQuote(outdir),
    shQuote(db_path),
    iterations,
    threshold_pct,
    threads
  )
  message("[CPDB-PY] ", cmd)
  status <- system(cmd)
  if (!identical(status, 0L)) stop("CPDB python runner failed with exit status ", status)
}

# Postprocess CPDB outputs to a compact long table:
# expects CPDB outputs: pvalues.txt and means.txt (or significant_means.txt depending version)
read_cpdb_outputs <- function(outdir) {
  # CPDB v5 writes timestamped files like:
  # statistical_analysis_means_*.txt, statistical_analysis_pvalues_*.txt, statistical_analysis_significant_means_*.txt

  pick_latest <- function(pattern) {
    f <- list.files(outdir, pattern = pattern, full.names = TRUE)
    if (!length(f)) return(NA_character_)
    f[which.max(file.info(f)$mtime)]
  }

  means_file <- pick_latest("^statistical_analysis_means_.*\\.txt$")
  pvals_file <- pick_latest("^statistical_analysis_pvalues_.*\\.txt$")

  # fallback: significant_means (some workflows only keep these)
  if (is.na(means_file)) {
    means_file <- pick_latest("^statistical_analysis_significant_means_.*\\.txt$")
  }

  if (is.na(means_file) || is.na(pvals_file)) {
    stop("Could not find CPDB outputs in: ", outdir,
         "\nFound files:\n", paste(basename(list.files(outdir)), collapse = "\n"))
  }

  means <- readr::read_tsv(means_file, show_col_types = FALSE)
  pvals <- readr::read_tsv(pvals_file, show_col_types = FALSE)

  list(means = means, pvals = pvals, means_file = means_file, pvals_file = pvals_file)
}

# Convert CPDB wide interaction table to long per pair
# CPDB pair columns look like "CellTypeA|CellTypeB"
wide_to_long_pairs <- function(means, pvals) {

  # pair columns in CPDB are typically "CellTypeA|CellTypeB"
  pair_cols <- grep("\\|", colnames(means), value = TRUE)

  if (!length(pair_cols)) {
    stop("No CPDB pair columns with '|' found. First columns are:\n",
         paste(colnames(means)[1:min(30, ncol(means))], collapse = ", "))
  }

  core_cols <- setdiff(colnames(means), pair_cols)

  m_long <- means %>%
    dplyr::select(dplyr::all_of(core_cols), dplyr::all_of(pair_cols)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(pair_cols), names_to = "pair", values_to = "mean")

  # pvals should have the same pair columns; if not, intersect
  pair_cols_p <- grep("\\|", colnames(pvals), value = TRUE)
  core_cols_p <- setdiff(colnames(pvals), pair_cols_p)

  p_long <- pvals %>%
    dplyr::select(dplyr::all_of(intersect(core_cols, core_cols_p)), dplyr::all_of(intersect(pair_cols, pair_cols_p))) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(intersect(pair_cols, pair_cols_p)), names_to = "pair", values_to = "pval")

  out <- m_long %>%
    dplyr::left_join(p_long, by = c(intersect(core_cols, core_cols_p), "pair")) %>%
    dplyr::mutate(
      pair = as.character(pair),
      sender = stringr::str_split_fixed(pair, "\\|", 2)[, 1],
      receiver = stringr::str_split_fixed(pair, "\\|", 2)[, 2]
    )

  out
}

# -------------------- PREPARE SEURAT SUBSETS --------------------
Idents(obj) <- celltype_col

obj_A <- subset_by_filters_safe(obj, condition_A_filter)
obj_B <- subset_by_filters_safe(obj, condition_B_filter)

obj_A <- min_cells_filter(obj_A, celltype_col, unique(c(senders, receivers)), min_cells_per_group)
obj_B <- min_cells_filter(obj_B, celltype_col, unique(c(senders, receivers)), min_cells_per_group)

# downsample for CPDB feasibility
obj_A <- downsample_by_celltype(obj_A, celltype_col, cap_cells_per_celltype)
obj_B <- downsample_by_celltype(obj_B, celltype_col, cap_cells_per_celltype)

# extract expression (log1p) and make pseudo-counts
obj_A <- min_cells_filter(obj_A, celltype_col, unique(c(senders, receivers)), min_cells_per_group)
obj_B <- min_cells_filter(obj_B, celltype_col, unique(c(senders, receivers)), min_cells_per_group)

# downsample for CPDB feasibility (DO THIS HERE)
obj_A <- downsample_by_celltype(obj_A, celltype_col, cap_cells_per_celltype)
obj_B <- downsample_by_celltype(obj_B, celltype_col, cap_cells_per_celltype)

# now extract expression
expr_A_log <- get_mat(obj_A)
expr_B_log <- get_mat(obj_B)

# intersect genes to keep consistent universe
genes <- intersect(rownames(expr_A_log), rownames(expr_B_log))
expr_A_log <- expr_A_log[genes, , drop = FALSE]
expr_B_log <- expr_B_log[genes, , drop = FALSE]

expr_A_counts <- logdata_to_pseudocounts(expr_A_log)
expr_B_counts <- logdata_to_pseudocounts(expr_B_log)

meta_A <- obj_A[[]] %>% mutate(cell = rownames(.))
meta_B <- obj_B[[]] %>% mutate(cell = rownames(.))
suppressPackageStartupMessages({
  library(tidyr)
})
# -------------------- RUN CPDB PER CONDITION --------------------
out_main <- safe_dir(file.path(out_root, paste0(delivery_use, "_", condition_B_name, "_vs_", condition_A_name)))

run_one_condition <- function(expr_counts, meta_df, condition_name) {
  od_in  <- safe_dir(file.path(out_main, paste0("cpdb_input_", condition_name)))
  od_out <- safe_dir(file.path(out_main, paste0("cpdb_out_", condition_name)))

  write_cpdb_inputs(expr_counts, meta_df, celltype_col, od_in)
  run_cpdb(od_in, od_out, db_path = cpdb_db, threads = threads,
           iterations = iterations, threshold_pct = threshold_pct)

  res <- read_cpdb_outputs(od_out)
  long <- wide_to_long_pairs(res$means, res$pvals) %>%
    mutate(Condition = condition_name)

  write_csv(long, file.path(out_main, paste0("cpdb_pairs_long_", condition_name, ".csv")))
  invisible(list(long = long, outdir = od_out))
}

res_A <- run_one_condition(expr_A_counts, meta_A, condition_A_name)
res_B <- run_one_condition(expr_B_counts, meta_B, condition_B_name)

pairs_A <- res_A$long
pairs_B <- res_B$long

# -------------------- DIFF SUMMARY (PTL - Control) --------------------
# merge A/B by interaction + sender/receiver pair and compute delta mean + delta significance
key_cols <- c("interacting_pair","gene_a","gene_b","partner_a","partner_b","annotation","classification","pair","sender","receiver")

merged <- pairs_B %>%
  select(all_of(key_cols), mean_B = mean, pval_B = pval) %>%
  full_join(
    pairs_A %>% select(all_of(key_cols), mean_A = mean, pval_A = pval),
    by = key_cols
  ) %>%
  mutate(mean_A = ifelse(is.na(mean_A), 0, mean_A),
         mean_B = ifelse(is.na(mean_B), 0, mean_B),
         pval_A = ifelse(is.na(pval_A), 1, pval_A),
         pval_B = ifelse(is.na(pval_B), 1, pval_B),
         delta_mean = mean_B - mean_A,
         sig_A = pval_A < 0.05,
         sig_B = pval_B < 0.05,
         sig_class = case_when(
           sig_B & !sig_A ~ "PTL_only",
           sig_A & !sig_B ~ "Control_only",
           sig_A & sig_B  ~ "Both",
           TRUE ~ "Neither"
         )) %>%
  arrange(desc(delta_mean))

write_csv(merged, file.path(out_main, "cpdb_pairs_diff_PTL_minus_Control.csv"))

# -------------------- PER RECEIVER OUTPUTS (mirror NicheNet) --------------------
out_per_receiver <- safe_dir(file.path(out_main, "per_receiver"))

for (recv in receivers) {
  od <- safe_dir(file.path(out_per_receiver, paste0("receiver_", gsub("[^A-Za-z0-9]+","_", recv))))
  tab <- merged %>%
    filter(receiver == recv, sender %in% senders) %>%
    arrange(desc(delta_mean))
  write_csv(tab, file.path(od, "cpdb_pairs_diff.csv"))

  # Top interactions in PTL for that receiver (rank by mean_B with pval_B < 0.05)
  top_ptl <- merged %>%
    filter(receiver == recv, sender %in% senders, pval_B < 0.05) %>%
    arrange(desc(mean_B)) %>%
    slice_head(n = 200)
  write_csv(top_ptl, file.path(od, "top_interactions_PTL.csv"))

  # Top interactions in Control
  top_ctl <- merged %>%
    filter(receiver == recv, sender %in% senders, pval_A < 0.05) %>%
    arrange(desc(mean_A)) %>%
    slice_head(n = 200)
  write_csv(top_ctl, file.path(od, "top_interactions_Control.csv"))
}

saveRDS(list(
  pairs_A = pairs_A,
  pairs_B = pairs_B,
  merged_diff = merged,
  params = list(
    condition_A_filter = condition_A_filter,
    condition_B_filter = condition_B_filter,
    senders = senders, receivers = receivers,
    min_cells_per_group = min_cells_per_group,
    cap_cells_per_celltype = cap_cells_per_celltype,
    iterations = iterations, threshold_pct = threshold_pct
  )
), file.path(out_main, "cellphonedb_results.rds"))

message("DONE. CellPhoneDB outputs under: ", normalizePath(out_main))



























suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr)
})

pick_latest <- function(outdir, pattern) {
  f <- list.files(outdir, pattern = pattern, full.names = TRUE)
  if (!length(f)) return(NA_character_)
  f[which.max(file.info(f)$mtime)]
}

read_cpdb_v5 <- function(outdir) {
  means_file <- pick_latest(outdir, "^statistical_analysis_means_.*\\.txt$")
  pvals_file <- pick_latest(outdir, "^statistical_analysis_pvalues_.*\\.txt$")
  decon_file <- pick_latest(outdir, "^statistical_analysis_deconvoluted_.*\\.txt$")
  deconp_file <- pick_latest(outdir, "^statistical_analysis_deconvoluted_percents_.*\\.txt$")
  sigmeans_file <- pick_latest(outdir, "^statistical_analysis_significant_means_.*\\.txt$")

  stopifnot(!is.na(means_file), !is.na(pvals_file))

  list(
    means = read_tsv(means_file, show_col_types = FALSE),
    pvals  = read_tsv(pvals_file,  show_col_types = FALSE),
    decon  = if (!is.na(decon_file)) read_tsv(decon_file, show_col_types = FALSE) else NULL,
    decon_perc = if (!is.na(deconp_file)) read_tsv(deconp_file, show_col_types = FALSE) else NULL,
    sigmeans = if (!is.na(sigmeans_file)) read_tsv(sigmeans_file, show_col_types = FALSE) else NULL,
    files = c(means=means_file, pvals=pvals_file, decon=decon_file, decon_perc=deconp_file, sigmeans=sigmeans_file)
  )
}

wide_to_long_pairs_v5 <- function(means, pvals) {
  pair_cols <- grep("\\|", colnames(means), value = TRUE)
  core_cols <- setdiff(colnames(means), pair_cols)

  m_long <- means %>%
    select(all_of(core_cols), all_of(pair_cols)) %>%
    pivot_longer(all_of(pair_cols), names_to="pair", values_to="mean")

  pair_cols_p <- grep("\\|", colnames(pvals), value = TRUE)
  core_cols_p <- setdiff(colnames(pvals), pair_cols_p)

  p_long <- pvals %>%
    select(all_of(intersect(core_cols, core_cols_p)), all_of(intersect(pair_cols, pair_cols_p))) %>%
    pivot_longer(all_of(intersect(pair_cols, pair_cols_p)), names_to="pair", values_to="pval")

  out <- m_long %>%
    left_join(p_long, by=c(intersect(core_cols, core_cols_p), "pair")) %>%
    mutate(
      sender = str_split_fixed(pair, "\\|", 2)[,1],
      receiver = str_split_fixed(pair, "\\|", 2)[,2]
    )
  out
}

# paths (your run directories)
od_A <- file.path(out_main, paste0("cpdb_out_", condition_A_name))
od_B <- file.path(out_main, paste0("cpdb_out_", condition_B_name))

cpA <- read_cpdb_v5(od_A)
cpB <- read_cpdb_v5(od_B)

pairs_A <- wide_to_long_pairs_v5(cpA$means, cpA$pvals) %>% mutate(Condition = condition_A_name)
pairs_B <- wide_to_long_pairs_v5(cpB$means, cpB$pvals) %>% mutate(Condition = condition_B_name)

# merge for PTL-Control delta
key_cols <- c("id_cp_interaction","interacting_pair","partner_a","partner_b","gene_a","gene_b",
              "annotation_strategy","is_integrin","directionality","classification","pair","sender","receiver")

merged <- pairs_B %>%
  select(any_of(key_cols), mean_B = mean, pval_B = pval) %>%
  full_join(pairs_A %>% select(any_of(key_cols), mean_A = mean, pval_A = pval),
            by = intersect(key_cols, colnames(pairs_B))) %>%
  mutate(
    mean_A = ifelse(is.na(mean_A), 0, mean_A),
    mean_B = ifelse(is.na(mean_B), 0, mean_B),
    pval_A = ifelse(is.na(pval_A), 1, pval_A),
    pval_B = ifelse(is.na(pval_B), 1, pval_B),
    delta_mean = mean_B - mean_A,
    sig_A = pval_A < 0.05,
    sig_B = pval_B < 0.05
  )

# col_map: named by Seurat levels with spaces, e.g. "NK Cell"
# make a cpdb-compatible map with underscores
col_map_cpdb <- setNames(
  col_map,
  gsub("[^A-Za-z0-9_]+", "_", names(col_map))
)

# quick check
stopifnot(all(unique(pairs_A$sender) %in% names(col_map_cpdb)))





library(ggplot2)

plot_receiver_heatmap <- function(merged, receiver_cpdb, top_n = 40, only_sig_in_B = TRUE) {
  df <- merged %>%
    filter(receiver == receiver_cpdb) %>%
    { if (only_sig_in_B) filter(., pval_B < 0.05) else . } %>%
    arrange(desc(delta_mean)) %>%
    slice_head(n = top_n) %>%
    mutate(
      interaction = paste0(interacting_pair),
      sender = factor(sender, levels = names(col_map_cpdb)),
      interaction = factor(interaction, levels = rev(unique(interaction)))
    )

  ggplot(df, aes(x = sender, y = interaction, fill = delta_mean)) +
    geom_tile() +
    geom_point(aes(size = -log10(pval_B)), alpha = 0.7) +
    scale_size_continuous(range = c(0.5, 4)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Sender", y = "Interaction", fill = "PTL - Control", size = "-log10(pval PTL)")
}

p1 <- plot_receiver_heatmap(merged, receiver_cpdb="Extravillous_Trophoblast", top_n=50)
ggsave(file.path(out_main, "viz_EVT_delta_heatmap_top50.png"), p1, width=10, height=12, dpi=300, bg="white")





plot_receiver_edges <- function(merged, receiver_cpdb, top_n = 30) {
  df <- merged %>%
    filter(receiver == receiver_cpdb, pval_B < 0.05) %>%
    arrange(desc(delta_mean)) %>%
    slice_head(n = top_n) %>%
    mutate(sender = factor(sender, levels = names(col_map_cpdb)))

  ggplot(df, aes(x = sender, y = delta_mean)) +
    geom_col() +
    coord_flip() +
    theme_classic() +
    labs(x = "Sender", y = "PTL - Control mean", title = paste0(receiver_cpdb, ": top ", top_n, " Δ interactions (sig in PTL)"))
}

p2 <- plot_receiver_edges(merged, "Extravillous_Trophoblast", 25)
ggsave(file.path(out_main, "viz_EVT_top25_delta_bar.png"), p2, width=8, height=6, dpi=300, bg="white")





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





suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(ktplots)   # devtools::install_github("zktuong/ktplots", dependencies=TRUE)
})

# Make sure the celltype labels match CPDB columns (underscores)
comb_ds2$celltype_cpdb <- gsub("[^A-Za-z0-9_]+", "_", comb_ds2[[celltype_col]][,1])

sce <- as.SingleCellExperiment(comb_ds2, assay = "RNA")  # you created RNA assay already
colData(sce)$celltype_cpdb <- comb_ds2$celltype_cpdb

means_A <- cpA$means
pvals_A <- cpA$pvals
decon_A <- cpA$decon

means_A_df <- as.data.frame(means_A)
pvals_A_df <- as.data.frame(pvals_A)
decon_A_df <- as.data.frame(decon_A)

p_dot <- plot_cpdb2(
  scdata = sce,
  cell_type1 = "Extravillous_Trophoblast",
  cell_type2 = "Macrophage",
  celltype_key = "celltype_cpdb",
  means = means_A_df,
  pvals = pvals_A_df,
  deconvoluted = decon_A_df,
  keep_significant_only = TRUE,
  standard_scale = TRUE,
  frac = 0.10,
  remove_self = TRUE,
  node_group_colors = col_map_cpdb
)

ggsave(file.path(out_main, "ktplots_dot_EVT_to_Macrophage_Control.png"),
       p_dot, width=10, height=6, dpi=300, bg="white")


p_hm <- plot_cpdb_heatmap(
  scdata = sce,
  celltype_key = "celltype_cpdb",
  means = means_A,
  pvals = pvals_A
)
ggsave(file.path(out_main, "ktplots_heatmap_Control.png"),
       p_hm, width=12, height=10, dpi=300, bg="white")

# choose a small, interpretable set (example: Notch + PGF)
keep_interactions <- c("JAG1_NOTCH2", "JAG1_NOTCH3", "PGF_FLT1", "PGF_NRP1")

means_B <- cpB$means
pvals_B <- cpB$pvals
decon_B <- cpB$decon

means_B_df <- as.data.frame(means_B)
pvals_B_df <- as.data.frame(pvals_B)
decon_B_df <- as.data.frame(decon_B)

p_chord <- plot_cpdb4(
  scdata = sce,
  cell_type1 = "Extravillous_Trophoblast",
  cell_type2 = "Macrophage",
  celltype_key = "celltype_cpdb",
  means = means_B_df,            # PTL
  pvals = pvals_B_df,
  deconvoluted = decon_B_df,
  interaction = keep_interactions,
  keep_significant_only = TRUE,
  grid_colors = col_map_cpdb,
  alpha = 0.6
)

ggsave(file.path(out_main, "ktplots_chord_EVT_Macrophage_PTL_select.png"),
       p_chord, width=10, height=8, dpi=300, bg="white")

decp <- cpA$decon_perc
# Example: plot percent expression for a handful of genes
genes_oi <- c("PGF","FLT1","JAG1","NOTCH2","NOTCH3")

dfp <- decp %>%
  filter(gene_name %in% genes_oi) %>%
  pivot:contentReference[oaicite:10]{index=10}, uniprot, is_complex, protein_name, complex_name, id_cp_interaction, gene),
               names_to = "celltype", values_to = "pct") %>%
  mutate(celltype = factor(celltype, levels = names(col_map_cpdb)))

p_qc <- ggplot(dfp, aes(x = celltype, y = pct, fill = gene_name)) +
  geom_col(position = "dodge") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x=NULL, y="Percent expressing (deconvoluted_percents)", fill="Gene")

ggsave(file.path(out_main, "qc_deconvoluted_percents_genes_oi.png"),
       p_qc, width=10, height=5, dpi=300, bg="white")






















suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

# robust: only columns that look like CPDB pair columns
pair_cols <- function(df) grep("\\|", colnames(df), value = TRUE)

# return top interactions for a specific directed pair (e.g., "A|B")
top_interactions_for_pair <- function(means, pvals, pair,
                                      p_cut = 0.05, mean_min = 0.0,
                                      n = 20) {

  stopifnot(pair %in% colnames(means), pair %in% colnames(pvals))

  out <- means %>%
    transmute(
      interacting_pair = interacting_pair,
      gene_a = gene_a, gene_b = gene_b,
      classification = classification,
      mean = .data[[pair]]
    ) %>%
    left_join(
      pvals %>% transmute(interacting_pair = interacting_pair, pval = .data[[pair]]),
      by = "interacting_pair"
    ) %>%
    filter(!is.na(mean), !is.na(pval), pval < p_cut, mean > mean_min) %>%
    arrange(desc(mean)) %>%
    slice_head(n = n)

  unique(out$interacting_pair)
}

# quick “does this pair have anything to plot?”
n_sig_for_pair <- function(means, pvals, pair, p_cut = 0.05) {
  sum(!is.na(means[[pair]]) & !is.na(pvals[[pair]]) & pvals[[pair]] < p_cut & means[[pair]] > 0)
}

keep_interactions <- keep_interactions_B  # e.g., PTL top interactions for pair1

if (length(keep_interactions) < 1) {
  message("No significant interactions for this pair/condition. Try another pair or lower strictness.")
} else {
  p_chord <- plot_cpdb4(
    scdata = sce,
    cell_type1 = "Extravillous_Trophoblast",
    cell_type2 = "Macrophage",
    celltype_key = "celltype_cpdb",
    means = means_B_df,
    pvals = pvals_B_df,
    deconvoluted = decon_B_df,
    interaction = keep_interactions,
    keep_significant_only = TRUE,
    grid_colors = col_map_cpdb,
    alpha = 0.6
  )
}



pairs <- pair_cols(pvals_B_df)

sig_counts <- sapply(pairs, function(pr) n_sig_for_pair(means_B_df, pvals_B_df, pr))
sort(sig_counts, decreasing = TRUE)[1:15]




pair_ok <- "Villous_Cytotrophoblast|T_Cell"

# auto-pick top interactions for this pair in PTL
keep_interactions <- top_interactions_for_pair(means_B_df, pvals_B_df, pair_ok, n = 25)

# dotplot (PTL)
p_dot_vct_t <- plot_cpdb2(
  scdata = sce,
  cell_type1 = "Villous_Cytotrophoblast",
  cell_type2 = "T_Cell",
  celltype_key = "celltype_cpdb",
  means = as.data.frame(means_B_df),
  pvals = as.data.frame(pvals_B_df),
  deconvoluted = as.data.frame(decon_B_df),
  keep_significant_only = TRUE,
  standard_scale = TRUE,
  frac = 0.10,
  remove_self = TRUE,
  node_group_colors = col_map_cpdb
)

ggsave(file.path(out_main, "ktplots_dot_VCT_TCell_PTL.png"),
       p_dot_vct_t, width = 11, height = 6, dpi = 300, bg = "white")

pair1 <- "Macrophage|Villous_Cytotrophoblast"
pair2 <- "Villous_Cytotrophoblast|Macrophage"

c(n_sig_for_pair(means_B_df, pvals_B_df, pair1),
  n_sig_for_pair(means_B_df, pvals_B_df, pair2))

pair_ok <- "Macrophage|Villous_Cytotrophoblast"  # pick the better one
keep_interactions <- top_interactions_for_pair(means_B_df, pvals_B_df, pair_ok, n = 12)

p_chord <- plot_cpdb4(
  scdata = sce,
  cell_type1 = "Macrophage",
  cell_type2 = "Villous_Cytotrophoblast",
  celltype_key = "celltype_cpdb",
  means = as.data.frame(means_B_df),
  pvals = as.data.frame(pvals_B_df),
  deconvoluted = as.data.frame(decon_B_df),
  interaction = keep_interactions,
  keep_significant_only = TRUE,
  grid_colors = col_map_cpdb,
  alpha = 0.6
)

ggsave(file.path(out_main, "ktplots_chord_Mac_VCT_PTL_top12.png"),
       p_chord, width = 10, height = 8, dpi = 300, bg = "white")
















library(ggplot2)

sig_matrix <- function(means, pvals, p_cut = 0.05) {
  pcs <- grep("\\|", colnames(means), value = TRUE)
  tibble(pair = pcs) %>%
    mutate(
      n_sig = sapply(pair, \(pr) n_sig_for_pair(means, pvals, pr, p_cut)),
      sender = str_split_fixed(pair, "\\|", 2)[,1],
      receiver = str_split_fixed(pair, "\\|", 2)[,2]
    )
}

smB <- sig_matrix(means_B_df, pvals_B_df)  # PTL
p_sig <- ggplot(smB, aes(x = sender, y = receiver, fill = n_sig)) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title="PTL: # significant interactions per pair", x=NULL, y=NULL, fill="# sig")

ggsave(file.path(out_main, "qc_sigcount_matrix_PTL.png"),
       p_sig, width=9, height=7, dpi=300, bg="white")