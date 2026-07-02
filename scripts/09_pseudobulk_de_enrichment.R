# ============================================================
# Fig. 3 pseudobulk DE and pathway enrichment
# Purpose:
#   Reproduce subject-aware pseudobulk high-vs-low MNP differential
#   expression and enrichment source tables from the public Seurat object.
#
# Inputs:
#   obj_public loaded by scripts/02_load_public_object.R
#
# Outputs:
#   outputs/tables/pseudobulk_de/
#   outputs/figures/pseudobulk_de/
#
# Notes:
#   This script consolidates the manuscript pseudobulk/clusterProfiler
#   workflows into a cleaner public-facing workflow. Original exploratory
#   scripts are archived in provenance/archived_scripts/.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(purrr)
  library(ggplot2)
})

required_pkgs <- c("limma", "matrixStats", "clusterProfiler", "org.Hs.eg.db", "msigdbr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  dir.create("outputs/tables/pseudobulk_de", recursive = TRUE, showWarnings = FALSE)
  write_csv(
    tibble(missing_required_package = missing_pkgs),
    "outputs/tables/pseudobulk_de/pseudobulk_de_missing_packages.csv"
  )
  message("Skipping pseudobulk DE/enrichment because packages are missing: ", paste(missing_pkgs, collapse = ", "))
  stop("Install missing Bioconductor/CRAN packages before running pseudobulk DE/enrichment.")
}

suppressPackageStartupMessages({
  library(limma)
  library(matrixStats)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
})

if (!exists("obj_public")) {
  if (exists("obj")) obj_public <- obj else stop("Expected obj_public from scripts/02_load_public_object.R")
}

obj <- obj_public

out_root <- "outputs/tables/pseudobulk_de"
fig_root <- "outputs/figures/pseudobulk_de"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_root, recursive = TRUE, showWarnings = FALSE)

# ------------------------- config -------------------------
annotation_col <- Sys.getenv("PSEUDOBULK_ANNOTATION_COL", unset = "predicted.celltype")
if (!annotation_col %in% colnames(obj@meta.data) && "nn.predicted.celltype" %in% colnames(obj@meta.data)) {
  annotation_col <- "nn.predicted.celltype"
}

q_cutoff_de <- as.numeric(Sys.getenv("PSEUDOBULK_Q_CUTOFF", unset = "0.10"))
l2fc_min <- as.numeric(Sys.getenv("PSEUDOBULK_L2FC_MIN", unset = "0.25"))
min_samples_per_group <- as.integer(Sys.getenv("PSEUDOBULK_MIN_SAMPLES", unset = "2"))
var_eps <- as.numeric(Sys.getenv("PSEUDOBULK_VAR_EPS", unset = "1e-8"))
enr_min_sig_genes <- as.integer(Sys.getenv("PSEUDOBULK_ENR_MIN_SIG_GENES", unset = "10"))
padj_cut_family <- as.numeric(Sys.getenv("PSEUDOBULK_PATHWAY_FAMILY_FDR", unset = "0.05"))

Y_GENES <- c("DDX3Y","KDM5D","RPS4Y1","EIF1AY","UTY","USP9Y","ZFY","TXLNGY","TMSB4Y",
             "PRKY","PCDH11Y","NLGN4Y","AMELY")
lnc_pattern <- "^RP[0-9]+-|^AC[0-9]+\\.|^AL[0-9]+\\.|^LINC|^MIR|^SNHG|^HCG|^CT[CDE]-|^ENSG"
PATHOGEN_GENES <- c(
  "CMV-UL117","CMV-UL122","CMV-UL123","CMV-UL138","CMV-UL5","CMV-UL54","CMV-US33",
  "GBS-fabD","GBS-fbsB","GBS-hylB","GBS-lmb","GBS-mfd","GBS-pbsP","GBS-scpB",
  "HPV-E1","HPV-E2","HPV-E5","HPV-E6","HPV-L1","HPV-L2",
  "HSV-1-LAT","HSV-1-RL1","HSV-1-RL2","HSV-1-RS1","HSV-1-UL14","HSV-1-UL29","HSV-1-UL48","HSV-1-UL54","HSV-1-US1",
  "HSV-2-LAT","HSV-2-RS1","HSV-2-UL30","HSV-2-UL41","HSV-2-UL52","HSV-2-US1","HSV-2-US3",
  "SARS-CoV-2-M","SARS-CoV-2-N","SARS-CoV-2-ORF10","SARS-CoV-2-ORF1ab","SARS-CoV-2-ORF3a","SARS-CoV-2-ORF8","SARS-CoV-2-S",
  "Syphilis-dnaK","Syphilis-rplC","Syphilis-rplF","Syphilis-rpsE","Syphilis-rpsM","Syphilis-secY","Syphilis-tuf"
)

vecify_chr <- function(x) {
  if (is.list(x)) vapply(x, function(y) paste(as.character(y), collapse = ";"), character(1)) else as.character(x)
}

sanitize_filename <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)

safe_dup_cor <- function(expr_mat, design, block) {
  if (is.null(block)) return(0)
  tb <- table(block)
  if (all(tb <= 1)) return(0)
  res <- try(limma::duplicateCorrelation(expr_mat, design, block = block), silent = TRUE)
  if (inherits(res, "try-error")) return(0)
  res$consensus
}

run_go_bp <- function(genes_sym, universe_sym) {
  genes_sym <- unique(setdiff(genes_sym, Y_GENES))
  if (length(genes_sym) < enr_min_sig_genes) return(NULL)
  mm <- clusterProfiler::bitr(genes_sym, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) %>%
    distinct(SYMBOL, .keep_all = TRUE)
  mm_uni <- clusterProfiler::bitr(universe_sym, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) %>%
    distinct(SYMBOL, .keep_all = TRUE)
  if (nrow(mm) < enr_min_sig_genes || nrow(mm_uni) < 100) return(NULL)
  ego <- clusterProfiler::enrichGO(
    gene = mm$ENTREZID,
    universe = mm_uni$ENTREZID,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.10,
    readable = TRUE
  )
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)
  suppressWarnings(clusterProfiler::simplify(ego, cutoff = 0.6, by = "p.adjust", select_fun = min))
}

run_hallmark <- function(genes_sym, universe_sym) {
  genes_sym <- unique(setdiff(genes_sym, Y_GENES))
  if (length(genes_sym) < enr_min_sig_genes) return(NULL)
  mm <- clusterProfiler::bitr(genes_sym, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) %>%
    distinct(SYMBOL, .keep_all = TRUE)
  mm_uni <- clusterProfiler::bitr(universe_sym, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) %>%
    distinct(SYMBOL, .keep_all = TRUE)
  if (nrow(mm) < enr_min_sig_genes || nrow(mm_uni) < 100) return(NULL)
  H <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>%
    select(gs_name, entrez_gene) %>%
    distinct()
  enr <- clusterProfiler::enricher(mm$ENTREZID, TERM2GENE = H, universe = mm_uni$ENTREZID, pAdjustMethod = "BH")
  if (is.null(enr) || nrow(as.data.frame(enr)) == 0) return(NULL)
  enr
}

# ------------------------- object prep -------------------------
needed <- c(annotation_col, "MNPBurden", "Delivery_class", "study_id")
missing_cols <- setdiff(needed, colnames(obj@meta.data))
if (length(missing_cols) > 0) stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))

obj@meta.data$pseudobulk_celltype <- vecify_chr(obj@meta.data[[annotation_col]])
obj@meta.data$MNPBurden <- sub("^(High|Low).*", "\\1", vecify_chr(obj@meta.data$MNPBurden))
obj@meta.data$Delivery_class <- vecify_chr(obj@meta.data$Delivery_class)
obj@meta.data$study_id <- vecify_chr(obj@meta.data$study_id)

assay_use <- if ("SCT" %in% names(obj@assays)) "SCT" else DefaultAssay(obj)
DefaultAssay(obj) <- assay_use

# ------------------------- delivery-level pseudobulk -------------------------
run_delivery_pseudobulk <- function(delivery_to_run) {
  message("[pseudobulk] Running ", delivery_to_run)
  outdir <- file.path(out_root, paste0("pb_", delivery_to_run, "_MNPBurden_limma"))
  enrdir <- file.path(outdir, "cp_results")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(enrdir, recursive = TRUE, showWarnings = FALSE)

  obj_use <- subset(obj, subset = Delivery_class == delivery_to_run & MNPBurden %in% c("Low", "High"))
  if (ncol(obj_use) == 0) stop("No cells for delivery class: ", delivery_to_run)

  bulk <- Seurat::AggregateExpression(
    obj_use,
    slot = "data",
    group.by = c("pseudobulk_celltype", "MNPBurden", "study_id"),
    return.seurat = TRUE
  )
  DefaultAssay(bulk) <- if (assay_use %in% names(bulk@assays)) assay_use else DefaultAssay(bulk)

  pb_meta <- tibble(pb_id = colnames(bulk)) %>%
    bind_cols(bulk@meta.data[colnames(bulk), c("pseudobulk_celltype", "MNPBurden", "study_id"), drop = FALSE] %>% as_tibble()) %>%
    mutate(
      pseudobulk_celltype = vecify_chr(pseudobulk_celltype),
      MNPBurden = sub("^(High|Low).*", "\\1", vecify_chr(MNPBurden)),
      study_id = vecify_chr(study_id)
    )

  cell_counts <- obj_use@meta.data %>%
    as_tibble(rownames = "cell") %>%
    transmute(
      pseudobulk_celltype = vecify_chr(pseudobulk_celltype),
      MNPBurden = sub("^(High|Low).*", "\\1", vecify_chr(MNPBurden)),
      study_id = vecify_chr(study_id)
    ) %>%
    count(pseudobulk_celltype, MNPBurden, study_id, name = "n_cells")

  pb_meta <- left_join(pb_meta, cell_counts, by = c("pseudobulk_celltype", "MNPBurden", "study_id"))
  pb_meta$n_cells[is.na(pb_meta$n_cells)] <- 1L

  pb_counts <- pb_meta %>%
    count(pseudobulk_celltype, MNPBurden, name = "n_pseudobulk") %>%
    arrange(pseudobulk_celltype, MNPBurden)
  write_tsv(pb_counts, file.path(outdir, "pseudobulk_sample_counts_by_celltype.tsv"))

  E <- as.matrix(GetAssayData(bulk, slot = "data"))
  keep_genes <- rownames(E)[!grepl(lnc_pattern, rownames(E), perl = TRUE) & !(rownames(E) %in% PATHOGEN_GENES)]
  E <- E[keep_genes, , drop = FALSE]
  stopifnot(ncol(E) == nrow(pb_meta))

  run_one_celltype <- function(ct) {
    sel <- which(pb_meta$pseudobulk_celltype == ct)
    if (length(sel) < 2 * min_samples_per_group) return(NULL)
    meta_ct <- pb_meta[sel, , drop = FALSE]
    E_ct <- E[, sel, drop = FALSE]
    grp_counts <- table(meta_ct$MNPBurden)
    if (!all(c("Low", "High") %in% names(grp_counts))) return(NULL)
    if (any(grp_counts[c("Low", "High")] < min_samples_per_group)) return(NULL)

    v <- matrixStats::rowVars(E_ct)
    keep <- which(is.finite(v) & v > var_eps)
    if (length(keep) < 50) return(NULL)
    E_ct <- E_ct[keep, , drop = FALSE]

    meta_ct$MNPBurden <- factor(meta_ct$MNPBurden, levels = c("Low", "High"))
    design <- model.matrix(~ 0 + MNPBurden, data = meta_ct)
    w <- sqrt(meta_ct$n_cells)
    corval <- safe_dup_cor(E_ct, design, block = meta_ct$study_id)

    fit0 <- limma::lmFit(
      E_ct,
      design,
      block = if (corval > 0) meta_ct$study_id else NULL,
      correlation = if (corval > 0) corval else NULL,
      weights = w
    )
    cont <- limma::makeContrasts(High_vs_Low = MNPBurdenHigh - MNPBurdenLow, levels = design)
    fit2 <- limma::eBayes(limma::contrasts.fit(fit0, cont), trend = TRUE)
    tt <- limma::topTable(fit2, number = Inf, sort.by = "P") %>%
      rownames_to_column("gene") %>%
      mutate(
        mean_High = rowMeans(E_ct[gene, meta_ct$MNPBurden == "High", drop = FALSE]),
        mean_Low = rowMeans(E_ct[gene, meta_ct$MNPBurden == "Low", drop = FALSE]),
        sig_q10 = adj.P.Val < q_cutoff_de & abs(logFC) >= l2fc_min,
        sig_q05 = adj.P.Val < 0.05 & abs(logFC) >= l2fc_min,
        direction = if_else(logFC > 0, "UP_in_High", "DOWN_in_High"),
        predicted.celltype = ct,
        Delivery_class = delivery_to_run,
        n_samples_High = as.integer(grp_counts["High"]),
        n_samples_Low = as.integer(grp_counts["Low"]),
        used_duplicate_correlation = corval > 0,
        duplicate_correlation = corval
      )

    fp <- file.path(outdir, sprintf("pb_DE_%s.tsv", sanitize_filename(ct)))
    write_tsv(tt, fp)

    tibble(
      Delivery_class = delivery_to_run,
      predicted.celltype = ct,
      n_samples_High = as.integer(grp_counts["High"]),
      n_samples_Low = as.integer(grp_counts["Low"]),
      n_sig_q10 = sum(tt$sig_q10, na.rm = TRUE),
      n_sig_q05 = sum(tt$sig_q05, na.rm = TRUE),
      n_up_q10 = sum(tt$sig_q10 & tt$logFC > 0, na.rm = TRUE),
      n_down_q10 = sum(tt$sig_q10 & tt$logFC < 0, na.rm = TRUE),
      used_duplicate_correlation = corval > 0
    )
  }

  cts <- sort(unique(pb_meta$pseudobulk_celltype))
  summary_tbl <- bind_rows(lapply(cts, run_one_celltype))
  write_tsv(summary_tbl, file.path(outdir, "pb_DE_summary_by_celltype.tsv"))

  all_de_files <- list.files(outdir, pattern = "^pb_DE_.*\\.tsv$", full.names = TRUE)
  all_de <- bind_rows(lapply(all_de_files, read_tsv, show_col_types = FALSE))
  if (nrow(all_de) > 0) write_tsv(all_de, file.path(outdir, "pb_DE_ALL.tsv"))

  # Enrichment
  cts_run <- summary_tbl %>%
    mutate(total_sig = coalesce(n_sig_q10, 0L)) %>%
    filter(total_sig >= enr_min_sig_genes) %>%
    pull(predicted.celltype)

  enr_all <- list()
  for (ct in cts_run) {
    de_fp <- file.path(outdir, sprintf("pb_DE_%s.tsv", sanitize_filename(ct)))
    de <- read_tsv(de_fp, show_col_types = FALSE)
    universe <- unique(de$gene)
    up <- de %>% filter(sig_q10, logFC > 0) %>% pull(gene) %>% setdiff(Y_GENES)
    down <- de %>% filter(sig_q10, logFC <= 0) %>% pull(gene) %>% setdiff(Y_GENES)

    analyses <- list(
      GO_BP_UP = run_go_bp(up, universe),
      GO_BP_DOWN = run_go_bp(down, universe),
      MSIG_HALLMARK_UP = run_hallmark(up, universe),
      MSIG_HALLMARK_DOWN = run_hallmark(down, universe)
    )

    for (nm in names(analyses)) {
      enr <- analyses[[nm]]
      if (is.null(enr)) next
      df <- as.data.frame(enr)
      if (!nrow(df)) next
      fp <- file.path(enrdir, sprintf("%s_%s.csv", nm, sanitize_filename(ct)))
      write_csv(df, fp)
      family <- ifelse(startsWith(nm, "GO_BP"), "GO_BP", "HALLMARK")
      direction <- ifelse(grepl("_UP$|_UP_", nm), "UP", "DOWN")
      enr_all[[paste(delivery_to_run, ct, nm, sep = "__")]] <- df %>%
        mutate(
          Delivery_class = delivery_to_run,
          predicted.celltype = ct,
          pathway_family = family,
          direction = direction,
          source_file = basename(fp)
        )
    }
  }

  enr_all <- bind_rows(enr_all)
  if (nrow(enr_all) > 0) {
    enr_all <- enr_all %>%
      group_by(Delivery_class, pathway_family, predicted.celltype, direction) %>%
      mutate(
        p_adj_family = p.adjust(pvalue, method = "BH"),
        sig_family = p_adj_family <= padj_cut_family
      ) %>%
      ungroup()
    write_tsv(enr_all, file.path(enrdir, "enrichment_all.familyFDR.q05.tsv"))
  }

  list(summary = summary_tbl, de = all_de, enrichment = enr_all)
}

results <- list(
  Term = run_delivery_pseudobulk("Term"),
  Preterm = run_delivery_pseudobulk("Preterm")
)

# ------------------------- combined source tables -------------------------
all_de <- bind_rows(lapply(results, `[[`, "de"))
all_summary <- bind_rows(lapply(results, `[[`, "summary"))
all_enr <- bind_rows(lapply(results, `[[`, "enrichment"))

write_tsv(all_summary, file.path(out_root, "Fig3_pseudobulk_DE_summary_by_celltype_all_delivery.tsv"))
write_tsv(all_de, file.path(out_root, "Fig3_pseudobulk_DE_ALL_delivery.tsv"))
if (nrow(all_enr) > 0) write_tsv(all_enr, file.path(out_root, "Fig3_pathway_enrichment_ALL_delivery.tsv"))

# Fig. 3a heatmap source: top DE genes by absolute effect among significant hits,
# plus fallback top genes so the source table is never empty.
heatmap_source <- all_de %>%
  filter(!is.na(logFC)) %>%
  group_by(Delivery_class, predicted.celltype) %>%
  arrange(adj.P.Val, desc(abs(logFC)), .by_group = TRUE) %>%
  slice_head(n = 25) %>%
  ungroup() %>%
  select(Delivery_class, predicted.celltype, gene, logFC, AveExpr, t, P.Value, adj.P.Val, sig_q10, direction,
         n_samples_High, n_samples_Low)
write_tsv(heatmap_source, file.path(out_root, "Fig3a_pseudobulk_logFC_heatmap_source.tsv"))

if (nrow(heatmap_source) > 0) {
  plot_df <- heatmap_source %>%
    mutate(
      panel = paste(Delivery_class, predicted.celltype, sep = " | "),
      gene = factor(gene, levels = rev(unique(gene[order(abs(logFC), decreasing = TRUE)])))
    )

  p_heat <- ggplot(plot_df, aes(x = panel, y = gene, fill = logFC)) +
    geom_tile(color = "white", linewidth = 0.1) +
    scale_fill_gradient2(name = "log2FC\nHigh vs Low") +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 8) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.ticks = element_blank(),
      axis.line = element_blank()
    )
  ggsave(file.path(fig_root, "Fig3a_pseudobulk_logFC_heatmap_draft.png"), p_heat, width = 10.5, height = 8.0, dpi = 450)
}

if (nrow(all_enr) > 0) {
  enr_plot_source <- all_enr %>%
    filter(!is.na(p_adj_family)) %>%
    group_by(Delivery_class, predicted.celltype, direction, pathway_family) %>%
    arrange(p_adj_family, .by_group = TRUE) %>%
    slice_head(n = 8) %>%
    ungroup() %>%
    mutate(
      neglog10_family_fdr = -log10(pmax(p_adj_family, .Machine$double.xmin)),
      Description = str_wrap(Description, width = 45)
    )
  write_tsv(enr_plot_source, file.path(out_root, "Fig3bc_pathway_enrichment_plot_source.tsv"))

  p_enr <- ggplot(enr_plot_source, aes(x = neglog10_family_fdr, y = reorder(Description, neglog10_family_fdr), size = Count, color = direction)) +
    geom_point(alpha = 0.9) +
    facet_grid(Delivery_class ~ predicted.celltype, scales = "free_y", space = "free_y") +
    labs(x = "-log10 family-wise FDR", y = NULL, color = "Direction", size = "Genes") +
    theme_classic(base_size = 8) +
    theme(strip.text = element_text(face = "bold"))
  ggsave(file.path(fig_root, "Fig3bc_pathway_enrichment_draft.png"), p_enr, width = 12, height = 8, dpi = 450)
}

message("Pseudobulk DE/enrichment workflow complete. Outputs written to: ", out_root)
