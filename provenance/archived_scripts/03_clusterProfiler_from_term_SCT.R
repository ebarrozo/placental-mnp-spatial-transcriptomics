# 03_clusterProflier_from_term_SCT.R
# ======================= TERM placentae pipeline → family-wise FDR 0.05 =======================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr); library(tidyr); library(purrr); library(stringr)
  library(readr); library(data.table)
  library(limma); library(matrixStats)
  library(clusterProfiler); library(enrichplot)
  library(org.Hs.eg.db); library(msigdbr)
  library(ggplot2)
})

# ------------------------- CONFIG -------------------------
DELIVERY_TO_RUN        <- "Term"   # <- key change vs PTB run
OUTDIR_DE              <- "pb_Term_MNPBurden_limma_SCT"
OUTDIR_ENR             <- file.path(OUTDIR_DE, "cp_results")
dir.create(OUTDIR_DE,  showWarnings = FALSE, recursive = TRUE)
dir.create(OUTDIR_ENR, showWarnings = FALSE, recursive = TRUE)

Q_CUTOFF_DE            <- 0.10     # DE FDR
L2FC_MIN               <- 0.25
MIN_SAMPLES_PER_GROUP  <- 2
VAR_EPS                <- 1e-8
ENR_MIN_SIG_GENES      <- 10       # run ORA only if ≥ this many sig genes (q<0.10)
PADJ_CUT_FAMILY        <- 0.05     # family-wise FDR threshold (requested)

# Filters
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

# --------------------- UTILITIES (robust) ---------------------
vecify_chr <- function(x) {
  if (is.list(x)) vapply(x, function(y) paste(as.character(y), collapse = ";"), character(1))
  else as.character(x)
}

safe_dup_cor <- function(expr_mat, design, block) {
  if (is.null(block)) return(0)
  tb <- table(block)
  if (all(tb <= 1)) return(0)
  res <- try(duplicateCorrelation(expr_mat, design, block = block), silent = TRUE)
  if (inherits(res, "try-error")) return(0)
  res$consensus
}

run_go_bp <- function(genes_sym, universe_sym) {
  if (length(genes_sym) < 10) return(NULL)
  mm     <- clusterProfiler::bitr(genes_sym,    fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db) %>% dplyr::distinct(SYMBOL, .keep_all = TRUE)
  mm_uni <- clusterProfiler::bitr(universe_sym, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db) %>% dplyr::distinct(SYMBOL, .keep_all = TRUE)
  if (nrow(mm) < 10 || nrow(mm_uni) < 100) return(NULL)
  ego <- enrichGO(
    gene          = mm$ENTREZID,
    universe      = mm_uni$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.10,
    readable      = TRUE
  )
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)
  suppressWarnings(simplify(ego, cutoff = 0.6, by = "p.adjust", select_fun = min))
}

run_hallmark <- function(genes_sym, universe_sym) {
  if (length(genes_sym) < 10) return(NULL)
  mm     <- clusterProfiler::bitr(genes_sym,    fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db) %>% dplyr::distinct(SYMBOL, .keep_all = TRUE)
  mm_uni <- clusterProfiler::bitr(universe_sym, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db) %>% dplyr::distinct(SYMBOL, .keep_all = TRUE)
  if (nrow(mm) < 10 || nrow(mm_uni) < 100) return(NULL)
  H <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)
  enricher(mm$ENTREZID, TERM2GENE = H, universe = mm_uni$ENTREZID, pAdjustMethod = "BH")
}

# --------------------- 1) SUBSET → TERM & PSEUDOBULK ---------------------
stopifnot(exists("obj"), !is.null(obj@meta.data))

md <- obj@meta.data
if (!"Delivery_class" %in% colnames(md)) stop("Delivery_class not found in obj@meta.data")
if (!"MNPBurden" %in% colnames(md) && "MNPBurden_Delivery_class" %in% colnames(md)) {
  md$MNPBurden <- sub("^(High|Low).*", "\\1", md$MNPBurden_Delivery_class)
}
obj@meta.data <- md

obj_term <- subset(obj, subset = Delivery_class == DELIVERY_TO_RUN)
if (ncol(obj_term) == 0) stop("No cells with Delivery_class == '", DELIVERY_TO_RUN, "'.")

DefaultAssay(obj_term) <- if ("SCT" %in% names(obj_term@assays)) "SCT" else DefaultAssay(obj_term)

bulk <- AggregateExpression(
  obj_term, slot = "data",
  group.by = c("predicted.celltype","MNPBurden","study_id"),
  return.seurat = TRUE
)
DefaultAssay(bulk) <- if ("SCT" %in% names(bulk@assays)) "SCT" else DefaultAssay(bulk)

# pseudobulk meta
pb_meta <- tibble(pb_id = colnames(bulk)) %>%
  dplyr::bind_cols(bulk@meta.data[colnames(bulk), c("predicted.celltype","MNPBurden","study_id"), drop = FALSE] %>% tibble::as_tibble()) %>%
  dplyr::mutate(
    predicted.celltype = vecify_chr(predicted.celltype),
    MNPBurden          = sub("^(High|Low).*", "\\1", vecify_chr(MNPBurden)),
    study_id           = vecify_chr(study_id)
  )

# contributing cell counts per pseudo-sample (weights)
cell_counts <- obj_term@meta.data %>%
  tibble::as_tibble(rownames = "cell") %>%
  dplyr::mutate(
    predicted.celltype = vecify_chr(.data$predicted.celltype),
    MNPBurden          = sub("^(High|Low).*", "\\1", vecify_chr(.data$MNPBurden)),
    study_id           = vecify_chr(.data$study_id)
  ) %>%
  dplyr::count(predicted.celltype, MNPBurden, study_id, name = "n_cells")

pb_meta <- dplyr::left_join(pb_meta, cell_counts,
                            by = c("predicted.celltype","MNPBurden","study_id"))

# write QC counts
pb_meta %>%
  dplyr::group_by(predicted.celltype, MNPBurden, .drop = FALSE) %>%
  dplyr::summarise(n_pseudobulk = dplyr::n(), .groups = "drop") %>%
  dplyr::arrange(predicted.celltype, MNPBurden) %>%
  readr::write_tsv(file.path(OUTDIR_DE, "pseudobulk_sample_counts_by_celltype.tsv"))

E <- as.matrix(GetAssayData(bulk, slot = "data"))
keep_genes <- rownames(E)[ !grepl(lnc_pattern, rownames(E), perl = TRUE) &
                           !(rownames(E) %in% PATHOGEN_GENES) ]
E <- E[keep_genes, , drop = FALSE]
stopifnot(ncol(E) == nrow(pb_meta))

# --------------------- 2) DE: limma-trend (High vs Low) ---------------------
run_one_celltype <- function(ct) {
  sel <- which(pb_meta$predicted.celltype == ct)
  if (length(sel) < 2*MIN_SAMPLES_PER_GROUP) return(NULL)

  meta_ct <- pb_meta[sel, , drop = FALSE]
  E_ct    <- E[, sel, drop = FALSE]

  grp_counts <- table(meta_ct$MNPBurden)
  if (!all(c("Low","High") %in% names(grp_counts))) return(NULL)
  if (any(grp_counts[c("Low","High")] < MIN_SAMPLES_PER_GROUP)) return(NULL)

  v <- rowVars(E_ct)
  keep <- which(is.finite(v) & v > VAR_EPS)
  if (length(keep) < 50) return(NULL)
  E_ct <- E_ct[keep, , drop = FALSE]

  meta_ct$MNPBurden <- factor(meta_ct$MNPBurden, levels = c("Low","High"))
  design <- model.matrix(~ 0 + MNPBurden, data = meta_ct)

  w <- meta_ct$n_cells; if (is.null(w) || all(is.na(w))) w <- NULL else { w[is.na(w)] <- median(w, na.rm = TRUE); w <- sqrt(w) }
  corval <- safe_dup_cor(E_ct, design, block = meta_ct$study_id)

  fit0 <- lmFit(E_ct, design,
                block       = if (corval > 0) meta_ct$study_id else NULL,
                correlation = if (corval > 0) corval else NULL,
                weights     = w)
  cont <- makeContrasts(High_vs_Low = MNPBurdenHigh - MNPBurdenLow, levels = design)
  fit1 <- contrasts.fit(fit0, cont)
  fit2 <- eBayes(fit1, trend = TRUE)

  tt <- topTable(fit2, number = Inf, sort.by = "P")
  tt$gene <- rownames(tt)
  m_hi <- rowMeans(E_ct[, meta_ct$MNPBurden == "High", drop = FALSE])
  m_lo <- rowMeans(E_ct[, meta_ct$MNPBurden == "Low",  drop = FALSE])
  tt$mean_High <- m_hi[tt$gene]; tt$mean_Low <- m_lo[tt$gene]

  tt <- tt %>%
    dplyr::mutate(
      sig_q10   = adj.P.Val < Q_CUTOFF_DE & abs(logFC) >= L2FC_MIN,
      sig_q05   = adj.P.Val < 0.05        & abs(logFC) >= L2FC_MIN,
      direction = ifelse(logFC > 0, "UP_in_High", "DOWN_in_High"),
      predicted.celltype = ct
    )

  fp <- file.path(OUTDIR_DE, sprintf("pb_DE_%s.tsv", gsub("[^A-Za-z0-9._-]", "_", ct)))
  readr::write_tsv(tt, fp)

  tibble(
    predicted.celltype = ct,
    n_samples_High = as.integer(grp_counts["High"]),
    n_samples_Low  = as.integer(grp_counts["Low"]),
    n_sig_q10      = sum(tt$sig_q10, na.rm = TRUE),
    n_sig_q05      = sum(tt$sig_q05, na.rm = TRUE)
  )
}

cts <- sort(unique(pb_meta$predicted.celltype))
summ <- dplyr::bind_rows(lapply(cts, run_one_celltype))
readr::write_tsv(summ, file.path(OUTDIR_DE, "pb_DE_summary_by_celltype.tsv"))

# --------------------- 3) Enrichment (ORA) per CT ---------------------
# choose CTs with enough signal
cts_run <- summ %>%
  dplyr::mutate(total_sig = dplyr::coalesce(n_sig_q10, 0L)) %>%
  dplyr::filter(total_sig >= ENR_MIN_SIG_GENES) %>%
  dplyr::arrange(dplyr::desc(total_sig)) %>%
  dplyr::pull(predicted.celltype)

if (length(cts_run) == 0) {
  message("No CTs with ≥ ", ENR_MIN_SIG_GENES, " sig genes at q<", Q_CUTOFF_DE, ". Family-wise step will be empty.")
} else {
  for (ct in cts_run) {
    de_fp <- file.path(OUTDIR_DE, sprintf("pb_DE_%s.tsv", gsub("[^A-Za-z0-9._-]", "_", ct)))
    de <- readr::read_tsv(de_fp, show_col_types = FALSE)
    universe <- unique(de$gene)
    up   <- de %>% dplyr::filter(sig_q10, logFC >  0) %>% dplyr::pull(gene) %>% setdiff(Y_GENES)
    down <- de %>% dplyr::filter(sig_q10, logFC <= 0) %>% dplyr::pull(gene) %>% setdiff(Y_GENES)

    # GO BP
    ego_up   <- run_go_bp(up,   universe)
    ego_down <- run_go_bp(down, universe)
    if (!is.null(ego_up))   readr::write_csv(as.data.frame(ego_up),   file.path(OUTDIR_ENR, sprintf("GO_BP_UP_%s.csv",   ct)))
    if (!is.null(ego_down)) readr::write_csv(as.data.frame(ego_down), file.path(OUTDIR_ENR, sprintf("GO_BP_DOWN_%s.csv", ct)))

    # Hallmark
    hm_up   <- run_hallmark(up,   universe)
    hm_down <- run_hallmark(down, universe)
    if (!is.null(hm_up) && nrow(as.data.frame(hm_up)) > 0)
      readr::write_csv(as.data.frame(hm_up),   file.path(OUTDIR_ENR, sprintf("MSIG_HALLMARK_UP_%s.csv",   ct)))
    if (!is.null(hm_down) && nrow(as.data.frame(hm_down)) > 0)
      readr::write_csv(as.data.frame(hm_down), file.path(OUTDIR_ENR, sprintf("MSIG_HALLMARK_DOWN_%s.csv", ct)))
  }
}

# --------------------- 4) Family-wise BH @ 0.05 across ORA outputs ---------------------
pats <- c("^GO_BP_UP_.*\\.csv$", "^GO_BP_DOWN_.*\\.csv$",
          "^MSIG_HALLMARK_UP_.*\\.csv$", "^MSIG_HALLMARK_DOWN_.*\\.csv$")
files <- unlist(lapply(pats, function(p) list.files(OUTDIR_ENR, pattern = p, full.names = TRUE)))

if (!length(files)) {
  message("No enrichment CSVs found in: ", OUTDIR_ENR)
} else {
  enr <- lapply(files, function(fp) {
    df <- suppressMessages(readr::read_csv(fp, show_col_types = FALSE))
    if (!nrow(df)) return(NULL)
    base <- tools::file_path_sans_ext(basename(fp))
    family    <- if (startsWith(base, "GO_BP")) "GO_BP" else "HALLMARK"
    direction <- if (grepl("_UP_", base)) "UP" else "DOWN"
    celltype  <- sub("^(GO_BP|MSIG_HALLMARK)_(UP|DOWN)_", "", base)
    df %>% dplyr::mutate(pathway_family = family, direction = direction, predicted.celltype = celltype)
  }) %>% dplyr::bind_rows()

  if (!nrow(enr)) {
    message("Parsed zero rows from enrichment CSVs.")
  } else {
    if (!("pvalue" %in% names(enr))) stop("ClusterProfiler CSVs need a 'pvalue' column.")
    enr_fw <- enr %>%
      dplyr::group_by(pathway_family, predicted.celltype, direction) %>%
      dplyr::mutate(p_adj_family = p.adjust(pvalue, method = "BH"),
                    sig_family   = p_adj_family <= PADJ_CUT_FAMILY) %>%
      dplyr::ungroup()

    out_fw <- file.path(OUTDIR_ENR, "enrichment_all.familyFDR.q05.tsv")
    readr::write_tsv(enr_fw, out_fw)
    message("Wrote family-wise table: ", out_fw)

    # console summary
    totals_by_ct <- enr_fw %>%
      dplyr::filter(sig_family) %>%
      dplyr::count(predicted.celltype, name = "n_sig_pathways") %>%
      dplyr::arrange(dplyr::desc(n_sig_pathways))

    print(totals_by_ct, n = Inf)
    total_sig <- sum(totals_by_ct$n_sig_pathways)
    cat(sprintf("\nTERM placentae — family-wise BH (q ≤ %.2f): %d significant pathways total across %d cell types.\n",
                PADJ_CUT_FAMILY, total_sig, nrow(totals_by_ct)))
  }
}

# -------------------------------- RUN --------------------------------
# Nothing else to call; everything runs on load.

# ================== TERM visualization to match PTB macrophage style ==================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(ggplot2)
  library(clusterProfiler); library(enrichplot); library(msigdbr); library(org.Hs.eg.db)
})

# --------- Paths (adjust if needed) ---------
OUTDIR_DE   <- "pb_Term_MNPBurden_limma_SCT"
OUTDIR_ENR  <- file.path(OUTDIR_DE, "cp_results")
FW_TABLE    <- file.path(OUTDIR_ENR, "enrichment_all.familyFDR.q05.tsv")  # family-wise FDR ≤ 0.05

dir.create(OUTDIR_ENR, showWarnings = FALSE, recursive = TRUE)

# --------- Helper: GeneRatio numeric ---------
.gene_ratio_num <- function(gr) {
  if (is.numeric(gr)) return(gr)
  as.numeric(vapply(strsplit(as.character(gr), "/"), function(v) {
    if (length(v)==2) suppressWarnings(as.numeric(v[1]) / as.numeric(v[2])) else NA_real_
  }, 0.0))
}

# --------- Preference list for NK tolerance / remodeling cnet focus ---------
PREF_TOLERANCE <- c(
  "collagen fibril organization","extracellular matrix organization","collagen metabolic process",
  "regulation of cell adhesion","negative regulation of immune effector process",
  "response to transforming growth factor beta","epithelial to mesenchymal transition",
  "myogenesis","hypoxia","xenobiotic metabolic process"
)

# --------- ORA helpers (same params as your pipeline) ---------
run_go_bp <- function(genes_sym, universe_sym) {
  if (length(genes_sym) < 10) return(NULL)
  mm     <- clusterProfiler::bitr(genes_sym,    fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db) %>% dplyr::distinct(SYMBOL, .keep_all = TRUE)
  mm_uni <- clusterProfiler::bitr(universe_sym, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db) %>% dplyr::distinct(SYMBOL, .keep_all = TRUE)
  if (nrow(mm) < 10 || nrow(mm_uni) < 100) return(NULL)
  ego <- enrichGO(
    gene          = mm$ENTREZID, universe = mm_uni$ENTREZID,
    OrgDb         = org.Hs.eg.db, keyType = "ENTREZID",
    ont           = "BP", pAdjustMethod = "BH",
    pvalueCutoff  = 0.05, qvalueCutoff = 0.10, readable = TRUE
  )
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)
  suppressWarnings(simplify(ego, cutoff = 0.6, by = "p.adjust", select_fun = min))
}
run_hallmark <- function(genes_sym, universe_sym) {
  if (length(genes_sym) < 10) return(NULL)
  mm     <- clusterProfiler::bitr(genes_sym,    fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db) %>% dplyr::distinct(SYMBOL, .keep_all = TRUE)
  mm_uni <- clusterProfiler::bitr(universe_sym, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db) %>% dplyr::distinct(SYMBOL, .keep_all = TRUE)
  if (nrow(mm) < 10 || nrow(mm_uni) < 100) return(NULL)
  H <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)
  enricher(mm$ENTREZID, TERM2GENE = H, universe = mm_uni$ENTREZID, pAdjustMethod = "BH")
}

# --------- Family-wise rings: fetch terms that pass q≤0.05 for a CT/family/direction ---------
get_fw_sig_terms <- function(ct, family, direction, fw_table = FW_TABLE) {
  if (!file.exists(fw_table)) return(character(0))
  fw <- suppressMessages(readr::read_tsv(fw_table, show_col_types = FALSE))
  req <- c("predicted.celltype","pathway_family","direction","Description","p_adj_family","sig_family")
  if (!all(req %in% colnames(fw))) return(character(0))
  fw %>%
    dplyr::filter(predicted.celltype == ct,
                  pathway_family == family,
                  direction == direction,
                  isTRUE(sig_family) | p_adj_family <= 0.05) %>%
    dplyr::pull(Description) %>% unique()
}

# --------- Dotplot (theme_classic) matching style; Top-10 highlighted; FW q≤0.05 ring ---------
dotplot_match_ptb <- function(ego, title, n_show = 30, top_n = 10, wrap = 45,
                              fw_terms = character(0),
                              width = 10, height = 9, outfile) {
  if (is.null(ego)) return(invisible(NULL))
  df <- as.data.frame(ego); if (!nrow(df)) return(invisible(NULL))
  df$GeneRatio_num <- if ("GeneRatio" %in% names(df)) .gene_ratio_num(df$GeneRatio) else NA_real_
  if (all(is.na(df$GeneRatio_num))) {
    df$GeneRatio_num <- if ("setSize" %in% names(df)) df$Count / pmax(df$setSize, df$Count)
                        else df$Count / max(df$Count, na.rm = TRUE)
  }
  df <- df[order(df$p.adjust, -df$GeneRatio_num), ]
  show <- head(df, min(n_show, nrow(df)))
  show$is_top <- FALSE; if (nrow(show)) show$is_top[seq_len(min(top_n, nrow(show)))] <- TRUE
  show$is_fw  <- show$Description %in% fw_terms
  show$Description <- factor(show$Description, levels = rev(show$Description))

  p <- ggplot(show, aes(GeneRatio_num, Description)) +
    # base points
    geom_point(aes(size = Count, fill = p.adjust), shape = 21, color = NA, alpha = 0.65) +
    # ring = family-wise q≤0.05
    geom_point(data = subset(show, is_fw),
               aes(size = Count, fill = p.adjust),
               shape = 21, color = "red3", stroke = 1.2, alpha = 1) +
    # bold outline = Top-10
    geom_point(data = subset(show, is_top),
               aes(size = Count, fill = p.adjust),
               shape = 21, color = "black", stroke = 1.1, alpha = 1) +
    scale_fill_continuous(name = "p.adjust", trans = "reverse") +
    scale_size_continuous(name = "Count") +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = wrap)) +
    labs(x = "GeneRatio", y = NULL, title = title) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.y = element_text(hjust = 1),
          legend.position = "bottom",
          plot.margin = margin(10, 20, 10, 140)) +
    coord_cartesian(clip = "off")

  ggsave(outfile, p, width = width, height = height, dpi = 300, bg = "white")
  invisible(p)
}

# --------- Subset enrichResult by description (for cnet) ---------
subset_enrich_by_description <- function(ego, terms) {
  if (is.null(ego)) return(NULL)
  df <- as.data.frame(ego); if (!nrow(df)) return(NULL)
  idx <- which(df$Description %in% terms)
  if (!length(idx)) return(NULL)
  ego2 <- ego; ego2@result <- ego@result[idx, , drop = FALSE]
  ego2
}

# --------- MAIN: build plots for CTs (NK Cell, T Cell if available) ---------
cts_to_plot <- c("NK Cell", "T Cell")

for (ct in cts_to_plot) {
  de_fp <- file.path(OUTDIR_DE, sprintf("pb_DE_%s.tsv", gsub("[^A-Za-z0-9._-]", "_", ct)))
  if (!file.exists(de_fp)) { message("No DE table for: ", ct); next }
  de <- suppressMessages(readr::read_tsv(de_fp, show_col_types = FALSE))
  if (!nrow(de)) { message("Empty DE table for: ", ct); next }

  # define universe + up/down gene sets (match your DE flags)
  universe <- unique(de$gene)
  up   <- de %>% dplyr::filter(sig_q10, logFC >  0) %>% dplyr::pull(gene)
  down <- de %>% dplyr::filter(sig_q10, logFC <= 0) %>% dplyr::pull(gene)

  # Run ORA fresh (ensures enrichResult objects for cnet)
  ego_up   <- run_go_bp(up,   universe)
  ego_down <- run_go_bp(down, universe)
  hm_up    <- run_hallmark(up,   universe)
  hm_down  <- run_hallmark(down, universe)

  # family-wise rings for this CT
  fw_go_up   <- get_fw_sig_terms(ct, "GO_BP",     "UP")
  fw_go_down <- get_fw_sig_terms(ct, "GO_BP",     "DOWN")
  fw_hm_up   <- get_fw_sig_terms(ct, "HALLMARK",  "UP")
  fw_hm_down <- get_fw_sig_terms(ct, "HALLMARK",  "DOWN")

  # ---- DOTPLOTS (match PTB macrophage style) ----
  if (!is.null(ego_up)) {
    dotplot_match_ptb(
      ego_up,
      title   = sprintf("GO BP (UP) — %s (TERM; Top-10 highlighted; red ring = family-wise q≤0.05)", ct),
      n_show  = 30, top_n = 10, wrap = 48,
      fw_terms = fw_go_up,
      outfile = file.path(OUTDIR_ENR, sprintf("TERM_%s_GO_BP_UP.matchPTB.png", gsub("[^A-Za-z0-9._-]", "_", ct)))
    )
  }
  if (!is.null(ego_down)) {
    dotplot_match_ptb(
      ego_down,
      title   = sprintf("GO BP (DOWN) — %s (TERM; Top-10 highlighted; red ring = family-wise q≤0.05)", ct),
      n_show  = 30, top_n = 10, wrap = 48,
      fw_terms = fw_go_down,
      outfile = file.path(OUTDIR_ENR, sprintf("TERM_%s_GO_BP_DOWN.matchPTB.png", gsub("[^A-Za-z0-9._-]", "_", ct)))
    )
  }
  if (!is.null(hm_up) && nrow(as.data.frame(hm_up)) > 0) {
    dotplot_match_ptb(
      hm_up,
      title   = sprintf("Hallmark (UP) — %s (TERM; Top-10 highlighted; red ring = family-wise q≤0.05)", ct),
      n_show  = 20, top_n = 10, wrap = 48,
      fw_terms = fw_hm_up,
      outfile = file.path(OUTDIR_ENR, sprintf("TERM_%s_HALLMARK_UP.matchPTB.png", gsub("[^A-Za-z0-9._-]", "_", ct)))
    )
  }
  if (!is.null(hm_down) && nrow(as.data.frame(hm_down)) > 0) {
    dotplot_match_ptb(
      hm_down,
      title   = sprintf("Hallmark (DOWN) — %s (TERM; Top-10 highlighted; red ring = family-wise q≤0.05)", ct),
      n_show  = 20, top_n = 10, wrap = 48,
      fw_terms = fw_hm_down,
      outfile = file.path(OUTDIR_ENR, sprintf("TERM_%s_HALLMARK_DOWN.matchPTB.png", gsub("[^A-Za-z0-9._-]", "_", ct)))
    )
  }

  # ---- CNETPLOTS (tolerance/remodeling preferred terms) ----
  if (!is.null(ego_up)) {
    # choose present preferred terms; fallback to top few by p.adjust
    ego_df <- as.data.frame(ego_up)
    cats <- intersect(PREF_TOLERANCE, ego_df$Description)
    if (length(cats) < 3) cats <- head(ego_df$Description[order(ego_df$p.adjust)], 3)

    ego_up_sub <- subset_enrich_by_description(ego_up, cats)
    if (!is.null(ego_up_sub)) {
      fc_named <- setNames(de$logFC, de$gene)  # SYMBOL -> logFC
      ggsave(
        file.path(OUTDIR_ENR, sprintf("TERM_%s_GO_BP_UP.cnet.matchPTB.png", gsub("[^A-Za-z0-9._-]", "_", ct))),
        enrichplot::cnetplot(ego_up_sub, circular = TRUE, colorEdge = TRUE,
                             foldChange = fc_named,
                             showCategory = nrow(as.data.frame(ego_up_sub))),
        width = 10, height = 9, dpi = 300, bg = "white"
      )
      writeLines(cats, file.path(OUTDIR_ENR, sprintf("TERM_%s_GO_BP_UP.cnet_terms.txt", gsub("[^A-Za-z0-9._-]", "_", ct))))
    }
  }
  if (!is.null(ego_down)) {
    ego_df <- as.data.frame(ego_down)
    cats <- intersect(PREF_TOLERANCE, ego_df$Description)
    if (length(cats) < 3) cats <- head(ego_df$Description[order(ego_df$p.adjust)], 3)
    ego_dn_sub <- subset_enrich_by_description(ego_down, cats)
    if (!is.null(ego_dn_sub)) {
      fc_named <- setNames(de$logFC, de$gene)
      ggsave(
        file.path(OUTDIR_ENR, sprintf("TERM_%s_GO_BP_DOWN.cnet.matchPTB.png", gsub("[^A-Za-z0-9._-]", "_", ct))),
        enrichplot::cnetplot(ego_dn_sub, circular = TRUE, colorEdge = TRUE,
                             foldChange = fc_named,
                             showCategory = nrow(as.data.frame(ego_dn_sub))),
        width = 18, height = 5, dpi = 300, bg = "white"
      )
      writeLines(cats, file.path(OUTDIR_ENR, sprintf("TERM_%s_GO_BP_DOWN.cnet_terms.txt", gsub("[^A-Za-z0-9._-]", "_", ct))))
    }
  }
}

message("TERM visualization complete → ", OUTDIR_ENR,
        " (files named TERM_<CT>_*_matchPTB*.png and cnet_terms.txt)")



# ================== TERM — T Cell visualization (match PTB macrophage style) ==================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(ggplot2)
  library(clusterProfiler); library(enrichplot); library(msigdbr); library(org.Hs.eg.db)
})

OUTDIR_DE   <- "pb_Term_MNPBurden_limma_SCT"
OUTDIR_ENR  <- file.path(OUTDIR_DE, "cp_results")
FW_TABLE    <- file.path(OUTDIR_ENR, "enrichment_all.familyFDR.q05.tsv")  # for family-wise rings
dir.create(OUTDIR_ENR, showWarnings = FALSE, recursive = TRUE)

.ct <- "T Cell"  # <— target cell type

# ---------- helpers ----------
.gene_ratio_num <- function(gr) {
  if (is.numeric(gr)) return(gr)
  as.numeric(vapply(strsplit(as.character(gr), "/"), function(v) {
    if (length(v)==2) suppressWarnings(as.numeric(v[1]) / as.numeric(v[2])) else NA_real_
  }, 0.0))
}

run_go_bp <- function(genes_sym, universe_sym) {
  if (length(genes_sym) < 5) return(NULL)           # allow small sets; may return 0 rows
  mm     <- clusterProfiler::bitr(genes_sym,    "SYMBOL","ENTREZID", org.Hs.eg.db) %>% dplyr::distinct(SYMBOL, .keep_all = TRUE)
  mm_uni <- clusterProfiler::bitr(universe_sym, "SYMBOL","ENTREZID", org.Hs.eg.db) %>% dplyr::distinct(SYMBOL, .keep_all = TRUE)
  if (nrow(mm) < 5 || nrow(mm_uni) < 50) return(NULL)
  ego <- enrichGO(
    gene          = mm$ENTREZID, universe = mm_uni$ENTREZID,
    OrgDb         = org.Hs.eg.db, keyType = "ENTREZID",
    ont           = "BP", pAdjustMethod = "BH",
    pvalueCutoff  = 0.05, qvalueCutoff = 0.10, readable = TRUE
  )
  if (is.null(ego) || !nrow(as.data.frame(ego))) return(NULL)
  suppressWarnings(simplify(ego, cutoff = 0.6, by = "p.adjust", select_fun = min))
}

run_hallmark <- function(genes_sym, universe_sym) {
  if (length(genes_sym) < 5) return(NULL)
  mm     <- clusterProfiler::bitr(genes_sym,    "SYMBOL","ENTREZID", org.Hs.eg.db) %>% dplyr::distinct(SYMBOL, .keep_all = TRUE)
  mm_uni <- clusterProfiler::bitr(universe_sym, "SYMBOL","ENTREZID", org.Hs.eg.db) %>% dplyr::distinct(SYMBOL, .keep_all = TRUE)
  if (nrow(mm) < 5 || nrow(mm_uni) < 50) return(NULL)
  H <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)
  enricher(mm$ENTREZID, TERM2GENE = H, universe = mm_uni$ENTREZID, pAdjustMethod = "BH")
}

get_fw_sig_terms <- function(ct, family, direction, fw_table = FW_TABLE) {
  if (!file.exists(fw_table)) return(character(0))
  fw <- suppressMessages(readr::read_tsv(fw_table, show_col_types = FALSE))
  req <- c("predicted.celltype","pathway_family","direction","Description","p_adj_family","sig_family")
  if (!all(req %in% colnames(fw))) return(character(0))
  fw %>%
    dplyr::filter(predicted.celltype == ct,
                  pathway_family == family,
                  direction == direction,
                  isTRUE(sig_family) | p_adj_family <= 0.05) %>%
    dplyr::pull(Description) %>% unique()
}

dotplot_match_ptb <- function(ego, title, n_show = 30, top_n = 10, wrap = 48,
                              fw_terms = character(0),
                              width = 18, height = 5, outfile) {
  if (is.null(ego)) return(invisible(NULL))
  df <- as.data.frame(ego); if (!nrow(df)) return(invisible(NULL))
  df$GeneRatio_num <- if ("GeneRatio" %in% names(df)) .gene_ratio_num(df$GeneRatio) else NA_real_
  if (all(is.na(df$GeneRatio_num))) {
    df$GeneRatio_num <- if ("setSize" %in% names(df)) df$Count / pmax(df$setSize, df$Count)
                        else df$Count / max(df$Count, na.rm = TRUE)
  }
  df <- df[order(df$p.adjust, -df$GeneRatio_num), ]
  show <- head(df, min(n_show, nrow(df)))
  show$is_top <- FALSE; if (nrow(show)) show$is_top[seq_len(min(top_n, nrow(show)))] <- TRUE
  show$is_fw  <- show$Description %in% fw_terms
  show$Description <- factor(show$Description, levels = rev(show$Description))

  p <- ggplot(show, aes(GeneRatio_num, Description)) +
    geom_point(aes(size = Count, fill = p.adjust), shape = 21, color = NA, alpha = 0.65) +
    geom_point(data = subset(show, is_fw),
               aes(size = Count, fill = p.adjust),
               shape = 21, color = "red3", stroke = 1.2, alpha = 1) +
    geom_point(data = subset(show, is_top),
               aes(size = Count, fill = p.adjust),
               shape = 21, color = "black", stroke = 1.1, alpha = 1) +
    scale_fill_continuous(name = "p.adjust", trans = "reverse") +
    scale_size_continuous(name = "Count") +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = wrap)) +
    labs(x = "GeneRatio", y = NULL, title = title) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.y = element_text(hjust = 1),
          legend.position = "bottom",
          plot.margin = margin(10, 20, 10, 140)) +
    coord_cartesian(clip = "off")

  ggsave(outfile, p, width = width, height = height, dpi = 300, bg = "white")
  invisible(p)
}

subset_enrich_by_description <- function(ego, terms) {
  if (is.null(ego)) return(NULL)
  df <- as.data.frame(ego); if (!nrow(df)) return(NULL)
  idx <- which(df$Description %in% terms)
  if (!length(idx)) return(NULL)
  ego2 <- ego; ego2@result <- ego@result[idx, , drop = FALSE]
  ego2
}

# Preference list for T-cell tolerance/exhaustion/remodeling cnet focus
PREF_T_TOL <- c(
  "negative regulation of immune effector process","negative regulation of T cell activation",
  "regulation of T cell differentiation","T cell migration","response to transforming growth factor beta",
  "epithelial to mesenchymal transition","hypoxia","HALLMARK_TGF_BETA_SIGNALING","HALLMARK_IL2_STAT5_SIGNALING"
)

# ---------- load DE and run enrichment ----------
de_fp <- file.path(OUTDIR_DE, sprintf("pb_DE_%s.tsv", gsub("[^A-Za-z0-9._-]", "_", .ct)))
if (!file.exists(de_fp)) stop("DE table not found for T Cell: ", de_fp)
de <- suppressMessages(readr::read_tsv(de_fp, show_col_types = FALSE))
if (!nrow(de)) stop("Empty DE table for T Cell.")

universe <- unique(de$gene)
up   <- de %>% dplyr::filter(sig_q10, logFC >  0) %>% dplyr::pull(gene)
down <- de %>% dplyr::filter(sig_q10, logFC <= 0) %>% dplyr::pull(gene)

ego_up   <- run_go_bp(up,   universe)
ego_down <- run_go_bp(down, universe)
hm_up    <- run_hallmark(up,   universe)
hm_down  <- run_hallmark(down, universe)

# family-wise significant terms to ring (often none for T cells at q≤0.05)
fw_go_up   <- get_fw_sig_terms(.ct, "GO_BP",    "UP")
fw_go_down <- get_fw_sig_terms(.ct, "GO_BP",    "DOWN")
fw_hm_up   <- get_fw_sig_terms(.ct, "HALLMARK", "UP")
fw_hm_down <- get_fw_sig_terms(.ct, "HALLMARK", "DOWN")

# ---------- DOTPLOTS ----------
made_any <- FALSE

if (!is.null(ego_up)) {
  dotplot_match_ptb(
    ego_up,
    title   = sprintf("GO BP (UP) — %s (TERM; Top-10 highlighted; red ring = family-wise q≤0.05)", .ct),
    n_show  = 30, top_n = 10, wrap = 48,
    fw_terms = fw_go_up,
    outfile = file.path(OUTDIR_ENR, sprintf("TERM_%s_GO_BP_UP.matchPTB.png", gsub("[^A-Za-z0-9._-]", "_", .ct)))
  )
  made_any <- TRUE
}
if (!is.null(ego_down)) {
  dotplot_match_ptb(
    ego_down,
    title   = sprintf("GO BP (DOWN) — %s (TERM; Top-10 highlighted; red ring = family-wise q≤0.05)", .ct),
    n_show  = 30, top_n = 10, wrap = 48,
    fw_terms = fw_go_down,
    outfile = file.path(OUTDIR_ENR, sprintf("TERM_%s_GO_BP_DOWN.matchPTB.png", gsub("[^A-Za-z0-9._-]", "_", .ct)))
  )
  made_any <- TRUE
}
if (!is.null(hm_up) && nrow(as.data.frame(hm_up)) > 0) {
  dotplot_match_ptb(
    hm_up,
    title   = sprintf("Hallmark (UP) — %s (TERM; Top-10 highlighted; red ring = family-wise q≤0.05)", .ct),
    n_show  = 20, top_n = 10, wrap = 48,
    fw_terms = fw_hm_up,
    outfile = file.path(OUTDIR_ENR, sprintf("TERM_%s_HALLMARK_UP.matchPTB.png", gsub("[^A-Za-z0-9._-]", "_", .ct)))
  )
  made_any <- TRUE
}
if (!is.null(hm_down) && nrow(as.data.frame(hm_down)) > 0) {
  dotplot_match_ptb(
    hm_down,
    title   = sprintf("Hallmark (DOWN) — %s (TERM; Top-10 highlighted; red ring = family-wise q≤0.05)", .ct),
    n_show  = 20, top_n = 10, wrap = 48,
    fw_terms = fw_hm_down,
    outfile = file.path(OUTDIR_ENR, sprintf("TERM_%s_HALLMARK_DOWN.matchPTB.png", gsub("[^A-Za-z0-9._-]", "_", .ct)))
  )
  made_any <- TRUE
}

# Fallback: if up/down are too small and none plotted, build an ALL-genes GO BP panel
if (!made_any) {
  all_sig <- de %>% dplyr::filter(sig_q10) %>% dplyr::pull(gene)
  ego_all <- run_go_bp(all_sig, universe)
  if (!is.null(ego_all)) {
    dotplot_match_ptb(
      ego_all,
      title   = sprintf("GO BP (ALL sig genes) — %s (TERM; exploratory)", .ct),
      n_show  = 20, top_n = 10, wrap = 48,
      fw_terms = character(0),
      outfile = file.path(OUTDIR_ENR, sprintf("TERM_%s_GO_BP_ALLsig.matchPTB.png", gsub("[^A-Za-z0-9._-]", "_", .ct)))
    )
    made_any <- TRUE
  } else {
    message("T Cell: insufficient signal for ORA (even with ALL sig genes).")
  }
}

# ---------- CNET (optional; tolerance/remodeling targets) ----------
if (!is.null(ego_up)) {
  ego_df <- as.data.frame(ego_up)
  cats <- intersect(PREF_T_TOL, ego_df$Description)
  if (length(cats) < 3) cats <- head(ego_df$Description[order(ego_df$p.adjust)], 3)
  ego_up_sub <- subset_enrich_by_description(ego_up, cats)
  if (!is.null(ego_up_sub)) {
    fc_named <- setNames(de$logFC, de$gene)
    ggsave(
      file.path(OUTDIR_ENR, sprintf("TERM_%s_GO_BP_UP.cnet.matchPTB.png", gsub("[^A-Za-z0-9._-]", "_", .ct))),
      enrichplot::cnetplot(ego_up_sub, circular = TRUE, colorEdge = TRUE,
                           foldChange = fc_named,
                           showCategory = nrow(as.data.frame(ego_up_sub))),
      width = 10, height = 9, dpi = 300, bg = "white"
    )
    writeLines(cats, file.path(OUTDIR_ENR, sprintf("TERM_%s_GO_BP_UP.cnet_terms.txt", gsub("[^A-Za-z0-9._-]", "_", .ct))))
  }
}

message("T Cell visualization complete -> ", OUTDIR_ENR)
# ================== TERM — Combined GO+Hallmark plots for NK & T cells ==================
suppressPackageStartupMessages({
    library(dplyr); library(readr); library(stringr); library(ggplot2)
})

OUTDIR_DE   <- "pb_Term_MNPBurden_limma_SCT"
OUTDIR_ENR  <- file.path(OUTDIR_DE, "cp_results")
FW_TABLE    <- file.path(OUTDIR_ENR, "enrichment_all.familyFDR.q05.tsv")  # family-wise BH table
dir.create(OUTDIR_ENR, showWarnings = FALSE, recursive = TRUE)

# -- helpers -------------------------------------------------------------------
gene_ratio_num <- function(gr, Count = NULL, setSize = NULL) {
    if ("GeneRatio" %in% names(gr)) gr <- gr$GeneRatio
    if (is.numeric(gr)) return(gr)
    v <- strsplit(as.character(gr), "/")
    out <- suppressWarnings(as.numeric(vapply(v, function(x) if (length(x)==2) as.numeric(x[1])/as.numeric(x[2]) else NA_real_, 0.0)))
    if (all(is.na(out))) {
        if (!is.null(setSize) && !is.null(Count)) out <- Count / pmax(setSize, Count)
        else if (!is.null(Count)) out <- Count / max(Count, na.rm = TRUE)
    }
    out
}

reorder_within <- function(x, by, within, fun = mean, sep = "___") {
    x2 <- paste(x, within, sep = sep); stats::reorder(x2, by, FUN = fun)
}
scale_y_reordered <- function(..., sep = "___") {
    ggplot2::scale_y_discrete(labels = function(x) gsub(paste0(sep, ".*$"), "", x), ...)
}

# -- core plotting fn ----------------------------------------------------------
plot_combined_ct <- function(ct, q_main = 0.05, q_fallback = 0.3, n_per_family = 10,
                             wrap = 48, width = 18, height = 5) {
    if (!file.exists(FW_TABLE)) stop("Family-wise table not found: ", FW_TABLE)
    fw <- suppressMessages(readr::read_tsv(FW_TABLE, show_col_types = FALSE))
    req <- c("predicted.celltype","pathway_family","direction","Description","p_adj_family","sig_family","Count")
    if (!all(req %in% names(fw))) stop("Missing columns in family-wise table: ", paste(setdiff(req, names(fw)), collapse=", "))
    
    df_ct <- fw %>% dplyr::filter(predicted.celltype == ct) %>% dplyr::mutate(
        GeneRatio_num = gene_ratio_num(., Count = Count, setSize = if ("setSize" %in% names(.)) .$setSize else NULL)
    )
    
    # choose significance set
    sel <- df_ct %>% dplyr::filter(isTRUE(sig_family) | p_adj_family <= q_main)
    label_q <- sprintf("family-wise FDR \u2264 %.2f", q_main)
    
    # fallback if empty
    exploratory <- FALSE
    if (!nrow(sel)) {
        sel <- df_ct %>% dplyr::filter(is.finite(p_adj_family), p_adj_family <= q_fallback)
        exploratory <- TRUE
        label_q <- sprintf("exploratory (family-wise FDR \u2264 %.2f; none at \u2264 %.2f)", q_fallback, q_main)
        if (!nrow(sel)) {
            message("No terms for ", ct, " even at q≤", q_fallback, ". Skipping.")
            return(invisible(NULL))
        }
    }
    
    # Top-N per family
    top_by_family <- sel %>%
        dplyr::arrange(p_adj_family, dplyr::desc(GeneRatio_num)) %>%
        dplyr::group_by(pathway_family) %>%
        dplyr::slice_head(n = n_per_family) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(y_lab = reorder_within(Description, -p_adj_family, pathway_family))
    
    # Plot
    p <- ggplot(top_by_family, aes(x = GeneRatio_num, y = y_lab)) +
        geom_point(aes(size = Count, fill = p_adj_family, shape = direction),
                   color = "black", stroke = 0.6) +
        facet_wrap(~ pathway_family, scales = "free_y") +
        scale_y_reordered() +
        scale_fill_continuous(name = "FDR (family)", trans = "reverse") +
        scale_size_continuous(name = "Count") +
        scale_shape_manual(values = c(21, 24)) +
        labs(x = "GeneRatio", y = NULL,
             title = sprintf("%s — Top 10 per pathway family (%s)", ct, label_q)) +
        theme_classic(base_size = 12) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              axis.text.y = element_text(hjust = 1),
              legend.position = "right",
              plot.margin = margin(10, 20, 10, 140)) +
        coord_cartesian(clip = "off")
    
    stub <- sprintf("TERM_%s.Combined_GO_BP_HALLMARK.Top10", gsub("[^A-Za-z0-9._-]", "_", ct))
    if (exploratory) stub <- paste0(stub, ".exploratory_q", gsub("\\.","",sprintf("%.2f", q_fallback)))
    png_out <- file.path(OUTDIR_ENR, paste0(stub, ".png"))
    pdf_out <- file.path(OUTDIR_ENR, paste0(stub, ".pdf"))
    
    ggsave(png_out, p, width = width, height = height, dpi = 300, bg = "white")
    ggsave(pdf_out, p, width = width, height = height, device = cairo_pdf, bg = "white")
    
    message("Saved: ", png_out, " and PDF")
    invisible(p)
}

# ------------------------- RUN for NK Cell and T Cell -------------------------
plot_combined_ct("NK Cell", q_main = 0.05, q_fallback = 0.10, n_per_family = 10)
plot_combined_ct("T Cell",  q_main = 0.05, q_fallback = 0.3, n_per_family = 10)





 "Delivery_class","MNPBurden","MNPBurden_Delivery_class","GA_group",
  "platform", "Preeclampsia_bin","Infections_any","Ethnicity_Race","Infant gender",
  "predicted.celltype","sample_id_std","study_id",