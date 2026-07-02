# pb_preterm_MNP_SCT_DE_enrich_dropin.R
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr); library(tidyr); library(purrr); library(stringr)
  library(readr); library(data.table)
  library(limma); library(matrixStats)
  library(clusterProfiler); library(enrichplot)
  library(org.Hs.eg.db); library(msigdbr)
})

# ====================== CONFIG ======================
OUTDIR_DE   <- "pb_Preterm_MNPBurden_limma_SCT"
OUTDIR_ENR  <- file.path(OUTDIR_DE, "cp_results")
dir.create(OUTDIR_DE,  showWarnings = FALSE, recursive = TRUE)
dir.create(OUTDIR_ENR, showWarnings = FALSE, recursive = TRUE)

Q_CUTOFF             <- 0.10   # FDR threshold for DE
L2FC_MIN             <- 0.25   # |log2FC| minimum for DE
MIN_SAMPLES_PER_GROUP<- 2      # require >=2 pseudo-samples/group
VAR_EPS              <- 1e-8   # min variance to keep gene in a CT
RUN_ENRICHMENT       <- TRUE   # set FALSE to stop after DE summary
ENR_MIN_SIG_GENES    <- 10     # per CT, min # of sig genes to run ORA

# immunomod categories for cnet (edit as needed)
IMMUNE_CATEGORIES <- c(
  "cytokine-mediated signaling pathway",
  "regulation of innate immune response",
  "inflammatory response",
  "leukocyte mediated immunity",
  "response to interferon-gamma", "response to interferon-alpha",
  "NF-kappaB signaling", "MAPK cascade",
  "response to oxidative stress", "reactive oxygen species metabolic process",
  "xenobiotic metabolic process", "cellular response to toxic substance",
  "phagocytosis", "autophagy", "hypoxia"
)

# genes to exclude from enrichment (Y-chromosome)
Y_GENES <- c("DDX3Y","KDM5D","RPS4Y1","EIF1AY","UTY","USP9Y","ZFY","TXLNGY","TMSB4Y",
             "PRKY","PCDH11Y","NLGN4Y","AMELY")

# lnc/ENSG & pathogen filters (already applied upstream but keep here too)
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

# ====================== META + SUBSET ======================
stopifnot(exists("obj"), !is.null(obj@meta.data))
md <- obj@meta.data

# Delivery_class (Preterm/Term); if absent, derive from MNPBurden_Delivery_class
if (!"Delivery_class" %in% colnames(md)) {
  if ("MNPBurden_Delivery_class" %in% colnames(md)) {
    md$Delivery_class <- sub("^.*-(Preterm|Term)$", "\\1", md$MNPBurden_Delivery_class)
  } else stop("Delivery_class not found; add it or provide MNPBurden_Delivery_class.")
}

# Normalize MNPBurden to High/Low
if (!"MNPBurden" %in% colnames(md)) {
  if ("MNPBurden_Delivery_class" %in% colnames(md)) {
    md$MNPBurden <- sub("^(High|Low)-.*$", "\\1", md$MNPBurden_Delivery_class)
  } else stop("MNPBurden missing; add it or MNPBurden_Delivery_class.")
} else {
  md$MNPBurden <- sub("^(High|Low).*", "\\1", md$MNPBurden)
}

obj@meta.data <- md
req_cols <- c("predicted.celltype","MNPBurden","Delivery_class","study_id")
miss <- setdiff(req_cols, colnames(obj@meta.data))
if (length(miss)) stop("Missing meta columns: ", paste(miss, collapse = ", "))

# Subset to Preterm
obj_ptb <- subset(obj, subset = Delivery_class == "Preterm")
if (ncol(obj_ptb) == 0) stop("No cells with Delivery_class == 'Preterm'.")

# ====================== PSEUDOBULK (SCT/data) ======================
DefaultAssay(obj_ptb) <- if ("SCT" %in% names(obj_ptb@assays)) "SCT" else DefaultAssay(obj_ptb)
bulk <- AggregateExpression(
  obj_ptb, slot = "data",
  group.by = c("predicted.celltype","MNPBurden","study_id"),
  return.seurat = TRUE
)
DefaultAssay(bulk) <- if ("SCT" %in% names(bulk@assays)) "SCT" else DefaultAssay(bulk)

# Build pseudobulk meta
# --- make sure pb_meta join keys are atomic (not list-cols) ---
pb_meta <- pb_meta %>%
  dplyr::mutate(
    predicted.celltype = as.character(predicted.celltype),
    MNPBurden          = as.character(MNPBurden),
    study_id           = as.character(study_id)
  )

# --- robust helper: coerce any list-col to a scalar string ---
vecify_chr <- function(x) {
  if (is.list(x)) vapply(x, function(y) paste(as.character(y), collapse = ";"), character(1))
  else as.character(x)
}

# --- build cell_counts safely from obj_ptb meta ---
cell_counts <- obj_ptb@meta.data %>%
  tibble::as_tibble(rownames = "cell") %>%
  dplyr::mutate(
    predicted.celltype = vecify_chr(.data$predicted.celltype),
    MNPBurden          = vecify_chr(.data$MNPBurden),
    study_id           = vecify_chr(.data$study_id)
  ) %>%
  dplyr::mutate(
    # normalize burden labels to High/Low only
    MNPBurden = sub("^(High|Low).*", "\\1", MNPBurden)
  ) %>%
  dplyr::count(predicted.celltype, MNPBurden, study_id, name = "n_cells")

# --- now the join will work ---
pb_meta <- pb_meta %>%
  dplyr::left_join(cell_counts,
                   by = c("predicted.celltype", "MNPBurden", "study_id"))
# explicit namespace (in case plyr is loaded somewhere)
avail_tbl <- dplyr::summarise(
  dplyr::group_by(pb_meta,
                  predicted.celltype = vecify_chr(predicted.celltype),
                  MNPBurden          = vecify_chr(MNPBurden),
                  .drop = FALSE),
  n_pseudobulk = dplyr::n(),
  .groups = "drop"
) %>%
  dplyr::arrange(predicted.celltype, MNPBurden)

# OR base-R fallback:
tmp <- with(pb_meta, table(vecify_chr(predicted.celltype), vecify_chr(MNPBurden)))
avail_tbl_base <- as.data.frame(tmp)
colnames(avail_tbl_base) <- c("predicted.celltype","MNPBurden","n_pseudobulk")
readr::write_tsv(avail_tbl_base, file.path(OUTDIR_DE, "pseudobulk_sample_counts_by_celltype.tsv"))


# Expression matrix (genes x pseudo-samples)
E <- as.matrix(GetAssayData(bulk, slot = "data"))
stopifnot(ncol(E) == nrow(pb_meta))

# Drop lnc/ENSG + pathogen genes
keep_genes <- rownames(E)[ !grepl(lnc_pattern, rownames(E), perl = TRUE) &
                           !(rownames(E) %in% PATHOGEN_GENES) ]
E <- E[keep_genes, , drop = FALSE]

# ====================== DE per CT (limma-trend on log means) ======================
safe_dup_cor <- function(expr_mat, design, block) {
  if (is.null(block)) return(0)
  tb <- table(block)
  if (all(tb <= 1)) return(0)
  res <- try(duplicateCorrelation(expr_mat, design, block = block), silent = TRUE)
  if (inherits(res, "try-error")) return(0)
  res$consensus
}

run_one_celltype <- function(ct) {
  sel <- which(pb_meta$predicted.celltype == ct)
  if (length(sel) < 2*MIN_SAMPLES_PER_GROUP) return(NULL)

  meta_ct <- pb_meta[sel, , drop = FALSE]
  E_ct    <- E[, sel, drop = FALSE]

  grp_counts <- table(meta_ct$MNPBurden)
  if (!all(c("Low","High") %in% names(grp_counts))) return(NULL)
  if (any(grp_counts[c("Low","High")] < MIN_SAMPLES_PER_GROUP)) return(NULL)

  # Filter zero/near-zero variance genes
  v <- rowVars(E_ct)
  keep <- which(is.finite(v) & v > VAR_EPS)
  if (length(keep) < 50) return(NULL)
  E_ct <- E_ct[keep, , drop = FALSE]

  meta_ct$MNPBurden <- factor(meta_ct$MNPBurden, levels = c("Low","High"))
  design <- model.matrix(~ 0 + MNPBurden, data = meta_ct) # Low, High

  # optional weights ~ sqrt(n_cells)
  w <- meta_ct$n_cells
  if (is.null(w) || all(is.na(w))) {
    weights_vec <- NULL
  } else {
    w[is.na(w)] <- median(w, na.rm = TRUE)
    weights_vec <- sqrt(w)
  }

  corval <- safe_dup_cor(E_ct, design, block = meta_ct$study_id)

  fit0 <- lmFit(E_ct, design,
                block       = if (corval > 0) meta_ct$study_id else NULL,
                correlation = if (corval > 0) corval else NULL,
                weights     = weights_vec)

  cont <- makeContrasts(High_vs_Low = MNPBurdenHigh - MNPBurdenLow, levels = design)
  fit1 <- contrasts.fit(fit0, cont)
  fit2 <- eBayes(fit1, trend = TRUE, robust = FALSE)

  tt <- topTable(fit2, number = Inf, sort.by = "P")
  tt$gene <- rownames(tt)

  m_hi <- rowMeans(E_ct[, meta_ct$MNPBurden == "High", drop = FALSE])
  m_lo <- rowMeans(E_ct[, meta_ct$MNPBurden == "Low",  drop = FALSE])
  tt$mean_High <- m_hi[tt$gene]
  tt$mean_Low  <- m_lo[tt$gene]

  tt <- tt %>%
    mutate(
      sig_q10   = adj.P.Val < Q_CUTOFF & abs(logFC) >= L2FC_MIN,
      sig_q05   = adj.P.Val < 0.05     & abs(logFC) >= L2FC_MIN,
      direction = ifelse(logFC > 0, "UP_in_High", "DOWN_in_High"),
      predicted.celltype = ct
    )

  fp <- file.path(OUTDIR_DE, sprintf("pb_DE_%s.tsv", gsub("[^A-Za-z0-9._-]", "_", ct)))
  write_tsv(tt, fp)

  write_tsv(tibble(gene = tt$gene[tt$sig_q10 & tt$logFC > 0]),
            file.path(OUTDIR_DE, sprintf("genes_UP_q10_%s.tsv",   gsub("[^A-Za-z0-9._-]", "_", ct))))
  write_tsv(tibble(gene = tt$gene[tt$sig_q10 & tt$logFC < 0]),
            file.path(OUTDIR_DE, sprintf("genes_DOWN_q10_%s.tsv", gsub("[^A-Za-z0-9._-]", "_", ct))))

  tibble(
    predicted.celltype = ct,
    n_samples_High = as.integer(grp_counts["High"]),
    n_samples_Low  = as.integer(grp_counts["Low"]),
    n_sig_q10      = sum(tt$sig_q10, na.rm = TRUE),
    n_sig_q05      = sum(tt$sig_q05, na.rm = TRUE),
    n_up_q10       = sum(tt$sig_q10 & tt$logFC > 0, na.rm = TRUE),
    n_down_q10     = sum(tt$sig_q10 & tt$logFC < 0, na.rm = TRUE),
    used_dupCor    = corval > 0
  )
}

cts <- sort(unique(pb_meta$predicted.celltype))
summ <- bind_rows(lapply(cts, run_one_celltype))
write_tsv(summ, file.path(OUTDIR_DE, "pb_DE_summary_by_celltype.tsv"))

# merged table for quick scan
all_tables <- list.files(OUTDIR_DE, pattern = "^pb_DE_.*\\.tsv$", full.names = TRUE)
if (length(all_tables)) {
  merged <- rbindlist(lapply(all_tables, fread), fill = TRUE)
  fwrite(merged, file.path(OUTDIR_DE, "pb_DE_ALL.tsv"), sep = "\t")
}

message("DE complete. Inspect:\n - ", file.path(OUTDIR_DE, "pseudobulk_sample_counts_by_celltype.tsv"),
        "\n - ", file.path(OUTDIR_DE, "pb_DE_summary_by_celltype.tsv"))

if (!RUN_ENRICHMENT) quit(save = "no")

# ====================== ENRICHMENT (ORA) ======================
# helpers
subset_enrich_by_description <- function(ego, terms_desc) {
  if (is.null(ego)) return(NULL)
  df <- as.data.frame(ego)
  idx <- which(tolower(df$Description) %in% tolower(terms_desc))
  if (!length(idx)) return(NULL)
  ego2 <- ego; ego2@result <- ego@result[idx, , drop = FALSE]; ego2
}

select_terms <- function(ego, pref = IMMUNE_CATEGORIES, n_fallback = 6) {
  if (is.null(ego)) return(character(0))
  df <- as.data.frame(ego); df$lo <- tolower(df$Description)
  pref_lo <- tolower(pref)
  hit <- df$Description[match(pref_lo, df$lo, nomatch = 0)]
  hit <- unique(hit[hit != ""])
  if (length(hit) > 0) return(head(hit, n_fallback))
  keep <- grepl("immune|inflamm|cytokine|chemokin|leukocyte|neutrophil|interferon|oxidative|ros|xenobiotic|toxic|nf-?kappa?b|mapk|autophag|phagocyt|hypoxia",
                df$lo, perl = TRUE)
  cand <- df$Description[keep]
  if (length(cand) >= n_fallback) return(head(cand, n_fallback))
  head(df$Description, n = min(n_fallback, nrow(df)))
}

run_go_bp <- function(genes_sym, universe_sym) {
  if (length(genes_sym) < 10) return(NULL)
  mm     <- bitr(genes_sym,   fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db) %>% distinct(SYMBOL, .keep_all = TRUE)
  mm_uni <- bitr(universe_sym,fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db) %>% distinct(SYMBOL, .keep_all = TRUE)
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
  mm     <- bitr(genes_sym,   fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db) %>% distinct(SYMBOL, .keep_all = TRUE)
  mm_uni <- bitr(universe_sym,fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db) %>% distinct(SYMBOL, .keep_all = TRUE)
  if (nrow(mm) < 10 || nrow(mm_uni) < 100) return(NULL)
  H <- msigdbr(species = "Homo sapiens", category = "H") %>% select(gs_name, entrez_gene)
  enricher(mm$ENTREZID, TERM2GENE = H, universe = mm_uni$ENTREZID, pAdjustMethod = "BH")
}

# choose CTs to enrich
summ <- readr::read_tsv(file.path(OUTDIR_DE, "pb_DE_summary_by_celltype.tsv"), show_col_types = FALSE)
cts_run <- summ %>%
  mutate(total_sig = coalesce(n_sig_q10, 0L)) %>%
  filter(total_sig >= ENR_MIN_SIG_GENES) %>%
  arrange(desc(total_sig)) %>%
  pull(predicted.celltype)

if (!length(cts_run)) {
  message("No CTs with ≥ ", ENR_MIN_SIG_GENES, " significant genes (q<", Q_CUTOFF, "). Consider lowering threshold or use ranked GSEA.")
  quit(save = "no")
}

for (ct in cts_run) {
  message("=== Enrichment: ", ct, " ===")
  de_fp <- file.path(OUTDIR_DE, sprintf("pb_DE_%s.tsv", gsub("[^A-Za-z0-9._-]", "_", ct)))
  if (!file.exists(de_fp)) next
  de <- readr::read_tsv(de_fp, show_col_types = FALSE)

  universe <- unique(de$gene)
  up   <- de %>% filter(sig_q10, logFC >  0) %>% pull(gene) %>% setdiff(Y_GENES)
  down <- de %>% filter(sig_q10, logFC <= 0) %>% pull(gene) %>% setdiff(Y_GENES)

  writeLines(up,   file.path(OUTDIR_ENR, sprintf("genes_UP_q10_noY_%s.txt",   gsub("[^A-Za-z0-9._-]", "_", ct))))
  writeLines(down, file.path(OUTDIR_ENR, sprintf("genes_DOWN_q10_noY_%s.txt", gsub("[^A-Za-z0-9._-]", "_", ct))))

  # GO-BP
  ego_up   <- run_go_bp(up,   universe)
  ego_down <- run_go_bp(down, universe)

  if (!is.null(ego_up)) {
    write.csv(as.data.frame(ego_up), file.path(OUTDIR_ENR, sprintf("GO_BP_UP_%s.csv", ct)), row.names = FALSE)
    ggsave(file.path(OUTDIR_ENR, sprintf("GO_BP_UP_%s.png", ct)),
           dotplot(ego_up, showCategory = 20) + ggtitle(sprintf("GO BP (UP) — %s", ct)),
           width = 9, height = 7, dpi = 300)

    cats <- select_terms(ego_up, IMMUNE_CATEGORIES, n_fallback = 6)
    ego_up_sub <- subset_enrich_by_description(ego_up, cats)
    if (!is.null(ego_up_sub)) {
      fc_named <- setNames(de$logFC, de$gene)
      ggsave(file.path(OUTDIR_ENR, sprintf("CNET_GO_BP_UP_%s.png", ct)),
             cnetplot(ego_up_sub, circular = TRUE, colorEdge = TRUE,
                      foldChange = fc_named, showCategory = nrow(as.data.frame(ego_up_sub))),
             width = 10, height = 9, dpi = 300)
      writeLines(as.character(cats), file.path(OUTDIR_ENR, sprintf("CNET_terms_used_UP_%s.txt", ct)))
    }
  }

  if (!is.null(ego_down)) {
    write.csv(as.data.frame(ego_down), file.path(OUTDIR_ENR, sprintf("GO_BP_DOWN_%s.csv", ct)), row.names = FALSE)
    ggsave(file.path(OUTDIR_ENR, sprintf("GO_BP_DOWN_%s.png", ct)),
           dotplot(ego_down, showCategory = 20) + ggtitle(sprintf("GO BP (DOWN) — %s", ct)),
           width = 9, height = 7, dpi = 300)

    cats <- select_terms(ego_down, IMMUNE_CATEGORIES, n_fallback = 6)
    ego_down_sub <- subset_enrich_by_description(ego_down, cats)
    if (!is.null(ego_down_sub)) {
      fc_named <- setNames(de$logFC, de$gene)
      ggsave(file.path(OUTDIR_ENR, sprintf("CNET_GO_BP_DOWN_%s.png", ct)),
             cnetplot(ego_down_sub, circular = TRUE, colorEdge = TRUE,
                      foldChange = fc_named, showCategory = nrow(as.data.frame(ego_down_sub))),
             width = 10, height = 9, dpi = 300)
      writeLines(as.character(cats), file.path(OUTDIR_ENR, sprintf("CNET_terms_used_DOWN_%s.txt", ct)))
    }
  }

  # Hallmark overview
  hm_up   <- run_hallmark(up,   universe)
  hm_down <- run_hallmark(down, universe)
  if (!is.null(hm_up) && nrow(as.data.frame(hm_up)) > 0) {
    write.csv(as.data.frame(hm_up), file.path(OUTDIR_ENR, sprintf("MSIG_HALLMARK_UP_%s.csv", ct)), row.names = FALSE)
    ggsave(file.path(OUTDIR_ENR, sprintf("MSIG_HALLMARK_UP_%s.png", ct)),
           dotplot(hm_up, showCategory = 15) + ggtitle(sprintf("Hallmark (UP) — %s", ct)),
           width = 9, height = 7, dpi = 300)
  }
  if (!is.null(hm_down) && nrow(as.data.frame(hm_down)) > 0) {
    write.csv(as.data.frame(hm_down), file.path(OUTDIR_ENR, sprintf("MSIG_HALLMARK_DOWN_%s.csv", ct)), row.names = FALSE)
    ggsave(file.path(OUTDIR_ENR, sprintf("MSIG_HALLMARK_DOWN_%s.png", ct)),
           dotplot(hm_down, showCategory = 15) + ggtitle(sprintf("Hallmark (DOWN) — %s", ct)),
           width = 9, height = 7, dpi = 300)
  }
}

message("Enrichment complete. Outputs in: ", OUTDIR_ENR)
# 02_cp_enrichment_from_pb_SCT_fix_select.R
suppressPackageStartupMessages({
  library(tidyverse)    # we'll fully-qualify: dplyr::, readr::, etc.
  library(data.table)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(msigdbr)
})

INDIR  <- "pb_Preterm_MNPBurden_limma_SCT"
OUTDIR <- file.path(INDIR, "cp_results")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# --- Immunomod categories for cnet focus ---
IMMUNE_CATEGORIES <- c(
  "cytokine-mediated signaling pathway",
  "regulation of innate immune response",
  "inflammatory response",
  "leukocyte mediated immunity",
  "response to interferon-gamma", "response to interferon-alpha",
  "NF-kappaB signaling", "MAPK cascade",
  "response to oxidative stress", "reactive oxygen species metabolic process",
  "xenobiotic metabolic process", "cellular response to toxic substance",
  "phagocytosis", "autophagy", "hypoxia"
)

# --- Drop likely sex-mix genes from enrichment sets ---
Y_GENES <- c("DDX3Y","KDM5D","RPS4Y1","EIF1AY","UTY","USP9Y","ZFY","TXLNGY","TMSB4Y",
             "PRKY","PCDH11Y","NLGN4Y","AMELY")

# ---------- helpers ----------
subset_enrich_by_description <- function(ego, terms_desc) {
  if (is.null(ego)) return(NULL)
  df  <- as.data.frame(ego)
  idx <- which(tolower(df$Description) %in% tolower(terms_desc))
  if (!length(idx)) return(NULL)
  ego2 <- ego
  ego2@result <- ego@result[idx, , drop = FALSE]
  ego2
}

select_terms <- function(ego, pref = IMMUNE_CATEGORIES, n_fallback = 6) {
  if (is.null(ego)) return(character(0))
  df <- as.data.frame(ego)
  df$lo <- tolower(df$Description)
  pref_lo <- tolower(pref)
  hit <- df$Description[match(pref_lo, df$lo, nomatch = 0)]
  hit <- unique(hit[hit != ""])
  if (length(hit) > 0) return(head(hit, n_fallback))
  # fallback to immune-ish terms
  keep <- grepl("immune|inflamm|cytokine|chemokin|leukocyte|neutrophil|interferon|oxidative|ros|xenobiotic|toxic|nf-?kappa?b|mapk|autophag|phagocyt|hypoxia",
                df$lo, perl = TRUE)
  cand <- df$Description[keep]
  if (length(cand) >= n_fallback) return(head(cand, n_fallback))
  head(df$Description, n = min(n_fallback, nrow(df)))
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

  H <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>%
       dplyr::select(gs_name, entrez_gene)
  enricher(mm$ENTREZID, TERM2GENE = H, universe = mm_uni$ENTREZID, pAdjustMethod = "BH")
}

# ---------- choose CTs to enrich ----------
summ_fp <- file.path(INDIR, "pb_DE_summary_by_celltype.tsv")
if (!file.exists(summ_fp)) stop("Summary file not found: ", summ_fp)
summ <- readr::read_tsv(summ_fp, show_col_types = FALSE)

ENR_MIN_SIG_GENES <- 10
cts_run <- summ %>%
  dplyr::mutate(total_sig = dplyr::coalesce(n_sig_q10, 0L)) %>%
  dplyr::filter(total_sig >= ENR_MIN_SIG_GENES) %>%
  dplyr::arrange(dplyr::desc(total_sig)) %>%
  dplyr::pull(predicted.celltype)

if (!length(cts_run)) {
  message("No cell types with ≥ ", ENR_MIN_SIG_GENES, " sig genes. Consider lowering threshold or using ranked GSEA.")
  quit(save = "no")
}

# ---------- run per CT ----------
for (ct in cts_run) {
  message("=== Enrichment: ", ct, " ===")
  de_fp <- file.path(INDIR, sprintf("pb_DE_%s.tsv", gsub("[^A-Za-z0-9._-]", "_", ct)))
  if (!file.exists(de_fp)) next
  de <- readr::read_tsv(de_fp, show_col_types = FALSE)

  universe <- unique(de$gene)
  up   <- de %>% dplyr::filter(sig_q10, logFC >  0) %>% dplyr::pull(gene) %>% setdiff(Y_GENES)
  down <- de %>% dplyr::filter(sig_q10, logFC <= 0) %>% dplyr::pull(gene) %>% setdiff(Y_GENES)

  # Keep a record of exactly what went into ORA
  writeLines(up,   file.path(OUTDIR, sprintf("genes_UP_q10_noY_%s.txt",   gsub("[^A-Za-z0-9._-]", "_", ct))))
  writeLines(down, file.path(OUTDIR, sprintf("genes_DOWN_q10_noY_%s.txt", gsub("[^A-Za-z0-9._-]", "_", ct))))

  # GO-BP
  ego_up   <- run_go_bp(up,   universe)
  ego_down <- run_go_bp(down, universe)

  if (!is.null(ego_up)) {
    readr::write_csv(as.data.frame(ego_up), file.path(OUTDIR, sprintf("GO_BP_UP_%s.csv", ct)))
    ggplot2::ggsave(file.path(OUTDIR, sprintf("GO_BP_UP_%s.png", ct)),
                    enrichplot::dotplot(ego_up, showCategory = 20) + ggplot2::ggtitle(sprintf("GO BP (UP) — %s", ct)),
                    width = 9, height = 7, dpi = 300)

    cats <- select_terms(ego_up, IMMUNE_CATEGORIES, n_fallback = 6)
    ego_up_sub <- subset_enrich_by_description(ego_up, cats)
    if (!is.null(ego_up_sub)) {
      fc_named <- setNames(de$logFC, de$gene)  # SYMBOL -> logFC
      ggplot2::ggsave(file.path(OUTDIR, sprintf("CNET_GO_BP_UP_%s.png", ct)),
                      enrichplot::cnetplot(ego_up_sub, circular = TRUE, colorEdge = TRUE,
                                           foldChange = fc_named,
                                           showCategory = nrow(as.data.frame(ego_up_sub))),
                      width = 10, height = 9, dpi = 300)
      writeLines(as.character(cats), file.path(OUTDIR, sprintf("CNET_terms_used_UP_%s.txt", ct)))
    }
  }

  if (!is.null(ego_down)) {
    readr::write_csv(as.data.frame(ego_down), file.path(OUTDIR, sprintf("GO_BP_DOWN_%s.csv", ct)))
    ggplot2::ggsave(file.path(OUTDIR, sprintf("GO_BP_DOWN_%s.png", ct)),
                    enrichplot::dotplot(ego_down, showCategory = 20) + ggplot2::ggtitle(sprintf("GO BP (DOWN) — %s", ct)),
                    width = 9, height = 7, dpi = 300)

    cats <- select_terms(ego_down, IMMUNE_CATEGORIES, n_fallback = 6)
    ego_down_sub <- subset_enrich_by_description(ego_down, cats)
    if (!is.null(ego_down_sub)) {
      fc_named <- setNames(de$logFC, de$gene)
      ggplot2::ggsave(file.path(OUTDIR, sprintf("CNET_GO_BP_DOWN_%s.png", ct)),
                      enrichplot::cnetplot(ego_down_sub, circular = TRUE, colorEdge = TRUE,
                                           foldChange = fc_named,
                                           showCategory = nrow(as.data.frame(ego_down_sub))),
                      width = 10, height = 9, dpi = 300)
      writeLines(as.character(cats), file.path(OUTDIR, sprintf("CNET_terms_used_DOWN_%s.txt", ct)))
    }
  }

  # Hallmark overview
  hm_up   <- run_hallmark(up,   universe)
  hm_down <- run_hallmark(down, universe)

  if (!is.null(hm_up) && nrow(as.data.frame(hm_up)) > 0) {
    readr::write_csv(as.data.frame(hm_up), file.path(OUTDIR, sprintf("MSIG_HALLMARK_UP_%s.csv", ct)))
    ggplot2::ggsave(file.path(OUTDIR, sprintf("MSIG_HALLMARK_UP_%s.png", ct)),
                    enrichplot::dotplot(hm_up, showCategory = 15) + ggplot2::ggtitle(sprintf("Hallmark (UP) — %s", ct)),
                    width = 9, height = 11, dpi = 300)
  }
  if (!is.null(hm_down) && nrow(as.data.frame(hm_down)) > 0) {
    readr::write_csv(as.data.frame(hm_down), file.path(OUTDIR, sprintf("MSIG_HALLMARK_DOWN_%s.csv", ct)))
    ggplot2::ggsave(file.path(OUTDIR, sprintf("MSIG_HALLMARK_DOWN_%s.png", ct)),
                    enrichplot::dotplot(hm_down, showCategory = 15) + ggplot2::ggtitle(sprintf("Hallmark (DOWN) — %s", ct)),
                    width = 9, height = 11, dpi = 300)
  }
}

message("Enrichment complete. Outputs in: ", OUTDIR)
# ===== Append-only: make readable copies of your dotplots =====
suppressPackageStartupMessages({ library(ggplot2); library(stringr) })

make_readable_dotplot <- function(ego, title,
                                  n = 20, wrap = 45, base_size = 12,
                                  left_margin_pts = 140,
                                  width = 10, height = 9, outfile) {
  if (is.null(ego)) return(invisible(NULL))
  df <- try(as.data.frame(ego), silent = TRUE)
  if (inherits(df, "try-error") || !nrow(df)) return(invisible(NULL))

  p <- enrichplot::dotplot(ego, showCategory = n, label_format = wrap) +
    ggtitle(title) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = wrap)) +
    theme_minimal(base_size = base_size) +
    theme(
      axis.text.y  = element_text(size = base_size, hjust = 1),
      axis.title.y = element_blank(),
      plot.title   = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      plot.margin  = margin(10, 20, 10, left_margin_pts)
    ) +
    coord_cartesian(clip = "off")

  ggsave(outfile, p, width = width, height = height, dpi = 300)
  invisible(p)
}

# Call this right after your existing ggsave() calls inside the `for (ct in cts_run)` loop:

# GO BP (UP) readable copy
if (exists("ego_up") && !is.null(ego_up) && nrow(as.data.frame(ego_up)) > 0) {
  make_readable_dotplot(
    ego_up,
    title   = sprintf("GO BP (UP) — %s", ct),
    n       = 20, wrap = 45, base_size = 12, left_margin_pts = 140,
    width   = 10, height = 9,
    outfile = file.path(OUTDIR, sprintf("GO_BP_UP_%s.readable.png", ct))
  )
}

# GO BP (DOWN) readable copy
if (exists("ego_down") && !is.null(ego_down) && nrow(as.data.frame(ego_down)) > 0) {
  make_readable_dotplot(
    ego_down,
    title   = sprintf("GO BP (DOWN) — %s", ct),
    n       = 20, wrap = 45, base_size = 12, left_margin_pts = 140,
    width   = 10, height = 9,
    outfile = file.path(OUTDIR, sprintf("GO_BP_DOWN_%s.readable.png", ct))
  )
}

# Hallmark (UP) readable copy
if (exists("hm_up") && !is.null(hm_up) && nrow(as.data.frame(hm_up)) > 0) {
  make_readable_dotplot(
    hm_up,
    title   = sprintf("Hallmark (UP) — %s", ct),
    n       = 15, wrap = 45, base_size = 12, left_margin_pts = 140,
    width   = 10, height = 8,
    outfile = file.path(OUTDIR, sprintf("MSIG_HALLMARK_UP_%s.readable.png", ct))
  )
}

# Hallmark (DOWN) readable copy
if (exists("hm_down") && !is.null(hm_down) && nrow(as.data.frame(hm_down)) > 0) {
  make_readable_dotplot(
    hm_down,
    title   = sprintf("Hallmark (DOWN) — %s", ct),
    n       = 15, wrap = 45, base_size = 12, left_margin_pts = 140,
    width   = 10, height = 8,
    outfile = file.path(OUTDIR, sprintf("MSIG_HALLMARK_DOWN_%s.readable.png", ct))
  )
}




## top 10 

# ======================= POST-HOC ENRICHMENT PLOTS =======================
suppressPackageStartupMessages({
    library(ggplot2)
    library(readr)
    library(stringr)
    library(tools)
})

# where your enrichment CSVs live (change if needed)
DEFAULT_OUTDIR <- file.path("pb_Preterm_MNPBurden_limma_SCT", "cp_results")

# robust parser for GeneRatio (handles "7/23" or numeric)
.gene_ratio_num <- function(gr) {
    if (is.numeric(gr)) return(gr)
    as.numeric(vapply(strsplit(as.character(gr), "/"), function(v) {
        if (length(v) != 2) return(NA_real_)
        suppressWarnings(as.numeric(v[1]) / as.numeric(v[2]))
    }, 0.0))
}

# plot from a data.frame that looks like as.data.frame(enrichResult)
.plot_df_highlight_classic <- function(df, title,
                                       n_show = 10, top_n = 10,
                                       wrap = 45, base_size = 12,
                                       left_margin_pts = 140,
                                       width = 10, height = 7, outfile) {
    if (is.null(df) || !nrow(df)) return(invisible(NULL))
    
    # Ensure required cols exist
    req <- c("Description","p.adjust","Count")
    if (!all(req %in% names(df))) {
        warning("Missing required columns in df: ", paste(setdiff(req, names(df)), collapse=", "))
        return(invisible(NULL))
    }
    
    # Compute numeric GeneRatio
    if ("GeneRatio" %in% names(df)) {
        df$GeneRatio_num <- .gene_ratio_num(df$GeneRatio)
    } else {
        df$GeneRatio_num <- NA_real_
    }
    # Fallbacks if GeneRatio missing or all NA
    if (all(is.na(df$GeneRatio_num))) {
        if ("setSize" %in% names(df)) {
            df$GeneRatio_num <- df$Count / pmax(df$setSize, df$Count)
        } else {
            df$GeneRatio_num <- df$Count / max(df$Count, na.rm = TRUE)
        }
    }
    
    # Keep complete cases
    df <- df[is.finite(df$p.adjust) & is.finite(df$GeneRatio_num), , drop = FALSE]
    if (!nrow(df)) return(invisible(NULL))
    
    # Order and slice
    df <- df[order(df$p.adjust, -df$GeneRatio_num), ]
    df_show <- head(df, n_show)
    
    # mark top_n within shown set
    top_n <- min(top_n, nrow(df_show))
    df_show$is_top <- FALSE
    if (top_n > 0) df_show$is_top[seq_len(top_n)] <- TRUE
    
    # y order: best at top
    df_show$Description <- factor(df_show$Description, levels = rev(df_show$Description))
    
    p <- ggplot(df_show, aes(x = GeneRatio_num, y = Description)) +
        geom_point(aes(size = Count, fill = p.adjust),
                   shape = 21, color = NA, alpha = 0.65) +
        geom_point(data = subset(df_show, is_top),
                   aes(size = Count, fill = p.adjust),
                   shape = 21, color = "black", stroke = 1.1, alpha = 1) +
        scale_fill_continuous(name = "p.adjust", trans = "reverse") +
        scale_size_continuous(name = "Count") +
        scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = wrap)) +
        labs(x = "GeneRatio", y = NULL, title = title) +
        theme_classic(base_size = base_size) +
        theme(
            plot.title   = element_text(hjust = 0.5, face = "bold"),
            axis.text.y  = element_text(size = base_size, hjust = 1),
            plot.margin  = margin(10, 20, 10, left_margin_pts)
        ) +
        coord_cartesian(clip = "off")
    
    ggsave(outfile, p, width = width, height = height, dpi = 300, bg = "white")
    invisible(p)
}

# Main: scan cp_results for CSVs and emit highlighted/classic PNGs
rebuild_enrichment_plots <- function(outdir = DEFAULT_OUTDIR,
                                     n_show = 10, top_n = 10, wrap = 45) {
    if (!dir.exists(outdir)) stop("Directory not found: ", outdir)
    
    # Patterns to find files you already saved
    patterns <- c("^GO_BP_UP_.*\\.csv$", "^GO_BP_DOWN_.*\\.csv$",
                  "^MSIG_HALLMARK_UP_.*\\.csv$", "^MSIG_HALLMARK_DOWN_.*\\.csv$")
    
    files <- unlist(lapply(patterns, function(p) list.files(outdir, pattern = p, full.names = TRUE)))
    if (!length(files)) {
        message("No enrichment CSVs found in: ", outdir)
        return(invisible(NULL))
    }
    
    for (fp in files) {
        df <- try(readr::read_csv(fp, show_col_types = FALSE), silent = TRUE)
        if (inherits(df, "try-error")) next
        
        base <- file_path_sans_ext(basename(fp))
        # Derive a nice title
        dir_tag <- if (grepl("^GO_BP_", base)) "GO BP" else "Hallmark"
        updown  <- if (grepl("_UP_", base)) "UP" else "DOWN"
        ct      <- sub("^(GO_BP|MSIG_HALLMARK)_(UP|DOWN)_", "", base)
        title   <- sprintf("%s (%s) — %s (Top 10 highlighted)", dir_tag, updown, ct)
        
        out_png <- file.path(outdir, paste0(base, ".posthoc.highlightTop10.classic.png"))
        .plot_df_highlight_classic(df, title, n_show = n_show, top_n = top_n,
                                   wrap = wrap, outfile = out_png)
    }
    
    message("Done. Wrote highlighted plots next to the CSVs in: ", outdir)
}

# ======================= USAGE =======================
# If your working directory is the project root and you used my earlier script:
# rebuild_enrichment_plots()
# Otherwise, point to your results folder explicitly, e.g.:
# rebuild_enrichment_plots(outdir = "/path/to/pb_Preterm_MNPBurden_limma_SCT/cp_results")
rebuild_enrichment_plots(outdir = "/home/ebarrozo/visiumhd_mnp/results/seurat_all_v1/seurat_all_ref_merged/de_MNPBurden_Delivery_class_bulk_sensitive_go/pb_Preterm_MNPBurden_limma_SCT/cp_results")


## now bin by family
# ============== Bin cell types into families + meta-FDR + Top10 plots ==============
suppressPackageStartupMessages({
    library(dplyr); library(readr); library(stringr); library(ggplot2); library(tools)
})

OUTDIR <- file.path("pb_Preterm_MNPBurden_limma_SCT","cp_results")
if (!dir.exists(OUTDIR)) stop("Directory not found: ", OUTDIR)
PADJ_CUT <- 0.05  # significance threshold on family-wise FDR

# ----- helpers -------------------------------------------------------------------
read_if <- function(fp) if (file.exists(fp)) suppressMessages(readr::read_csv(fp, show_col_types = FALSE)) else NULL

gene_ratio_num <- function(gr) {
    if (is.numeric(gr)) return(gr)
    as.numeric(vapply(strsplit(as.character(gr), "/"),
                      function(v) if (length(v)==2) as.numeric(v[1])/as.numeric(v[2]) else NA_real_,
                      0.0))
}

# cell-type -> family bins (edit if you want different bins)
map_cell_family <- function(ct) {
    if (ct %in% c("Macrophage","Hofbauer Cell")) return("Macrophage family")
    if (ct %in% c("T Cell","Treg"))              return("T cell family")
    return(ct)  # leave others as-is
}

# meta-analyze p-values via Fisher's method
fisher_meta <- function(pvals) {
    pvals <- pvals[is.finite(pvals) & pvals > 0 & pvals <= 1]
    if (!length(pvals)) return(NA_real_)
    stat <- -2 * sum(log(pvals))
    k <- length(pvals)
    pchisq(stat, df = 2 * k, lower.tail = FALSE)
}

# reorder within facet; and scale for y labels
reorder_within <- function(x, by, within, fun = mean, sep = "___") {
    x2 <- paste(x, within, sep = sep)
    stats::reorder(x2, by, FUN = fun)
}
scale_y_reordered <- function(..., sep = "___") {
    ggplot2::scale_y_discrete(labels = function(x) gsub(paste0(sep, ".*$"), "", x), ...)
}

# ----- load all enrichment CSVs ---------------------------------------------------
pats <- c("^GO_BP_UP_.*\\.csv$", "^GO_BP_DOWN_.*\\.csv$",
          "^MSIG_HALLMARK_UP_.*\\.csv$", "^MSIG_HALLMARK_DOWN_.*\\.csv$")
files <- unlist(lapply(pats, function(p) list.files(OUTDIR, pattern = p, full.names = TRUE)))
if (!length(files)) stop("No enrichment CSVs found in: ", OUTDIR)

enr <- lapply(files, function(fp) {
    df <- read_if(fp); if (is.null(df) || !nrow(df)) return(NULL)
    base <- file_path_sans_ext(basename(fp))
    pathway_family <- if (startsWith(base, "GO_BP")) "GO_BP" else "HALLMARK"
    direction <- if (grepl("_UP_", base)) "UP" else "DOWN"
    ct <- sub("^(GO_BP|MSIG_HALLMARK)_(UP|DOWN)_", "", base)
    
    df %>%
        mutate(
            pathway_family   = pathway_family,
            direction        = direction,
            predicted.celltype = ct
        )
})
enr <- bind_rows(enr)
if (!nrow(enr)) stop("Parsed zero rows from enrichment CSVs.")

# robust GeneRatio numeric
enr$GeneRatio_num <- if ("GeneRatio" %in% names(enr)) gene_ratio_num(enr$GeneRatio) else NA_real_
if (all(is.na(enr$GeneRatio_num))) {
    enr$GeneRatio_num <- if ("setSize" %in% names(enr)) enr$Count / pmax(enr$setSize, enr$Count)
    else enr$Count / max(enr$Count, na.rm = TRUE)
}

# ensure pvalue column exists
if (!("pvalue" %in% names(enr))) stop("Input CSVs are missing column 'pvalue' from clusterProfiler output.")

# ----- BIN: cell-type -> cell-family ---------------------------------------------
enr <- enr %>%
    mutate(cell_family = vapply(predicted.celltype, map_cell_family, character(1)))

# ----- META-ANALYSIS: combine within (pathway_family × direction × cell_family × term)
meta <- enr %>%
    group_by(pathway_family, direction, cell_family, Description) %>%
    summarise(
        p_meta       = fisher_meta(pvalue),
        Count_sum    = sum(Count, na.rm = TRUE),
        GeneRatio_max= max(GeneRatio_num, na.rm = TRUE),
        n_celltypes  = n_distinct(predicted.celltype),
        .groups = "drop"
    ) %>%
    filter(is.finite(p_meta))

# ----- FDR within (pathway_family × direction × cell_family)
meta <- meta %>%
    group_by(pathway_family, direction, cell_family) %>%
    mutate(p_adj_family = p.adjust(p_meta, method = "BH"),
           sig_family   = p_adj_family <= PADJ_CUT) %>%
    ungroup()

# write combined TSV
out_tsv <- file.path(OUTDIR, "enrichment_binned.cellfamily.tsv")
readr::write_tsv(meta, out_tsv)
message("Wrote: ", out_tsv)

# ----- PLOTS: Top-10 per pathway-family for Macrophage and T cell families -------
plot_top10_cellfamily <- function(df_cf, cf_label, out_stub) {
    df_cf <- df_cf %>% filter(sig_family)
    if (!nrow(df_cf)) { message("No significant pathways for ", cf_label); return(invisible(NULL)) }
    
    # pick top 10 per pathway_family by FDR then GeneRatio
    top_by_pf <- df_cf %>%
        arrange(p_adj_family, desc(GeneRatio_max)) %>%
        group_by(pathway_family) %>%
        slice_head(n = 10) %>%
        ungroup() %>%
        mutate(y_lab = reorder_within(Description, -p_adj_family, pathway_family))
    
    p <- ggplot(top_by_pf,
                aes(x = GeneRatio_max, y = y_lab)) +
        geom_point(aes(size = Count_sum, fill = p_adj_family, shape = direction),
                   color = "black", stroke = 0.6) +
        facet_wrap(~ pathway_family, scales = "free_y") +
        scale_y_reordered() +
        scale_fill_continuous(name = "FDR (family)", trans = "reverse") +
        scale_size_continuous(name = "Count (sum)") +
        scale_shape_manual(values = c(21, 24)) +
        labs(x = "GeneRatio (max across CTs in family)", y = NULL,
             title = sprintf("%s — Top 10 per pathway family (family-wise FDR ≤ %.2f)", cf_label, PADJ_CUT)) +
        theme_classic(base_size = 12) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              axis.text.y = element_text(hjust = 1),
              legend.position = "right",
              plot.margin = margin(10, 20, 10, 140)) +
        coord_cartesian(clip = "off")
    
    ggsave(file.path(OUTDIR, paste0(out_stub, ".Top10.byPathwayFamily.classic.png")),
           p, width = 18, height = 5, dpi = 300, bg = "white")
    #ggsave(file.path(OUTDIR, paste0(out_stub, ".Top10.byPathwayFamily.classic.pdf")),
    #       p, width = 12, height = 9, device = cairo_pdf, bg = "white")
    
    invisible(p)
}

# Macrophage family
mac_meta <- meta %>% filter(cell_family == "Macrophage family")
plot_top10_cellfamily(mac_meta, "Macrophage family", "Macrophage_family")

# T cell family
t_meta <- meta %>% filter(cell_family == "T cell family")
plot_top10_cellfamily(t_meta, "T cell family", "Tcell_family")

# ----- Totals: print overall + per-family ----------------------------------------
#totals_by_cf <- meta %>% filter(sig_family) %>%
#    count(cell_family, name = "n_sig_pathways") %>%
#    arrange(desc(n_sig_pathways))
#print(totals_by_cf, n = Inf)

total_sig <- sum(totals_by_cf$n_sig_pathways)
n_mac <- totals_by_cf$n_sig_pathways[totals_by_cf$cell_family == "Macrophage family"]; if (!length(n_mac)) n_mac <- 0L
n_t   <- totals_by_cf$n_sig_pathways[totals_by_cf$cell_family == "T cell family"];   if (!length(n_t)) n_t <- 0L

cat(sprintf("Family-binned meta-analysis (Fisher) with BH per pathway family × direction × cell family, q ≤ %.2f:\n", PADJ_CUT))
cat(sprintf("  Total significant pathways: %d\n", total_sig))
cat(sprintf("  Macrophage family: %d\n", n_mac))
cat(sprintf("  T cell family: %d\n", n_t))

suppressPackageStartupMessages({ library(dplyr); library(tibble) })

# keep your PADJ_CUT; if not defined, set it here
if (!exists("PADJ_CUT")) PADJ_CUT <- 0.10

# make sure meta is a plain tibble with atomic columns
meta_tbl <- meta %>%
  tibble::as_tibble() %>%
  dplyr::mutate(cell_family = as.character(cell_family))

# totals by cell-family using dplyr::count (avoids plyr collision)
totals_by_cf <- meta_tbl %>%
  dplyr::filter(sig_family) %>%
  dplyr::count(cell_family, name = "n_sig_pathways") %>%
  dplyr::arrange(dplyr::desc(n_sig_pathways))

print(totals_by_cf, n = Inf)

# overall + per-family counts
total_sig <- sum(totals_by_cf$n_sig_pathways)
n_mac <- totals_by_cf$n_sig_pathways[totals_by_cf$cell_family == "Macrophage family"]; if (!length(n_mac)) n_mac <- 0L
n_t   <- totals_by_cf$n_sig_pathways[totals_by_cf$cell_family == "T cell family"];   if (!length(n_t)) n_t <- 0L

cat(sprintf(
  "Family-binned meta-analysis (Fisher) with BH per pathway family × direction × cell family, q ≤ %.2f:\n",
  PADJ_CUT
))
cat(sprintf("  Total significant pathways: %d\n", total_sig))
cat(sprintf("  Macrophage family: %d\n", n_mac))
cat(sprintf("  T cell family: %d\n", n_t))

# one-liner you can paste into the text
cat(sprintf(
  "High MNP-burden differentially expressed genes were associated with %d significant pathways (family-wise FDR ≤ %.2f), highlighting altered immunoregulation in macrophage-lineage cells.\n",
  n_mac, PADJ_CUT
))
#### all celltypes by family

# ===================== Cell-type → Family table (standalone) =====================
suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(readr); library(stringr)
})

# Where your enrichment outputs live (adjust if needed)
OUTDIR <- file.path("pb_Preterm_MNPBurden_limma_SCT", "cp_results")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# --- helper: robust char coercion (avoids list-cols) ---
vecify_chr <- function(x) {
  if (is.list(x)) vapply(x, function(y) paste(as.character(y), collapse = ";"), character(1))
  else as.character(x)
}

# --- 1) Define cell-type → family bins (edit to taste) ---
map_cell_family <- function(ct) {
  ct <- as.character(ct)
  if (ct %in% c("Macrophage","Hofbauer Cell"))              return("Macrophage family")
  if (ct %in% c("T Cell","Treg"))                           return("T cell family")
  if (ct %in% c("B Cell","Plasma Cell"))                    return("B cell family")
  if (ct %in% c("NK Cell","Decidual NK Cell"))              return("NK cell family")
  if (ct %in% c("ILC3"))                                    return("Innate lymphoid family")
  if (ct %in% c("Syncytiotrophoblast","Villous Cytotrophoblast","Extravillous Trophoblast"))
                                                            return("Trophoblast family")
  if (ct %in% c("Placental Endothelial","Perivascular Cell")) return("Endothelial/perivascular family")
  if (ct %in% c("Fibroblast"))                               return("Stromal fibroblast family")
  if (grepl("lowconf|low_conf", ct, ignore.case = TRUE))     return("Unassigned/low confidence")
  return(ct)  # default: keep as-is
}

# --- 2) Pull all cell types from the Seurat object and count cells ---
stopifnot(exists("obj"), !is.null(obj@meta.data))
md <- obj@meta.data %>% tibble::as_tibble()
if (!"predicted.celltype" %in% colnames(md)) stop("predicted.celltype not found in obj@meta.data")

cell_counts <- md %>%
  dplyr::mutate(predicted.celltype = vecify_chr(.data$predicted.celltype)) %>%
  dplyr::count(predicted.celltype, name = "n_cells") %>%
  dplyr::mutate(cell_family = vapply(predicted.celltype, map_cell_family, character(1))) %>%
  dplyr::group_by(cell_family) %>%
  dplyr::mutate(family_total_cells = sum(n_cells)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(frac_of_all_cells = n_cells / sum(n_cells)) %>%
  dplyr::arrange(cell_family, dplyr::desc(n_cells))

# --- 3) If enrichment present, add # significant pathways per family ------------
# (Uses either in-memory 'meta' from my previous binning script, or reads file)
PADJ_CUT <- 0.10
sig_by_family <- NULL

if (exists("meta") && all(c("cell_family","sig_family") %in% colnames(meta))) {
  sig_by_family <- meta %>%
    dplyr::mutate(cell_family = vecify_chr(.data$cell_family)) %>%
    dplyr::filter(sig_family) %>%
    dplyr::count(cell_family, name = "n_sig_pathways_familyFDR")
} else {
  meta_fp <- file.path(OUTDIR, "enrichment_binned.cellfamily.tsv")
  if (file.exists(meta_fp)) {
    meta_df <- suppressMessages(readr::read_tsv(meta_fp, show_col_types = FALSE))
    if (all(c("cell_family","sig_family") %in% colnames(meta_df))) {
      sig_by_family <- meta_df %>%
        dplyr::mutate(cell_family = vecify_chr(.data$cell_family)) %>%
        dplyr::filter(sig_family) %>%
        dplyr::count(cell_family, name = "n_sig_pathways_familyFDR")
    }
  }
}

# --- 4) Build final table and write to disk --------------------------------------
ct_family_tbl <- cell_counts %>%
  dplyr::left_join(sig_by_family %||% tibble(cell_family = character(), n_sig_pathways_familyFDR = integer()),
                   by = "cell_family") %>%
  dplyr::mutate(
    frac_of_all_cells = round(frac_of_all_cells, 4),
    n_sig_pathways_familyFDR = dplyr::coalesce(n_sig_pathways_familyFDR, 0L)
  ) %>%
  dplyr::select(predicted.celltype, cell_family, n_cells, family_total_cells,
                frac_of_all_cells, n_sig_pathways_familyFDR)

readr::write_tsv(ct_family_tbl, file.path(OUTDIR, "celltype_family_table.tsv"))

# Print to console
print(ct_family_tbl, n = Inf)
cat("\nWrote: ", file.path(OUTDIR, "celltype_family_table.tsv"), "\n", sep = "")
## now in term

