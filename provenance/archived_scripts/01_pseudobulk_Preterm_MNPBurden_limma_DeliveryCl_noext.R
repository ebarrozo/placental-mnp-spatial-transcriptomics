# 01_pseudobulk_Preterm_MNPBurden_limma_DeliveryClass.R
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(readr); library(data.table)
  library(limma); library(Matrix)
})

# ----------- Outputs -----------
OUTDIR <- "pb_Preterm_MNPBurden_limma"
dir.create(OUTDIR, showWarnings = FALSE)

# ----------- Meta checks -----------
stopifnot(!is.null(obj@meta.data))
md <- obj@meta.data

# Ensure MNPBurden is present (derive from MNPBurden_Delivery_class if needed)
if (!"MNPBurden" %in% colnames(md)) {
  if (!"MNPBurden_Delivery_class" %in% colnames(md))
    stop("Need either MNPBurden or MNPBurden_Delivery_class in obj@meta.data")
  md$MNPBurden <- sub("^(High|Low)-.*$", "\\1", md$MNPBurden_Delivery_class)
}

# Use Delivery_class exactly as provided (values should be 'Preterm' or 'Term')
if (!"Delivery_class" %in% colnames(md)) {
  # fallback if user kept old combined column
  if ("MNPBurden_Delivery_class" %in% colnames(md)) {
    md$Delivery_class <- sub("^.*-(Preterm|Term)$", "\\1", md$MNPBurden_Delivery_class)
  } else {
    stop("Delivery_class not found in obj@meta.data")
  }
}

# Push back to object
obj@meta.data <- md

req_cols <- c("predicted.celltype","MNPBurden","Delivery_class","study_id")
missing <- setdiff(req_cols, colnames(obj@meta.data))
if (length(missing) > 0) stop("Missing required meta columns: ", paste(missing, collapse = ", "))

# ----------- Subset: PRETERM only -----------
obj_ptb <- subset(obj, subset = Delivery_class == "Preterm")
if (ncol(obj_ptb) == 0) stop("No cells with Delivery_class == 'Preterm'.")

# ----------- Pseudobulk aggregation (your spec) -----------
bulk <- AggregateExpression(
  obj_ptb,
  slot = "data",
  group.by = c("predicted.celltype", "MNPBurden", "study_id"),
  return.seurat = TRUE
)

assay_name <- if ("SCT" %in% names(bulk@assays)) "SCT" else if ("RNA" %in% names(bulk@assays)) "RNA" else names(bulk@assays)[1]
DefaultAssay(bulk) <- assay_name

# Derive metadata for pseudobulk columns
pb_cells <- colnames(bulk)
pb_meta <- tibble(pb_id = pb_cells)

# Use available meta on bulk if present, else parse pb_id safely
if (all(c("predicted.celltype","MNPBurden","study_id") %in% colnames(bulk@meta.data))) {
  pb_meta <- bind_cols(
    pb_meta,
    bulk@meta.data[pb_cells, c("predicted.celltype","MNPBurden","study_id"), drop = FALSE] %>% as_tibble()
  )
} else {
  # Parse by splitting from the end to tolerate underscores in study_id
  parsed <- strsplit(pb_meta$pb_id, "_")
  pb_meta <- pb_meta %>%
    mutate(
      n_tok = lengths(parsed),
      study_id = map2_chr(parsed, n_tok, ~ .x[.y]),
      MNPBurden = map2_chr(parsed, n_tok, ~ .x[.y - 1]),
      predicted.celltype = map_chr(parsed, ~ if (length(.x) > 2) paste(.x[1:(length(.x)-2)], collapse = "_") else .x[1])
    ) %>%
    select(-n_tok)
}

# Count contributing cells per pseudobulk sample (QC/weights if needed)
cell_counts <- obj_ptb@meta.data %>%
  as_tibble(rownames = "cell") %>%
  count(predicted.celltype, MNPBurden, study_id, name = "n_cells")

pb_meta <- pb_meta %>% left_join(cell_counts, by = c("predicted.celltype","MNPBurden","study_id"))

# QC table: how many pseudobulk samples per celltype x burden
avail_tbl <- pb_meta %>%
  count(predicted.celltype, MNPBurden, name = "n_pseudobulk") %>%
  arrange(predicted.celltype, MNPBurden)
write_tsv(avail_tbl, file.path(OUTDIR, "pseudobulk_sample_counts_by_celltype.tsv"))

# ----------- Extract expression matrix (log-normalized means) -----------
E <- GetAssayData(bulk, slot = "data")  # genes x pseudobulk samples
stopifnot(ncol(E) == nrow(pb_meta))

# ----------- Gene sanitization -----------
lnc_pattern <- "^RP[0-9]+-|^AC[0-9]+\\.|^AL[0-9]+\\.|^LINC|^MIR|^SNHG|^HCG|^CT[CDE]-|^ENSG"
pathogen.genes <- c(
  "CMV-UL117","CMV-UL122","CMV-UL123","CMV-UL138","CMV-UL5","CMV-UL54","CMV-US33",
  "GBS-fabD","GBS-fbsB","GBS-hylB","GBS-lmb","GBS-mfd","GBS-pbsP","GBS-scpB",
  "HPV-E1","HPV-E2","HPV-E5","HPV-E6","HPV-L1","HPV-L2",
  "HSV-1-LAT","HSV-1-RL1","HSV-1-RL2","HSV-1-RS1","HSV-1-UL14","HSV-1-UL29","HSV-1-UL48","HSV-1-UL54","HSV-1-US1",
  "HSV-2-LAT","HSV-2-RS1","HSV-2-UL30","HSV-2-UL41","HSV-2-UL52","HSV-2-US1","HSV-2-US3",
  "SARS-CoV-2-M","SARS-CoV-2-N","SARS-CoV-2-ORF10","SARS-CoV-2-ORF1ab","SARS-CoV-2-ORF3a","SARS-CoV-2-ORF8","SARS-CoV-2-S",
  "Syphilis-dnaK","Syphilis-rplC","Syphilis-rplF","Syphilis-rpsE","Syphilis-rpsM","Syphilis-secY","Syphilis-tuf"
)

keep_genes <- rownames(E)[ !grepl(lnc_pattern, rownames(E), perl = TRUE) &
                           !(rownames(E) %in% pathogen.genes) ]
E <- E[keep_genes, , drop = FALSE]

# ----------- DE params -----------
Q_CUTOFF <- 0.10     # FDR (BH) threshold
L2FC_MIN <- 0.25     # absolute log2 fold-change threshold
MIN_SAMPLES_PER_GROUP <- 3  # require ≥ this many pseudo-samples per group

# ----------- DE per cell type: High vs Low (limma + duplicateCorrelation) -----------
run_one_celltype <- function(ct) {
  idx <- which(pb_meta$predicted.celltype == ct)
  if (length(idx) < 4) return(NULL)

  meta_ct <- pb_meta[idx, , drop = FALSE]
  E_ct    <- E[, idx, drop = FALSE]

  # Need both groups with enough samples
  grp_counts <- table(meta_ct$MNPBurden)
  if (!all(c("High","Low") %in% names(grp_counts))) return(NULL)
  if (any(grp_counts[c("High","Low")] < MIN_SAMPLES_PER_GROUP)) return(NULL)

  meta_ct$MNPBurden <- factor(meta_ct$MNPBurden, levels = c("Low","High"))
  design <- model.matrix(~ 0 + MNPBurden, data = meta_ct)  # columns: MNPBurdenLow, MNPBurdenHigh

  # Within-study correlation
  corfit <- try(duplicateCorrelation(E_ct, design, block = meta_ct$study_id), silent = TRUE)
  corval <- if (inherits(corfit, "try-error")) 0 else corfit$consensus

  fit  <- lmFit(E_ct, design, block = meta_ct$study_id, correlation = corval)
  cont <- makeContrasts(High_vs_Low = MNPBurdenHigh - MNPBurdenLow, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, cont), trend = TRUE, robust = TRUE)

  tt <- topTable(fit2, number = Inf, sort.by = "P")
  tt$gene <- rownames(tt)

  # Group means (on log scale used in E_ct)
  m_hi <- rowMeans(E_ct[, meta_ct$MNPBurden == "High", drop = FALSE])
  m_lo <- rowMeans(E_ct[, meta_ct$MNPBurden == "Low",  drop = FALSE])
  tt$mean_High <- m_hi[tt$gene]
  tt$mean_Low  <- m_lo[tt$gene]

  # Flags
  tt <- tt %>%
    mutate(
      sig_q10 = adj.P.Val < Q_CUTOFF & abs(logFC) >= L2FC_MIN,
      sig_q05 = adj.P.Val < 0.05     & abs(logFC) >= L2FC_MIN,
      direction = ifelse(logFC > 0, "UP_in_High", "DOWN_in_High")
    )

  tt$predicted.celltype <- ct

  # Write per-celltype table
  out_fp <- file.path(OUTDIR, sprintf("pb_DE_%s.tsv", gsub("[^A-Za-z0-9._-]", "_", ct)))
  write_tsv(tt, out_fp)

  # Export gene lists for downstream enrichment (q<0.10)
  up   <- tt %>% filter(sig_q10, logFC >  0) %>% pull(gene)
  down <- tt %>% filter(sig_q10, logFC <= 0) %>% pull(gene)
  write_tsv(tibble(gene = up),   file.path(OUTDIR, sprintf("genes_UP_q10_%s.tsv",   gsub("[^A-Za-z0-9._-]", "_", ct))))
  write_tsv(tibble(gene = down), file.path(OUTDIR, sprintf("genes_DOWN_q10_%s.tsv", gsub("[^A-Za-z0-9._-]", "_", ct))))

  tibble(
    predicted.celltype = ct,
    n_samples_High = as.integer(grp_counts["High"]),
    n_samples_Low  = as.integer(grp_counts["Low"]),
    n_sig_q10      = sum(tt$sig_q10, na.rm = TRUE),
    n_sig_q05      = sum(tt$sig_q05, na.rm = TRUE),
    n_up_q10       = sum(tt$sig_q10 & tt$logFC > 0, na.rm = TRUE),
    n_down_q10     = sum(tt$sig_q10 & tt$logFC < 0, na.rm = TRUE),
    corfit = corval
  )
}

cts <- sort(unique(pb_meta$predicted.celltype))
summ_list <- lapply(cts, run_one_celltype)
summ <- bind_rows(summ_list)
write_tsv(summ, file.path(OUTDIR, "pb_DE_summary_by_celltype.tsv"))

# Also bind all per-celltype tables into one file for quick scanning
all_tables <- list.files(OUTDIR, pattern = "^pb_DE_.*\\.tsv$", full.names = TRUE)
if (length(all_tables) > 0) {
  merged <- rbindlist(lapply(all_tables, fread), fill = TRUE)
  fwrite(merged, file.path(OUTDIR, "pb_DE_ALL.tsv"), sep = "\t")
}

# Console preview: which cell types have enough signal for ORA right away?
if (nrow(summ)) {
  message("Top cell types by n_sig_q10 (q<0.10 & |log2FC|≥", L2FC_MIN, "):")
  print(summ %>% arrange(desc(n_sig_q10), desc(n_up_q10 + n_down_q10)) %>% head(10))
}

message("Done. Check: ",
        file.path(OUTDIR, "pseudobulk_sample_counts_by_celltype.tsv"), " and ",
        file.path(OUTDIR, "pb_DE_summary_by_celltype.tsv"))
