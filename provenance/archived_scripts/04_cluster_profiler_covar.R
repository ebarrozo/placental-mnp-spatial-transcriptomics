# 04_cluster_profiler_covar.R
# ===================== Re-aggregate with study_id (not sample_id_std) =====================
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tidyr); library(tibble); library(readr)
  library(limma);  library(matrixStats)
})

# ---- config (keep your existing values; set defaults if missing) ----
if (!exists("TARGET_DELIVERY")) TARGET_DELIVERY <- "Term"  # or "Preterm"
if (!exists("OUTDIR_DE")) OUTDIR_DE <- sprintf("pb_%s_MNPBurden_limma_SCT", TARGET_DELIVERY)
if (!exists("L2FC_MIN"))     L2FC_MIN    <- 0.25
if (!exists("Q_CUTOFF_DE"))  Q_CUTOFF_DE <- 0.10
if (!exists("MIN_SAMPLES"))   MIN_SAMPLES <- 2
if (!exists("VAR_EPS"))       VAR_EPS     <- 1e-8
dir.create(OUTDIR_DE, showWarnings = FALSE, recursive = TRUE)

# ---- helper ----
vecify_chr  <- function(x) if (is.list(x)) vapply(x, \(y) paste(as.character(y), collapse=";"), "") else as.character(x)
clean_piece <- function(x) gsub("\\|", "/", vecify_chr(x))  # forbid '|' inside pieces

# ---- 1) subset to delivery class; clean metadata we need ----
stopifnot(exists("obj"))
obj_use <- subset(obj, subset = Delivery_class == TARGET_DELIVERY)
DefaultAssay(obj_use) <- if ("SCT" %in% names(obj_use@assays)) "SCT" else DefaultAssay(obj_use)

md0 <- obj_use@meta.data %>%
  tibble::as_tibble(rownames = "cell") %>%
  dplyr::mutate(
    predicted.celltype = clean_piece(predicted.celltype),
    MNPBurden          = sub("^(High|Low).*", "\\1", clean_piece(MNPBurden)),
    study_id           = clean_piece(study_id),
    Hypertension_TRUE  = factor(clean_piece(Hypertension_TRUE), levels = c("No","Yes")),
    Prior_PTB_TRUE     = factor(clean_piece(Prior_PTB_TRUE),    levels = c("No","Yes")),
    Twins_TRUE         = factor(clean_piece(Twins_TRUE),        levels = c("No","Yes"))
  )

# keep the cleaned columns in the object (optional)
obj_use@meta.data[colnames(obj_use),
                  c("predicted.celltype","MNPBurden","study_id",
                    "Hypertension_TRUE","Prior_PTB_TRUE","Twins_TRUE")] <-
  md0 %>%
  dplyr::select(predicted.celltype, MNPBurden, study_id,
                Hypertension_TRUE, Prior_PTB_TRUE, Twins_TRUE) %>%
  as.data.frame()

# ---- 2) build a single key using STUDY_ID and aggregate on it ----
obj_use$pb_key <- with(obj_use@meta.data,
                       paste(predicted.celltype, MNPBurden, study_id, sep = "|"))

bulk <- AggregateExpression(
  obj_use, slot = "data",
  group.by = "pb_key",                 # <-- use the single key (with study_id)
  return.seurat = TRUE
)
DefaultAssay(bulk) <- if ("SCT" %in% names(bulk@assays)) "SCT" else DefaultAssay(bulk)

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(tibble); library(readr); library(Seurat) })

# canonize CT labels in cell-level metadata so counts can join later
unify_ct <- function(x) gsub("[-_]+","-", as.character(x))  # always use hyphen form

# 1) identities that Seurat used for columns (lossless truth)
pb_key_bulk <- as.character(Seurat::Idents(bulk))
stopifnot(length(pb_key_bulk) == ncol(bulk))

# 2) parse identities into columns
pb_meta <- tibble(pb_key = pb_key_bulk) %>%
  tidyr::separate(pb_key, into = c("predicted.celltype","MNPBurden","study_id"),
                  sep = "\\|", remove = FALSE)

# 3) donor-level covariates (distinct by study_id), then join
donor_cov <- obj_use@meta.data %>%
  as_tibble(rownames = "cell") %>%
  transmute(
    study_id          = as.character(study_id),
    Hypertension_TRUE = factor(as.character(Hypertension_TRUE), levels = c("No","Yes")),
    Prior_PTB_TRUE    = factor(as.character(Prior_PTB_TRUE),    levels = c("No","Yes")),
    Twins_TRUE        = factor(as.character(Twins_TRUE),        levels = c("No","Yes"))
  ) %>%
  distinct(study_id, .keep_all = TRUE)

pb_meta <- pb_meta %>%
  left_join(donor_cov, by = "study_id")

# sanity: any missing donor covariates?
if (anyNA(pb_meta$Hypertension_TRUE) | anyNA(pb_meta$Prior_PTB_TRUE) | anyNA(pb_meta$Twins_TRUE)) {
  na_rows <- which(is.na(pb_meta$Hypertension_TRUE) | is.na(pb_meta$Prior_PTB_TRUE) | is.na(pb_meta$Twins_TRUE))
  message("Warning: missing donor covariates for some study_id rows: ",
          paste(head(pb_meta$study_id[na_rows], 10), collapse = ", "),
          if (length(na_rows) > 10) " …" else "")
}

# 4) contributing cell counts as weights — group using the SAME label conventions
cell_counts <- obj_use@meta.data %>%
  as_tibble(rownames = "cell") %>%
  transmute(
    predicted.celltype = unify_ct(predicted.celltype),
    MNPBurden          = sub("^(High|Low).*", "\\1", as.character(MNPBurden)),
    study_id           = as.character(study_id)
  ) %>%
  count(predicted.celltype, MNPBurden, study_id, name = "n_cells")

pb_meta <- pb_meta %>%
  mutate(predicted.celltype = unify_ct(predicted.celltype)) %>%
  left_join(cell_counts, by = c("predicted.celltype","MNPBurden","study_id"))

# fill any missing n_cells with 1 (conservative) to keep weights finite later
pb_meta$n_cells[is.na(pb_meta$n_cells)] <- 1L

# 5) (optional) carry metadata on bulk for easy inspection
bulk@meta.data <- pb_meta %>% select(-pb_key) %>% as.data.frame()
rownames(bulk@meta.data) <- colnames(bulk)

# 6) expression matrix aligned 1:1 with pb_meta
E <- as.matrix(GetAssayData(bulk, slot = "data"))
stopifnot(ncol(E) == nrow(pb_meta))

# quick QC: replication per CT × burden
pb_meta %>%
  count(predicted.celltype, MNPBurden, name = "n_pseudobulk") %>%
  arrange(predicted.celltype, MNPBurden) %>%
  write_tsv(file.path(OUTDIR_DE, "pseudobulk_sample_counts_by_celltype.tsv"))


# ---- 5) DE per cell type (adjusting for Hypertension/Prior_PTB/Twins) -------------
run_one_celltype <- function(ct) {
  sel <- which(pb_meta$predicted.celltype == ct)
  if (length(sel) < 2*MIN_SAMPLES) return(NULL)

  meta_ct <- pb_meta[sel, , drop = FALSE]
  E_ct    <- E[, sel, drop = FALSE]

  # group sizes
  tbl <- table(meta_ct$MNPBurden)
  if (!all(c("Low","High") %in% names(tbl))) return(NULL)
  if (any(tbl[c("Low","High")] < MIN_SAMPLES)) return(NULL)

  # gene filter
  v <- matrixStats::rowVars(E_ct)
  keep <- which(is.finite(v) & v > VAR_EPS)
  if (length(keep) < 50) return(NULL)
  E_ct <- E_ct[keep, , drop = FALSE]

  # factors
  meta_ct$MNPBurden <- factor(meta_ct$MNPBurden, levels = c("Low","High"))

  # covariates: drop if constant in this CT
  covars <- c("Hypertension_TRUE","Prior_PTB_TRUE","Twins_TRUE")
  covars <- covars[vapply(meta_ct[covars], function(x) nlevels(factor(x)) > 1, logical(1))]

  # design and fit
  fml <- as.formula(paste0("~ 0 + MNPBurden",
                           if (length(covars)) paste0(" + ", paste(covars, collapse = " + ")) else ""))
  design <- model.matrix(fml, data = meta_ct)
  ne <- limma::nonEstimable(design)
  if (!is.null(ne) && length(ne)) design <- design[, !(colnames(design) %in% ne), drop = FALSE]

  w <- meta_ct$n_cells; w[is.na(w)] <- median(w, na.rm = TRUE); w <- sqrt(w)

  fit0 <- lmFit(E_ct, design, weights = w)
  contr <- makeContrasts(High_vs_Low = MNPBurdenHigh - MNPBurdenLow, levels = design)
  fit1 <- contrasts.fit(fit0, contr)
  fit2 <- eBayes(fit1, trend = TRUE)

  tt <- topTable(fit2, number = Inf, sort.by = "P")
  tt$gene <- rownames(tt)
  # means
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

  out_fp <- file.path(OUTDIR_DE, sprintf("pb_DE_%s.tsv", gsub("[^A-Za-z0-9._-]", "_", ct)))
  readr::write_tsv(tt, out_fp)

  data.frame(
    predicted.celltype = ct,
    n_samples_High = as.integer(tbl["High"]),
    n_samples_Low  = as.integer(tbl["Low"]),
    covars_in_model = paste(covars, collapse = ";"),
    n_sig_q10      = sum(tt$sig_q10, na.rm = TRUE),
    n_sig_q05      = sum(tt$sig_q05, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

cts  <- sort(unique(pb_meta$predicted.celltype))
summ <- dplyr::bind_rows(lapply(cts, run_one_celltype))
readr::write_tsv(summ, file.path(OUTDIR_DE, "pb_DE_summary_by_celltype.tsv"))

print(summ %>% dplyr::arrange(dplyr::desc(n_sig_q10)), n = Inf)
message("Adjusted DE complete (grouped by study_id) → ", OUTDIR_DE)
