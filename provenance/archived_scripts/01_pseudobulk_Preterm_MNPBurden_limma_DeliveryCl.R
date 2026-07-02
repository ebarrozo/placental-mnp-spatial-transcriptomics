# 01d_pseudobulk_Preterm_MNPBurden_limma_SCT_data.R
suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr); library(tidyr); library(purrr); library(stringr)
    library(readr); library(data.table)
    library(limma)
})

OUTDIR <- "pb_Preterm_MNPBurden_limma_SCT"
dir.create(OUTDIR, showWarnings = FALSE)

stopifnot(!is.null(obj@meta.data))
md <- obj@meta.data

# --- Use Delivery_class exactly (Preterm/Term), subset Preterm ---
if (!"Delivery_class" %in% colnames(md)) {
    if ("MNPBurden_Delivery_class" %in% colnames(md)) {
        md$Delivery_class <- sub("^.*-(Preterm|Term)$", "\\1", md$MNPBurden_Delivery_class)
    } else stop("Delivery_class not found.")
}

# --- Normalize MNPBurden to High/Low only ---
if (!"MNPBurden" %in% colnames(md)) {
    if ("MNPBurden_Delivery_class" %in% colnames(md)) {
        md$MNPBurden <- sub("^(High|Low)-.*$", "\\1", md$MNPBurden_Delivery_class)
    } else stop("Need MNPBurden or MNPBurden_Delivery_class.")
} else {
    md$MNPBurden <- sub("^(High|Low).*", "\\1", md$MNPBurden)
}

# Push back
obj@meta.data <- md

req <- c("predicted.celltype","MNPBurden","Delivery_class","study_id")
miss <- setdiff(req, colnames(obj@meta.data))
if (length(miss)) stop("Missing meta cols: ", paste(miss, collapse=", "))

# --- Subset to Preterm ---
obj_ptb <- subset(obj, subset = Delivery_class == "Preterm")
if (ncol(obj_ptb) == 0) stop("No Preterm cells.")

# --- AggregateExpression on SCT/data ---
DefaultAssay(obj_ptb) <- if ("SCT" %in% names(obj_ptb@assays)) "SCT" else DefaultAssay(obj_ptb)
bulk <- AggregateExpression(
    obj_ptb, slot = "data",
    group.by = c("predicted.celltype", "MNPBurden", "study_id"),
    return.seurat = TRUE
)

assay_name <- if ("SCT" %in% names(bulk@assays)) "SCT" else DefaultAssay(bulk)
DefaultAssay(bulk) <- assay_name

# --- Build pseudobulk meta for aggregated columns ---
pb_cells <- colnames(bulk)
pb_meta <- tibble(pb_id = pb_cells)

if (all(c("predicted.celltype","MNPBurden","study_id") %in% colnames(bulk@meta.data))) {
    pb_meta <- bind_cols(pb_meta, bulk@meta.data[pb_cells, c("predicted.celltype","MNPBurden","study_id"), drop=FALSE] %>% as_tibble())
} else {
    # Parse "predicted.celltype_MNPBurden_study_id" but allow underscores in CT
    parts <- strsplit(pb_meta$pb_id, "_")
    pb_meta <- pb_meta %>%
        mutate(
            n = lengths(parts),
            study_id = map2_chr(parts, n, ~ .x[.y]),
            MNPBurden = map2_chr(parts, n, ~ .x[.y - 1]),
            predicted.celltype = map_chr(parts, ~ if (length(.x) > 2) paste(.x[1:(length(.x)-2)], collapse="_") else .x[1])
        ) %>% select(-n)
}

# Normalize any lingering 'High-Preterm' etc.
pb_meta$MNPBurden <- sub("^(High|Low).*", "\\1", pb_meta$MNPBurden)

# --- Contributing cell counts per pseudo-sample (weights) ---
cell_counts <- obj_ptb@meta.data %>%
    as_tibble(rownames = "cell") %>%
    mutate(MNPBurden = sub("^(High|Low).*", "\\1", MNPBurden)) %>%
    count(predicted.celltype, MNPBurden, study_id, name = "n_cells")

pb_meta <- pb_meta %>% left_join(cell_counts, by = c("predicted.celltype","MNPBurden","study_id"))
if (!"n_cells" %in% colnames(pb_meta)) pb_meta$n_cells <- NA_integer_

# QC table
avail_tbl <- pb_meta %>%
    count(predicted.celltype, MNPBurden, name = "n_pseudobulk") %>%
    arrange(predicted.celltype, MNPBurden)
write_tsv(avail_tbl, file.path(OUTDIR, "pseudobulk_sample_counts_by_celltype.tsv"))

# --- Expression matrix (aggregated log-normalized means) ---
E <- as.matrix(GetAssayData(bulk, slot = "data"))  # genes x pseudo-samples
stopifnot(ncol(E) == nrow(pb_meta))

# --- Gene sanitization (drop lncRNA/ENSG/pathogens) ---
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
keep <- rownames(E)[ !grepl(lnc_pattern, rownames(E), perl=TRUE) & !(rownames(E) %in% pathogen.genes) ]
E <- E[keep, , drop=FALSE]

# --- DE params (tweakable) ---
Q_CUTOFF <- 0.10
L2FC_MIN <- 0.25
MIN_SAMPLES_PER_GROUP <- 2
VAR_EPS <- 1e-8   # minimum variance to keep a gene within a CT

safe_dup_cor <- function(expr_mat, design, block) {
    ok <- FALSE; corval <- 0
    if (!is.null(block)) {
        tb <- table(block)
        if (all(tb <= 1)) return(0)  # all blocks size 1
        res <- try(duplicateCorrelation(expr_mat, design, block = block), silent = TRUE)
        if (!inherits(res, "try-error")) {
            ok <- TRUE; corval <- res$consensus
        }
    }
    corval
}

run_one_celltype <- function(ct) {
    sel <- which(pb_meta$predicted.celltype == ct)
    if (length(sel) < 2*MIN_SAMPLES_PER_GROUP) return(NULL)
    
    meta_ct <- pb_meta[sel, , drop=FALSE]
    E_ct    <- E[, sel, drop=FALSE]
    
    # Both groups present with enough samples
    grp_counts <- table(meta_ct$MNPBurden)
    if (!all(c("Low","High") %in% names(grp_counts))) return(NULL)
    if (any(grp_counts[c("Low","High")] < MIN_SAMPLES_PER_GROUP)) return(NULL)
    
    # Filter genes with (near-)zero variance across pseudo-samples
    v <- matrixStats::rowVars(E_ct)
    keep_g <- which(is.finite(v) & v > VAR_EPS)
    if (length(keep_g) < 50) return(NULL)  # skip if too few genes left
    E_ct <- E_ct[keep_g, , drop=FALSE]
    
    # Design: 0+group
    meta_ct$MNPBurden <- factor(meta_ct$MNPBurden, levels=c("Low","High"))
    design <- model.matrix(~ 0 + MNPBurden, data = meta_ct)  # columns: Low, High
    
    # Optional weights by sample precision (more cells -> higher weight)
    w <- meta_ct$n_cells
    if (all(is.na(w))) {
        weights_vec <- NULL
    } else {
        w[is.na(w)] <- median(w, na.rm=TRUE)
        # monotone transform; sqrt is conservative
        weights_vec <- sqrt(w)
    }
    
    # Duplicate correlation only if blocks have >1 samples
    corval <- safe_dup_cor(E_ct, design, block = meta_ct$study_id)
    
    fit0 <- lmFit(E_ct, design, block = if (corval>0) meta_ct$study_id else NULL,
                  correlation = if (corval>0) corval else NULL,
                  weights = weights_vec)
    
    cont <- makeContrasts(High_vs_Low = MNPBurdenHigh - MNPBurdenLow, levels = design)
    fit1 <- contrasts.fit(fit0, cont)
    
    # robust=FALSE avoids fitFDistRobustly failures on low variance sets
    fit2 <- eBayes(fit1, trend = TRUE, robust = FALSE)
    
    tt <- topTable(fit2, number = Inf, sort.by = "P")
    tt$gene <- rownames(tt)
    
    # Add group means (on aggregated log scale)
    m_hi <- rowMeans(E_ct[, meta_ct$MNPBurden=="High", drop=FALSE])
    m_lo <- rowMeans(E_ct[, meta_ct$MNPBurden=="Low",  drop=FALSE])
    tt$mean_High <- m_hi[tt$gene]
    tt$mean_Low  <- m_lo[tt$gene]
    
    tt <- tt %>%
        mutate(
            sig_q10   = adj.P.Val < Q_CUTOFF & abs(logFC) >= L2FC_MIN,
            sig_q05   = adj.P.Val < 0.05     & abs(logFC) >= L2FC_MIN,
            direction = ifelse(logFC > 0, "UP_in_High", "DOWN_in_High"),
            predicted.celltype = ct
        )
    
    out_fp <- file.path(OUTDIR, sprintf("pb_DE_%s.tsv", gsub("[^A-Za-z0-9._-]", "_", ct)))
    write_tsv(tt, out_fp)
    
    write_tsv(tibble(gene = tt$gene[tt$sig_q10 & tt$logFC > 0]),
              file.path(OUTDIR, sprintf("genes_UP_q10_%s.tsv",   gsub("[^A-Za-z0-9._-]", "_", ct))))
    write_tsv(tibble(gene = tt$gene[tt$sig_q10 & tt$logFC < 0]),
              file.path(OUTDIR, sprintf("genes_DOWN_q10_%s.tsv", gsub("[^A-Za-z0-9._-]", "_", ct))))
    
    tibble(
        predicted.celltype = ct,
        n_samples_High = as.integer(grp_counts["High"]),
        n_samples_Low  = as.integer(grp_counts["Low"]),
        n_sig_q10      = sum(tt$sig_q10, na.rm=TRUE),
        n_sig_q05      = sum(tt$sig_q05, na.rm=TRUE),
        n_up_q10       = sum(tt$sig_q10 & tt$logFC > 0, na.rm=TRUE),
        n_down_q10     = sum(tt$sig_q10 & tt$logFC < 0, na.rm=TRUE),
        used_dupCor    = corval > 0
    )
}

cts <- sort(unique(pb_meta$predicted.celltype))
summ <- bind_rows(lapply(cts, run_one_celltype))

write_tsv(summ, file.path(OUTDIR, "pb_DE_summary_by_celltype.tsv"))

# Merge all CT tables for quick scan
all_tables <- list.files(OUTDIR, pattern = "^pb_DE_.*\\.tsv$", full.names = TRUE)
if (length(all_tables)) {
    merged <- rbindlist(lapply(all_tables, fread), fill = TRUE)
    fwrite(merged, file.path(OUTDIR, "pb_DE_ALL.tsv"), sep = "\t")
}

message("Done. Inspect:\n  - ", file.path(OUTDIR, "pseudobulk_sample_counts_by_celltype.tsv"),
        "\n  - ", file.path(OUTDIR, "pb_DE_summary_by_celltype.tsv"))
