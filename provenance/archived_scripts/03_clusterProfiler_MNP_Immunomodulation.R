# 03_clusterProfiler_MNP_Immunomodulation.R
suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(stringr)
  library(clusterProfiler); library(enrichplot)
  library(org.Hs.eg.db); library(msigdbr)
  # Optional, if you want Reactome:
  # library(ReactomePA)
})

# ---------- I/O ----------
IN_DE  <- "de_markers.with_minpct.tsv"  # from your previous step
OUTDIR <- "cp_MNP_immunomod"
dir.create(OUTDIR, showWarnings = FALSE)

# ---------- Parameters ----------
q_cutoff        <- 0.05          # FDR threshold for DE -> enrichment genes
logfc_min       <- 0.25          # min |log2FC| for DE -> enrichment genes
min_genes_enrich<- 10            # skip enrich if fewer genes after filtering
simplify_cutoff <- 0.6           # GO term redundancy pruning
top_show        <- 20            # top terms in dotplots
IMMUNE_CATEGORIES <- c(          # curate cnet “modules” you care about
  "regulation of innate immune response", "cytokine-mediated signaling pathway",
  "leukocyte mediated immunity", "neutrophil activation",
  "inflammatory response", "response to interferon-gamma", "response to interferon-alpha",
  "NF-kappaB signaling", "MAPK cascade",
  "response to oxidative stress", "reactive oxygen species metabolic process",
  "xenobiotic metabolic process", "cellular response to toxic substance",
  "response to hypoxia", "autophagy", "phagocytosis"
)

# Pathogen/viral features to exclude (matches your earlier DE prefiltering)
pathogen.genes <- c("CMV-UL117","CMV-UL122","CMV-UL123","CMV-UL138","CMV-UL5","CMV-UL54","CMV-US33",
                    "GBS-fabD","GBS-fbsB","GBS-hylB","GBS-lmb","GBS-mfd","GBS-pbsP","GBS-scpB",
                    "HPV-E1","HPV-E2","HPV-E5","HPV-E6","HPV-L1","HPV-L2",
                    "HSV-1-LAT","HSV-1-RL1","HSV-1-RL2","HSV-1-RS1","HSV-1-UL14","HSV-1-UL29","HSV-1-UL48","HSV-1-UL54","HSV-1-US1",
                    "HSV-2-LAT","HSV-2-RS1","HSV-2-UL30","HSV-2-UL41","HSV-2-UL52","HSV-2-US1","HSV-2-US3",
                    "SARS-CoV-2-M","SARS-CoV-2-N","SARS-CoV-2-ORF10","SARS-CoV-2-ORF1ab","SARS-CoV-2-ORF3a","SARS-CoV-2-ORF8","SARS-CoV-2-S",
                    "Syphilis-dnaK","Syphilis-rplC","Syphilis-rplF","Syphilis-rpsE","Syphilis-rpsM","Syphilis-secY","Syphilis-tuf")

lnc_pattern <- "^RP[0-9]+-|^AC[0-9]+\\.|^AL[0-9]+\\.|^LINC|^MIR|^SNHG|^HCG|^CT[CDE]-|^ENSG"

# ---------- Helpers ----------
sanitize_genes <- function(sym) {
  sym <- sym[!is.na(sym)]
  sym <- sym[!grepl(lnc_pattern, sym, perl=TRUE)]
  sym <- setdiff(sym, pathogen.genes)
  unique(sym)
}

bitr_drop_dups <- function(genes_symbol) {
  # Map SYMBOL -> ENTREZ, drop SYMBOLs with no map
  mm <- bitr(genes_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  mm <- distinct(mm, SYMBOL, .keep_all = TRUE)
  mm
}

run_go_bp <- function(entrez, universe_entrez) {
  if (length(entrez) < min_genes_enrich) return(NULL)
  ego <- enrichGO(
    gene          = entrez,
    universe      = universe_entrez,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.10,
    readable      = TRUE
  )
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)
  suppressWarnings(simplify(ego, cutoff = simplify_cutoff, by = "p.adjust", select_fun = min))
}

run_hallmark <- function(entrez, universe_entrez) {
  if (length(entrez) < min_genes_enrich) return(NULL)
  msig_h <- msigdbr(species = "Homo sapiens", category = "H") %>%
    select(gs_name, entrez_gene)
  TERM2GENE <- msig_h %>% distinct(gs_name, entrez_gene)
  enricher(entrez, TERM2GENE = TERM2GENE, universe = universe_entrez, pAdjustMethod = "BH")
}

# Optional Reactome
# run_reactome <- function(entrez, universe_entrez) {
#   if (length(entrez) < min_genes_enrich) return(NULL)
#   enrichPathway(gene = entrez, universe = universe_entrez, pAdjustMethod = "BH", readable = TRUE)
# }

# Pick immunomod-focused categories present in the object
select_immuno_terms <- function(egores, preferred = IMMUNE_CATEGORIES, n_fallback = 6) {
  if (is.null(egores)) return(character(0))
  df <- as.data.frame(egores)
  # exact/loose matching
  df$term_lower <- tolower(df$Description)
  pref_lower <- tolower(preferred)
  hit <- df$Description[match(pref_lower, df$term_lower, nomatch = 0)]
  hit <- hit[hit != ""] %>% unique()
  if (length(hit) > 0) return(hit)
  # fallback: pick top inflammation/immune/oxidative/toxic terms
  keep <- grepl("immune|inflamm|cytokine|chemokin|leukocyte|neutrophil|interferon|oxidative|ros|xenobiotic|toxic|nf-?kappa?b|mapk|autophag|phagocyt|hypoxia",
                df$term_lower, perl = TRUE)
  cand <- df$Description[keep]
  if (length(cand) >= n_fallback) return(cand[seq_len(n_fallback)])
  head(df$Description, n = min(n_fallback, nrow(df)))
}

# ---------- Load DE ----------
de <- readr::read_tsv(IN_DE, show_col_types = FALSE)
stopifnot(all(c("gene","cluster","avg_log2FC","p_val","p_val_adj") %in% colnames(de)))

# Universe: all tested symbols after sanitization
universe_symbols <- sanitize_genes(unique(de$gene))
mm_universe <- bitr_drop_dups(universe_symbols)
universe_entrez <- mm_universe$ENTREZID

# Iterate per class (cluster = MNPBurden_Delivery_class)
classes <- unique(as.character(de$cluster))
message("Found classes: ", paste(classes, collapse = ", "))

for (cl in classes) {
  message("=== Running enrichment for class: ", cl, " ===")
  dsub <- de %>% filter(cluster == cl)

  # Keep FDR-sig + effect size
  dsub_clean <- dsub %>%
    filter(!is.na(gene)) %>%
    mutate(gene = as.character(gene)) %>%
    filter(!grepl(lnc_pattern, gene, perl=TRUE)) %>%
    filter(!(gene %in% pathogen.genes))

  up   <- dsub_clean %>% filter(p_val_adj < q_cutoff, avg_log2FC >=  logfc_min)
  down <- dsub_clean %>% filter(p_val_adj < q_cutoff, avg_log2FC <= -logfc_min)

  # Save the filtered DE used for enrichment
  readr::write_tsv(up,   file.path(OUTDIR, sprintf("DE_forEnrich_UP_%s.tsv",   cl)))
  readr::write_tsv(down, file.path(OUTDIR, sprintf("DE_forEnrich_DOWN_%s.tsv", cl)))

  # Map to ENTREZ
  up_map   <- bitr_drop_dups(up$gene)
  down_map <- bitr_drop_dups(down$gene)

  entrez_up   <- up_map$ENTREZID
  entrez_down <- down_map$ENTREZID

  # --- GO:BP
  ego_up   <- run_go_bp(entrez_up,   universe_entrez)
  ego_down <- run_go_bp(entrez_down, universe_entrez)

  if (!is.null(ego_up)) {
    # tables
    write.csv(as.data.frame(ego_up), file.path(OUTDIR, sprintf("GO_BP_UP_%s.csv", cl)), row.names = FALSE)
    # dot
    p_bp_up <- dotplot(ego_up, showCategory = top_show) + ggtitle(sprintf("GO BP (UP) — %s", cl))
    ggsave(file.path(OUTDIR, sprintf("GO_BP_UP_%s.png", cl)), p_bp_up, width = 9, height = 7, dpi = 300)

    # cnet: pick immunomod categories or fallback best
    cats <- select_immuno_terms(ego_up, IMMUNE_CATEGORIES, n_fallback = 6)
    if (length(cats) > 0) {
      # foldChange vector (SYMBOL -> log2FC) for edge coloring
      fc_named <- setNames(up$avg_log2FC, up$gene)
      # restrict ego to chosen categories
      ego_up_sel <- ego_up[cats]
      # cnet
      p_cnet <- cnetplot(ego_up_sel, circular = TRUE, colorEdge = TRUE,
                         foldChange = fc_named, showCategory = cats)
      ggsave(file.path(OUTDIR, sprintf("CNET_GO_BP_UP_%s.png", cl)), p_cnet, width = 10, height = 9, dpi = 300)
    }
  }

  if (!is.null(ego_down)) {
    write.csv(as.data.frame(ego_down), file.path(OUTDIR, sprintf("GO_BP_DOWN_%s.csv", cl)), row.names = FALSE)
    p_bp_dn <- dotplot(ego_down, showCategory = top_show) + ggtitle(sprintf("GO BP (DOWN) — %s", cl))
    ggsave(file.path(OUTDIR, sprintf("GO_BP_DOWN_%s.png", cl)), p_bp_dn, width = 9, height = 7, dpi = 300)

    cats <- select_immuno_terms(ego_down, IMMUNE_CATEGORIES, n_fallback = 6)
    if (length(cats) > 0) {
      fc_named <- setNames(down$avg_log2FC, down$gene)
      ego_dn_sel <- ego_down[cats]
      p_cnet <- cnetplot(ego_dn_sel, circular = TRUE, colorEdge = TRUE,
                         foldChange = fc_named, showCategory = cats)
      ggsave(file.path(OUTDIR, sprintf("CNET_GO_BP_DOWN_%s.png", cl)), p_cnet, width = 10, height = 9, dpi = 300)
    }
  }

  # --- Hallmark MSigDB (immunomod signal panels)
  hm_up   <- run_hallmark(entrez_up,   universe_entrez)
  hm_down <- run_hallmark(entrez_down, universe_entrez)

  if (!is.null(hm_up) && nrow(as.data.frame(hm_up)) > 0) {
    write.csv(as.data.frame(hm_up), file.path(OUTDIR, sprintf("MSIG_HALLMARK_UP_%s.csv", cl)), row.names = FALSE)
    p_hm_up <- dotplot(hm_up, showCategory = 15) + ggtitle(sprintf("Hallmark (UP) — %s", cl))
    ggsave(file.path(OUTDIR, sprintf("MSIG_HALLMARK_UP_%s.png", cl)), p_hm_up, width = 9, height = 7, dpi=300)
  }
  if (!is.null(hm_down) && nrow(as.data.frame(hm_down)) > 0) {
    write.csv(as.data.frame(hm_down), file.path(OUTDIR, sprintf("MSIG_HALLMARK_DOWN_%s.csv", cl)), row.names = FALSE)
    p_hm_dn <- dotplot(hm_down, showCategory = 15) + ggtitle(sprintf("Hallmark (DOWN) — %s", cl))
    ggsave(file.path(OUTDIR, sprintf("MSIG_HALLMARK_DOWN_%s.png", cl)), p_hm_dn, width = 9, height = 7, dpi=300)
  }

  # --- Optional: Reactome pathways (lots of xenobiotic/toxic & innate immunity)
  # if (length(entrez_up) >= min_genes_enrich) {
  #   rx_up <- run_reactome(entrez_up, universe_entrez)
  #   if (!is.null(rx_up)) {
  #     write.csv(as.data.frame(rx_up), file.path(OUTDIR, sprintf("REACTOME_UP_%s.csv", cl)), row.names = FALSE)
  #     ggsave(file.path(OUTDIR, sprintf("REACTOME_UP_%s.png", cl)),
  #            dotplot(rx_up, showCategory = 20) + ggtitle(sprintf("Reactome (UP) — %s", cl)),
  #            width = 9, height = 7, dpi = 300)
  #   }
  # }
}

# Export session info (reproducibility)
sink(file.path(OUTDIR, "sessionInfo.txt")); print(sessionInfo()); sink()
message("Done. Outputs in: ", OUTDIR)
