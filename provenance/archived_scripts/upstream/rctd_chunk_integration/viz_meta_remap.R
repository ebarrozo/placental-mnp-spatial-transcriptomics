## viz_meta_remap.R

set.seed(1) # set seed for reproducibility
options(future.globals.maxSize = 400e9) ## 300 gb ram (out of 440 on server 3; adjust if on different server)

setwd("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations")
suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(tibble)
})

# ---- helpers (same as before) ----
parse_platform <- function(x) case_when(
  str_ends(x, regex("_xen$", ignore_case=TRUE)) ~ "Xenium",
  str_ends(x, regex("_hd$",  ignore_case=TRUE)) ~ "VisiumHD",
  str_ends(x, regex("_v1$",  ignore_case=TRUE)) ~ "VisiumV1",
  TRUE ~ NA_character_
)
strip_platform_suffix <- function(x) str_remove(x, "_(xen|hd|v1)$")
infer_study_id <- function(sample_id) sub("([A-Za-z]+\\d+).*", "\\1", sample_id)

# ---- 0) load/build chunk manifest 'man' (unchanged) ----
chunk_dir <- "SCT_chunks"
man_file  <- file.path(chunk_dir, "chunk_manifest.csv")
if (!exists("man")) {
  if (file.exists(man_file)) {
    man <- read.csv(man_file, stringsAsFactors = FALSE)
  } else {
    fns <- list.files(chunk_dir, "^sctrans_.*_chunk_\\d{3}_SCTonly\\.rds$", full.names = TRUE)
    stopifnot(length(fns) > 0)
    parse_row <- function(fp){
      bn <- basename(fp)
      m  <- regexec("^sctrans_(.+)_chunk_(\\d{3})_SCTonly\\.rds$", bn)
      mm <- regmatches(bn, m)[[1]]
      o  <- readRDS(fp); nc <- ncol(o); rm(o); gc()
      tibble(sample = mm[2], chunk = as.integer(mm[3]), cells = nc, file = fp)
    }
    man <- lapply(fns, parse_row) |> bind_rows() |> arrange(sample, chunk)
    write.csv(man, man_file, row.names = FALSE)
  }
}

# ---- 1) robust parse from cell names (works WITH or WITHOUT platform suffix) ----
cell_df <- tibble(cell = colnames(viz)) %>%
  mutate(
    sample = sub("_c\\d{3}_.+$", "", cell),                 # everything before _cNNN_
    chunk  = suppressWarnings(as.integer(sub("^.*_c(\\d{3})_.+$", "\\1", cell)))
  )

# quick sanity: any unparsable cells?
bad <- filter(cell_df, is.na(chunk) | !str_detect(cell, "_c\\d{3}_"))
if (nrow(bad)) write.csv(bad, "cells_unparsable_regex.csv", row.names = FALSE)

# count cells per (sample, chunk) present in viz
viz_chunks <- cell_df %>%
  count(sample, chunk, name = "cells_in_viz") %>%
  arrange(sample, chunk)

# ---- 2) join to manifest (this is your “chunks in viz” counterpart) ----
chunks_in_viz <- viz_chunks %>%
  left_join(man, by = c("sample","chunk")) %>%
  rename(cells_chunk_manifest = cells, file_chunk = file) %>%
  mutate(
    frac_in_viz = ifelse(!is.na(cells_chunk_manifest) & cells_chunk_manifest > 0,
                         pmin(1, cells_in_viz / cells_chunk_manifest), NA_real_)
  )

# ---- 3) add platform + standardized IDs (handles samples w/o suffix) ----
manifest_samples <- man %>%
  group_by(sample) %>%
  summarise(n_chunks = n(), cells_chunks = sum(cells), .groups="drop") %>%
  mutate(platform_guess = parse_platform(sample),
         sample_core    = strip_platform_suffix(sample),
         study_id_guess = infer_study_id(sample_core))

manual_map <- c(
  "HP1" = "Xenium", "HP13a" = "VisiumV1", "HT2b" = "Xenium",
  "HT3" = "Xenium", "LP2" = "Xenium", "placental_ref_Jul_2025" = "RNAref"
)
suffix_map <- c(Xenium="xen", VisiumHD="hd", VisiumV1="v1", RNAref=NA)

manifest_fixed_local <- manifest_samples %>%
  mutate(
    platform = case_when(
      !is.na(platform_guess) ~ platform_guess,
      sample %in% names(manual_map) ~ unname(manual_map[sample]),
      TRUE ~ "VisiumV1"  # default for remaining NA per your rule
    ),
    seurat_id = case_when(
      platform == "RNAref" ~ sample,
      str_detect(sample, "(_xen|_hd|_v1)$") ~ sample,
      TRUE ~ paste0(sample_core, "_", suffix_map[platform])
    ),
    study_id = study_id_guess
  ) %>%
  select(sample, platform, seurat_id, study_id)

chunks_in_viz <- chunks_in_viz %>%
  left_join(manifest_fixed_local, by = "sample") %>%
  relocate(study_id, platform, seurat_id, .before = sample) %>%
  arrange(study_id, platform, sample, chunk)

# ---- 4) per-sample roll-up (only chunks that landed in viz) ----
chunks_in_viz_by_sample <- chunks_in_viz %>%
  group_by(study_id, platform, seurat_id, sample) %>%
  summarise(
    cells_in_viz_sum         = sum(cells_in_viz, na.rm = TRUE),
    cells_chunk_manifest_sum = sum(cells_chunk_manifest, na.rm = TRUE),
    n_chunks_in_viz          = dplyr::n(),
    .groups = "drop"
  ) %>%
  mutate(frac_in_viz_vs_manifest = ifelse(cells_chunk_manifest_sum > 0,
                                          pmin(1, cells_in_viz_sum / cells_chunk_manifest_sum), NA_real_)) %>%
  arrange(study_id, platform, seurat_id)

# ---- 5) sanity + write out ----
tot_viz <- ncol(viz)
tot_from_chunks <- sum(viz_chunks$cells_in_viz)
cat("Cells in viz:", tot_viz, "\nCells accounted by (sample,chunk):", tot_from_chunks, "\nDelta:", tot_viz - tot_from_chunks, "\n")

write.csv(chunks_in_viz,           "chunks_in_viz_detailed.csv", row.names = FALSE)
write.csv(chunks_in_viz_by_sample, "chunks_in_viz_by_sample.csv", row.names = FALSE)



###### 


library(Seurat)
library(dplyr)

## 0) sanity: make sure the column exists and has no NAs
if (!"sample_id_std" %in% colnames(viz@meta.data)) stop("sample_id_std not found in viz@meta.data")

# fill NAs with a sentinel (Seurat Idents cannot be NA)
viz$sample_id_std[is.na(viz$sample_id_std)] <- "Unknown_sample"
viz$predicted.celltype[is.na(viz$predicted.celltype)] <- "RCTD_lowconf"

## 1) Choose useful level orders (edit if you want a custom order)
# samples: alphabetical by StudyID then platform, then sample_id_std
#ord_samples <- viz@meta.data %>%
#  tibble::rownames_to_column("cell") %>%
#  distinct(sample_id_std, study_id, platform) %>%
#  arrange(study_id, factor(platform, levels = c("Xenium","VisiumHD","VisiumV1","RNAref")), sample_id_std) %>%
#  pull(sample_id_std)

# platforms
ord_platform <- c("Xenium","VisiumHD","VisiumV1","RNAref")

# cell types (use your palette order; edit if needed)
ord_celltypes <- c(
  "Syncytiotrophoblast","Villous Cytotrophoblast","Extravillous Trophoblast",
  "Fibroblast","Perivascular Cell","Placental Endothelial",
  "Hofbauer Cell","Macrophage","NK Cell","Decidual NK Cell","ILC3",
  "Plasma Cell","B Cell","T Cell",
  "Treg", "RCTD_lowconf"
)

## 2) Cast to factors with desired levels
viz$orig.ident           <- factor(viz$orig.ident)                                # 28 levels, as before
# viz$sample_id_std        <- factor(viz$sample_id_std, levels = unique(ord_samples))
viz$platform             <- factor(viz$platform,     levels = ord_platform)
if ("predicted.celltype" %in% colnames(viz@meta.data)) {
  # keep unknowns at the end if any
  known <- intersect(ord_celltypes, unique(viz$predicted.celltype))
  other <- setdiff(unique(viz$predicted.celltype), known)
  viz$predicted.celltype <- factor(viz$predicted.celltype, levels = c(known, sort(other)))
}

## 3) Set identities to sample_id_std
Idents(viz) <- "sample_id_std"   # this works even if the column was character; now it's factor
levels(viz)                       # inspect identity levels

## 4) quick checks
stopifnot(all(!is.na(Idents(viz))))             # no NA identities
table(Idents(viz))[1:5]                         # first few counts


# If you want identities to be cell types instead:
# Idents(viz) <- "predicted.celltype"

# If you want to relevel after the fact (e.g., put a specific sample first):
# Idents(viz) <- factor(Idents(viz), levels = c("HT2a_xen", setdiff(levels(viz), "HT2a_xen")))



suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(tibble)
})

# --- 0) load cohort CSV once ---
CSV_PATH <- "/path/to/local_home/visiumhd_mnp/docs/spatial_cohort_annotation_spatialMedian_quickannot_final.csv"
ext <- read_csv(CSV_PATH, show_col_types = FALSE)

stopifnot("StudyID" %in% names(ext))
stopifnot("study_id" %in% colnames(viz@meta.data))  # you created this earlier

# --- 1) throw away platform-specific ID columns so we get pure StudyID metadata ---
drop_cols <- c("XeniumID_seq","VisiumV1_seq","VisiumHD_seq",
               "XenID_seurat","VisV1_seurat","VisHD_seurat")
keep_cols <- setdiff(names(ext), intersect(drop_cols, names(ext)))
ext_study <- ext %>% select(all_of(keep_cols))

# --- 2) collapse to ONE row per StudyID (no duplicates) ---
collapse_char <- function(x){
  u <- unique(na.omit(as.character(x)))
  if (length(u) == 0) NA_character_ else paste(u, collapse = " | ")
}
collapse_num <- function(x){
  ux <- na.omit(x)
  if (!length(ux)) NA_real_ else ux[1]  # pick first non-NA; swap to median(ux) if you prefer
}

study_meta <- ext_study %>%
  group_by(StudyID) %>%
  summarise(
    across(where(is.numeric), collapse_num),
    across(!where(is.numeric), collapse_char),
    .groups = "drop"
  ) %>%
  rename(study_id = StudyID)

# sanity
stopifnot(nrow(study_meta) == nrow(distinct(study_meta, study_id)))

# --- 3) build a per-cell table: cell, study_id -> study-level metadata ---
cell_by_study <- tibble(
  cell     = rownames(viz@meta.data),
  study_id = viz$study_id
) %>%
  left_join(study_meta, by = "study_id")

# --- 4) write new columns back into viz (don’t overwrite existing core fields) ---
protect <- c("cell","study_id","platform","sample_id_std","predicted.celltype","orig.ident")
new_cols <- setdiff(names(cell_by_study), protect)
for (nm in new_cols) viz[[nm]] <- cell_by_study[[nm]]

# --- 5) convert added character columns to factors (bounded #levels to avoid huge factors) ---
to_factor_if_ok <- function(v, max_levels = 2000) {
  if (is.character(v)) {
    lv <- unique(na.omit(v))
    if (length(lv) <= max_levels) return(factor(v, levels = sort(lv)))
  }
  v
}
for (nm in new_cols) viz[[nm]] <- to_factor_if_ok(viz[[nm]])

# (re)factor core fields too (if you haven’t already)
viz$sample_id_std    <- factor(viz$sample_id_std)                           # identities can use this
viz$platform         <- factor(viz$platform, levels = c("Xenium","VisiumHD","VisiumV1","RNAref"))
if ("predicted.celltype" %in% colnames(viz@meta.data))
  viz$predicted.celltype <- factor(viz$predicted.celltype)

# --- 6) (optional) set identities to sample or study ---
# Idents(viz) <- "sample_id_std"
# Idents(viz) <- "study_id"

# --- 7) quick diagnostics ---
cat("study_meta rows:", nrow(study_meta), "\n")
cat("cells with missing StudyID metadata:", sum(is.na(viz[[new_cols[1]]])), "\n")  # crude check on the first new column




library(ggplot2)
library(cowplot)
library(patchwork)
library(ggrepel)
library(ggpubr)
library(ggplot2)


mypal3 <-get_palette("igv",40)
mypal3 <-get_palette("igv",40)

Idents(viz) <- "predicted.celltype"
p2<-DimPlot(viz, reduction = "umap", group.by = "predicted.celltype", alpha = 0.1, label = FALSE,raster=T, cols = mypal3) + ggtitle("Predicted Cell Types (3,969,532 Cells)") + NoAxes()
p2
ggsave("viz_predicted.predicted.celltype.png", plot = p2, device = "png", width = 8, height = 5, units = "in", dpi = 1200)


## rename NA to RCTD_lowconf for lowconfidence not mapped 
Idents(viz) <- "predicted.celltype"
levels(viz)         # 
levels(viz)         # 

new.metadata <- c(
  "Syncytiotrophoblast","Villous Cytotrophoblast","Extravillous Trophoblast",
  "Fibroblast","Perivascular Cell","Placental Endothelial",
  "Hofbauer Cell","Macrophage","NK Cell","Decidual NK Cell","ILC3",
  "Plasma Cell","B Cell","T Cell",
  "Treg", "RCTD_lowconf"
)

names(new.metadata) <- levels(viz)
viz <- RenameIdents(viz, new.metadata)
viz$predicted.celltype <- Idents(viz)
Idents(viz) <- "predicted.celltype"



# Save for reuse
saveRDS(viz, file = file.path(OUTDIR_MAPPED, "viz_re-metadata.rds"))

rm(list = ls()[!ls() %in% c("viz", "mypal3")])
gc()





library(ggplot2)
library(cowplot)
library(patchwork)
library(ggrepel)
library(ggpubr)
library(ggplot2)


mypal3 <-get_palette("igv",40)

Idents(viz) <- "predicted.celltype"
p2<-DimPlot(viz, reduction = "umap", group.by = "predicted.celltype", alpha = 0.1, label = FALSE,raster=T, cols = mypal3) + ggtitle("Predicted Cell Types (3,969,532 Cells)") + NoAxes()
p2
ggsave("viz_predicted.predicted.celltype.png", plot = p2, device = "png", width = 8, height = 5, units = "in", dpi = 1200)

viz <- readRDS("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/SCT_chunks_mapped/viz_re-metadata.rds") ## this is  5 gb



setwd("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations")


viz@meta.data[["MNPBurden_Delivery_class"]]
viz@meta.data[["MNPBurden"]]
mypal3 <-get_palette("igv",40)

Idents(viz) <- "MNPBurden_Delivery_class"
p2<-DimPlot(viz, reduction = "umap", group.by = "predicted.celltype", alpha = 0.1, label = FALSE,raster=T, cols = mypal3) + ggtitle("Predicted Cell Types (3,969,532 Cells)") + NoAxes()
p2
ggsave("viz_predicted.predicted.celltype.png", plot = p2, device = "png", width = 8, height = 5, units = "in", dpi = 1200)












## I need counts to subset
Idents(viz) <- "MNPBurden_Delivery_class"
levels(viz)
## remove samples w/o MNP data
ncol(viz) # 3969532
viz.mnp <- subset(x = viz, idents = c("High_Term", "High_Preterm", "Low_Term", "Low_Preterm"))
ncol(viz.mnp)
gc()





