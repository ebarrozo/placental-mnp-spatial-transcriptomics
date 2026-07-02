# maternal_fetal_origin.R
# 4.20.2026
# Enrico Barrozo, PhD

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(stringr)
})

# -------------------- USER CONFIG --------------------
setwd("/home/ebarrozo/visiumhd_mnp/results/audra")

obj_path <- "/media/jochum00/Aagaard_Raid3/ebarrozo/seurat_all_ref_merged/with_annotations/SCT_chunks_mapped/audra/mnp_obj_w_nn_predicted_celltype_fullcells_with_RNAassay.rds"

celltype_col <- "nn.predicted.celltype"   # change if needed
out_name <- "mnp_obj_fullcells_maternal_fetal_rule_based.rds"
out_path <- file.path(getwd(), out_name)

# set TRUE only for known male fetuses / male samples
male_fetus <- TRUE

# -------------------- LOAD --------------------
obj <- readRDS(obj_path)

sort(unique(as.character(obj$nn.predicted.celltype)))

if (!celltype_col %in% colnames(obj[[]])) {
  stop(paste("Metadata column not found:", celltype_col))
}

if (!"RNA" %in% names(obj@assays)) {
  stop("RNA assay not found in object.")
}

DefaultAssay(obj) <- "RNA"

## remove RNA reference data 
Idents(obj) <- "platform"
ncol(obj)
levels(obj)
#DefaultAssay(obj) <- "SCT"
obj <- subset(x = obj, idents = c("Xenium","VisiumHD","VisiumV1"))
ncol(obj)

# -------------------- STANDARDIZE CELLTYPE LABELS --------------------
ct <- as.character(obj[[celltype_col]][,1])
ct_clean <- str_squish(ct)

obj$celltype_clean <- ct_clean

# -------------------- INITIALIZE OUTPUT COLUMNS --------------------
obj$origin_rule_based   <- "ambiguous"
obj$origin_confidence   <- "low"
obj$origin_reason       <- "no_rule_matched"
obj$origin_sex_supported <- NA

# -------------------- RULE SETS --------------------
# adjust these patterns to match your actual labels
fetal_patterns <- c(
  "Syncytiotrophoblast",
  "Extravillous Trophoblast", 
  "Villous Cytotrophoblast", 
  "Hofbauer Cell", "Placental Endothelial", "Perivascular Cell"
)

maternal_patterns <- c(
  "Decidual NK Cell",
  "B Cell"
)

# helper
match_any_pattern <- function(x, patterns) {
  out <- rep(FALSE, length(x))
  for (p in patterns) {
    out <- out | str_detect(str_to_lower(x), fixed(str_to_lower(p)))
  }
  out
}

is_fetal_ct <- match_any_pattern(obj$celltype_clean, fetal_patterns)
is_maternal_ct <- match_any_pattern(obj$celltype_clean, maternal_patterns)

# assign high-confidence labels
obj$origin_rule_based[is_fetal_ct] <- "fetal"
obj$origin_confidence[is_fetal_ct] <- "high"
obj$origin_reason[is_fetal_ct] <- "fetal_lineage_annotation"

obj$origin_rule_based[is_maternal_ct] <- "maternal"
obj$origin_confidence[is_maternal_ct] <- "high"
obj$origin_reason[is_maternal_ct] <- "maternal_decidual_annotation"

# if any labels matched both, mark for review
both <- is_fetal_ct & is_maternal_ct
if (any(both)) {
  obj$origin_rule_based[both] <- "ambiguous"
  obj$origin_confidence[both] <- "low"
  obj$origin_reason[both] <- "conflicting_annotation_patterns"
}

# -------------------- SEX-LINKED SUPPORT --------------------
xist_genes <- intersect(c("XIST"), rownames(obj))
y_genes <- intersect(c("RPS4Y1", "KDM5D", "DDX3Y", "EIF1AY", "UTY", "ZFY"), rownames(obj))

if (length(xist_genes) > 0) {
  obj <- AddModuleScore(
    object = obj,
    features = list(xist_genes),
    name = "XIST_score",
    assay = "RNA"
  )
  obj$XIST_score <- obj$XIST_score1
}

if (male_fetus && length(y_genes) > 0) {
  obj <- AddModuleScore(
    object = obj,
    features = list(y_genes),
    name = "Y_score",
    assay = "RNA"
  )
  obj$Y_score <- obj$Y_score1
}

# -------------------- OPTIONAL SUPPORT / UPGRADE RULES --------------------
# only support calls; do not indiscriminately force every ambiguous cell

if (male_fetus && "Y_score" %in% colnames(obj[[]])) {

  # strong Y signal supports fetal origin
  y_hi <- obj$Y_score > 0.25

  obj$origin_sex_supported[y_hi & obj$origin_rule_based == "fetal"] <- "supported_by_Y"
  obj$origin_sex_supported[y_hi & obj$origin_rule_based == "maternal"] <- "discordant_Y_in_maternal_call"
  obj$origin_sex_supported[y_hi & obj$origin_rule_based == "ambiguous"] <- "Y_supports_fetal"

  # cautiously upgrade ambiguous cells with strong Y signal
  upgrade_to_fetal <- y_hi & obj$origin_rule_based == "ambiguous"
  obj$origin_rule_based[upgrade_to_fetal] <- "fetal"
  obj$origin_confidence[upgrade_to_fetal] <- "medium"
  obj$origin_reason[upgrade_to_fetal] <- "male_fetus_Y_signal"
}

if ("XIST_score" %in% colnames(obj[[]])) {

  xist_hi <- obj$XIST_score > 0.25

  obj$origin_sex_supported[xist_hi & obj$origin_rule_based == "maternal" &
                             is.na(obj$origin_sex_supported)] <- "supported_by_XIST"

  # Do NOT automatically relabel ambiguous cells as maternal from XIST alone
  # XIST is less specific than Y in this setting
}

# -------------------- OPTIONAL: KEEP GENERIC IMMUNE/STROMAL CELLS AMBIGUOUS --------------------
# This prevents overcalling maternal/fetal when labels are broad.
generic_ambiguous_patterns <- c(
  "Macrophage",
  "NK Cell",
  "T Cell", "Treg",
  "Plasma Cell",
  "B cell",
  "Perivascular Cell",
  "Fibroblast",
  "RCTD_lowconf"
)

generic_ambiguous <- match_any_pattern(obj$celltype_clean, generic_ambiguous_patterns)

# preserve high-confidence decidual / trophoblast assignments
generic_to_reset <- generic_ambiguous &
  !(obj$origin_reason %in% c("fetal_lineage_annotation", "maternal_decidual_annotation"))

obj$origin_rule_based[generic_to_reset & obj$origin_reason == "no_rule_matched"] <- "ambiguous"
obj$origin_confidence[generic_to_reset & obj$origin_reason == "no_rule_matched"] <- "low"
obj$origin_reason[generic_to_reset & obj$origin_reason == "no_rule_matched"] <- "generic_celltype_ambiguous"

# -------------------- SUMMARIES --------------------
cat("\nOrigin counts:\n")
print(table(obj$origin_rule_based, useNA = "ifany"))

cat("\nConfidence counts:\n")
print(table(obj$origin_confidence, useNA = "ifany"))

cat("\nReason counts:\n")
print(sort(table(obj$origin_reason, useNA = "ifany"), decreasing = TRUE))

cat("\nCross-tab by cell type:\n")
print(table(obj$origin_rule_based, obj$celltype_clean, useNA = "ifany"))

t1<- table(obj$origin_rule_based, obj$nn.predicted.celltype)
write.csv(t1, "obj-origin_rule_based_nn.predicted.celltype.csv", row.names = T)
t1<- table(obj$origin_confidence, obj$nn.predicted.celltype)
write.csv(t1, "origin_confidence_nn.predicted.celltype.csv", row.names = T)
t1<- table(obj$origin_rule_based, obj$sample_id_std)
write.csv(t1, "obj-origin_rule_based_sample_id_std.csv", row.names = T)


# -------------------- SAVE --------------------
out_name <- "mnp_obj_fullcells_maternal_fetal_rule_based.rds"
saveRDS(obj, out_path)
cat("\nSaved:\n", out_path, "\n")








obj <- readRDS("/home/ebarrozo/visiumhd_mnp/results/audra/mnp_obj_fullcells_maternal_fetal_rule_based.rds")
DefaultAssay(obj) <- "RNA"

table(obj$origin_rule_based)
table(obj$origin_confidence)
table(obj$origin_rule_based, obj$nn.predicted.celltype)

library(ggplot2)
library(cowplot)
library(patchwork)
library(ggrepel)
library(ggpubr)
library(ggplot2)
mypal3 <-get_palette("igv",40)

Idents(obj) <- "origin_rule_based"
p2<-DimPlot(obj, reduction = "umap", group.by = "origin_rule_based", alpha = 0.1, label = FALSE,raster=T, cols = mypal3) + ggtitle("Predicted Cell Types (3,969,532 Cells)") + NoAxes()
p2
ggsave("obj_group.by_origin_rule_based.png", plot = p2, device = "png", width = 8, height = 5, units = "in", dpi = 1200)


Idents(obj) <- "origin_rule_based"
p2<-DimPlot(obj, reduction = "umap", split.by = "origin_rule_based", alpha = 0.1, label = FALSE,raster=T, cols = mypal3) + ggtitle("Predicted Cell Types (3,969,532 Cells) | 1,037 maternal origin") + NoAxes()
p2
ggsave("obj_split.by_origin_rule_based.png", plot = p2, device = "png", width = 18, height = 5, units = "in", dpi = 1200)


Idents(obj) <- "nn.predicted.celltype"
p2<-DimPlot(obj, reduction = "umap", group.by = "nn.predicted.celltype", alpha = 0.1, label = FALSE,raster=T, cols = mypal3) + ggtitle("Predicted Cell Types (3,969,532 Cells)") + NoAxes()
p2
ggsave("obj_nn.predicted.celltype.png", plot = p2, device = "png", width = 8, height = 5, units = "in", dpi = 1200)



DimPlot(obj, group.by = "origin_rule_based", reduction = "umap")
FeaturePlot(obj, features = c("XIST", "RPS4Y1", "KDM5D"), reduction = "umap")