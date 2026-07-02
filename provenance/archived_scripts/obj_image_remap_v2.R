## obj_image_remap_v2.R
## server3 
## Enrico Barrozo, Ph.D.
## Assistant Professor, Department of OB/GYN, Baylor College of Medicine and Texas Children’s Hospital
# 10.7.25

## just build niches and export for Xenium Explorer

setwd("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations")
obj <-readRDS("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/SCT_chunks_mapped/final_sct_atlas.rds") # 7.5 gb object with SCT data, all 3,969,532 cells, and umap
suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix); library(data.table)
})

##########################################################
### HP13b fov S23b
HP13b_xen <-readRDS("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/seurat_HP13b_xen_withPredCols.rds") 
HP13b_xen <- BuildNicheAssay(object = HP13b_xen, fov="S23b", group.by = "predicted.celltype")
Idents(HP13b_xen) <- "niches"

# rename NAs in RCTD predicted celltypes to low confidence mapped
HP13b_xen$predicted.celltype[is.na(obj$predicted.celltype)] <- "RCTD_lowconf"

# Export predicted for Xenium Explorer
HP13b_xen <- RenameCells(HP13b_xen, new.names = gsub("_[0-9]+$", "", Cells(HP13b_xen)))
# save the annotations for uploading into Xenium Explorer
Idents(HP13b_xen) <- "predicted.celltype"
predicted.celltype = as.data.frame(Idents(HP13b_xen))
predicted.celltype$cell_id = (rownames(predicted.celltype))
predicted.celltype$group = (predicted.celltype[,1])
predicted.celltype = predicted.celltype[,2:3]
write.csv(predicted.celltype, "HP13b_xen_predicted.celltype_forXeniumExplorer.csv", quote=FALSE, row.names = F)
HP13b_xen$predicted.celltype[is.na(obj$predicted.celltype)] <- "RCTD_lowconf"


saveRDS(HP13b_xen, "HP13b_xen_niches.rds")
#################################################################################
### HP2_xen fov S31542
HP2_xen <-readRDS("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/seurat_HP2_xen_withPredCols.rds") 
HP2_xen <- BuildNicheAssay(object = HP2_xen, fov="S31542", group.by = "predicted.celltype")
Idents(HP2_xen) <- "niches"

# rename NAs in RCTD predicted celltypes to low confidence mapped
HP2_xen$predicted.celltype[is.na(obj$predicted.celltype)] <- "RCTD_lowconf"

# Export predicted for Xenium Explorer
HP2_xen <- RenameCells(HP2_xen, new.names = gsub("_[0-9]+$", "", Cells(HP2_xen)))
# save the annotations for uploading into Xenium Explorer
Idents(HP2_xen) <- "predicted.celltype"
predicted.celltype = as.data.frame(Idents(HP2_xen))
predicted.celltype$cell_id = (rownames(predicted.celltype))
predicted.celltype$group = (predicted.celltype[,1])
predicted.celltype = predicted.celltype[,2:3]
write.csv(predicted.celltype, "HP2_xen_predicted.celltype_forXeniumExplorer.csv", quote=FALSE, row.names = F)
Idents(HP2_xen) <- "niches"
niches = as.data.frame(Idents(HP2_xen))
niches$cell_id = (rownames(niches))
niches$group = (niches[,1])
niches = niches[,2:3]
write.csv(niches, "HP2_xen_niches_forXeniumExplorer.csv", quote=FALSE, row.names = F)

saveRDS(HP2_xen, "HP2_xen_niches.rds")
#################################################################################
### HT1_xen fov S11231
HT1_xen <-readRDS("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/seurat_HT1_xen_withPredCols.rds") 
HT1_xen <- BuildNicheAssay(object = HT1_xen, fov="S11231", group.by = "predicted.celltype")
Idents(HT1_xen) <- "niches"

# rename NAs in RCTD predicted celltypes to low confidence mapped
HT1_xen$predicted.celltype[is.na(obj$predicted.celltype)] <- "RCTD_lowconf"

# Export predicted for Xenium Explorer
HT1_xen <- RenameCells(HT1_xen, new.names = gsub("_[0-9]+$", "", Cells(HT1_xen)))
# save the annotations for uploading into Xenium Explorer
Idents(HT1_xen) <- "predicted.celltype"
predicted.celltype = as.data.frame(Idents(HT1_xen))
predicted.celltype$cell_id = (rownames(predicted.celltype))
predicted.celltype$group = (predicted.celltype[,1])
predicted.celltype = predicted.celltype[,2:3]
write.csv(predicted.celltype, "HT1_xen_predicted.celltype_forXeniumExplorer.csv", quote=FALSE, row.names = F)
Idents(HT1_xen) <- "niches"
niches = as.data.frame(Idents(HT1_xen))
niches$cell_id = (rownames(niches))
niches$group = (niches[,1])
niches = niches[,2:3]
write.csv(niches, "HT1_xen_niches_forXeniumExplorer.csv", quote=FALSE, row.names = F)

saveRDS(HT1_xen, "HT1_xen_niches.rds")
#################################################################################
### HT2a_xen fov S13164a
HT2a_xen <-readRDS("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/seurat_HT2a_xen_withPredCols.rds") 
HT2a_xen <- BuildNicheAssay(object = HT2a_xen, fov="S13164a", group.by = "predicted.celltype")
Idents(HT2a_xen) <- "niches"

# rename NAs in RCTD predicted celltypes to low confidence mapped
HT2a_xen$predicted.celltype[is.na(obj$predicted.celltype)] <- "RCTD_lowconf"

# Export predicted for Xenium Explorer
HT2a_xen <- RenameCells(HT2a_xen, new.names = gsub("_[0-9]+$", "", Cells(HT2a_xen)))
# save the annotations for uploading into Xenium Explorer
Idents(HT2a_xen) <- "predicted.celltype"
predicted.celltype = as.data.frame(Idents(HT2a_xen))
predicted.celltype$cell_id = (rownames(predicted.celltype))
predicted.celltype$group = (predicted.celltype[,1])
predicted.celltype = predicted.celltype[,2:3]
write.csv(predicted.celltype, "HT2a_xen_predicted.celltype_forXeniumExplorer.csv", quote=FALSE, row.names = F)
Idents(HT2a_xen) <- "niches"
niches = as.data.frame(Idents(HT2a_xen))
niches$cell_id = (rownames(niches))
niches$group = (niches[,1])
niches = niches[,2:3]
write.csv(niches, "HT2a_xen_niches_forXeniumExplorer.csv", quote=FALSE, row.names = F)

saveRDS(HT2a_xen, "HT2a_xen_niches.rds")
#################################################################################
### HT3_xen fov S3065
HT3_xen <-readRDS("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/seurat_HT3_xen_withPredCols.rds") 
HT3_xen <- BuildNicheAssay(object = HT3_xen, fov="S3065", group.by = "predicted.celltype")
Idents(HT3_xen) <- "niches"

# rename NAs in RCTD predicted celltypes to low confidence mapped
HT3_xen$predicted.celltype[is.na(obj$predicted.celltype)] <- "RCTD_lowconf"

# Export predicted for Xenium Explorer
HT3_xen <- RenameCells(HT3_xen, new.names = gsub("_[0-9]+$", "", Cells(HT3_xen)))
# save the annotations for uploading into Xenium Explorer
Idents(HT3_xen) <- "predicted.celltype"
predicted.celltype = as.data.frame(Idents(HT3_xen))
predicted.celltype$cell_id = (rownames(predicted.celltype))
predicted.celltype$group = (predicted.celltype[,1])
predicted.celltype = predicted.celltype[,2:3]
write.csv(predicted.celltype, "HT3_xen_predicted.celltype_forXeniumExplorer.csv", quote=FALSE, row.names = F)
Idents(HT3_xen) <- "niches"
niches = as.data.frame(Idents(HT3_xen))
niches$cell_id = (rownames(niches))
niches$group = (niches[,1])
niches = niches[,2:3]
write.csv(niches, "HT3_xen_niches_forXeniumExplorer.csv", quote=FALSE, row.names = F)

saveRDS(HT3_xen, "HT3_xen_niches.rds")
#################################################################################
### HP1 fov S9380
HP1 <-readRDS("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/seurat_HP1_withPredCols.rds") 
HP1 <- BuildNicheAssay(object = HP1, fov="S9380", group.by = "predicted.celltype")
Idents(HP1) <- "niches"

# rename NAs in RCTD predicted celltypes to low confidence mapped
HP1$predicted.celltype[is.na(obj$predicted.celltype)] <- "RCTD_lowconf"

# Export predicted for Xenium Explorer
HP1 <- RenameCells(HP1, new.names = gsub("_[0-9]+$", "", Cells(HP1)))
# save the annotations for uploading into Xenium Explorer
Idents(HP1) <- "predicted.celltype"
predicted.celltype = as.data.frame(Idents(HP1))
predicted.celltype$cell_id = (rownames(predicted.celltype))
predicted.celltype$group = (predicted.celltype[,1])
predicted.celltype = predicted.celltype[,2:3]
write.csv(predicted.celltype, "HP1_predicted.celltype_forXeniumExplorer.csv", quote=FALSE, row.names = F)
Idents(HP1) <- "niches"
niches = as.data.frame(Idents(HP1))
niches$cell_id = (rownames(niches))
niches$group = (niches[,1])
niches = niches[,2:3]
write.csv(niches, "HP1_niches_forXeniumExplorer.csv", quote=FALSE, row.names = F)

saveRDS(HP1, "HP1_niches.rds")
#################################################################################
### HT2b fov S13164b
HT2b <-readRDS("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/seurat_HT2b_withPredCols.rds") 
HT2b <- BuildNicheAssay(object = HT2b, fov="S13164b", group.by = "predicted.celltype")
Idents(HT2b) <- "niches"

# rename NAs in RCTD predicted celltypes to low confidence mapped
HT2b$predicted.celltype[is.na(obj$predicted.celltype)] <- "RCTD_lowconf"

# Export predicted for Xenium Explorer
HT2b <- RenameCells(HT2b, new.names = gsub("_[0-9]+$", "", Cells(HT2b)))
# save the annotations for uploading into Xenium Explorer
Idents(HT2b) <- "predicted.celltype"
predicted.celltype = as.data.frame(Idents(HT2b))
predicted.celltype$cell_id = (rownames(predicted.celltype))
predicted.celltype$group = (predicted.celltype[,1])
predicted.celltype = predicted.celltype[,2:3]
write.csv(predicted.celltype, "HT2b_predicted.celltype_forXeniumExplorer.csv", quote=FALSE, row.names = F)
Idents(HT2b) <- "niches"
niches = as.data.frame(Idents(HT2b))
niches$cell_id = (rownames(niches))
niches$group = (niches[,1])
niches = niches[,2:3]
write.csv(niches, "HT2b_niches_forXeniumExplorer.csv", quote=FALSE, row.names = F)

saveRDS(HT2b, "HT2b_niches.rds")
#################################################################################
### LP2 fov S10963
LP2 <-readRDS("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/seurat_LP2_withPredCols.rds") 
LP2 <- BuildNicheAssay(object = LP2, fov="S10963", group.by = "predicted.celltype")
Idents(LP2) <- "niches"

# rename NAs in RCTD predicted celltypes to low confidence mapped
LP2$predicted.celltype[is.na(obj$predicted.celltype)] <- "RCTD_lowconf"

# Export predicted for Xenium Explorer
LP2 <- RenameCells(LP2, new.names = gsub("_[0-9]+$", "", Cells(LP2)))
# save the annotations for uploading into Xenium Explorer
Idents(LP2) <- "predicted.celltype"
predicted.celltype = as.data.frame(Idents(LP2))
predicted.celltype$cell_id = (rownames(predicted.celltype))
predicted.celltype$group = (predicted.celltype[,1])
predicted.celltype = predicted.celltype[,2:3]
write.csv(predicted.celltype, "LP2_predicted.celltype_forXeniumExplorer.csv", quote=FALSE, row.names = F)
Idents(LP2) <- "niches"
niches = as.data.frame(Idents(LP2))
niches$cell_id = (rownames(niches))
niches$group = (niches[,1])
niches = niches[,2:3]
write.csv(niches, "LP2_niches_forXeniumExplorer.csv", quote=FALSE, row.names = F)

saveRDS(LP2, "LP2_niches.rds")
#################################################################################
### LT2 fov S30578
LT2 <-readRDS("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/seurat_LT2_withPredCols.rds") 
LT2 <- BuildNicheAssay(object = LT2, fov="S30578", group.by = "predicted.celltype")
Idents(LT2) <- "niches"

# rename NAs in RCTD predicted celltypes to low confidence mapped
LT2$predicted.celltype[is.na(obj$predicted.celltype)] <- "RCTD_lowconf"

# Export predicted for Xenium Explorer
LT2 <- RenameCells(LT2, new.names = gsub("_[0-9]+$", "", Cells(LT2)))
# save the annotations for uploading into Xenium Explorer
Idents(LT2) <- "predicted.celltype"
predicted.celltype = as.data.frame(Idents(LT2))
predicted.celltype$cell_id = (rownames(predicted.celltype))
predicted.celltype$group = (predicted.celltype[,1])
predicted.celltype = predicted.celltype[,2:3]
write.csv(predicted.celltype, "LT2_predicted.celltype_forXeniumExplorer.csv", quote=FALSE, row.names = F)
Idents(LT2) <- "niches"
niches = as.data.frame(Idents(LT2))
niches$cell_id = (rownames(niches))
niches$group = (niches[,1])
niches = niches[,2:3]
write.csv(niches, "LT2_niches_forXeniumExplorer.csv", quote=FALSE, row.names = F)

saveRDS(LT2, "LT2_niches.rds")
#################################################################################
#################################################################################
### NC4_xen fov S26
NC4_xen <-readRDS("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/seurat_NC4_xen_withPredCols.rds") 
NC4_xen <- BuildNicheAssay(object = NC4_xen, fov="S26", group.by = "predicted.celltype")
Idents(NC4_xen) <- "niches"

# rename NAs in RCTD predicted celltypes to low confidence mapped
NC4_xen$predicted.celltype[is.na(obj$predicted.celltype)] <- "RCTD_lowconf"

# Export predicted for Xenium Explorer
NC4_xen <- RenameCells(NC4_xen, new.names = gsub("_[0-9]+$", "", Cells(NC4_xen)))
# save the annotations for uploading into Xenium Explorer
Idents(NC4_xen) <- "predicted.celltype"
predicted.celltype = as.data.frame(Idents(NC4_xen))
predicted.celltype$cell_id = (rownames(predicted.celltype))
predicted.celltype$group = (predicted.celltype[,1])
predicted.celltype = predicted.celltype[,2:3]
write.csv(predicted.celltype, "NC4_xen_predicted.celltype_forXeniumExplorer.csv", quote=FALSE, row.names = F)
Idents(NC4_xen) <- "niches"
niches = as.data.frame(Idents(NC4_xen))
niches$cell_id = (rownames(niches))
niches$group = (niches[,1])
niches = niches[,2:3]
write.csv(niches, "NC4_xen_niches_forXeniumExplorer.csv", quote=FALSE, row.names = F)

saveRDS(NC4_xen, "NC4_xen_niches.rds")
#################################################################################
### LT1_xen fov S3065
LT1_xen <-readRDS("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/seurat_LT1_xen_withPredCols.rds") 
LT1_xen <- BuildNicheAssay(object = LT1_xen, fov="S31175", group.by = "predicted.celltype")
Idents(LT1_xen) <- "niches"

# rename NAs in RCTD predicted celltypes to low confidence mapped
LT1_xen$predicted.celltype[is.na(obj$predicted.celltype)] <- "RCTD_lowconf"

# Export predicted for Xenium Explorer
LT1_xen <- RenameCells(LT1_xen, new.names = gsub("_[0-9]+$", "", Cells(LT1_xen)))
# save the annotations for uploading into Xenium Explorer
Idents(LT1_xen) <- "predicted.celltype"
predicted.celltype = as.data.frame(Idents(LT1_xen))
predicted.celltype$cell_id = (rownames(predicted.celltype))
predicted.celltype$group = (predicted.celltype[,1])
predicted.celltype = predicted.celltype[,2:3]
write.csv(predicted.celltype, "LT1_xen_predicted.celltype_forXeniumExplorer.csv", quote=FALSE, row.names = F)
Idents(LT1_xen) <- "niches"
niches = as.data.frame(Idents(LT1_xen))
niches$cell_id = (rownames(niches))
niches$group = (niches[,1])
niches = niches[,2:3]
write.csv(niches, "LT1_xen_niches_forXeniumExplorer.csv", quote=FALSE, row.names = F)

saveRDS(LT1_xen, "LT1_xen_niches.rds")
#################################################################################
### LP1_xen fov S3065
LP1_xen <-readRDS("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/seurat_LP1_xen_withPredCols.rds") 
LP1_xen <- BuildNicheAssay(object = LP1_xen, fov="S12129", group.by = "predicted.celltype")
Idents(LP1_xen) <- "niches"

# rename NAs in RCTD predicted celltypes to low confidence mapped
LP1_xen$predicted.celltype[is.na(obj$predicted.celltype)] <- "RCTD_lowconf"

# Export predicted for Xenium Explorer
LP1_xen <- RenameCells(LP1_xen, new.names = gsub("_[0-9]+$", "", Cells(LP1_xen)))
# save the annotations for uploading into Xenium Explorer
Idents(LP1_xen) <- "predicted.celltype"
predicted.celltype = as.data.frame(Idents(LP1_xen))
predicted.celltype$cell_id = (rownames(predicted.celltype))
predicted.celltype$group = (predicted.celltype[,1])
predicted.celltype = predicted.celltype[,2:3]
write.csv(predicted.celltype, "LP1_xen_predicted.celltype_forXeniumExplorer.csv", quote=FALSE, row.names = F)
# save the annotations for uploading into Xenium Explorer
Idents(LP1_xen) <- "niches"
niches = as.data.frame(Idents(LP1_xen))
niches$cell_id = (rownames(niches))
niches$group = (niches[,1])
niches = niches[,2:3]
write.csv(niches, "LP1_xen_niches_forXeniumExplorer.csv", quote=FALSE, row.names = F)

saveRDS(LP1_xen, "LP1_xen_niches.rds")













### Now visumHD samples? 
not run
#################################################################################
### HP2_hd fov S3065
HP2_hd <-readRDS("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/seurat_HP2_hd_withPredCols.rds") 
HP2_hd <- BuildNicheAssay(object = HP2_hd, fov="S3065", group.by = "predicted.celltype")
Idents(HP2_hd) <- "niches"

# rename NAs in RCTD predicted celltypes to low confidence mapped
HP2_hd$predicted.celltype[is.na(obj$predicted.celltype)] <- "RCTD_lowconf"

# Export predicted for Xenium Explorer
HP2_hd <- RenameCells(HP2_hd, new.names = gsub("_[0-9]+$", "", Cells(HP2_hd)))
# save the annotations for uploading into Xenium Explorer
Idents(HP2_hd) <- "predicted.celltype"
predicted.celltype = as.data.frame(Idents(HP2_hd))
predicted.celltype$cell_id = (rownames(predicted.celltype))
predicted.celltype$group = (predicted.celltype[,1])
predicted.celltype = predicted.celltype[,2:3]
write.csv(predicted.celltype, "HP2_hd_predicted.celltype_forXeniumExplorer.csv", quote=FALSE, row.names = F)

saveRDS(HP2_hd, "HP2_hd_niches.rds")
#################################################################################
#################################################################################
### HT2_hd fov S3065
HT2_hd <-readRDS("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/seurat_HT2_hd_withPredCols.rds") 
HT2_hd <- BuildNicheAssay(object = HT2_hd, fov="S3065", group.by = "predicted.celltype")
Idents(HT2_hd) <- "niches"

# rename NAs in RCTD predicted celltypes to low confidence mapped
HT2_hd$predicted.celltype[is.na(obj$predicted.celltype)] <- "RCTD_lowconf"

# Export predicted for Xenium Explorer
HT2_hd <- RenameCells(HT2_hd, new.names = gsub("_[0-9]+$", "", Cells(HT2_hd)))
# save the annotations for uploading into Xenium Explorer
Idents(HT2_hd) <- "predicted.celltype"
predicted.celltype = as.data.frame(Idents(HT2_hd))
predicted.celltype$cell_id = (rownames(predicted.celltype))
predicted.celltype$group = (predicted.celltype[,1])
predicted.celltype = predicted.celltype[,2:3]
write.csv(predicted.celltype, "HT2_hd_predicted.celltype_forXeniumExplorer.csv", quote=FALSE, row.names = F)

saveRDS(HT2_hd, "HT2_hd_niches.rds")
#################################################################################
### HP2_hd fov S3065
HP2_hd <-readRDS("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/seurat_HP2_hd_withPredCols.rds") 
HP2_hd <- BuildNicheAssay(object = HP2_hd, fov="S3065", group.by = "predicted.celltype")
Idents(HP2_hd) <- "niches"

# rename NAs in RCTD predicted celltypes to low confidence mapped
HP2_hd$predicted.celltype[is.na(obj$predicted.celltype)] <- "RCTD_lowconf"

# Export predicted for Xenium Explorer
HP2_hd <- RenameCells(HP2_hd, new.names = gsub("_[0-9]+$", "", Cells(HP2_hd)))
# save the annotations for uploading into Xenium Explorer
Idents(HP2_hd) <- "predicted.celltype"
predicted.celltype = as.data.frame(Idents(HP2_hd))
predicted.celltype$cell_id = (rownames(predicted.celltype))
predicted.celltype$group = (predicted.celltype[,1])
predicted.celltype = predicted.celltype[,2:3]
write.csv(predicted.celltype, "HP2_hd_predicted.celltype_forXeniumExplorer.csv", quote=FALSE, row.names = F)

saveRDS(HP2_hd, "HP2_hd_niches.rds")
#################################################################################
### HT3 fov S3065
HT3 <-readRDS("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/seurat_HT3_withPredCols.rds") 
HT3 <- BuildNicheAssay(object = HT3, fov="S3065", group.by = "predicted.celltype")
Idents(HT3) <- "niches"

# rename NAs in RCTD predicted celltypes to low confidence mapped
HT3$predicted.celltype[is.na(obj$predicted.celltype)] <- "RCTD_lowconf"

# Export predicted for Xenium Explorer
HT3 <- RenameCells(HT3, new.names = gsub("_[0-9]+$", "", Cells(HT3)))
# save the annotations for uploading into Xenium Explorer
Idents(HT3) <- "predicted.celltype"
predicted.celltype = as.data.frame(Idents(HT3))
predicted.celltype$cell_id = (rownames(predicted.celltype))
predicted.celltype$group = (predicted.celltype[,1])
predicted.celltype = predicted.celltype[,2:3]
write.csv(predicted.celltype, "HT3_predicted.celltype_forXeniumExplorer.csv", quote=FALSE, row.names = F)
Idents(HT3) <- "niches"
niches = as.data.frame(Idents(HT3))
niches$cell_id = (rownames(niches))
niches$group = (niches[,1])
niches = niches[,2:3]
write.csv(niches, "HT3_niches_forXeniumExplorer.csv", quote=FALSE, row.names = F)

saveRDS(HT3, "HT3_niches.rds")

#################################################################################
### HT3 fov HT3.008um ## this is a visiumHD obj
HT3 <-readRDS("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/seurat_HT3_withPredCols.rds") 
HT3 <- BuildNicheAssay(object = HT3, fov="HT3.008um", group.by = "predicted.celltype")
Idents(HT3) <- "niches"

# rename NAs in RCTD predicted celltypes to low confidence mapped
HT3$predicted.celltype[is.na(obj$predicted.celltype)] <- "RCTD_lowconf"

# Export predicted for Xenium Explorer
HT3 <- RenameCells(HT3, new.names = gsub("_[0-9]+$", "", Cells(HT3)))
# save the annotations for uploading into Xenium Explorer
Idents(HT3) <- "predicted.celltype"
predicted.celltype = as.data.frame(Idents(HT3))
predicted.celltype$cell_id = (rownames(predicted.celltype))
predicted.celltype$group = (predicted.celltype[,1])
predicted.celltype = predicted.celltype[,2:3]
write.csv(predicted.celltype, "HT3_predicted.celltype_forXeniumExplorer.csv", quote=FALSE, row.names = F)
Idents(HT3) <- "niches"
niches = as.data.frame(Idents(HT3))
niches$cell_id = (rownames(niches))
niches$group = (niches[,1])
niches = niches[,2:3]
write.csv(niches, "HT3_niches_forXeniumExplorer.csv", quote=FALSE, row.names = F)

saveRDS(HT3, "HT3_niches.rds")
#################################################################################
### HT2_hd fov HT2_hd.008um ## this is a visiumHD obj
HT2_hd <-readRDS("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/seurat_HT2_hd_withPredCols.rds") 
HT2_hd <- BuildNicheAssay(object = HT2_hd, fov="HT2_hd.008um", group.by = "predicted.celltype")
Idents(HT2_hd) <- "niches"

# rename NAs in RCTD predicted celltypes to low confidence mapped
HT2_hd$predicted.celltype[is.na(obj$predicted.celltype)] <- "RCTD_lowconf"

# Export predicted for Xenium Explorer
HT2_hd <- RenameCells(HT2_hd, new.names = gsub("_[0-9]+$", "", Cells(HT2_hd)))
# save the annotations for uploading into Xenium Explorer
Idents(HT2_hd) <- "predicted.celltype"
predicted.celltype = as.data.frame(Idents(HT2_hd))
predicted.celltype$cell_id = (rownames(predicted.celltype))
predicted.celltype$group = (predicted.celltype[,1])
predicted.celltype = predicted.celltype[,2:3]
write.csv(predicted.celltype, "HT2_hd_predicted.celltype_forXeniumExplorer.csv", quote=FALSE, row.names = F)
Idents(HT2_hd) <- "niches"
niches = as.data.frame(Idents(HT2_hd))
niches$cell_id = (rownames(niches))
niches$group = (niches[,1])
niches = niches[,2:3]
write.csv(niches, "HT2_hd_niches_forXeniumExplorer.csv", quote=FALSE, row.names = F)

saveRDS(HT2_hd, "HT2_hd_niches.rds")
#################################################################################
### HP2_hd fov HP2_hd.008um ## this is a visiumHD obj
HP2_hd <-readRDS("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/seurat_HP2_hd_withPredCols.rds") 
HP2_hd <- BuildNicheAssay(object = HP2_hd, fov="HP2.008um", group.by = "predicted.celltype")
Idents(HP2_hd) <- "niches"

# rename NAs in RCTD predicted celltypes to low confidence mapped
HP2_hd$predicted.celltype[is.na(obj$predicted.celltype)] <- "RCTD_lowconf"

# Export predicted for Xenium Explorer
HP2_hd <- RenameCells(HP2_hd, new.names = gsub("_[0-9]+$", "", Cells(HP2_hd)))
# save the annotations for uploading into Xenium Explorer
Idents(HP2_hd) <- "predicted.celltype"
predicted.celltype = as.data.frame(Idents(HP2_hd))
predicted.celltype$cell_id = (rownames(predicted.celltype))
predicted.celltype$group = (predicted.celltype[,1])
predicted.celltype = predicted.celltype[,2:3]
write.csv(predicted.celltype, "HP2_hd_predicted.celltype_forXeniumExplorer.csv", quote=FALSE, row.names = F)
Idents(HP2_hd) <- "niches"
niches = as.data.frame(Idents(HP2_hd))
niches$cell_id = (rownames(niches))
niches$group = (niches[,1])
niches = niches[,2:3]
write.csv(niches, "HP2_hd_niches_forXeniumExplorer.csv", quote=FALSE, row.names = F)

saveRDS(HP2_hd, "HP2_hd_niches.rds")
#################################################################################
### LP1_hd fov LP1_hd.008um ## this is a visiumHD obj
LP1_hd <-readRDS("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/seurat_LP1_hd_withPredCols.rds") 
LP1_hd <- BuildNicheAssay(object = LP1_hd, fov="LP1_hd.008um", group.by = "predicted.celltype")
Idents(LP1_hd) <- "niches"

# rename NAs in RCTD predicted celltypes to low confidence mapped
LP1_hd$predicted.celltype[is.na(LP1_hd$predicted.celltype)] <- "RCTD_lowconf"

# Export predicted for Xenium Explorer
LP1_hd <- RenameCells(LP1_hd, new.names = gsub("_[0-9]+$", "", Cells(LP1_hd)))
# save the annotations for uploading into Xenium Explorer
Idents(LP1_hd) <- "predicted.celltype"
predicted.celltype = as.data.frame(Idents(LP1_hd))
predicted.celltype$cell_id = (rownames(predicted.celltype))
predicted.celltype$group = (predicted.celltype[,1])
predicted.celltype = predicted.celltype[,2:3]
write.csv(predicted.celltype, "LP1_hd_predicted.celltype_forXeniumExplorer.csv", quote=FALSE, row.names = F)
Idents(LP1_hd) <- "niches"
niches = as.data.frame(Idents(LP1_hd))
niches$cell_id = (rownames(niches))
niches$group = (niches[,1])
niches = niches[,2:3]
#setwd("/path/to/local_home/visiumhd_mnp/results/seurat_all_v1")
write.csv(niches, "LP1_hd_niches_forXeniumExplorer.csv", quote=FALSE, row.names = F)

saveRDS(LP1_hd, "LP1_hd_niches.rds")

# annotations are here: /path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations

## visiumhd outs directories are here on SERVER 2: 
# /path/to/local_storage/raid2/project/visiumhd_mnp/customref
# S3057 is HT3
# S11732 is LP1_hd
# S13164 is HT2_hd
# S29033 is HP2_hd









# make_spatial_dimplots_by_niche.R
# Author: Orion for Enrico Barrozo, Ph.D. (BCM/TCH)
# Date: 2025-10-13
# Purpose: Load every Seurat object saved as '*_niches.rds' and export
#          PNG SpatialDimPlot overlays for each FOV, colored by:
#            (1) 'predicted.celltype' (using Enrico's order + palette)
#            (2) 'niches' (from BuildNicheAssay)
# Also emits:
#            (3) predicted.celltype (no lowconf)
#            (4) predicted.celltype (T cell / Treg only)
# Output:  spatial_figs/predicted.celltype/*.png,
#          spatial_figs/predicted.celltype_noLowConf/*.png,
#          spatial_figs/predicted.celltype_TcellTreg/*.png,
#          spatial_figs/niches/*.png

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(methods)
})

# ===================== PATHS =====================
BASE <- "/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations"
setwd(BASE)

OUT_DIR         <- file.path(BASE, "spatial_figs")
OUT_PRED_DIR    <- file.path(OUT_DIR, "predicted.celltype")
OUT_PRED_NOLOW_DIR <- file.path(OUT_DIR, "predicted.celltype_noLowConf")
OUT_PRED_TSUB_DIR  <- file.path(OUT_DIR, "predicted.celltype_TcellTreg")
OUT_NICHES_DIR  <- file.path(OUT_DIR, "niches")
for (d in c(OUT_DIR, OUT_PRED_DIR, OUT_PRED_NOLOW_DIR, OUT_PRED_TSUB_DIR, OUT_NICHES_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ===================== TOGGLES =====================
W_IN     <- 6        # width in inches
H_IN     <- 6        # height in inches
DPI      <- 300
PT_SIZE  <- 1.2      # point size for manual scatter; SpatialDimPlot uses pt.size.factor
SHOW_LEG <- TRUE     # legend visibility
BG_COL   <- "white"  # png background

# ===================== ORDER + PALETTE (from Enrico) =====================
ord_celltypes <- c(
  "Syncytiotrophoblast","Villous Cytotrophoblast","Extravillous Trophoblast",
  "Fibroblast","Perivascular Cell","Placental Endothelial",
  "Hofbauer Cell","Macrophage","NK Cell","Decidual NK Cell","ILC3",
  "Plasma Cell","B Cell","T Cell","Treg","RCTD_lowconf"
)

ct_colors <- c(
  "Syncytiotrophoblast"     = "#5050FF",
  "Villous Cytotrophoblast" = "#CE3D32",
  "Extravillous Trophoblast"= "#749B58",
  "Fibroblast"              = "#F0E685",
  "Perivascular Cell"       = "#466983",
  "Placental Endothelial"   = "#BA6338",
  "Hofbauer Cell"           = "#5DB1DD",
  "Macrophage"              = "#802268",
  "NK Cell"                 = "#6BD76B",
  "Decidual NK Cell"        = "#D595A7",
  "ILC3"                    = "#924822",
  "Plasma Cell"             = "#837B8D",
  "B Cell"                  = "#C75127",
  "T Cell"                  = "#D58F5C",
  "Treg"                    = "#7A65A5",
  "RCTD_lowconf"            = "#BFBFBF"
)

# ===================== HELPERS =====================
# Cells in a given FOV across backends
cells_by_fov <- function(obj, fov) {
  cells <- tryCatch(Seurat::CellsByImage(obj, fov), error = function(e) character(0))
  cells <- intersect(colnames(obj), cells)
  if (length(cells) > 0) return(cells)

  img <- obj@images[[fov]]

  # Visium-style: coordinates rownames are barcodes
  if ("coordinates" %in% slotNames(img)) {
    coords <- try(slot(img, "coordinates"), silent = TRUE)
    if (!inherits(coords, "try-error") && !is.null(coords) && nrow(coords) > 0) {
      cells2 <- intersect(colnames(obj), rownames(coords))
      if (length(cells2) > 0) return(cells2)
    }
  }

  # Xenium: boundaries$centroids@cells
  bnd <- try(img@boundaries[["centroids"]], silent = TRUE)
  if (!inherits(bnd, "try-error") && !is.null(bnd)) {
    cvec <- try(bnd@cells, silent = TRUE)
    if (!inherits(cvec, "try-error") && length(cvec) > 0) {
      cells3 <- intersect(colnames(obj), cvec)
      if (length(cells3) > 0) return(cells3)
    }
  }

  character(0)
}

# Coordinates dataframe (x,y,cell) for a FOV; robust across Xenium/Visium
get_coords_by_fov <- function(obj, fov) {
  # 1) Seurat helper
  coords <- tryCatch(Seurat::GetTissueCoordinates(obj, image = fov), error = function(e) NULL)
  if (!is.null(coords) && nrow(coords) > 0) {
    df <- as.data.frame(coords)
    if (!"cell" %in% names(df)) {
      if (!is.null(rownames(df))) df$cell <- rownames(df)
      else if ("barcode" %in% names(df)) df$cell <- df$barcode
      else if ("Cell" %in% names(df)) df$cell <- df$Cell
    }
    if (!all(c("x","y") %in% names(df))) {
      if (all(c("imagecol","imagerow") %in% names(df))) {
        df$x <- df$imagecol; df$y <- df$imagerow
      } else {
        num_cols <- names(df)[vapply(df, is.numeric, TRUE)]
        if (length(num_cols) >= 2) { df$x <- df[[num_cols[1]]]; df$y <- df[[num_cols[2]]] }
      }
    }
    if (all(c("x","y","cell") %in% names(df))) return(df[, c("x","y","cell")])
  }

  # 2) Visium fallback: coordinates slot
  img <- obj@images[[fov]]
  if ("coordinates" %in% slotNames(img)) {
    coords <- try(slot(img, "coordinates"), silent = TRUE)
    if (!inherits(coords, "try-error") && !is.null(coords) && nrow(coords) > 0) {
      df <- as.data.frame(coords)
      df$cell <- rownames(df)
      if (all(c("imagecol","imagerow") %in% names(df))) {
        df$x <- df$imagecol; df$y <- df$imagerow
      } else {
        num_cols <- names(df)[vapply(df, is.numeric, TRUE)]
        df$x <- df[[num_cols[1]]]; df$y <- df[[num_cols[2]]]
      }
      return(df[, c("x","y","cell")])
    }
  }

  # 3) Xenium fallback: centroids from boundaries
  bnd <- try(img@boundaries[["centroids"]], silent = TRUE)
  if (!inherits(bnd, "try-error") && !is.null(bnd)) {
    bx <- try(slot(bnd, "cells"), silent = TRUE)
    cr <- NULL
    if ("coords" %in% slotNames(bnd)) {
      cr <- as.data.frame(slot(bnd, "coords"))
    } else if ("data" %in% slotNames(bnd)) {
      tmp <- slot(bnd, "data")
      num_cols <- names(tmp)[vapply(tmp, is.numeric, TRUE)]
      if (length(num_cols) >= 2) cr <- as.data.frame(tmp[, num_cols[1:2], drop = FALSE])
    }
    if (!inherits(bx, "try-error") && !is.null(cr)) {
      if (nrow(cr) == length(bx)) {
        df <- cr; df$cell <- as.character(bx)
      } else {
        rn <- rownames(cr)
        if (!is.null(rn)) {
          keep <- intersect(rn, as.character(bx))
          cr2 <- cr[keep, , drop = FALSE]; cr2$cell <- keep; df <- cr2
        } else {
          df <- NULL
        }
      }
      if (!is.null(df)) {
        names(df)[1:2] <- c("x","y")
        return(df[, c("x","y","cell")])
      }
    }
  }

  stop(sprintf("No coordinates available for FOV '%s'", fov))
}

# Build colors aligned to factor levels actually present
make_order_and_colors <- function(vec_chr) {
  present <- unique(vec_chr)
  extras  <- setdiff(present, ord_celltypes)
  levs    <- c(intersect(ord_celltypes, present), sort(extras))
  cols    <- c(unname(ct_colors[intersect(ord_celltypes, present)]),
               rep("gray80", length(extras)))
  list(levels = levs, colors = cols)
}

save_gg <- function(plot, filename, w = W_IN, h = H_IN, dpi = DPI, bg = BG_COL) {
  ggsave(filename, plot, width = w, height = h, dpi = dpi, limitsize = FALSE, bg = bg)
}

# Safer save that retries with no background image; on failure, fallback to manual scatter
safe_ggsave_spatial <- function(build_plot_fn, manual_plot_fn, outfile, w = W_IN, h = H_IN, dpi = DPI, bg = BG_COL) {
  ok <- TRUE
  p1 <- build_plot_fn(image_alpha = NA)
  tryCatch({ save_gg(p1, outfile, w, h, dpi, bg) }, error = function(e) { ok <<- FALSE })
  if (ok) return(invisible(TRUE))

  ok2 <- TRUE
  p2 <- build_plot_fn(image_alpha = 0)
  tryCatch({ save_gg(p2, outfile, w, h, dpi, bg) }, error = function(e) { ok2 <<- FALSE })
  if (ok2) return(invisible(TRUE))

  p3 <- manual_plot_fn()
  save_gg(p3, outfile, w, h, dpi, bg)
  invisible(TRUE)
}

# Manual no-background scatter using per-FOV coordinates (bypasses geom_spatial)
manual_spatial_scatter <- function(obj, fov, group_vec, cols_map, title = NULL, cells_keep = NULL) {
  df <- get_coords_by_fov(obj, fov)
  df$cell <- as.character(if (!is.null(df$cell)) df$cell else rownames(df))
  keep <- intersect(df$cell, colnames(obj))
  if (!is.null(cells_keep)) keep <- intersect(keep, cells_keep)
  df <- df[df$cell %in% keep, , drop = FALSE]
  if (nrow(df) == 0) stop(sprintf("No overlapping cells for FOV '%s'", fov))

  g <- as.character(group_vec[df$cell])
  g[is.na(g) | g == ""] <- "RCTD_lowconf"
  df$group <- factor(g, levels = names(cols_map))

  missing_lvls <- setdiff(levels(df$group), names(cols_map))
  if (length(missing_lvls) > 0) cols_map <- c(cols_map, setNames(rep("gray80", length(missing_lvls)), missing_lvls))

  p <- ggplot(df, aes(x = x, y = y, color = group)) +
    geom_point(size = PT_SIZE) +
    scale_color_manual(values = cols_map, drop = FALSE, na.translate = FALSE) +
    coord_equal() + theme_void() +
    theme(legend.position = if (SHOW_LEG) "right" else "none",
          legend.key.height = grid::unit(5, "pt"),
          legend.key.width  = grid::unit(5, "pt"),
          legend.text = element_text(size = 7),
          legend.title = element_blank())
  if (!is.null(title)) p <- p + ggtitle(title) + theme(plot.title = element_text(size = 10, face = "bold"))
  p
}

# ===================== CORE =====================
spatial_save_pair <- function(obj, sample_id, out_pred, out_niche) {
  # FOVs present?
  fovs <- Images(obj); if (length(fovs) == 0) fovs <- names(obj@images)
  if (length(fovs) == 0) stop("No images/FOVs found in object: ", sample_id)

  # predicted.celltype present?
  if (!"predicted.celltype" %in% colnames(obj@meta.data))
    stop("Object ", sample_id, " lacks 'predicted.celltype' in meta.data")

  # normalize + palette (incl. _hd objects)
  pct_raw <- as.character(obj$predicted.celltype)
  pct_raw[is.na(pct_raw) | pct_raw == ""] <- "RCTD_lowconf"
  ord <- make_order_and_colors(pct_raw)
  obj$predicted.celltype <- factor(pct_raw, levels = ord$levels)
  cols_use <- setNames(ord$colors, ord$levels)

  has_niches <- "niches" %in% colnames(obj@meta.data)

  for (fov in fovs) {
    fn_safe <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)
    fov_tag <- fn_safe(fov)

    # cells for this FOV (backend-aware)
    cells_in_fov <- cells_by_fov(obj, fov)
    if (length(cells_in_fov) == 0) {
      warning(sprintf("%s: no cells mapped to FOV '%s'; skipping.", sample_id, fov))
      next
    }

    ## ---------- predicted.celltype (FULL) ----------
    Idents(obj) <- "predicted.celltype"

    build_pred_plot <- function(image_alpha = NA) {
      SpatialDimPlot(
        obj, images = fov, cells = cells_in_fov,
        pt.size.factor = PT_SIZE, stroke = NA, combine = TRUE, label = FALSE,
        cols = cols_use, image.alpha = image_alpha
      ) +
        ggtitle(paste0(sample_id, " • ", fov, " • predicted.celltype")) +
        theme(plot.title = element_text(size = 10, face = "bold"),
              legend.position = if (SHOW_LEG) "right" else "none",
              legend.key.height = grid::unit(5, "pt"),
              legend.key.width  = grid::unit(5, "pt"),
              legend.text = element_text(size = 7),
              legend.title = element_blank())
    }

    manual_pred_plot <- function() {
      manual_spatial_scatter(
        obj, fov,
        group_vec = obj$predicted.celltype,
        cols_map  = cols_use,
        title     = paste0(sample_id, " • ", fov, " • predicted.celltype")
      )
    }

    out1 <- file.path(OUT_PRED_DIR, paste0(sample_id, "_", fov_tag, "_predicted.celltype.png"))
    safe_ggsave_spatial(build_pred_plot, manual_pred_plot, out1)

    ## ---------- predicted.celltype (NO RCTD_lowconf) ----------
    cells_no_low <- intersect(
      cells_in_fov,
      WhichCells(obj, expression = !is.na(predicted.celltype) & predicted.celltype != "RCTD_lowconf")
    )
    if (length(cells_no_low) > 0) {
      cols_no_low <- cols_use[names(cols_use) != "RCTD_lowconf"]

      build_pred_plot_noLow <- function(image_alpha = NA) {
        SpatialDimPlot(
          obj, images = fov, cells = cells_no_low,
          pt.size.factor = PT_SIZE, stroke = NA, combine = TRUE, label = FALSE,
          cols = cols_no_low, image.alpha = image_alpha
        ) +
          ggtitle(paste0(sample_id, " • ", fov, " • predicted.celltype (no lowconf)")) +
          theme(plot.title = element_text(size = 10, face = "bold"),
                legend.position = if (SHOW_LEG) "right" else "none",
                legend.key.height = grid::unit(5, "pt"),
                legend.key.width  = grid::unit(5, "pt"),
                legend.text = element_text(size = 7),
                legend.title = element_blank())
      }

      manual_pred_plot_noLow <- function() {
        manual_spatial_scatter(
          obj, fov,
          group_vec = obj$predicted.celltype,
          cols_map  = cols_no_low,
          title     = paste0(sample_id, " • ", fov, " • predicted.celltype (no lowconf)"),
          cells_keep = cells_no_low
        )
      }

      out1b <- file.path(OUT_PRED_NOLOW_DIR, paste0(sample_id, "_", fov_tag, "_predicted.celltype_noLowConf.png"))
      safe_ggsave_spatial(build_pred_plot_noLow, manual_pred_plot_noLow, out1b)
    } else {
      message(sprintf("%s/%s: no non-lowconf cells; skipping noLowConf figure.", sample_id, fov))
    }

    ## ---------- predicted.celltype (T Cell + Treg only) ----------
    Idents(obj) <- "predicted.celltype"
    cells_tsub <- intersect(cells_in_fov, WhichCells(obj, idents = c("T Cell", "Treg")))
    if (length(cells_tsub) > 0) {
      t_levels_present <- intersect(c("T Cell","Treg"), as.character(unique(obj$predicted.celltype[cells_tsub])))
      cols_t <- cols_use[t_levels_present]
      if (length(cols_t) == 0) cols_t <- cols_use[c("T Cell","Treg")]

      build_pred_plot_tsub <- function(image_alpha = NA) {
        SpatialDimPlot(
          obj, images = fov, cells = cells_tsub,
          pt.size.factor = PT_SIZE, stroke = NA, combine = TRUE, label = FALSE,
          cols = cols_t, image.alpha = image_alpha
        ) +
          ggtitle(paste0(sample_id, " • ", fov, " • T cell / Treg")) +
          theme(plot.title = element_text(size = 10, face = "bold"),
                legend.position = if (SHOW_LEG) "right" else "none",
                legend.key.height = grid::unit(5, "pt"),
                legend.key.width  = grid::unit(5, "pt"),
                legend.text = element_text(size = 7),
                legend.title = element_blank())
      }

      manual_pred_plot_tsub <- function() {
        manual_spatial_scatter(
          obj, fov,
          group_vec = obj$predicted.celltype,
          cols_map  = cols_t,
          title     = paste0(sample_id, " • ", fov, " • T cell / Treg"),
          cells_keep = cells_tsub
        )
      }

      out1c <- file.path(OUT_PRED_TSUB_DIR, paste0(sample_id, "_", fov_tag, "_predicted.celltype_TcellTreg.png"))
      safe_ggsave_spatial(build_pred_plot_tsub, manual_pred_plot_tsub, out1c)
    } else {
      message(sprintf("%s/%s: no T Cell or Treg cells present; skipping Tcell/Treg figure.", sample_id, fov))
    }

    ## ---------- niches ----------
    if (has_niches) {
      obj$niches <- factor(obj$niches)
      Idents(obj) <- "niches"
      niche_lvls <- levels(obj$niches)
      niche_cols <- setNames(scales::hue_pal()(length(niche_lvls)), niche_lvls)

      build_niche_plot <- function(image_alpha = NA) {
        SpatialDimPlot(
          obj, images = fov, cells = cells_in_fov,
          pt.size.factor = PT_SIZE, stroke = NA, combine = TRUE, label = FALSE,
          image.alpha = image_alpha
        ) +
          ggtitle(paste0(sample_id, " • ", fov, " • niches")) +
          theme(plot.title = element_text(size = 10, face = "bold"),
                legend.position = if (SHOW_LEG) "right" else "none",
                legend.key.height = grid::unit(5, "pt"),
                legend.key.width  = grid::unit(5, "pt"),
                legend.text = element_text(size = 7),
                legend.title = element_blank())
      }

      manual_niche_plot <- function() {
        manual_spatial_scatter(
          obj, fov,
          group_vec = obj$niches,
          cols_map  = niche_cols,
          title     = paste0(sample_id, " • ", fov, " • niches")
        )
      }

      out2 <- file.path(OUT_NICHES_DIR, paste0(sample_id, "_", fov_tag, "_niches.png"))
      safe_ggsave_spatial(build_niche_plot, manual_niche_plot, out2)
    } else {
      warning("Object ", sample_id, " has no 'niches'; skipping niches plots")
    }
  } # end for fov
} # end function

# ===================== MAIN =====================
files_all <- sort(unique(list.files(recursive = FALSE)))
files <- files_all[endsWith(files_all, "_niches.rds")]
if (length(files) == 0) stop("No *_niches.rds files found in ", BASE)

manifest_list <- vector("list", length(files))
for (i in seq_along(files)) {
  rds <- files[i]
  sample_id <- sub("_niches.rds", "", basename(rds), fixed = TRUE)
  message(">>> Processing: ", sample_id)

  obj <- readRDS(rds)
  try(spatial_save_pair(obj, sample_id, OUT_PRED_DIR, OUT_NICHES_DIR), silent = FALSE)

  fovs <- Images(obj); if (length(fovs) == 0) fovs <- names(obj@images)
  manifest_list[[i]] <- data.frame(
    file = rds, sample = sample_id, fov = paste(fovs, collapse = ";"),
    n_cells = ncol(obj), stringsAsFactors = FALSE
  )

  rm(obj); invisible(gc())
}

manifest <- do.call(rbind, manifest_list)
write.csv(manifest, file.path(OUT_DIR, "spatial_plot_manifest.csv"), row.names = FALSE)
message("Done. PNGs saved to:\n - ", OUT_PRED_DIR, "\n - ", OUT_PRED_NOLOW_DIR, "\n - ", OUT_PRED_TSUB_DIR, "\n - ", OUT_NICHES_DIR)




## now redo _hd samples

# make_spatial_dimplots_by_niche.R
# Author: Orion for Enrico Barrozo, Ph.D. (BCM/TCH)
# Date: 2025-10-13
# Purpose: Load every Seurat object saved as '*_niches.rds' and export
#          PNG SpatialDimPlot overlays for each FOV, colored by:
#            (1) 'predicted.celltype' (using Enrico's order + palette)
#            (2) 'niches' (from BuildNicheAssay)
# Also emits:
#            (3) predicted.celltype (no lowconf)
#            (4) predicted.celltype (T cell / Treg only)
# Output:  spatial_figs/predicted.celltype/*.png,
#          spatial_figs/predicted.celltype_noLowConf/*.png,
#          spatial_figs/predicted.celltype_TcellTreg/*.png,
#          spatial_figs/niches/*.png

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(methods)
})

# ===================== PATHS =====================
BASE <- "/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations"
setwd(BASE)

OUT_DIR         <- file.path(BASE, "spatial_figs")
OUT_PRED_DIR    <- file.path(OUT_DIR, "predicted.celltype")
OUT_PRED_NOLOW_DIR <- file.path(OUT_DIR, "predicted.celltype_noLowConf")
OUT_PRED_TSUB_DIR  <- file.path(OUT_DIR, "predicted.celltype_TcellTreg")
OUT_NICHES_DIR  <- file.path(OUT_DIR, "niches")
for (d in c(OUT_DIR, OUT_PRED_DIR, OUT_PRED_NOLOW_DIR, OUT_PRED_TSUB_DIR, OUT_NICHES_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ===================== TOGGLES =====================
W_IN     <- 6        # width in inches
H_IN     <- 6        # height in inches
DPI      <- 300
PT_SIZE  <- 1.2      # point size for manual scatter; SpatialDimPlot uses pt.size.factor
SHOW_LEG <- TRUE     # legend visibility
BG_COL   <- "white"  # png background

# ===================== ORDER + PALETTE (from Enrico) =====================
ord_celltypes <- c(
  "Syncytiotrophoblast","Villous Cytotrophoblast","Extravillous Trophoblast",
  "Fibroblast","Perivascular Cell","Placental Endothelial",
  "Hofbauer Cell","Macrophage","NK Cell","Decidual NK Cell","ILC3",
  "Plasma Cell","B Cell","T Cell","Treg","RCTD_lowconf"
)

ct_colors <- c(
  "Syncytiotrophoblast"     = "#5050FF",
  "Villous Cytotrophoblast" = "#CE3D32",
  "Extravillous Trophoblast"= "#749B58",
  "Fibroblast"              = "#F0E685",
  "Perivascular Cell"       = "#466983",
  "Placental Endothelial"   = "#BA6338",
  "Hofbauer Cell"           = "#5DB1DD",
  "Macrophage"              = "#802268",
  "NK Cell"                 = "#6BD76B",
  "Decidual NK Cell"        = "#D595A7",
  "ILC3"                    = "#924822",
  "Plasma Cell"             = "#837B8D",
  "B Cell"                  = "#C75127",
  "T Cell"                  = "#D58F5C",
  "Treg"                    = "#7A65A5",
  "RCTD_lowconf"            = "#BFBFBF"
)

# ===================== HELPERS =====================
# Subset object to a specific cell vector (robust) and keep only the given FOV image if present
restrict_to_cells <- function(obj, cells, fov) {
  obj2 <- tryCatch(Seurat::subset(obj, cells = cells), error = function(e) obj[, cells, drop = FALSE])
  # keep only target image when possible
  if (!is.null(obj2@images[[fov]])) obj2@images <- obj2@images[fov]
  obj2
}

# Cells in a given FOV across backends
cells_by_fov <- function(obj, fov) {
  cells <- tryCatch(Seurat::CellsByImage(obj, fov), error = function(e) character(0))
  cells <- intersect(colnames(obj), cells)
  if (length(cells) > 0) return(cells)

  img <- obj@images[[fov]]

  # Visium-style: coordinates rownames are barcodes
  if ("coordinates" %in% slotNames(img)) {
    coords <- try(slot(img, "coordinates"), silent = TRUE)
    if (!inherits(coords, "try-error") && !is.null(coords) && nrow(coords) > 0) {
      cells2 <- intersect(colnames(obj), rownames(coords))
      if (length(cells2) > 0) return(cells2)
    }
  }

  # Xenium: boundaries$centroids@cells
  bnd <- try(img@boundaries[["centroids"]], silent = TRUE)
  if (!inherits(bnd, "try-error") && !is.null(bnd)) {
    cvec <- try(bnd@cells, silent = TRUE)
    if (!inherits(cvec, "try-error") && length(cvec) > 0) {
      cells3 <- intersect(colnames(obj), cvec)
      if (length(cells3) > 0) return(cells3)
    }
  }

  character(0)
}

# Coordinates dataframe (x,y,cell) for a FOV; robust across Xenium/Visium
get_coords_by_fov <- function(obj, fov) {
  # 1) Seurat helper
  coords <- tryCatch(Seurat::GetTissueCoordinates(obj, image = fov), error = function(e) NULL)
  if (!is.null(coords) && nrow(coords) > 0) {
    df <- as.data.frame(coords)
    if (!"cell" %in% names(df)) {
      if (!is.null(rownames(df))) df$cell <- rownames(df)
      else if ("barcode" %in% names(df)) df$cell <- df$barcode
      else if ("Cell" %in% names(df)) df$cell <- df$Cell
    }
    if (!all(c("x","y") %in% names(df))) {
      if (all(c("imagecol","imagerow") %in% names(df))) {
        df$x <- df$imagecol; df$y <- df$imagerow
      } else {
        num_cols <- names(df)[vapply(df, is.numeric, TRUE)]
        if (length(num_cols) >= 2) { df$x <- df[[num_cols[1]]]; df$y <- df[[num_cols[2]]] }
      }
    }
    if (all(c("x","y","cell") %in% names(df))) return(df[, c("x","y","cell")])
  }

  # 2) Visium fallback: coordinates slot
  img <- obj@images[[fov]]
  if ("coordinates" %in% slotNames(img)) {
    coords <- try(slot(img, "coordinates"), silent = TRUE)
    if (!inherits(coords, "try-error") && !is.null(coords) && nrow(coords) > 0) {
      df <- as.data.frame(coords)
      df$cell <- rownames(df)
      if (all(c("imagecol","imagerow") %in% names(df))) {
        df$x <- df$imagecol; df$y <- df$imagerow
      } else {
        num_cols <- names(df)[vapply(df, is.numeric, TRUE)]
        df$x <- df[[num_cols[1]]]; df$y <- df[[num_cols[2]]]
      }
      return(df[, c("x","y","cell")])
    }
  }

  # 3) Xenium fallback: centroids from boundaries
  bnd <- try(img@boundaries[["centroids"]], silent = TRUE)
  if (!inherits(bnd, "try-error") && !is.null(bnd)) {
    bx <- try(slot(bnd, "cells"), silent = TRUE)
    cr <- NULL
    if ("coords" %in% slotNames(bnd)) {
      cr <- as.data.frame(slot(bnd, "coords"))
    } else if ("data" %in% slotNames(bnd)) {
      tmp <- slot(bnd, "data")
      num_cols <- names(tmp)[vapply(tmp, is.numeric, TRUE)]
      if (length(num_cols) >= 2) cr <- as.data.frame(tmp[, num_cols[1:2], drop = FALSE])
    }
    if (!inherits(bx, "try-error") && !is.null(cr)) {
      if (nrow(cr) == length(bx)) {
        df <- cr; df$cell <- as.character(bx)
      } else {
        rn <- rownames(cr)
        if (!is.null(rn)) {
          keep <- intersect(rn, as.character(bx))
          cr2 <- cr[keep, , drop = FALSE]; cr2$cell <- keep; df <- cr2
        } else {
          df <- NULL
        }
      }
      if (!is.null(df)) {
        names(df)[1:2] <- c("x","y")
        return(df[, c("x","y","cell")])
      }
    }
  }

  stop(sprintf("No coordinates available for FOV '%s'", fov))
}

# Build colors aligned to factor levels actually present
make_order_and_colors <- function(vec_chr) {
  present <- unique(vec_chr)
  extras  <- setdiff(present, ord_celltypes)
  levs    <- c(intersect(ord_celltypes, present), sort(extras))
  cols    <- c(unname(ct_colors[intersect(ord_celltypes, present)]),
               rep("gray80", length(extras)))
  list(levels = levs, colors = cols)
}

save_gg <- function(plot, filename, w = W_IN, h = H_IN, dpi = DPI, bg = BG_COL) {
  ggsave(filename, plot, width = w, height = h, dpi = dpi, limitsize = FALSE, bg = bg)
}

# Safer save that retries with no background image; on failure, fallback to manual scatter
safe_ggsave_spatial <- function(build_plot_fn, manual_plot_fn, outfile, w = W_IN, h = H_IN, dpi = DPI, bg = BG_COL) {
  ok <- TRUE
  p1 <- build_plot_fn(image_alpha = NA)
  tryCatch({ save_gg(p1, outfile, w, h, dpi, bg) }, error = function(e) { ok <<- FALSE })
  if (ok) return(invisible(TRUE))

  ok2 <- TRUE
  p2 <- build_plot_fn(image_alpha = 0)
  tryCatch({ save_gg(p2, outfile, w, h, dpi, bg) }, error = function(e) { ok2 <<- FALSE })
  if (ok2) return(invisible(TRUE))

  p3 <- manual_plot_fn()
  save_gg(p3, outfile, w, h, dpi, bg)
  invisible(TRUE)
}

# Manual no-background scatter using per-FOV coordinates (bypasses geom_spatial)
manual_spatial_scatter <- function(obj, fov, group_vec, cols_map, title = NULL, cells_keep = NULL) {
  df <- get_coords_by_fov(obj, fov)
  df$cell <- as.character(if (!is.null(df$cell)) df$cell else rownames(df))
  keep <- intersect(df$cell, colnames(obj))
  if (!is.null(cells_keep)) keep <- intersect(keep, cells_keep)
  df <- df[df$cell %in% keep, , drop = FALSE]
  if (nrow(df) == 0) stop(sprintf("No overlapping cells for FOV '%s'", fov))

  g <- as.character(group_vec[df$cell])
  g[is.na(g) | g == ""] <- "RCTD_lowconf"
  df$group <- factor(g, levels = names(cols_map))

  missing_lvls <- setdiff(levels(df$group), names(cols_map))
  if (length(missing_lvls) > 0) cols_map <- c(cols_map, setNames(rep("gray80", length(missing_lvls)), missing_lvls))

  p <- ggplot(df, aes(x = x, y = y, color = group)) +
    geom_point(size = PT_SIZE) +
    scale_color_manual(values = cols_map, drop = FALSE, na.translate = FALSE) +
    coord_equal() + theme_void() +
    theme(legend.position = if (SHOW_LEG) "right" else "none",
          legend.key.height = grid::unit(5, "pt"),
          legend.key.width  = grid::unit(5, "pt"),
          legend.text = element_text(size = 7),
          legend.title = element_blank())
  if (!is.null(title)) p <- p + ggtitle(title) + theme(plot.title = element_text(size = 10, face = "bold"))
  p
}

# ===================== CORE =====================
spatial_save_pair <- function(obj, sample_id, out_pred, out_niche) {
  # FOVs present?
  fovs <- Images(obj); if (length(fovs) == 0) fovs <- names(obj@images)
  if (length(fovs) == 0) stop("No images/FOVs found in object: ", sample_id)

  # predicted.celltype present?
  if (!"predicted.celltype" %in% colnames(obj@meta.data))
    stop("Object ", sample_id, " lacks 'predicted.celltype' in meta.data")

  # normalize + palette (incl. _hd objects)
  pct_raw <- as.character(obj$predicted.celltype)
  pct_raw[is.na(pct_raw) | pct_raw == ""] <- "RCTD_lowconf"
  ord <- make_order_and_colors(pct_raw)
  obj$predicted.celltype <- factor(pct_raw, levels = ord$levels)
  cols_use <- setNames(ord$colors, ord$levels)

  has_niches <- "niches" %in% colnames(obj@meta.data)

  for (fov in fovs) {
    fn_safe <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)
    fov_tag <- fn_safe(fov)

    # cells for this FOV (backend-aware)
    cells_in_fov <- cells_by_fov(obj, fov)
    if (length(cells_in_fov) == 0) {
      warning(sprintf("%s: no cells mapped to FOV '%s'; skipping.", sample_id, fov))
      next
    }

    ## ---------- predicted.celltype (FULL) ----------
    Idents(obj) <- "predicted.celltype"

    build_pred_plot <- function(image_alpha = NA) {
      obj_fov <- restrict_to_cells(obj, cells_in_fov, fov)
      Idents(obj_fov) <- "predicted.celltype"
      SpatialDimPlot(
        obj_fov, images = fov,
        pt.size.factor = PT_SIZE, stroke = NA, combine = TRUE, label = FALSE,
        cols = cols_use, image.alpha = image_alpha
      ) +
        ggtitle(paste0(sample_id, " • ", fov, " • predicted.celltype")) +
        theme(plot.title = element_text(size = 10, face = "bold"),
              legend.position = if (SHOW_LEG) "right" else "none",
              legend.key.height = grid::unit(5, "pt"),
              legend.key.width  = grid::unit(5, "pt"),
              legend.text = element_text(size = 7),
              legend.title = element_blank())
    }

    manual_pred_plot <- function() {
      manual_spatial_scatter(
        obj, fov,
        group_vec = obj$predicted.celltype,
        cols_map  = cols_use,
        title     = paste0(sample_id, " • ", fov, " • predicted.celltype")
      )
    }

    out1 <- file.path(OUT_PRED_DIR, paste0(sample_id, "_", fov_tag, "_predicted.celltype.png"))
    safe_ggsave_spatial(build_pred_plot, manual_pred_plot, out1)

    ## ---------- predicted.celltype (NO RCTD_lowconf) ----------
    cells_no_low <- intersect(
      cells_in_fov,
      WhichCells(obj, expression = !is.na(predicted.celltype) & predicted.celltype != "RCTD_lowconf")
    )
    if (length(cells_no_low) > 0) {
      cols_no_low <- cols_use[names(cols_use) != "RCTD_lowconf"]

      build_pred_plot_noLow <- function(image_alpha = NA) {
        obj_no_low <- restrict_to_cells(obj, cells_no_low, fov)
        Idents(obj_no_low) <- "predicted.celltype"
        SpatialDimPlot(
          obj_no_low, images = fov,
          pt.size.factor = PT_SIZE, stroke = NA, combine = TRUE, label = FALSE,
          cols = cols_no_low, image.alpha = image_alpha
        ) +
          ggtitle(paste0(sample_id, " • ", fov, " • predicted.celltype (no lowconf)")) +
          theme(plot.title = element_text(size = 10, face = "bold"),
                legend.position = if (SHOW_LEG) "right" else "none",
                legend.key.height = grid::unit(5, "pt"),
                legend.key.width  = grid::unit(5, "pt"),
                legend.text = element_text(size = 7),
                legend.title = element_blank())
      }

      manual_pred_plot_noLow <- function() {
        manual_spatial_scatter(
          obj, fov,
          group_vec = obj$predicted.celltype,
          cols_map  = cols_no_low,
          title     = paste0(sample_id, " • ", fov, " • predicted.celltype (no lowconf)"),
          cells_keep = cells_no_low
        )
      }

      out1b <- file.path(OUT_PRED_NOLOW_DIR, paste0(sample_id, "_", fov_tag, "_predicted.celltype_noLowConf.png"))
      safe_ggsave_spatial(build_pred_plot_noLow, manual_pred_plot_noLow, out1b)
    } else {
      message(sprintf("%s/%s: no non-lowconf cells; skipping noLowConf figure.", sample_id, fov))
    }

    ## ---------- predicted.celltype (T Cell + Treg only) ----------
    Idents(obj) <- "predicted.celltype"
    cells_tsub <- intersect(cells_in_fov, WhichCells(obj, idents = c("T Cell", "Treg")))
    if (length(cells_tsub) > 0) {
      t_levels_present <- intersect(c("T Cell","Treg"), as.character(unique(obj$predicted.celltype[cells_tsub])))
      cols_t <- cols_use[t_levels_present]
      if (length(cols_t) == 0) cols_t <- cols_use[c("T Cell","Treg")]

      build_pred_plot_tsub <- function(image_alpha = NA) {
        obj_tsub <- restrict_to_cells(obj, cells_tsub, fov)
        Idents(obj_tsub) <- "predicted.celltype"
        SpatialDimPlot(
          obj_tsub, images = fov,
          pt.size.factor = PT_SIZE, stroke = NA, combine = TRUE, label = FALSE,
          cols = cols_t, image.alpha = image_alpha
        ) +
          ggtitle(paste0(sample_id, " • ", fov, " • T cell / Treg")) +
          theme(plot.title = element_text(size = 10, face = "bold"),
                legend.position = if (SHOW_LEG) "right" else "none",
                legend.key.height = grid::unit(5, "pt"),
                legend.key.width  = grid::unit(5, "pt"),
                legend.text = element_text(size = 7),
                legend.title = element_blank())
      }

      manual_pred_plot_tsub <- function() {
        manual_spatial_scatter(
          obj, fov,
          group_vec = obj$predicted.celltype,
          cols_map  = cols_t,
          title     = paste0(sample_id, " • ", fov, " • T cell / Treg"),
          cells_keep = cells_tsub
        )
      }

      out1c <- file.path(OUT_PRED_TSUB_DIR, paste0(sample_id, "_", fov_tag, "_predicted.celltype_TcellTreg.png"))
      safe_ggsave_spatial(build_pred_plot_tsub, manual_pred_plot_tsub, out1c)
    } else {
      message(sprintf("%s/%s: no T Cell or Treg cells present; skipping Tcell/Treg figure.", sample_id, fov))
    }

    ## ---------- niches ----------
    if (has_niches) {
      obj$niches <- factor(obj$niches)
      Idents(obj) <- "niches"
      niche_lvls <- levels(obj$niches)
      niche_cols <- setNames(scales::hue_pal()(length(niche_lvls)), niche_lvls)

      build_niche_plot <- function(image_alpha = NA) {
        obj_fov_n <- restrict_to_cells(obj, cells_in_fov, fov)
        Idents(obj_fov_n) <- "niches"
        SpatialDimPlot(
          obj_fov_n, images = fov,
          pt.size.factor = PT_SIZE, stroke = NA, combine = TRUE, label = FALSE,
          image.alpha = image_alpha
        ) +
          ggtitle(paste0(sample_id, " • ", fov, " • niches")) +
          theme(plot.title = element_text(size = 10, face = "bold"),
                legend.position = if (SHOW_LEG) "right" else "none",
                legend.key.height = grid::unit(5, "pt"),
                legend.key.width  = grid::unit(5, "pt"),
                legend.text = element_text(size = 7),
                legend.title = element_blank())
      }

      manual_niche_plot <- function() {
        manual_spatial_scatter(
          obj, fov,
          group_vec = obj$niches,
          cols_map  = niche_cols,
          title     = paste0(sample_id, " • ", fov, " • niches")
        )
      }

      out2 <- file.path(OUT_NICHES_DIR, paste0(sample_id, "_", fov_tag, "_niches.png"))
      safe_ggsave_spatial(build_niche_plot, manual_niche_plot, out2)
    } else {
      warning("Object ", sample_id, " has no 'niches'; skipping niches plots")
    }
  } # end for fov
} # end function

# ===================== MAIN =====================
files_all <- sort(unique(list.files(recursive = FALSE)))
files <- files_all[endsWith(files_all, "_niches.rds")]
if (length(files) == 0) stop("No *_niches.rds files found in ", BASE)

manifest_list <- vector("list", length(files))
for (i in seq_along(files)) {
  rds <- files[i]
  sample_id <- sub("_niches.rds", "", basename(rds), fixed = TRUE)
  message(">>> Processing: ", sample_id)

  obj <- readRDS(rds)
  try(spatial_save_pair(obj, sample_id, OUT_PRED_DIR, OUT_NICHES_DIR), silent = FALSE)

  fovs <- Images(obj); if (length(fovs) == 0) fovs <- names(obj@images)
  manifest_list[[i]] <- data.frame(
    file = rds, sample = sample_id, fov = paste(fovs, collapse = ";"),
    n_cells = ncol(obj), stringsAsFactors = FALSE
  )

  rm(obj); invisible(gc())
}

manifest <- do.call(rbind, manifest_list)
write.csv(manifest, file.path(OUT_DIR, "spatial_plot_manifest.csv"), row.names = FALSE)
message("Done. PNGs saved to:
 - ", OUT_PRED_DIR, "
 - ", OUT_PRED_NOLOW_DIR, "
 - ", OUT_PRED_TSUB_DIR, "
 - ", OUT_NICHES_DIR)

# make_spatial_dimplots_by_niche.R
# Author: Orion for Enrico Barrozo, Ph.D. (BCM/TCH)
# Date: 2025-10-13
# Purpose: Load every Seurat object saved as '*_niches.rds' and export
#          PNG SpatialDimPlot overlays for each FOV, colored by:
#            (1) 'predicted.celltype' (using Enrico's order + palette)
#            (2) 'niches' (from BuildNicheAssay)
# Also emits:
#            (3) predicted.celltype (no lowconf)
#            (4) predicted.celltype (T cell / Treg only)
# Output:  spatial_figs/predicted.celltype/*.png,
#          spatial_figs/predicted.celltype_noLowConf/*.png,
#          spatial_figs/predicted.celltype_TcellTreg/*.png,
#          spatial_figs/niches/*.png

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(methods)
})

# ===================== PATHS =====================
BASE <- "/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations"
setwd(BASE)

OUT_DIR         <- file.path(BASE, "spatial_figs")
OUT_PRED_DIR    <- file.path(OUT_DIR, "predicted.celltype")
OUT_PRED_NOLOW_DIR <- file.path(OUT_DIR, "predicted.celltype_noLowConf")
OUT_PRED_TSUB_DIR  <- file.path(OUT_DIR, "predicted.celltype_TcellTreg")
OUT_NICHES_DIR  <- file.path(OUT_DIR, "niches")
for (d in c(OUT_DIR, OUT_PRED_DIR, OUT_PRED_NOLOW_DIR, OUT_PRED_TSUB_DIR, OUT_NICHES_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ===================== HD SAMPLE FILTER =====================
#HD_SAMPLES <- c("HT3", "LP1_hd", "HT2_hd", "HP2_hd")
HD_SAMPLES <- c("LP1_hd")
HD_FOV_MAP <- list(
  LP1_hd = "LP1_hd.008um"
)
#HD_FOV_MAP <- list(
 # HT3    = "HT3.008um",
 # LP1_hd = "LP1_hd.008um",
 # HT2_hd = "HT2_hd.008um",
 # HP2_hd = "HP2.008um"
#)

# ===================== TOGGLES =====================
W_IN     <- 6        # width in inches
H_IN     <- 6        # height in inches
DPI      <- 300
PT_SIZE  <- 1.2      # point size for manual scatter; SpatialDimPlot uses pt.size.factor
SHOW_LEG <- TRUE     # legend visibility
BG_COL   <- "white"  # png background

# ===================== ORDER + PALETTE (from Enrico) =====================
ord_celltypes <- c(
  "Syncytiotrophoblast","Villous Cytotrophoblast","Extravillous Trophoblast",
  "Fibroblast","Perivascular Cell","Placental Endothelial",
  "Hofbauer Cell","Macrophage","NK Cell","Decidual NK Cell","ILC3",
  "Plasma Cell","B Cell","T Cell","Treg","RCTD_lowconf"
)

ct_colors <- c(
  "Syncytiotrophoblast"     = "#5050FF",
  "Villous Cytotrophoblast" = "#CE3D32",
  "Extravillous Trophoblast"= "#749B58",
  "Fibroblast"              = "#F0E685",
  "Perivascular Cell"       = "#466983",
  "Placental Endothelial"   = "#BA6338",
  "Hofbauer Cell"           = "#5DB1DD",
  "Macrophage"              = "#802268",
  "NK Cell"                 = "#6BD76B",
  "Decidual NK Cell"        = "#D595A7",
  "ILC3"                    = "#924822",
  "Plasma Cell"             = "#837B8D",
  "B Cell"                  = "#C75127",
  "T Cell"                  = "#D58F5C",
  "Treg"                    = "#7A65A5",
  "RCTD_lowconf"            = "#BFBFBF"
)

# ===================== HELPERS =====================
# Subset object to a specific cell vector (robust) and keep only the given FOV image if present
restrict_to_cells <- function(obj, cells, fov) {
  obj2 <- tryCatch(Seurat::subset(obj, cells = cells), error = function(e) obj[, cells, drop = FALSE])
  # keep only target image when possible
  if (!is.null(obj2@images[[fov]])) obj2@images <- obj2@images[fov]
  obj2
}

# Cells in a given FOV across backends
cells_by_fov <- function(obj, fov) {
  cells <- tryCatch(Seurat::CellsByImage(obj, fov), error = function(e) character(0))
  cells <- intersect(colnames(obj), cells)
  if (length(cells) > 0) return(cells)

  img <- obj@images[[fov]]

  # Visium-style: coordinates rownames are barcodes
  if ("coordinates" %in% slotNames(img)) {
    coords <- try(slot(img, "coordinates"), silent = TRUE)
    if (!inherits(coords, "try-error") && !is.null(coords) && nrow(coords) > 0) {
      cells2 <- intersect(colnames(obj), rownames(coords))
      if (length(cells2) > 0) return(cells2)
    }
  }

  # Xenium: boundaries$centroids@cells
  bnd <- try(img@boundaries[["centroids"]], silent = TRUE)
  if (!inherits(bnd, "try-error") && !is.null(bnd)) {
    cvec <- try(bnd@cells, silent = TRUE)
    if (!inherits(cvec, "try-error") && length(cvec) > 0) {
      cells3 <- intersect(colnames(obj), cvec)
      if (length(cells3) > 0) return(cells3)
    }
  }

  character(0)
}

# Coordinates dataframe (x,y,cell) for a FOV; robust across Xenium/Visium
get_coords_by_fov <- function(obj, fov) {
  # 1) Seurat helper
  coords <- tryCatch(Seurat::GetTissueCoordinates(obj, image = fov), error = function(e) NULL)
  if (!is.null(coords) && nrow(coords) > 0) {
    df <- as.data.frame(coords)
    if (!"cell" %in% names(df)) {
      if (!is.null(rownames(df))) df$cell <- rownames(df)
      else if ("barcode" %in% names(df)) df$cell <- df$barcode
      else if ("Cell" %in% names(df)) df$cell <- df$Cell
    }
    if (!all(c("x","y") %in% names(df))) {
      if (all(c("imagecol","imagerow") %in% names(df))) {
        df$x <- df$imagecol; df$y <- df$imagerow
      } else {
        num_cols <- names(df)[vapply(df, is.numeric, TRUE)]
        if (length(num_cols) >= 2) { df$x <- df[[num_cols[1]]]; df$y <- df[[num_cols[2]]] }
      }
    }
    if (all(c("x","y","cell") %in% names(df))) return(df[, c("x","y","cell")])
  }

  # 2) Visium fallback: coordinates slot
  img <- obj@images[[fov]]
  if ("coordinates" %in% slotNames(img)) {
    coords <- try(slot(img, "coordinates"), silent = TRUE)
    if (!inherits(coords, "try-error") && !is.null(coords) && nrow(coords) > 0) {
      df <- as.data.frame(coords)
      df$cell <- rownames(df)
      if (all(c("imagecol","imagerow") %in% names(df))) {
        df$x <- df$imagecol; df$y <- df$imagerow
      } else {
        num_cols <- names(df)[vapply(df, is.numeric, TRUE)]
        df$x <- df[[num_cols[1]]]; df$y <- df[[num_cols[2]]]
      }
      return(df[, c("x","y","cell")])
    }
  }

  # 3) Xenium fallback: centroids from boundaries
  bnd <- try(img@boundaries[["centroids"]], silent = TRUE)
  if (!inherits(bnd, "try-error") && !is.null(bnd)) {
    bx <- try(slot(bnd, "cells"), silent = TRUE)
    cr <- NULL
    if ("coords" %in% slotNames(bnd)) {
      cr <- as.data.frame(slot(bnd, "coords"))
    } else if ("data" %in% slotNames(bnd)) {
      tmp <- slot(bnd, "data")
      num_cols <- names(tmp)[vapply(tmp, is.numeric, TRUE)]
      if (length(num_cols) >= 2) cr <- as.data.frame(tmp[, num_cols[1:2], drop = FALSE])
    }
    if (!inherits(bx, "try-error") && !is.null(cr)) {
      if (nrow(cr) == length(bx)) {
        df <- cr; df$cell <- as.character(bx)
      } else {
        rn <- rownames(cr)
        if (!is.null(rn)) {
          keep <- intersect(rn, as.character(bx))
          cr2 <- cr[keep, , drop = FALSE]; cr2$cell <- keep; df <- cr2
        } else {
          df <- NULL
        }
      }
      if (!is.null(df)) {
        names(df)[1:2] <- c("x","y")
        return(df[, c("x","y","cell")])
      }
    }
  }

  stop(sprintf("No coordinates available for FOV '%s'", fov))
}

# Build colors aligned to factor levels actually present
make_order_and_colors <- function(vec_chr) {
  present <- unique(vec_chr)
  extras  <- setdiff(present, ord_celltypes)
  levs    <- c(intersect(ord_celltypes, present), sort(extras))
  cols    <- c(unname(ct_colors[intersect(ord_celltypes, present)]),
               rep("gray80", length(extras)))
  list(levels = levs, colors = cols)
}

save_gg <- function(plot, filename, w = W_IN, h = H_IN, dpi = DPI, bg = BG_COL) {
  ggsave(filename, plot, width = w, height = h, dpi = dpi, limitsize = FALSE, bg = bg)
}

# Safer save that retries with no background image; on failure, fallback to manual scatter
safe_ggsave_spatial <- function(build_plot_fn, manual_plot_fn, outfile, w = W_IN, h = H_IN, dpi = DPI, bg = BG_COL) {
  ok <- TRUE
  p1 <- build_plot_fn(image_alpha = NA)
  tryCatch({ save_gg(p1, outfile, w, h, dpi, bg) }, error = function(e) { ok <<- FALSE })
  if (ok) return(invisible(TRUE))

  ok2 <- TRUE
  p2 <- build_plot_fn(image_alpha = 0)
  tryCatch({ save_gg(p2, outfile, w, h, dpi, bg) }, error = function(e) { ok2 <<- FALSE })
  if (ok2) return(invisible(TRUE))

  p3 <- manual_plot_fn()
  save_gg(p3, outfile, w, h, dpi, bg)
  invisible(TRUE)
}

# Manual no-background scatter using per-FOV coordinates (bypasses geom_spatial)
manual_spatial_scatter <- function(obj, fov, group_vec, cols_map, title = NULL, cells_keep = NULL) {
  df <- get_coords_by_fov(obj, fov)
  df$cell <- as.character(if (!is.null(df$cell)) df$cell else rownames(df))
  keep <- intersect(df$cell, colnames(obj))
  if (!is.null(cells_keep)) keep <- intersect(keep, cells_keep)
  df <- df[df$cell %in% keep, , drop = FALSE]
  if (nrow(df) == 0) stop(sprintf("No overlapping cells for FOV '%s'", fov))

  g <- as.character(group_vec[df$cell])
  g[is.na(g) | g == ""] <- "RCTD_lowconf"
  df$group <- factor(g, levels = names(cols_map))

  missing_lvls <- setdiff(levels(df$group), names(cols_map))
  if (length(missing_lvls) > 0) cols_map <- c(cols_map, setNames(rep("gray80", length(missing_lvls)), missing_lvls))

  p <- ggplot(df, aes(x = x, y = y, color = group)) +
    geom_point(size = PT_SIZE) +
    scale_color_manual(values = cols_map, drop = FALSE, na.translate = FALSE) +
    coord_equal() + theme_void() +
    theme(legend.position = if (SHOW_LEG) "right" else "none",
          legend.key.height = grid::unit(5, "pt"),
          legend.key.width  = grid::unit(5, "pt"),
          legend.text = element_text(size = 7),
          legend.title = element_blank())
  if (!is.null(title)) p <- p + ggtitle(title) + theme(plot.title = element_text(size = 10, face = "bold"))
  p
}

# ===================== CORE =====================
spatial_save_pair <- function(obj, sample_id, out_pred, out_niche, fov_whitelist = NULL) {
  # FOVs present?
  fovs <- Images(obj); if (length(fovs) == 0) fovs <- names(obj@images)
  # Optional: restrict to specific FOV names (e.g., HD .008um)
  if (!is.null(fov_whitelist)) {
    fovs <- intersect(fovs, as.character(fov_whitelist))
  }
  if (length(fovs) == 0) stop("No images/FOVs found in object: ", sample_id)

  # predicted.celltype present?
  if (!"predicted.celltype" %in% colnames(obj@meta.data))
    stop("Object ", sample_id, " lacks 'predicted.celltype' in meta.data")

  # normalize + palette (incl. _hd objects)
  pct_raw <- as.character(obj$predicted.celltype)
  pct_raw[is.na(pct_raw) | pct_raw == ""] <- "RCTD_lowconf"
  ord <- make_order_and_colors(pct_raw)
  obj$predicted.celltype <- factor(pct_raw, levels = ord$levels)
  cols_use <- setNames(ord$colors, ord$levels)

  has_niches <- "niches" %in% colnames(obj@meta.data)

  for (fov in fovs) {
    fn_safe <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)
    fov_tag <- fn_safe(fov)

    # cells for this FOV (backend-aware)
    cells_in_fov <- cells_by_fov(obj, fov)
    if (length(cells_in_fov) == 0) {
      warning(sprintf("%s: no cells mapped to FOV '%s'; skipping.", sample_id, fov))
      next
    }

    ## ---------- predicted.celltype (FULL) ----------
    Idents(obj) <- "predicted.celltype"

    build_pred_plot <- function(image_alpha = NA) {
      obj_fov <- restrict_to_cells(obj, cells_in_fov, fov)
      Idents(obj_fov) <- "predicted.celltype"
      SpatialDimPlot(
        obj_fov, images = fov,
        pt.size.factor = PT_SIZE, stroke = NA, combine = TRUE, label = FALSE,
        cols = cols_use, image.alpha = image_alpha
      ) +
        ggtitle(paste0(sample_id, " • ", fov, " • predicted.celltype")) +
        theme(plot.title = element_text(size = 10, face = "bold"),
              legend.position = if (SHOW_LEG) "right" else "none",
              legend.key.height = grid::unit(5, "pt"),
              legend.key.width  = grid::unit(5, "pt"),
              legend.text = element_text(size = 7),
              legend.title = element_blank())
    }

    manual_pred_plot <- function() {
      manual_spatial_scatter(
        obj, fov,
        group_vec = obj$predicted.celltype,
        cols_map  = cols_use,
        title     = paste0(sample_id, " • ", fov, " • predicted.celltype")
      )
    }

    out1 <- file.path(OUT_PRED_DIR, paste0(sample_id, "_", fov_tag, "_predicted.celltype.png"))
    safe_ggsave_spatial(build_pred_plot, manual_pred_plot, out1)

    ## ---------- predicted.celltype (NO RCTD_lowconf) ----------
    cells_no_low <- intersect(
      cells_in_fov,
      WhichCells(obj, expression = !is.na(predicted.celltype) & predicted.celltype != "RCTD_lowconf")
    )
    if (length(cells_no_low) > 0) {
      cols_no_low <- cols_use[names(cols_use) != "RCTD_lowconf"]

      build_pred_plot_noLow <- function(image_alpha = NA) {
        obj_no_low <- restrict_to_cells(obj, cells_no_low, fov)
        Idents(obj_no_low) <- "predicted.celltype"
        SpatialDimPlot(
          obj_no_low, images = fov,
          pt.size.factor = PT_SIZE, stroke = NA, combine = TRUE, label = FALSE,
          cols = cols_no_low, image.alpha = image_alpha
        ) +
          ggtitle(paste0(sample_id, " • ", fov, " • predicted.celltype (no lowconf)")) +
          theme(plot.title = element_text(size = 10, face = "bold"),
                legend.position = if (SHOW_LEG) "right" else "none",
                legend.key.height = grid::unit(5, "pt"),
                legend.key.width  = grid::unit(5, "pt"),
                legend.text = element_text(size = 7),
                legend.title = element_blank())
      }

      manual_pred_plot_noLow <- function() {
        manual_spatial_scatter(
          obj, fov,
          group_vec = obj$predicted.celltype,
          cols_map  = cols_no_low,
          title     = paste0(sample_id, " • ", fov, " • predicted.celltype (no lowconf)"),
          cells_keep = cells_no_low
        )
      }

      out1b <- file.path(OUT_PRED_NOLOW_DIR, paste0(sample_id, "_", fov_tag, "_predicted.celltype_noLowConf.png"))
      safe_ggsave_spatial(build_pred_plot_noLow, manual_pred_plot_noLow, out1b)
    } else {
      message(sprintf("%s/%s: no non-lowconf cells; skipping noLowConf figure.", sample_id, fov))
    }

    ## ---------- predicted.celltype (T Cell + Treg only) ----------
    Idents(obj) <- "predicted.celltype"
    cells_tsub <- intersect(cells_in_fov, WhichCells(obj, idents = c("T Cell", "Treg")))
    if (length(cells_tsub) > 0) {
      t_levels_present <- intersect(c("T Cell","Treg"), as.character(unique(obj$predicted.celltype[cells_tsub])))
      cols_t <- cols_use[t_levels_present]
      if (length(cols_t) == 0) cols_t <- cols_use[c("T Cell","Treg")]

      build_pred_plot_tsub <- function(image_alpha = NA) {
        obj_tsub <- restrict_to_cells(obj, cells_tsub, fov)
        Idents(obj_tsub) <- "predicted.celltype"
        SpatialDimPlot(
          obj_tsub, images = fov,
          pt.size.factor = PT_SIZE, stroke = NA, combine = TRUE, label = FALSE,
          cols = cols_t, image.alpha = image_alpha
        ) +
          ggtitle(paste0(sample_id, " • ", fov, " • T cell / Treg")) +
          theme(plot.title = element_text(size = 10, face = "bold"),
                legend.position = if (SHOW_LEG) "right" else "none",
                legend.key.height = grid::unit(5, "pt"),
                legend.key.width  = grid::unit(5, "pt"),
                legend.text = element_text(size = 7),
                legend.title = element_blank())
      }

      manual_pred_plot_tsub <- function() {
        manual_spatial_scatter(
          obj, fov,
          group_vec = obj$predicted.celltype,
          cols_map  = cols_t,
          title     = paste0(sample_id, " • ", fov, " • T cell / Treg"),
          cells_keep = cells_tsub
        )
      }

      out1c <- file.path(OUT_PRED_TSUB_DIR, paste0(sample_id, "_", fov_tag, "_predicted.celltype_TcellTreg.png"))
      safe_ggsave_spatial(build_pred_plot_tsub, manual_pred_plot_tsub, out1c)
    } else {
      message(sprintf("%s/%s: no T Cell or Treg cells present; skipping Tcell/Treg figure.", sample_id, fov))
    }

    ## ---------- niches ----------
    if (has_niches) {
      obj$niches <- factor(obj$niches)
      Idents(obj) <- "niches"
      niche_lvls <- levels(obj$niches)
      niche_cols <- setNames(scales::hue_pal()(length(niche_lvls)), niche_lvls)

      build_niche_plot <- function(image_alpha = NA) {
        obj_fov_n <- restrict_to_cells(obj, cells_in_fov, fov)
        Idents(obj_fov_n) <- "niches"
        SpatialDimPlot(
          obj_fov_n, images = fov,
          pt.size.factor = PT_SIZE, stroke = NA, combine = TRUE, label = FALSE,
          image.alpha = image_alpha
        ) +
          ggtitle(paste0(sample_id, " • ", fov, " • niches")) +
          theme(plot.title = element_text(size = 10, face = "bold"),
                legend.position = if (SHOW_LEG) "right" else "none",
                legend.key.height = grid::unit(5, "pt"),
                legend.key.width  = grid::unit(5, "pt"),
                legend.text = element_text(size = 7),
                legend.title = element_blank())
      }

      manual_niche_plot <- function() {
        manual_spatial_scatter(
          obj, fov,
          group_vec = obj$niches,
          cols_map  = niche_cols,
          title     = paste0(sample_id, " • ", fov, " • niches")
        )
      }

      out2 <- file.path(OUT_NICHES_DIR, paste0(sample_id, "_", fov_tag, "_niches.png"))
      safe_ggsave_spatial(build_niche_plot, manual_niche_plot, out2)
    } else {
      warning("Object ", sample_id, " has no 'niches'; skipping niches plots")
    }
  } # end for fov
} # end function

# ===================== MAIN =====================
files_all <- sort(unique(list.files(recursive = FALSE)))
files <- files_all[endsWith(files_all, "_niches.rds")]
# Only run on HD_SAMPLES
files <- files[sub("_niches.rds", "", basename(files), fixed = TRUE) %in% HD_SAMPLES]
if (length(files) == 0) stop("No *_niches.rds files found in ", BASE)

manifest_list <- vector("list", length(files))
for (i in seq_along(files)) {
  rds <- files[i]
  sample_id <- sub("_niches.rds", "", basename(rds), fixed = TRUE)
  message(">>> Processing: ", sample_id)

  obj <- readRDS(rds)
  forced_fov <- HD_FOV_MAP[[sample_id]]
  try(spatial_save_pair(obj, sample_id, OUT_PRED_DIR, OUT_NICHES_DIR, fov_whitelist = forced_fov), silent = FALSE)

  fovs <- Images(obj); if (length(fovs) == 0) fovs <- names(obj@images)
  manifest_list[[i]] <- data.frame(
    file = rds, sample = sample_id, fov = paste(fovs, collapse = ";"),
    n_cells = ncol(obj), stringsAsFactors = FALSE
  )

  rm(obj); invisible(gc())
}

manifest <- do.call(rbind, manifest_list)
write.csv(manifest, file.path(OUT_DIR, "spatial_plot_manifest.csv"), row.names = FALSE)
message("Done. (HD-only) PNGs saved to:
 - ", OUT_PRED_DIR, "
 - ", OUT_PRED_NOLOW_DIR, "
 - ", OUT_PRED_TSUB_DIR, "
 - ", OUT_NICHES_DIR)



HP2_hd <-readRDS("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations/HP2_hd_niches.rds") 
Idents(HP2_hd) <- "predicted.celltype"
HP2_hd$predicted.celltype[is.na(HP2_hd$predicted.celltype)] <- "RCTD_lowconf"

# object_subset <- subset(HP2_hd, cells = Cells(HP2_hd[["HT3.008um"]]), downsample = 10000)

Idents(HP2_hd) <- "predicted.celltype"
levels(HP2_hd)


> head(HP2_hd@images[["HP2.008um"]]@boundaries[["centroids"]]@cells, n=5)
[1] "s_008um_00269_00526-1" "s_008um_00547_00611-1" "s_008um_00693_00628-1"
[4] "s_008um_00260_00253-1" "s_008um_00433_00599-1"
> head(colnames(HP2_hd),n=5)
[1] "s_008um_00269_00526-1" "s_008um_00547_00611-1" "s_008um_00693_00628-1"
[4] "s_008um_00260_00253-1" "s_008um_00433_00599-1"
head(HP2_hd@images[["HP2.008um"]]@boundaries[["centroids"]]@coords,n=5)
> levels(HP2_hd)
 [1] "RCTD_lowconf"             "Syncytiotrophoblast"      "Villous Cytotrophoblast" 
 [4] "Extravillous Trophoblast" "T Cell"                   "NK Cell"                 
 [7] "Macrophage"               "Hofbauer Cell"            "Fibroblast"              
[10] "Plasma Cell"              "Perivascular Cell"        "B Cell"                  
[13] "Placental Endothelial"    "Decidual NK Cell"         "ILC3"
head(HP13b_xen@assays[["Xenium.500"]]@layers[["counts"]],n=5)
p <- SpatialDimPlot(HP2_hd,group.by = "predicted.celltype") 
ggsave("HP13b_xen_SpatialDimPlot_predicted.celltype.pdf", plot = p, device = "pdf", width = 10, height = 8, units = "in")


Idents(obj) <- "predicted.celltype"

Idents(obj) <- "predicted.celltype"
celltypes <- levels(obj)
celltypes.2 <- # remove the RCTD_lowconf from the list of celltypes 
obj <- subset(x = obj, idents = c(celltypes.2))

HP13b_xen@images[["S23b"]]@boundaries[["centroids"]]@cells
