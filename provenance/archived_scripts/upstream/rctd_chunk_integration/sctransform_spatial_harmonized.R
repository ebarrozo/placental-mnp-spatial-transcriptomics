## sctransform_spatial_harmonized.R

set.seed(1) # set seed for reproducibility
options(future.globals.maxSize = 400e9) ## 300 gb ram (out of 440 on server 3; adjust if on different server)

setwd("/path/to/local_storage/raid3/project/seurat_all_ref_merged/with_annotations")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(stringr)
})

## =========================
## Config
## =========================
set.seed(1)
options(Seurat.object.validate = FALSE)   # quiet FOV/Seurat warnings

# Chunking / memory
CHUNK_SIZE      <- 30000
KEEP_ALL_GENES  <- TRUE       # SCT model keeps all genes (ensures features exist in models)
OVERWRITE_CHUNK <- FALSE      # set TRUE to force re-run per-chunk

# Assay preference to pick source counts
ASSAY_PREF      <- c("Spatial.055um","Spatial.008um","Xenium.5k","Xenium.500","RNA","SCT")

# Paths
IN_PATTERN      <- "^seurat_.*_withPredCols\\.rds$"
OUTDIR_CHUNKS   <- "SCT_chunks"
MANIFEST_CSV    <- file.path(OUTDIR_CHUNKS, "chunk_manifest.csv")

# RCTD reference (assay = RNA), include in pipeline as another “sample”
REF_PATH        <- "/path/to/local_storage/raid2/project/spatial_MNP/placental_ref_Jul_2025.rds"
REF_ID          <- sub("\\.rds$","", basename(REF_PATH))   # "placental_ref_Jul_2025"

# Integration / mapping
OUTDIR_MAPPED   <- "SCT_chunks_mapped"
NFEATURES_INTEG <- 3000
DIMS            <- 1:30
EXPORT_EMBEDDINGS <- TRUE     # write one big table with UMAP + meta for all chunks
EMBED_OUT       <- file.path(OUTDIR_MAPPED, "all_embeddings.tsv.gz")

dir.create(OUTDIR_CHUNKS, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTDIR_MAPPED, showWarnings = FALSE, recursive = TRUE)

## =========================
## Helpers
## =========================
sample_id_from_file <- function(f) {
  bn <- basename(f)
  bn <- sub("^seurat_", "", bn)
  sub("_withPredCols\\.rds$", "", bn)
}
pick_source_assay <- function(obj) {
  cand <- intersect(ASSAY_PREF, names(obj@assays))
  if (!length(cand)) stop("No usable assay found")
  cand[1]
}
# Visium/RNA are UMI; Xenium not
should_correct_umi <- function(assay_name) assay_name %in% c("Spatial.055um","Spatial.008um","RNA")

subset_positive_counts_if_present <- function(obj, assay_name) {
  colname <- paste0("nCount_", assay_name)
  if (colname %in% colnames(obj@meta.data)) {
    keep <- which(obj@meta.data[[colname]] > 0)
    if (length(keep) > 0 && length(keep) < ncol(obj)) {
      obj <- subset(obj, cells = colnames(obj)[keep])
    }
  }
  obj
}
strip_to_SCT_only <- function(obj) {
  DefaultAssay(obj) <- "SCT"
  obj@assays     <- obj@assays["SCT"]
  obj@graphs     <- list()
  obj@reductions <- list()
  obj@images     <- list()   # drop Visium/Xenium images
  obj
}
prefix_cells <- function(obj, sid, chunk_i) {
  new_names <- paste0(sid, "_c", sprintf("%03d", chunk_i), "_", colnames(obj))
  # ensure uniqueness just in case
  if (anyDuplicated(new_names)) new_names <- make.unique(new_names, sep = "_")
  colnames(obj) <- new_names
  obj
}

sctransform_chunk <- function(chunk_obj, src, do_correct_umi, keep_all_genes = TRUE) {
  DefaultAssay(chunk_obj) <- src
  all_genes <- rownames(chunk_obj)
  chunk_obj <- SCTransform(
    object                = chunk_obj,
    assay                 = src,
    new.assay.name        = "SCT",
    vst.flavor            = "v2",
    residual.features     = all_genes,
    return.only.var.genes = !keep_all_genes,
    do.correct.umi        = do_correct_umi,
    conserve.memory       = TRUE,
    verbose               = TRUE
  )
  DefaultAssay(chunk_obj) <- "SCT"
  strip_to_SCT_only(chunk_obj)
}

## =========================
## 1) CHUNK + SCT: samples
## =========================
in_files <- list.files(pattern = IN_PATTERN, full.names = TRUE)
stopifnot(length(in_files) > 0)

manifest <- list()

process_one_object <- function(sid, obj, force_src_assay = NULL) {
  # big memory saver: drop images immediately
  obj@images <- list()

  src <- if (is.null(force_src_assay)) pick_source_assay(obj) else force_src_assay
  obj <- subset_positive_counts_if_present(obj, src)

  n_cells <- ncol(obj)
  if (n_cells == 0) {
    cat(sprintf("[%s] No cells after filtering; skip.\n", sid))
    return(invisible(NULL))
  }
  idx <- seq_len(n_cells)
  chunks <- if (n_cells <= CHUNK_SIZE) list(idx) else split(idx, ceiling(idx / CHUNK_SIZE))

  cat(sprintf("[%s] Cells: %d | Chunks: %d | Source: %s | UMI-correct: %s\n",
              sid, n_cells, length(chunks), src, ifelse(should_correct_umi(src), "TRUE","FALSE")))

  for (i in seq_along(chunks)) {
    out_file <- file.path(OUTDIR_CHUNKS, sprintf("sctrans_%s_chunk_%03d_SCTonly.rds", sid, i))
    if (!OVERWRITE_CHUNK && file.exists(out_file)) {
      cat(sprintf("  - [%s] chunk %03d exists; skip\n", sid, i))
      tmp <- readRDS(out_file); nc <- ncol(tmp); rm(tmp); gc()
      manifest[[length(manifest)+1]] <- data.frame(sample=sid, chunk=i, cells=nc, file=out_file)
      next
    }
    chunk_cells <- colnames(obj)[chunks[[i]]]
    chunk_obj   <- subset(obj, cells = chunk_cells)
    chunk_obj   <- prefix_cells(chunk_obj, sid, i)

    cat(sprintf("  - [%s] chunk %03d: %d cells | %d genes\n", sid, i, ncol(chunk_obj), nrow(chunk_obj)))

    do_umi <- should_correct_umi(src)
    ok <- TRUE; fit <- NULL; err <- NULL
    fit <- tryCatch(
      sctransform_chunk(chunk_obj, src, do_umi, KEEP_ALL_GENES),
      error = function(e) { ok <<- FALSE; err <<- e$message; NULL }
    )
    if (!ok) {
      message(sprintf("    retry SCT with do.correct.umi=FALSE due to: %s", err))
      fit <- sctransform_chunk(chunk_obj, src, FALSE, KEEP_ALL_GENES)
    }

    saveRDS(fit, out_file)
    cat(sprintf("    saved -> %s\n", out_file))
    manifest[[length(manifest)+1]] <- data.frame(sample=sid, chunk=i, cells=ncol(fit), file=out_file)

    rm(chunk_obj, fit); gc()
  }
}

# — samples
for (f in in_files) {
  sid <- sample_id_from_file(f)
  cat(sprintf("\n[%s] Loading %s\n", sid, f)); flush.console()
  obj <- readRDS(f)
  process_one_object(sid, obj)   # auto pick source assay
  rm(obj); gc()
}

## =========================
## 2) CHUNK + SCT: RCTD reference (assay = RNA)
## =========================
if (file.exists(REF_PATH)) {
  cat(sprintf("\n[REF %s] Loading %s\n", REF_ID, REF_PATH))
  ref_obj <- readRDS(REF_PATH)
  # force source assay = RNA for the reference
  if (!"RNA" %in% names(ref_obj@assays)) stop("Reference lacks RNA assay.")
  process_one_object(REF_ID, ref_obj, force_src_assay = "RNA")
  rm(ref_obj); gc()
} else {
  warning(sprintf("Reference file not found: %s (skipping inclusion)", REF_PATH))
}

## =========================
## 3) Write manifest
## =========================
## --- paths (works from with_annotations/ or with_annotations/SCT_chunks) ---
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(stringr)
})

## ---------- Paths / constants ----------
HERE <- getwd()
SCT_DIR <- if (basename(HERE) == "SCT_chunks") "." else "SCT_chunks"
OUTDIR_MAPPED <- if (SCT_DIR == ".") "../SCT_chunks_mapped" else "SCT_chunks_mapped"
dir.create(OUTDIR_MAPPED, showWarnings = FALSE, recursive = TRUE)

MANIFEST_CSV   <- file.path(SCT_DIR, "chunk_manifest.csv")
REF_ID_PATTERN <- "^placental_ref_Jul_2025$"  # keep this RNA reference in the integration ref
NFEATURES_INTEG <- 3000
DIMS <- 1:30

## make anchors lighter/safer
MAX_REF_CELLS_PER_OBJ <- 10000  # downsample per ref chunk for anchors
N_REF_DATASETS        <- 4      # include RNA ref + top-3 largest non-Xenium
K_ANCHOR              <- 5
K_FILTER              <- 150

options(Seurat.object.validate = FALSE)
set.seed(1)

## ---------- (1) Manifest from on-disk chunks ----------
if (!file.exists(MANIFEST_CSV)) {
  fns <- list.files(SCT_DIR, pattern = "^sctrans_.*_chunk_\\d+_SCTonly\\.rds$", full.names = TRUE)
  stopifnot(length(fns) > 0)

  parse_sc <- function(bn) {
    m <- regexec("^sctrans_(.+)_chunk_(\\d+)_SCTonly\\.rds$", bn)
    p <- regmatches(bn, m)[[1]]
    list(sample = p[2], chunk = as.integer(p[3]))
  }
  man <- lapply(fns, function(fp) {
    sc <- parse_sc(basename(fp))
    o  <- readRDS(fp); nc <- ncol(o); rm(o); gc()
    data.frame(sample = sc$sample, chunk = sc$chunk, cells = nc, file = fp, stringsAsFactors = FALSE)
  }) |> bind_rows() |> arrange(sample, chunk)

  write.csv(man, MANIFEST_CSV, row.names = FALSE)
  cat("Manifest built ->", MANIFEST_CSV, "\n")
} else {
  man <- read.csv(MANIFEST_CSV, stringsAsFactors = FALSE)
}

## ---------- (2) Pick reference set (non-Xenium + RNA ref), one chunk/sample ----------
man$platform <- ifelse(grepl("xen", man$sample, ignore.case = TRUE), "Xenium", "Visium_RNA")

cand <- man |>
  filter(platform != "Xenium" | grepl(REF_ID_PATTERN, sample)) |>
  group_by(sample) |>
  slice_max(order_by = cells, n = 1, with_ties = FALSE) |>
  ungroup()

stopifnot(nrow(cand) > 0)

# always include the RNA ref + top non-Xenium by size
ref_ids_force <- which(grepl(REF_ID_PATTERN, cand$sample))
ref_order     <- cand |> arrange(desc(cells))
ref_pick_idx  <- unique(c(ref_ids_force, seq_len(min(nrow(ref_order), N_REF_DATASETS))))
ref_rows      <- ref_order[ref_pick_idx, , drop = FALSE]

cat("Reference datasets:\n")
print(ref_rows[, c("sample","cells")], row.names = FALSE)

ref_list <- lapply(ref_rows$file, readRDS)
names(ref_list) <- paste0(ref_rows$sample, "_c", sprintf("%03d", ref_rows$chunk))

## ---------- (3) Downsample each ref chunk (≤ MAX_REF_CELLS_PER_OBJ) ----------
downsample_obj <- function(o, max_n = 10000) {
  DefaultAssay(o) <- "SCT"
  if (ncol(o) > max_n) {
    keep <- sample(colnames(o), max_n)
    o <- subset(o, cells = keep)
  }
  DefaultAssay(o) <- "SCT"
  o
}
ref_small <- lapply(ref_list, downsample_obj, max_n = MAX_REF_CELLS_PER_OBJ)

## ---------- (4) Robust shared features; pre-fill SCT scale layer ----------
genes_data_by <- lapply(ref_small, function(o) rownames(LayerData(o[["SCT"]], layer = "data")))
common_genes  <- Reduce(intersect, genes_data_by)
cat("Common genes across reference (SCT data):", length(common_genes), "\n")

features_ref <- SelectIntegrationFeatures(ref_small, nfeatures = NFEATURES_INTEG)
features_ref <- intersect(features_ref, common_genes)

if (length(features_ref) < 2000) {
  warning(sprintf("Anchor feature set small (n=%d). Falling back to full intersection (n=%d).",
                  length(features_ref), length(common_genes)))
  features_ref <- common_genes
}
if (length(features_ref) < 200) {
  stop(sprintf("Too few shared genes (n=%d). Harmonize gene names and retry.", length(features_ref)))
}

for (i in seq_along(ref_small)) {
  ref_small[[i]] <- ScaleData(ref_small[[i]], features = features_ref, verbose = FALSE)  # write SCT scale layer rows
}

## ---------- (5) Reference-based SCT integration (no all-pairs) ----------
ref_small <- PrepSCTIntegration(ref_small, anchor.features = features_ref, assay = "SCT")

# treat these as the reference datasets (reference-based integration avoids O(n^2) pairs)
ref_indices <- seq_along(ref_small)

ref_anchors <- FindIntegrationAnchors(
  object.list          = ref_small,
  normalization.method = "SCT",
  anchor.features      = features_ref,
  dims                 = DIMS,
  k.anchor             = K_ANCHOR,
  k.filter             = K_FILTER,
  reference            = ref_indices
)

ref_integrated <- IntegrateData(anchorset = ref_anchors, normalization.method = "SCT")
DefaultAssay(ref_integrated) <- "integrated"
ref_integrated <- ScaleData(ref_integrated, verbose = FALSE)
ref_integrated <- RunPCA(ref_integrated, dims = DIMS, verbose = FALSE)
ref_integrated <- RunUMAP(ref_integrated, dims = DIMS, return.model = TRUE, verbose = FALSE)

ref_path <- file.path(OUTDIR_MAPPED, "integrated_reference_SCT.rds")
saveRDS(ref_integrated, ref_path)
cat("Reference saved ->", ref_path, "\n")




## ---------- (6) Map ALL chunks (incl. Xenium) onto the reference ----------
## ---------- (6) Map ALL chunks (incl. Xenium) onto the reference ----------
all_files <- man$file
dir.create(OUTDIR_MAPPED, showWarnings = FALSE, recursive = TRUE)

# recover features_ref if not in env (best-effort fallback)
if (!exists("features_ref") || !length(features_ref)) {
  features_ref <- tryCatch(
    rownames(LayerData(ref_integrated[["integrated"]], layer = "scale.data")),
    error = function(...) character(0)
  )
  if (!length(features_ref)) {
    features_ref <- rownames(ref_integrated)  # coarse fallback
  }
}

for (ff in all_files) {
  outfile <- file.path(OUTDIR_MAPPED, sub("\\.rds$", "_mapped.rds", basename(ff)))
  if (file.exists(outfile)) {
    message("Skipping existing mapped: ", basename(outfile))
    next
  }

  q <- readRDS(ff)
  DefaultAssay(q) <- "SCT"

  # ensure the query has scale rows for the anchor features (avoid any re-normalization)
  feats_q <- intersect(features_ref, rownames(q))
  if (length(feats_q) < 200) {
    # tighter fallback: use SCT data-layer genes
    feats_q <- intersect(features_ref, rownames(LayerData(q[["SCT"]], layer = "data")))
  }
  if (length(feats_q) < 200) {
    warning("Very small overlap of features for ", basename(ff), " (n=", length(feats_q), ").")
  }
  q <- ScaleData(q, features = feats_q, verbose = FALSE)

  # anchors: critical -> recompute.residuals = FALSE (we stripped counts assays)
  anchors_q <- FindTransferAnchors(
    reference            = ref_integrated,
    query                = q,
    normalization.method = "SCT",
    reference.reduction  = "pca",
    dims                 = DIMS,
    features             = feats_q,
    recompute.residuals  = FALSE
  )

  # MapQuery: DO NOT pass recompute.residuals here in v5
  q_mapped <- MapQuery(
    anchorset            = anchors_q,
    query                = q,
    reference            = ref_integrated,
    refdata              = NULL,      # add a named list if you want to transfer labels
    reference.reduction  = "pca",
    reduction.model      = "umap",
    verbose              = FALSE
  )

  saveRDS(q_mapped, outfile)
  message("Saved mapped -> ", outfile)

  rm(q, q_mapped); gc()
}
message("Mapped chunks saved -> ", OUTDIR_MAPPED)




suppressPackageStartupMessages(library(Seurat))
# optional: make sure we’re fully single-threaded
suppressWarnings({
  if (requireNamespace("future", quietly = TRUE)) {
    future::plan(future::sequential)
    options(future.rng.onMisuse = "ignore")
  }
})

## ---------- (6) Map ALL chunks (incl. Xenium) to the reference via IntegrateEmbeddings + ProjectUMAP ----------
all_files <- man$file
dir.create(OUTDIR_MAPPED, showWarnings = FALSE, recursive = TRUE)

# recover anchor features from the reference if 'features_ref' isn’t in scope
if (!exists("features_ref") || !length(features_ref)) {
  features_ref <- tryCatch(
    rownames(LayerData(ref_integrated[["integrated"]], layer = "scale.data")),
    error = function(...) character(0)
  )
  if (!length(features_ref)) stop("Could not recover features_ref; re-run the reference build block to define it.")
}

# sanity: saved UMAP model must be present
stopifnot(!is.null(ref_integrated[["umap"]]) && !is.null(ref_integrated[["umap"]]@misc$model))

DIMS        <- 1:30
K_WEIGHT    <- 50   # lighter than default; reduces memory
N_FAIL      <- 0
fail_log    <- file.path(OUTDIR_MAPPED, "mapping_failures.log")

for (ff in all_files) {
  outfile <- file.path(OUTDIR_MAPPED, sub("\\.rds$", "_mapped.rds", basename(ff)))
  if (file.exists(outfile)) { message("Skipping existing mapped: ", basename(outfile)); next }

  ok <- TRUE
  msg <- NULL
  q_mapped <- NULL

  tryCatch({
    q <- readRDS(ff)
    DefaultAssay(q) <- "SCT"

    # prefill scale rows so Seurat never tries to recompute residuals
    feats_q <- intersect(features_ref, rownames(q))
    if (length(feats_q) < 200) feats_q <- intersect(features_ref, rownames(LayerData(q[["SCT"]], layer = "data")))
    q <- ScaleData(q, features = feats_q, verbose = FALSE)

    # anchors: critical -> recompute.residuals = FALSE (we stripped counts)
    anchors_q <- FindTransferAnchors(
      reference            = ref_integrated,
      query                = q,
      normalization.method = "SCT",
      reference.reduction  = "pca",
      dims                 = DIMS,
      features             = feats_q,
      recompute.residuals  = FALSE
    )

    # 1) project into reference PCA
    q <- IntegrateEmbeddings(
      anchorset           = anchors_q,
      reference           = ref_integrated,
      new.reduction.name  = "ref.pca",
      reductions          = "pca",
      dims                = DIMS,
      k.weight            = K_WEIGHT
    )

    # 2) project UMAP using the saved model (no neighbors built)
    q <- ProjectUMAP(
      object     = q,
      umap.model = ref_integrated[["umap"]]@misc$model,
      reduction  = "ref.pca"
    )

    q_mapped <- q
  }, error = function(e) {
    ok <<- FALSE
    msg <<- conditionMessage(e)
  })

  if (!ok) {
    N_FAIL <- N_FAIL + 1
    cat(sprintf("%s\tFAILED\t%s\n", basename(ff), msg), file = fail_log, append = TRUE)
    message("Failed on ", basename(ff), " -> ", msg)
    next
  }

  saveRDS(q_mapped, outfile)
  message("Saved mapped -> ", outfile)
  rm(q_mapped); gc()
}

if (N_FAIL > 0) {
  message("Completed with ", N_FAIL, " failures. See: ", fail_log)
} else {
  message("All chunks mapped -> ", OUTDIR_MAPPED)
}


### debug 

suppressPackageStartupMessages(library(Seurat))

## Assumes you already have:
##   - ref_integrated (RunUMAP(..., return.model=TRUE) was done earlier)
##   - features_ref   (the features used to build the reference)
##   - DIMS <- 1:30   (or your chosen dims)

stopifnot(!is.null(ref_integrated[["pca"]]))
DIMS       <- DIMS[DIMS <= ncol(Embeddings(ref_integrated, "pca"))]
K_WEIGHT   <- 50

# ---- ensure we have a usable UMAP model (with nn_index); rebuild if needed ----
ensure_umap_model <- function() {
    um_model <- ref_integrated[["umap"]]@misc$model
    # heuristic: uwot models typically have $nn_index; if missing, rebuild
    has_nn <- is.list(um_model) && !is.null(um_model$nn_index)
    if (!has_nn) {
        message("Rebuilding UMAP model from reference PCA (will also update ref UMAP coords)")
        stopifnot(requireNamespace("uwot", quietly = TRUE))
        Xref <- Embeddings(ref_integrated, "pca")[, DIMS, drop = FALSE]
        fit  <- uwot::umap(Xref, n_neighbors = 30, min_dist = 0.3,
                           metric = "euclidean", ret_model = TRUE, verbose = TRUE)
        # update ref UMAP reduction to keep everything consistent
        emb  <- as.matrix(fit$embedding)
        rownames(emb) <- rownames(Xref)
        colnames(emb) <- paste0("UMAP_", seq_len(ncol(emb)))
        ref_integrated[["umap"]] <- CreateDimReducObject(embeddings = emb, key = "UMAP_", assay = DefaultAssay(ref_integrated))
        ref_integrated[["umap"]]@misc$model <- fit
        um_model <- fit
    }
    um_model
}

# ---- UMAP projection robustly (prefers Seurat::ProjectUMAP in 'query' mode; else uwot) ----
project_umap_robust <- function(q, dims_use) {
    um_model <- ensure_umap_model()
    Xq <- Embeddings(q, "ref.pca")[, dims_use, drop = FALSE]
    
    # 1) Try Seurat's matrix-mode ProjectUMAP (some builds expect 'query=')
    out <- try(Seurat::ProjectUMAP(umap.model = um_model, query = as.matrix(Xq)), silent = TRUE)
    if (!inherits(out, "try-error")) {
        # returns a matrix of coords; attach as a reduction
        UM <- as.matrix(out)
        colnames(UM) <- if (is.null(colnames(UM))) paste0("UMAP_", seq_len(ncol(UM))) else colnames(UM)
        rownames(UM) <- colnames(q)
        q[["umap"]] <- CreateDimReducObject(embeddings = UM, key = "UMAP_", assay = DefaultAssay(q))
        return(q)
    }
    
    # 2) uwot fallback
    stopifnot(requireNamespace("uwot", quietly = TRUE))
    UM <- uwot::umap_transform(as.matrix(Xq), model = um_model)
    UM <- as.matrix(UM)
    colnames(UM) <- paste0("UMAP_", seq_len(ncol(UM)))
    rownames(UM) <- colnames(q)
    q[["umap"]] <- CreateDimReducObject(embeddings = UM, key = "UMAP_", assay = DefaultAssay(q))
    q
}

# ---- direct PCA projection (Plan B; no anchors). always coerce to base matrices ----
project_pca_direct <- function(q, dims_use) {
    L <- as.matrix(Loadings(ref_integrated, reduction = "pca"))       # features x PCs
    feats <- intersect(rownames(L), rownames(q))
    if (length(feats) < 50) stop("Too few overlapping genes with ref PCA (n=", length(feats), ").")
    X <- as.matrix(LayerData(q[["SCT"]], layer = "data")[feats, , drop = FALSE]) # genes x cells
    PC <- t(X) %*% L[feats, dims_use, drop = FALSE]                                # cells x dims
    PC <- as.matrix(PC)
    colnames(PC) <- paste0("refPC_", seq_along(dims_use))
    rownames(PC) <- colnames(q)
    q[["ref.pca"]] <- CreateDimReducObject(embeddings = PC, key = "refPC_", assay = DefaultAssay(q))
    q
}

# ---- signature-robust IntegrateEmbeddings (Plan A) ----
integrate_embeddings_try <- function(anchors_q, q, red_use, dims_use) {
    # try new API (with query/ref) then old; try 'reduction' then 'reductions'
    atts <- list(
        list(args = list(anchorset=anchors_q, reference=ref_integrated, query=q,
                         new.reduction.name="ref.pca", reduction=red_use,  dims=dims_use, k.weight=K_WEIGHT)),
        list(args = list(anchorset=anchors_q, reference=ref_integrated, query=q,
                         new.reduction.name="ref.pca", reductions=red_use, dims=dims_use, k.weight=K_WEIGHT)),
        list(args = list(anchorset=anchors_q, new.reduction.name="ref.pca",
                         reduction=red_use, dims=dims_use, k.weight=K_WEIGHT)),
        list(args = list(anchorset=anchors_q, new.reduction.name="ref.pca",
                         reductions=red_use, dims=dims_use, k.weight=K_WEIGHT))
    )
    last <- NULL
    for (a in atts) {
        res <- try(do.call(IntegrateEmbeddings, a$args), silent = TRUE)
        if (!inherits(res, "try-error")) return(res)
        last <- attr(res, "condition")$message
    }
    stop("IntegrateEmbeddings failed across signatures: ", last)
}

# ---- single-file mapper with automatic fallback ----
map_single <- function(ff) {
    message(">>> mapping single chunk: ", ff)
    q <- readRDS(ff)
    DefaultAssay(q) <- "SCT"
    
    # features present on the query; prefill scale rows to avoid residual recompute
    feats_q <- intersect(features_ref, rownames(q))
    if (length(feats_q) < 200)
        feats_q <- intersect(features_ref, rownames(LayerData(q[["SCT"]], layer = "data")))
    if (length(feats_q) < 50) stop("Too few overlapping features on query (n=", length(feats_q), ").")
    q <- ScaleData(q, features = feats_q, verbose = FALSE)
    
    # --- Plan A: anchors → integrate → project UMAP
    didA <- FALSE
    try({
        anchors_q <- FindTransferAnchors(
            reference            = ref_integrated,
            query                = q,
            normalization.method = "SCT",
            reference.reduction  = "pca",
            dims                 = DIMS,
            features             = feats_q,
            recompute.residuals  = FALSE
        )
        # get reduction name from anchorset, if present
        reds <- try(names(anchors_q@reductions), silent = TRUE)
        if (!inherits(reds, "try-error") && length(reds) > 0) {
            red_use <- if ("pcaproject" %in% reds) "pcaproject" else if ("pcaproject.l2" %in% reds) "pcaproject.l2" else reds[1]
            q <- integrate_embeddings_try(anchors_q, q, red_use, DIMS)
            q <- project_umap_robust(q, DIMS)
            didA <- TRUE
        }
    }, silent = TRUE)
    
    # --- Plan B: direct PCA projection → project UMAP
    if (!didA || !"umap" %in% Reductions(q)) {
        message("Plan A unavailable/failed -> using direct PCA + model projection (Plan B)")
        q <- project_pca_direct(q, DIMS)
        q <- project_umap_robust(q, DIMS)
    }
    
    message("Reductions: ", paste(Reductions(q), collapse = ", "))
    print(utils::head(Embeddings(q, "umap")))
    invisible(q)
}
ff <- "SCT_chunks/sctrans_HP2_hd_chunk_015_SCTonly.rds"
q_m <- map_single(ff)

# if it looks good:
saveRDS(q_m, file.path("SCT_chunks_mapped", sub("\\.rds$", "_mapped.rds", basename(ff))))
 

#### now loop this mmffff




suppressPackageStartupMessages(library(Seurat))
suppressWarnings({
  if (requireNamespace("future", quietly = TRUE)) {
    future::plan(future::sequential)  # keep it stable
    options(future.rng.onMisuse = "ignore")
  }
})
options(Seurat.object.validate = FALSE)

# ---------- ensure we have map_single() defined (the robust version you tested) ----------
if (!exists("map_single")) {
  stop("map_single() not found in this session. Re-run the single-file block we used, then re-run this loop.")
}

# ---------- helpers for repair (only used if any mapped object lacks UMAP) ----------
ensure_umap_model <- function() {
  um_model <- ref_integrated[["umap"]]@misc$model
  has_nn <- is.list(um_model) && !is.null(um_model$nn_index)
  if (!has_nn) {
    message("Rebuilding UMAP model from reference PCA (keeps ref/query consistent)")
    stopifnot(requireNamespace("uwot", quietly = TRUE))
    DIMS_USE <- DIMS[DIMS <= ncol(Embeddings(ref_integrated, "pca"))]
    Xref <- Embeddings(ref_integrated, "pca")[, DIMS_USE, drop = FALSE]
    fit  <- uwot::umap(Xref, n_neighbors = 30, min_dist = 0.3,
                       metric = "euclidean", ret_model = TRUE, verbose = TRUE)
    emb  <- as.matrix(fit$embedding)
    rownames(emb) <- rownames(Xref)
    colnames(emb) <- paste0("UMAP_", seq_len(ncol(emb)))
    ref_integrated[["umap"]] <- CreateDimReducObject(embeddings = emb, key = "UMAP_", assay = DefaultAssay(ref_integrated))
    ref_integrated[["umap"]]@misc$model <- fit
    um_model <- fit
  }
  um_model
}

repair_umap_if_missing <- function(mf) {
  o <- readRDS(mf)
  if ("umap" %in% Reductions(o)) return(invisible(TRUE))
  if (!"ref.pca" %in% Reductions(o)) return(invisible(FALSE))
  um_model <- ensure_umap_model()
  # prefer matrix-mode ProjectUMAP (query=)
  ok <- TRUE; err <- NULL
  o2 <- try(Seurat::ProjectUMAP(umap.model = um_model,
                                query = as.matrix(Embeddings(o, "ref.pca")[, DIMS, drop = FALSE])),
            silent = TRUE)
  if (inherits(o2, "try-error")) {
    ok <- FALSE; err <- attr(o2, "condition")$message
  }
  if (!ok) {
    # uwot fallback
    stopifnot(requireNamespace("uwot", quietly = TRUE))
    UM <- uwot::umap_transform(as.matrix(Embeddings(o, "ref.pca")[, DIMS, drop = FALSE]),
                               model = um_model)
    UM <- as.matrix(UM)
    colnames(UM) <- paste0("UMAP_", seq_len(ncol(UM)))
    rownames(UM) <- colnames(o)
    o[["umap"]] <- CreateDimReducObject(embeddings = UM, key = "UMAP_", assay = DefaultAssay(o))
  } else {
    UM <- as.matrix(o2)
    colnames(UM) <- if (is.null(colnames(UM))) paste0("UMAP_", seq_len(ncol(UM))) else colnames(UM)
    rownames(UM) <- colnames(o)
    o[["umap"]] <- CreateDimReducObject(embeddings = UM, key = "UMAP_", assay = DefaultAssay(o))
  }
  saveRDS(o, mf)
  message("Repaired UMAP -> ", basename(mf))
  invisible(TRUE)
}

# ---------- build targets (skip already-mapped) ----------
dir.create(OUTDIR_MAPPED, showWarnings = FALSE, recursive = TRUE)
mapped_now <- list.files(OUTDIR_MAPPED, pattern = "_mapped\\.rds$", full.names = FALSE)
mapped_keys <- sub("_mapped\\.rds$", ".rds", mapped_now)           # back to chunk base names

targets <- man$file[!(basename(man$file) %in% mapped_keys)]
cat("To map: ", length(targets), " of ", nrow(man), " total chunks\n", sep = "")

# ---------- failure log ----------
fail_log <- file.path(OUTDIR_MAPPED, "mapping_failures.log")
if (file.exists(fail_log)) file.remove(fail_log)

# ---------- loop ----------
n_ok <- 0L; n_fail <- 0L
for (ff in targets) {
  out <- file.path(OUTDIR_MAPPED, sub("\\.rds$", "_mapped.rds", basename(ff)))
  ok <- TRUE; err <- NULL
  q_m <- NULL

  cat(">> ", basename(ff), "\n", sep = "")
  tryCatch({
    q_m <- map_single(ff)                               # uses the robust Plan A/Plan B you tested
  }, error = function(e) { ok <<- FALSE; err <<- conditionMessage(e) })

  if (!ok) {
    n_fail <- n_fail + 1L
    cat(sprintf("%s\tFAILED\t%s\n", basename(ff), err), file = fail_log, append = TRUE)
    message("   failed -> ", err)
    next
  }

  # ensure UMAP present (it should be), else repair
  if (!("umap" %in% Reductions(q_m)) && "ref.pca" %in% Reductions(q_m)) {
    invisible(repair_umap_if_missing(out))  # uses ref model or rebuilds; saves inside
  } else {
    saveRDS(q_m, out)
  }
  n_ok <- n_ok + 1L
  rm(q_m); gc()
}

# ---------- optional: repair pass for already-mapped-but-missing-umap ----------
mapped_files <- list.files(OUTDIR_MAPPED, pattern = "_mapped\\.rds$", full.names = TRUE)
if (length(mapped_files)) {
  needs_repair <- vapply(mapped_files, function(mf) {
    o <- readRDS(mf); out <- !("umap" %in% Reductions(o)) && "ref.pca" %in% Reductions(o); rm(o); out
  }, logical(1))
  if (any(needs_repair)) {
    message("Repairing UMAP on ", sum(needs_repair), " mapped files missing UMAP")
    for (mf in mapped_files[needs_repair]) repair_umap_if_missing(mf)
  }
}

# ---------- summary ----------
mapped_files <- list.files(OUTDIR_MAPPED, pattern = "_mapped\\.rds$", full.names = TRUE)
cat("\nSummary:\n",
    "  mapped:   ", length(mapped_files), " files\n",
    "  new ok:   ", n_ok, "\n",
    "  new fail: ", n_fail, "\n",
    if (file.exists(fail_log)) paste0("  failures: ", fail_log, "\n") else "",
    sep = "")






suppressPackageStartupMessages(library(Seurat))

## Assumes in env:
##   OUTDIR_MAPPED, ref_integrated (with pca; UMAP may or may not have a transformable model)
##   DIMS (e.g., 1:30). If not set, we’ll derive it.
if (!exists("DIMS")) DIMS <- 1:ncol(Embeddings(ref_integrated, "pca"))
DIMS <- DIMS[DIMS <= ncol(Embeddings(ref_integrated, "pca"))]

# --- (1) Safe saver to avoid serialization crashes ---
safe_save_rds <- function(obj, path) {
  # scrub misc/tool slots that sometimes carry unserializable stuff
  if ("umap" %in% Reductions(obj))    obj[["umap"]]@misc    <- list()
  if ("ref.pca" %in% Reductions(obj)) obj[["ref.pca"]]@misc <- list()
  obj@tools <- list()

  # coerce meta to base types
  md <- obj@meta.data
  for (i in seq_len(ncol(md))) {
    xi <- md[[i]]
    if (is.factor(xi)) md[[i]] <- as.character(xi)
    else if (is.list(xi)) md[[i]] <- vapply(xi, function(z) if (length(z)) as.character(z)[1] else NA_character_, "")
  }
  obj@meta.data <- md

  # ensure reductions are base matrices
  for (rn in intersect(c("ref.pca","umap"), Reductions(obj))) {
    emb <- Embeddings(obj, rn)
    obj[[rn]] <- CreateDimReducObject(embeddings = as.matrix(emb), key = Key(obj[[rn]]), assay = DefaultAssay(obj))
  }

  # diet to SCT + reductions only (memory-lean)
  obj_small <- DietSeurat(
    obj,
    assays    = "SCT",
    dimreducs = intersect(c("ref.pca","umap"), Reductions(obj)),
    graphs    = NULL, counts = FALSE, data = TRUE, scale.data = FALSE
  )
  saveRDS(obj_small, path)
}

# --- (2) Ensure we have a transformable UMAP model (with nn_index); rebuild from ref PCA if needed ---
ensure_umap_model <- function() {
  um_model <- tryCatch(ref_integrated[["umap"]]@misc$model, error = function(e) NULL)
  ok <- is.list(um_model) && !is.null(um_model$nn_index)
  if (!ok) {
    message("Rebuilding UMAP model from reference PCA (consistent for all projections)")
    stopifnot(requireNamespace("uwot", quietly = TRUE))
    Xref <- Embeddings(ref_integrated, "pca")[, DIMS, drop = FALSE]
    fit  <- uwot::umap(
      Xref, n_neighbors = 30, min_dist = 0.3, metric = "euclidean",
      ret_model = TRUE, verbose = TRUE
    )
    emb  <- as.matrix(fit$embedding)
    rownames(emb) <- rownames(Xref)
    colnames(emb) <- paste0("UMAP_", seq_len(ncol(emb)))
    ref_integrated[["umap"]] <- CreateDimReducObject(embeddings = emb, key = "UMAP_", assay = DefaultAssay(ref_integrated))
    ref_integrated[["umap"]]@misc$model <- fit
    um_model <- fit
  }
  um_model
}

# --- (3) If ref.pca is missing on a mapped object, build it via direct projection ---
project_pca_direct <- function(o) {
  L <- as.matrix(Loadings(ref_integrated, reduction = "pca"))   # features x PCs
  feats <- intersect(rownames(L), rownames(o))
  if (length(feats) < 50) stop("Too few overlapping genes with ref PCA (n=", length(feats), ").")
  X <- as.matrix(LayerData(o[["SCT"]], layer = "data")[feats, , drop = FALSE])  # genes x cells
  PC <- t(X) %*% L[feats, DIMS, drop = FALSE]                                    # cells x dims
  colnames(PC) <- paste0("refPC_", seq_along(DIMS))
  rownames(PC) <- colnames(o)
  o[["ref.pca"]] <- CreateDimReducObject(embeddings = as.matrix(PC), key = "refPC_", assay = DefaultAssay(o))
  o
}

# --- (4) Project UMAP for an object that has ref.pca (matrix-mode first; uwot fallback) ---
project_umap_for_obj <- function(o) {
  um_model <- ensure_umap_model()
  Xq <- Embeddings(o, "ref.pca")[, DIMS, drop = FALSE]
  # Seurat matrix-mode (some builds want 'query=')
  out <- try(Seurat::ProjectUMAP(umap.model = um_model, query = as.matrix(Xq)), silent = TRUE)
  if (!inherits(out, "try-error")) {
    UM <- as.matrix(out)
    if (is.null(colnames(UM))) colnames(UM) <- paste0("UMAP_", seq_len(ncol(UM)))
    rownames(UM) <- colnames(o)
    o[["umap"]] <- CreateDimReducObject(embeddings = UM, key = "UMAP_", assay = DefaultAssay(o))
    return(o)
  }
  # uwot fallback
  stopifnot(requireNamespace("uwot", quietly = TRUE))
  UM <- uwot::umap_transform(as.matrix(Xq), model = um_model)
  UM <- as.matrix(UM)
  colnames(UM) <- paste0("UMAP_", seq_len(ncol(UM)))
  rownames(UM) <- colnames(o)
  o[["umap"]] <- CreateDimReducObject(embeddings = UM, key = "UMAP_", assay = DefaultAssay(o))
  o
}

# --- (5) Repair pass over already-mapped files ---
repair_all_mapped <- function(dir = OUTDIR_MAPPED) {
  files <- list.files(dir, pattern = "_mapped\\.rds$", full.names = TRUE)
  if (!length(files)) { message("No mapped files found in: ", dir); return(invisible(NULL)) }
  flog  <- file.path(dir, "mapping_repairs.log"); if (file.exists(flog)) file.remove(flog)

  n_checked <- 0L; n_repaired <- 0L; n_failed <- 0L
  for (fp in files) {
    n_checked <- n_checked + 1L
    ok <- TRUE; err <- NULL
    obj <- NULL
    tryCatch({
      o <- readRDS(fp)
      # if missing ref.pca, build via direct projection
      if (!"ref.pca" %in% Reductions(o)) o <- project_pca_direct(o)
      # if missing umap, project now
      if (!"umap" %in% Reductions(o)) {
        o <- project_umap_for_obj(o)
        n_repaired <- n_repaired + 1L
      }
      safe_save_rds(o, fp)
      rm(o)
    }, error = function(e) { ok <<- FALSE; err <<- conditionMessage(e) })

    if (!ok) {
      n_failed <- n_failed + 1L
      cat(sprintf("%s\tFAILED\t%s\n", basename(fp), err), file = flog, append = TRUE)
      message("Repair failed for ", basename(fp), " -> ", err)
    } else {
      message("OK -> ", basename(fp))
    }
    gc()
  }
  message(sprintf("Repair summary: checked=%d | repaired=%d | failed=%d", n_checked, n_repaired, n_failed))
  if (n_failed > 0) message("See log: ", flog)
  invisible(list(checked = n_checked, repaired = n_repaired, failed = n_failed))
}

# ---------------- one-liner to run the repair ----------------
repair_all_mapped()
mapped <- list.files(OUTDIR_MAPPED, "_mapped\\.rds$", full.names = TRUE)
all(has_umap <- vapply(mapped, function(f){ "umap" %in% Reductions(readRDS(f)) }, logical(1)))
table(has_umap)

# --- 0) which one failed? (from the repairs log) ---
fail_log <- file.path(OUTDIR_MAPPED, "mapping_repairs.log")
bad_map <- "sctrans_placental_ref_Jul_2025_chunk_006_SCTonly_mapped.rds"  # from your printout
# alternatively parse the log:
if (file.exists(fail_log)) {
  fl <- try(read.delim(fail_log, header = FALSE, sep = "\t", stringsAsFactors = FALSE), silent = TRUE)
  if (!inherits(fl, "try-error") && nrow(fl)) bad_map <- fl[[1]][1]
}
message("Targeting bad mapped file: ", bad_map)

# --- 1) map from mapped-filename -> original chunk path via manifest ---
get_chunk_src <- function(mapped_basename, manifest_df = man) {
  base_chunk <- sub("_mapped\\.rds$", ".rds", mapped_basename)
  hit <- manifest_df$file[basename(manifest_df$file) == base_chunk]
  if (!length(hit)) stop("Cannot find source chunk for ", mapped_basename, " in manifest.")
  normalizePath(hit[1], mustWork = TRUE)
}

# --- 2) redo that single chunk end-to-end (map_single + safe_save_rds) ---
redo_one_mapped <- function(mapped_basename) {
  src <- get_chunk_src(mapped_basename, man)
  out <- file.path(OUTDIR_MAPPED, sub("\\.rds$", "_mapped.rds", basename(src)))

  # remove the corrupted mapped file if it exists
  if (file.exists(out)) {
    message("Removing corrupted mapped file -> ", out)
    try(unlink(out), silent = TRUE)
  }

  message("Remapping from source -> ", src)
  q_m <- map_single(src)         # uses the robust Plan A/Plan B you validated
  safe_save_rds(q_m, out)        # diet + base matrices; avoids serialization crash
  message("OK: re-saved -> ", out)
  invisible(out)
}

# --- 3) run it ---
redo_one_mapped(bad_map)

# --- 4) quick sanity: count mapped vs manifest & list any still missing ---
mapped_files <- list.files(OUTDIR_MAPPED, pattern = "_mapped\\.rds$", full.names = FALSE)
mapped_keys  <- sub("_mapped\\.rds$", ".rds", mapped_files)
missing <- setdiff(basename(man$file), mapped_keys)
cat("\nSummary after redo:\n  mapped: ", length(mapped_files), " / ", nrow(man),
    "\n  missing: ", length(missing), "\n", sep = "")
if (length(missing)) print(utils::head(missing, 10))





# --- make sure DIMS is defined and valid ---
if (!exists("DIMS")) DIMS <- 1:ncol(Embeddings(ref_integrated, "pca"))
DIMS <- DIMS[DIMS <= ncol(Embeddings(ref_integrated, "pca"))]

# --- override: accept optional dims_use (and swallow extras) ---
project_pca_direct <- function(q, dims_use = DIMS, ...) {
  # loadings: features x PCs
  L <- as.matrix(Loadings(ref_integrated, reduction = "pca"))
  feats <- intersect(rownames(L), rownames(q))
  if (length(feats) < 50)
    stop("Too few overlapping genes with ref PCA (n=", length(feats), ").")

  # use SCT residuals; coerce to base matrix
  X  <- as.matrix(LayerData(q[["SCT"]], layer = "data")[feats, , drop = FALSE])  # genes x cells
  PC <- t(X) %*% L[feats, dims_use, drop = FALSE]                                # cells x dims
  PC <- as.matrix(PC)

  colnames(PC) <- paste0("refPC_", seq_along(dims_use))
  rownames(PC) <- colnames(q)
  q[["ref.pca"]] <- CreateDimReducObject(embeddings = PC, key = "refPC_", assay = DefaultAssay(q))
  q
}

redo_one_mapped(bad_map)









#### downstream



mapped <- list.files(OUTDIR_MAPPED, "_mapped\\.rds$", full.names=TRUE)
length(mapped)  # should equal nrow(man)
# confirm reductions exist
any(vapply(mapped, function(f){ !"umap" %in% Reductions(readRDS(f)) }, logical(1)))




library(data.table)
files <- list.files(OUTDIR_MAPPED, "_mapped\\.rds$", full.names=TRUE)
sink <- gzfile(file.path(OUTDIR_MAPPED,"all_embeddings.tsv.gz"), "wt")
wrote <- FALSE
for (fp in files) {
  o <- readRDS(fp)
  um <- Embeddings(o, "umap")
  rp <- Embeddings(o, "ref.pca")[, 1:30, drop=FALSE]
  dt <- as.data.table(um)[, cell := rownames(um)]
  dt <- cbind(dt, as.data.table(rp))
  dt[, sample := sub("^sctrans_(.+?)_chunk_\\d+_SCTonly.*$", "\\1", sub("\\.rds$","", basename(fp)))]
  keep <- intersect(c("predicted.celltype","predicted.celltype.score"), colnames(o@meta.data))
  if (length(keep)) dt <- cbind(dt, o@meta.data[rownames(um), keep, drop=FALSE])
  data.table::fwrite(dt, sink, sep="\t", col.names=!wrote)
  wrote <- TRUE
}
close(sink)
library(ClusterR)
X <- fread(file.path(OUTDIR_MAPPED,"all_embeddings.tsv.gz"), sep="\t")
pc_cols <- grep("^refPC_|^refPC_?\\d+$", names(X), value=TRUE)  # adjust to your column names
set.seed(1)
km <- Cluster_MiniBatchKmeans(as.matrix(X[, ..pc_cols]), clusters=50, batch_size=50000, num_init=5, max_iters=200)
X[, cluster_km := km$clusters]
fwrite(X, file.path(OUTDIR_MAPPED,"all_embeddings_with_km.tsv.gz"), sep="\t")
pred <- TransferData(anchorset = anchors_q, refdata = Idents(ref_integrated), dims = DIMS, k.weight = 50)
q <- AddMetaData(q, pred)
library(edgeR); library(Matrix)
meta <- fread(file.path(OUTDIR_MAPPED,"all_embeddings_with_km.tsv.gz"))
for (sid in unique(meta$sample)) {
  obj_path <- file.path("with_annotations", paste0("seurat_", sid, "_withPredCols.rds"))
  if (!file.exists(obj_path)) next
  obj <- readRDS(obj_path)
  src <- intersect(c("Spatial.055um","Spatial.008um","Xenium.5k","Xenium.500","RNA"), names(obj@assays))[1]
  DefaultAssay(obj) <- src
  # map subset
  m <- meta[sample==sid]
  m[, orig := sub("^[^_]+_c\\d+_", "", cell)]
  keep <- intersect(m$orig, colnames(obj))
  if (!length(keep)) next
  groups <- factor(m$predicted.celltype[match(keep, m$orig)])  # or cluster_km
  Y <- GetAssayData(obj, slot="counts", assay=src)[, keep, drop=FALSE]
  M <- as.matrix(Y)
  G <- model.matrix(~0 + groups); colnames(G) <- levels(groups)
  pb <- M %*% G  # genes x groups
  # store per sample
  saveRDS(pb, file.path("pseudobulk", paste0(sid,"_pb.rds")))
}
# merge sample-level pseudobulks across samples (align genes), make DGEList, TMM norm, design ~ condition + donor






