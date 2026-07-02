# Migration notes for `meta_update_2026-06-12.R`

The original working script was notebook-like and contained metadata updating, module scoring, composition summaries, robustness, polymer-specific analyses, Xenium neighborhood analyses, VisiumHD co-enrichment, annotation QC, and CellChat plotting in one file.

For the public repository, the code has been split into runnable downstream scripts that start from the public processed Seurat object. The original file is archived in `provenance/archived_scripts/` for transparency.

## Sections migrated in this repo revision

- Metadata/public object build logic: `scripts/01_build_public_object.R`
- Fig. 2 composition summaries: `scripts/04_celltype_composition.R`
- Fig. 3 module-score summaries, within-delivery contrasts, LOO, and permutation: `scripts/05_module_scores.R`

## Sections still needing final manuscript-specific migration

- Polymer-specific module heatmap and LOO sensitivity
- Xenium kNN neighborhood analysis
- VisiumHD local module co-enrichment
- Supplementary annotation QC / maternal-origin / NK marker support
- CellChat alluvial plots
- NicheNet network plots
- Macrophage/Hofbauer state UMAP and state-fraction contrasts
- Spatial marker maps for HOXB7 and PECAM1

## Public reproducibility decision

The default public workflow should regenerate figure/table outputs from the deposited public Seurat object. Upstream high-memory preprocessing, RCTD batching, SCT conversion, and cross-platform integration are documented for provenance but are not required for downstream figure reproduction.


## 2026-07-02 update: macrophage/Hofbauer state analysis

Archived the original notebook-style script as `provenance/archived_scripts/save_audra_macrophage_obj.R` and migrated the runnable downstream workflow into `scripts/06_macrophage_states.R`.

The migrated script starts from `obj_public`, subsets predicted `Macrophage` and `Hofbauer Cell` observations in the MNP model subset, rebuilds a macrophage-only PCA/UMAP, calculates five macrophage-state programs, assigns each observation to the highest z-scored program, exports subject-level state fractions and statistics, runs `FindAllMarkers`, and writes Fig. 5-style panels.

The older M0/M1/M2-style `polarization` assignment is retained only as `polarization_qc` in a QC table. The manuscript-facing labels use state/program terminology, including `Activated remodeling`, to avoid overinterpreting canonical polarization.
