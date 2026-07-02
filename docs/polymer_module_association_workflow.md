# Exploratory polymer-specific module-association workflow

This workflow reproduces the Supplementary Fig. 4c-d style exploratory polymer-specific analyses from the public processed Seurat object.

The script is:

```r
scripts/07_polymer_associations.R
```

It expects `scripts/05_module_scores.R` to have created:

```text
outputs/tables/module_scores/subject_level_module_scores_platform_z.csv
```

If that table is missing, the script runs the module-score workflow first.

## Inputs

The script starts from the public Seurat object loaded as `obj` and uses public polymer burden metadata columns such as:

```text
pet_ug_per_g
pu_ug_per_g
pvc_ug_per_g
sbr_ug_per_g
abs_ug_per_g
```

It computes `log2(ug/g + 1)` polymer terms internally so that the public workflow does not depend on private intermediate columns.

## Outputs

Tables are written to:

```text
outputs/tables/polymer_associations/
```

Key outputs include:

```text
polymer_subject_level_missingness.csv
polymer_module_model_input_subject_level.csv
polymer_module_within_delivery_slopes.csv
polymer_module_leave_one_subject_out_iterations.csv
polymer_module_leave_one_subject_out_summary.csv
```

Figures are written to:

```text
outputs/figures/polymer_associations/
```

Key outputs include:

```text
SuppFig4c_polymer_module_within_delivery_slopes_heatmap.png/pdf
SuppFig4d_polymer_module_focused_scatter_QC.png/pdf
```

## Interpretation

These models are exploratory and were used to prioritize candidate polymer-module relationships for future validation. They should not be interpreted as definitive polymer-specific mechanisms because of nested-cohort size, polymer sparsity, platform heterogeneity, and polymer collinearity.
