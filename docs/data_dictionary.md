# Public object data dictionary

This file should define each public metadata column retained in `spatial_mnp_publication_processed_seurat_v1.rds`.

Minimum columns to document:

| Column | Level | Type | Description | Public/derived note |
|---|---|---|---|---|
| donor_id | subject | categorical | Coded subject identifier | Deidentified/coded |
| pregnancy_id | pregnancy | categorical | Coded pregnancy identifier | Deidentified/coded |
| sample_id_standardized | sample/run | categorical | Harmonized spatial sample identifier | Deidentified/coded |
| Delivery_class | subject | categorical | Term or Preterm delivery class | Derived from gestational age |
| MNPBurden | subject/sample | categorical | Low or High cumulative MNP burden | Threshold = 146.3 µg/g |
| cumulative_mnp_ug_per_g | sample | numeric | Cumulative MNP burden in µg/g tissue | Py-GC/MS-derived |
| celltype_display | observation | categorical | Final display cell-type annotation | Derived |
| platform_clean | run | categorical | Harmonized spatial platform label | Derived |
