# Processed expression: scope and observations

These are observations from a local, processed-data audit. This pull request does not include the matrices, cell barcodes, patient-level data, or the locally generated audit tables. The public repository does not currently include all of the source matrices needed to repeat the audit. Availability and file locations must be checked for each environment.

## Skin single-cell data

The [GSE163973 study](https://doi.org/10.1038/s41467-021-24110-y) describes three keloid and three normal-scar samples. The authors report 40,655 cells after filtering and state that predicted doublets were removed. The local audit found six processed count matrices with 53,148 barcodes, of which 40,655 matched supplied cell annotations. The remaining 12,493 barcodes had no supplied cell annotation. These unannotated barcodes are not automatically usable cells.

The annotated set includes 13,437 cells labelled as fibroblasts. It also contains other author-labelled cell types. These labels have not been independently validated here. The audit reconstructed 47 locally retained sample-by-cell-type pooled count profiles exactly from the single-cell counts. These pooled profiles are derived from the cells; they are not additional independent samples. The comparison group is normal scar, not unaffected skin.

The paper states that cells with more than 6,000 UMI counts were removed, but 29,084 of the 40,655 annotated cells exceed that count in the deposited matrices. All annotated cells have at most 6,000 detected features and at most 15% mitochondrial counts. The reason for the apparent discrepancy is unknown. Do not automatically apply a 6,000-UMI cutoff to these matrices. The study reports removing doublets, but the retained metadata inspected locally does not provide individual doublet scores.

The six sample IDs do not prove which keloid and scar samples came from the same person. Grouping `KF1` with `NF1`, and similarly for suffixes 2 and 3, is an exploratory safeguard; it does **not** prevent leakage if a patient's samples have different suffixes. The paper reports three patients per condition, but a sample-ID-to-patient crosswalk has not been verified. For strict cross-condition donor isolation, verify those relationships or keep all six samples in one group. Within-condition sample folds also rely on an unverified one-sample-per-patient assumption.

## Cardiac and melanoma data

- GSE144424 has 84 cardiac profiles from 15 recorded donors; the local audit identified 60 primary-eligible profiles. All conditions and repeated differentiations from a donor must stay together. The documented hypoxia exposure is six hours, so this cohort alone does not establish chronic cardiac hypoxia.
- GSE95280 recount3 contains 27 locally identified clean oxygen-comparison profiles across two lineage groups: 501MEL and the related IGR37/IGR39 lines. Its gene-level values are coverage, not counts. Two lineage groups do not support a three-way melanoma-only split.
- In a local integration check, the skin pooled, cardiac and melanoma feature axes shared 17,527 strictly mapped Ensembl IDs. This is an identity overlap, not a validated model panel or evidence that their measurement scales are equivalent.

## Interpretation limits

Cells from one tissue sample are not independent donors. Sample split assignments alone do not enforce a cell-type holdout. Cross-domain comparisons can change cell type, study, and assay together. Counts, coverage, TPM and microarray intensities must retain their measurement definitions; identical gene names and sample scaling do not make them interchangeable.

The local audit used processed matrices and metadata only. Its exact reconstruction and numerical observations are historical checks of those local artifacts, not a claim that this repository contains all inputs needed to reproduce them. No completion model was trained or evaluated in that audit.
