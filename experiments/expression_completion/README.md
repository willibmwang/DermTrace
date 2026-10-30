# Experimental gene-expression completion

This package prepares processed expression profiles for a research question: can a model use measured genes from a new sample to predict its deliberately hidden genes across donors and cell types? It is an experimental data interface, not a trained model or evidence of biological transfer. See [DATA_REVIEW.md](DATA_REVIEW.md) for the source audit and its limits.

A separate, locally constructed DermTrace classification release uses within-sample ranks calculated before genes are hidden. That release and its classification loader are not part of current public `main`. Those saved ranks cannot be used directly as quantitative inputs for this completion task: a hidden gene may have influenced the ranks of visible genes. This package uses a separate completion interface and follows the local classification loader's validation principles; it does not call that loader.

## What is included

- Adapters for specific local processed-source packages: GSE163973 skin cell counts and pooled counts, GSE144424 cardiac counts, and GSE95280 melanoma gene coverage. The interfaces check sample groups and retain assay units. They are not general readers for every upstream dataset.
- Exact, unambiguous Ensembl gene alignment, with a separate mask for genes missing from an assay. A measured zero is not treated as a missing gene.
- Seeded context, feedback and evaluation gene masks. Fixed anchor genes remain visible, and only those visible anchors set the sample's normalization factor. Anchor selection and any fitted references must use training groups only.
- Separate model-facing observation, feedback and final evaluation interfaces. These are Python API boundaries, not a security sandbox. A real evaluation must keep held-out values and scores outside the adaptation loop and prevent them from influencing model or preprocessing choices.
- Checks that distinguish an exact source-file checksum match from a specifically recorded annotation-only Parquet footer update. Neither historical source manifests nor source expression tables are rewritten.
- A generator for exploratory sample-group split files and a read-only quality-review script. Their source metadata and processed counts are supplied separately.

Cell, donor, lineage and study assignments are a separate part of the experiment. Split files alone cannot prove sample independence, enforce a cell-type holdout, or make different assays comparable. Skin donor relationships remain unresolved; see the [data review](DATA_REVIEW.md).

## Tests from a clean checkout

With Python 3.11 or later and `numpy`, `pandas`, `scipy` and `pyarrow` installed, run the synthetic interface and integrity tests from the repository root:

```sh
python -m unittest experiments.expression_completion.test_completion_loader experiments.expression_completion.test_source_integrity experiments.expression_completion.test_splits experiments.expression_completion.test_review_qc -v
```

These tests exercise masking, visible-anchor normalization, missing versus measured-zero values, partition checks, cell eligibility and source-integrity behavior with synthetic fixtures. They do not require the external expression archives or establish model performance.

## Optional processed data

The public checkout does **not** include every recovered source matrix, source manifest, cell crosswalk or annotation receipt needed by these adapters. A real-data run requires the corresponding processed-source packages and their provenance files to be supplied locally. The GSE163973 archive has `RAW` in its GEO filename but contains processed count matrices; no sequencing-read processing is needed. Keep original input files unchanged.

`python -m experiments.expression_completion.splits --help` lists the required skin, cardiac and melanoma metadata inputs and output directory. `python -m experiments.expression_completion.review_qc --help` lists the required processed archive, cell crosswalk, cell metadata and recovery receipt. The QC script also needs `threadpoolctl`. Its optional barcode-level output is local QC metadata; do not commit it. A generated split manifest records hashes under portable source names; the loader requires explicit paths to those metadata inputs when reading a split.

The previous local integration check used 60 cardiac profiles, 47 pooled skin profiles, two single cells and 27 melanoma profiles. It exercised code and masking; its 100-anchor panel was selected from just two cells and is not a research-ready panel. It was performed on a separate local data recovery, not on a clean checkout, and is not a completion-model result.
