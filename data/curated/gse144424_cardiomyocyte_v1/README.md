# GSE144424 cardiomyocyte expression package v1

Standalone exploratory data addition to DermTrace. Existing frozen datasets, models and benchmarks are unchanged. This package does not establish predictive improvement or a clinical claim.

## Material and units

The deposited matrix contains 19,081 Ensembl genes × 84 profiles from 15 human donors and 21 differentiation experiments. Donors 18511, 18852 and 18855 each have three experiments. Each experiment has four condition subcultures. The matrix contains unnormalized, WASP-filtered featureCounts counts aligned to hg19, with native source annotations preserved. It is not TPM, TMM-normalized or RUV-corrected expression.

| Condition | Exposure | Harvest after experiment start |
|---|---|---:|
| A | 10% oxygen reference | 6 h |
| B | 1% oxygen for 6 h | 6 h |
| C | 1% for 6 h, then 10% for 6 h | 12 h |
| D | 1% for 6 h, then 10% for 24 h | 30 h |

The source describes iPSC-derived cardiomyocyte cultures (reported purity 40–97%, median 81%). Condition subcultures came from the same differentiation. Replicate differentiations are nested within donor. Later time-matched normoxia controls are not described, limiting a causal interpretation of C/D versus A. This experiment has a 10% reference, unlike the Tsao 20% reference; do not pool absolute expression values across studies.

## QC and provenance decisions

All 84 deposited profiles remain present. `H18852_A` is flagged as the paper's expression outlier (RIN=6); `H18852_D` is flagged for undetermined RIN. An unflagged row is not a new claim that it passed an independent QC assessment.

`published_primary_eligible` selects experiment 2 for each of the three repeated donors and the sole experiment for each other donor: 60 profiles, 15 donors. The paper used this selection for its DE/eQTL analyses. Its differential model also used TMM/voom/limma and four RUV factors; those analyses are not reproduced by this ingestion script.

GEO metadata says the counts contain only autosomal genes, but the actual file contains 18,226 autosomal genes and 855 nonautosomal genes (chrX, chrY, chrM). All rows are retained and their chromosomes are explicit in `genes.csv`.

Eight marker Ensembl IDs are verified from the existing Tsao human annotation, using exact gene-symbol matches and requiring one unique ID per marker and presence in the deposited counts. The local mapping is used only for gene identity, not expression values or outcomes; its gene-type labels are not used for filtering.

## Files

- `counts.parquet`: genes × 84 sample IDs; index `Geneid`; native integer counts.
- `genes.csv`: the six original annotation columns (`Geneid`, `Chr`, `Start`, `End`, `Strand`, `Length`).
- `samples.csv`: GEO accession, native matrix ID, donor, differentiation, condition, timing, published eligibility and explicit QC flags. Count-library totals refer only to the deposited gene matrix.
- `marker_mapping.csv`: eight exact symbol–Ensembl mappings and their annotation source hash.
- `provenance.json`: source URLs/hashes, decisions and audit results.
- `manifest.json`: input/output hashes and validation summary.

## Source hashes

- `inputs/GSE144424_Counts_RNA_MCW_NEB.txt.gz` — SHA-256 `cad9ac4c6514550ea9bfb2b491cc2934f6952894d7cbd17338d5054d03da6f7c` — https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144424/suppl/GSE144424_Counts_RNA_MCW_NEB.txt.gz
- `inputs/GSE144424_family.soft.gz` — SHA-256 `156b46d9e11fc4a59219530808c617887ac622e711d137964e1ccb88b0766907` — https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144424/soft/GSE144424_family.soft.gz
- `inputs/GSE144424_series_matrix.txt.gz` — SHA-256 `40dda761c345cb757eaa4d2c0eefcebad948c134137ffdcdc0af855384c90ec4` — https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144424/matrix/GSE144424_series_matrix.txt.gz
- Local gene mapping: `data/processed/tsaolab_o2/tsaolab_o2_human_genes.tsv.gz` (path basis: `repository_root`) — SHA-256 `b463ad00c0f6d6b3de0ec78fc8de3c3df44cb7451a56dd6e8d2a1b3f2dd43ccb`.
- Design, pairing and published exclusions: https://elifesciences.org/articles/57345.

Mapping paths with basis `repository_root` are relative to the checkout root. Custom mappings outside the checkout are recorded by filename with basis `external_input_filename`; supply that file through `curate(..., mapping_path=...)` and verify its SHA-256 to reproduce the package.

From the repository root, use Python 3.10 with `keloidbench/requirements-revision.txt`, then run `python scripts/ingest_gse144424.py`. The script is local-only, performs no downloads, and writes only this package. Raw reads were not downloaded.
