# DermTrace

DermTrace extends [KeloidBench](keloidbench/), the original keloid-only pipeline, to a broader collection of skin and fibrosis expression data. This repository holds the shared processed data, the pipelines that build new datasets, and the archived KeloidBench project.

## Repository layout

| Path | Contents |
|---|---|
| [`data/processed/`](data/processed/) | Every processed dataset, listed below |
| [`data/raw/`](data/raw/) | Cohort registries and download manifests. Raw downloads are not committed. |
| [`data/schemas/`](data/schemas/) | Feature and target schemas |
| [`scripts/tsaolab_rnaseq/`](scripts/tsaolab_rnaseq/) | Pipeline for the Tsao lab RNA-seq |
| [`keloidbench/`](keloidbench/) | KeloidBench as published: NeurIPS GenAI4Health poster and paper, code, results, demo and website source. The site is archived at <https://willibmwang.github.io/DermTrace/keloidbench/>. |

## Processed data

All processed data is in [`data/processed/`](data/processed/):

| Folder | Contents | Built by |
|---|---|---|
| [`tsaolab_o2/`](data/processed/tsaolab_o2/) | Tsao lab oxygen-tension bulk RNA-seq, 60 samples ([details below](#tsao-lab-oxygen-tension-rna-seq)) | `scripts/tsaolab_rnaseq/` |
| [`public_keloid/`](data/processed/public_keloid/) | 15 public GEO keloid cohorts, 169 samples: 12 development and 3 lockbox cohorts | `keloidbench/scripts/preprocess_public_keloid.py` |
| [`bulk_rnaseq/`](data/processed/bulk_rnaseq/) | Bulk RNA-seq from GSE158395, GSE188952 and Sun et al. (*Burns*), 50 samples | `keloidbench/scripts/preprocess_bulk_rnaseq.py` |
| [`microarray/`](data/processed/microarray/) | Six GEO microarray studies (GSE7890, GSE92566, GSE90051, GSE44270, GSE3189, GSE145725), 154 samples, plus keloid signatures | `keloidbench/scripts/preprocess_microarray.py`, `build_keloid_signature.py` |
| [`scrna_spatial/`](data/processed/scrna_spatial/) | 230 profiles: 228 from three single-cell and spatial studies (GSE163973, GSE181297, GSE160536), plus 2 bulk RNA-seq profiles of sorted keloid fibroblasts (GSE175866) from the GSE163973 paper | `keloidbench/scripts/preprocess_scrna_spatial.py` |
| [`spheroid/`](data/processed/spheroid/) | Keloid spheroid tables (Choi et al.) and spheroid morphology sources | `keloidbench/scripts/data_extraction.py`, `build_morphology_table.py` |
| [`training/`](data/processed/training/) | The combined KeloidBench training corpus: 963 profiles, with splits and feature views | `keloidbench/scripts/build_training_corpus.py` |
| [`archive/`](data/processed/archive/) | Earlier versions removed from git on 2026-06-11 and 2026-07-09, restored from history | See its [README](data/processed/archive/README.md) |

The public-cohort, bulk, microarray and single-cell folders hold labels, metadata and summaries. Their large matrices (`*_expression_wide.parquet`, `*_X.npy`, `*_encoder_decoder_samples.jsonl`) were never committed. The current versions are being added from the Misha cluster.

### Tsao lab oxygen-tension RNA-seq

The processed Tsao lab bulk RNA-seq is in [`data/processed/tsaolab_o2/`](data/processed/tsaolab_o2/). The dataset has 60 samples (48 human, 12 mouse) covering the ATCC keloid fibroblast line, BJ, AC16, HL-1 and HEK293-like cells. Each submitted line was grown in four oxygen arms (20% O₂, 1% O₂, or switched between the two), in triplicate.

| File | Contents |
|---|---|
| `tsaolab_o2_<species>_tpm_symbol.tsv.gz` | Gene symbol × sample TPM, not logged. Start here. |
| `tsaolab_o2_<species>_counts_symbol.tsv.gz` | Gene symbol × sample salmon estimated counts |
| `tsaolab_o2_<species>_sample_metadata.csv` | One row per sample: cell line, oxygen arm, QC flags |
| `per_order/` | log2(TPM+1) tables for each sequencing order |
| `salmon_quant/` | Transcript-level salmon output, readable by tximport |
| `qc/` | Read-level and cell-identity QC |

```python
import pandas as pd

tsao = "data/processed/tsaolab_o2"
tpm = pd.read_csv(f"{tsao}/tsaolab_o2_human_tpm_symbol.tsv.gz", sep="\t", index_col="gene_symbol")
meta = pd.read_csv(f"{tsao}/tsaolab_o2_human_sample_metadata.csv", index_col="sample_id")
```

Read the [dataset README](data/processed/tsaolab_o2/README.md) before analysis. Three things differ from the submission forms:

- In order 0013831748, samples 13–24 (all 12 submitted as NIH3T3 mouse fibroblasts) are HEK293-like human cells.
- In the same order, 10 of the 12 BJ samples carry HEK293-like admixture.
- In order 1030390, samples 1–6 were submitted as the ATCC keloid line but are BJ by SNP genotype. Only samples 7–12 are keloid (ATCC CRL-1762, one donor), and they were only grown in the arms that start at 1% O₂.

The raw FASTQs (about 258 GB) are too large for GitHub. They are on the lab Google Drive in the folders `0013831748-2024 Raw_Data` and `1030390-2025 Raw_Data`; ask William for access. [`scripts/tsaolab_rnaseq/`](scripts/tsaolab_rnaseq/) rebuilds every file above from those FASTQs.

## KeloidBench

[`keloidbench/`](keloidbench/) is the original keloid-only project, kept as it was published. Run its code from that folder (`cd keloidbench`); it reads processed data from the shared `data/` folder here. See [keloidbench/README.md](keloidbench/README.md) for its results and how to reproduce them.
