# DermTrace

[![Paper](https://img.shields.io/badge/paper-NeurIPS%202026%20GenAI4Health-bf684f)](https://neurips.cc/virtual/2026/loc/sydney/177621)
[![Site](https://img.shields.io/badge/site-willibmwang.github.io%2FDermTrace-3d806e)](https://willibmwang.github.io/DermTrace/)
[![License: MIT](https://img.shields.io/badge/license-MIT-d7a34b)](LICENSE)

DermTrace is a shared home for processed skin and fibrosis transcriptomic data: the datasets, the pipelines that build them, and the analysis code that uses them. It grew out of [KeloidBench](keloidbench/), the keloid-only predecessor, which is kept here as published.

> **KeloidBench: Interpretable machine learning on fibroblast program features for keloid detection.**
> William Wang, Joshua Tam, Hensin Tsao.
> [NeurIPS 2026 Workshop on Generative AI for Health](https://neurips.cc/virtual/2026/loc/sydney/177621), Sydney.

**New here? Start with the [Quickstart](QUICKSTART.md)** — it takes you from `git clone` to a loaded expression matrix in about five minutes.

| | |
|---|---|
| 🚀 **[Quickstart](QUICKSTART.md)** | Load data, reproduce KeloidBench, rebuild the RNA-seq pipeline |
| 🧬 **[Data catalog](#data)** | What is in the repository, and what is not |
| 🔬 **[Tsao lab O₂ RNA-seq](data/processed/tsaolab_o2/)** | The in-house dataset: 60 samples, four oxygen arms |
| 📄 **[NeurIPS paper](https://neurips.cc/virtual/2026/loc/sydney/177621)** | KeloidBench at GenAI4Health 2026 — also as [PDF](keloidbench/KeloidBench_NeurIPS_GenAI4Health_2026_.pdf) and [poster](keloidbench/KeloidBench%20Poster_WW.pdf) |

## Repository layout

| Path | Contents |
|---|---|
| [`data/processed/`](data/processed/) | Every processed dataset — see [Data](#data) |
| [`data/raw/`](data/raw/) | Cohort registries and download manifests. Raw downloads are not committed. |
| [`data/schemas/`](data/schemas/) | Feature and target schemas |
| [`scripts/tsaolab_rnaseq/`](scripts/tsaolab_rnaseq/) | Pipeline that builds the Tsao lab RNA-seq dataset from FASTQs |
| [`keloidbench/`](keloidbench/) | KeloidBench as published: poster, paper, code, results, demo and website source |
| [`website/`](website/) | Source of the project site at <https://willibmwang.github.io/DermTrace/> |

## Data

### What you get when you clone

Not every dataset ships its expression matrix. Large matrices were either removed from tracking or never committed, so check this table before planning an analysis.

| Folder | Expression matrices | What is committed |
|---|---|---|
| [`tsaolab_o2/`](data/processed/tsaolab_o2/) | ✅ **Yes** | Full TPM and count matrices, sample metadata, QC, salmon output |
| [`training/`](data/processed/training/) | ✅ **Yes** (feature views) | 18 feature-view parquets, profile manifest, splits, coverage audits |
| [`public_keloid/`](data/processed/public_keloid/) | ❌ No | Encoder labels, sample metadata, gene vocabulary, cohort roles, summary |
| [`bulk_rnaseq/`](data/processed/bulk_rnaseq/) | ❌ No | Encoder labels and a summary JSON |
| [`microarray/`](data/processed/microarray/) | ❌ No | Encoder labels, summary JSON, and one signature summary |
| [`scrna_spatial/`](data/processed/scrna_spatial/) | ❌ No | Encoder labels, summary JSON, metadata inventories |
| [`spheroid/`](data/processed/spheroid/) | ❌ No | Three extraction and morphology summaries |
| [`archive/`](data/processed/archive/) | ✅ Historical | Matrices and tables removed from tracking, restored from git history |

Two different things are going on:

- The **public-cohort** matrices (`public_keloid_*_expression_wide.parquet`, `*_X.npy`, `*_encoder_decoder_samples.jsonl`) were **never committed** and are gitignored.
- The **bulk, microarray, single-cell and spheroid** matrices and tables *were* committed until 2026-07-09, when commit `e625f2f` stopped tracking them. Those versions are preserved under [`data/processed/archive/2026-07-09/`](data/processed/archive/2026-07-09/) — including the keloid signature tables and all thirteen spheroid tables.

Anything not committed is rebuilt by the script in the **Built by** column below.

### Dataset catalog

| Folder | Contents | Built by |
|---|---|---|
| [`tsaolab_o2/`](data/processed/tsaolab_o2/) | Tsao lab oxygen-tension bulk RNA-seq, 60 samples ([details below](#tsao-lab-oxygen-tension-rna-seq)) | [`scripts/tsaolab_rnaseq/`](scripts/tsaolab_rnaseq/) |
| [`public_keloid/`](data/processed/public_keloid/) | 169 samples across 16 public GEO keloid accessions: 123 development, 46 lockbox | `keloidbench/scripts/preprocess_public_keloid.py` |
| [`bulk_rnaseq/`](data/processed/bulk_rnaseq/) | Bulk RNA-seq, 50 samples: GSE158395 (13), GSE188952 (12), Sun et al. *Burns* (25) | `keloidbench/scripts/preprocess_bulk_rnaseq.py` |
| [`microarray/`](data/processed/microarray/) | Six GEO microarray studies, 154 samples: GSE7890 (19), GSE92566 (7), GSE90051 (7), GSE44270 (32), GSE3189 (70), GSE145725 (19) | `keloidbench/scripts/preprocess_microarray.py`, `build_keloid_signature.py` |
| [`scrna_spatial/`](data/processed/scrna_spatial/) | 230 profiles: 228 pseudobulk from GSE163973 (220), GSE181297 (7) and GSE160536 (1), plus 2 bulk profiles of sorted keloid fibroblasts (GSE175866) | `keloidbench/scripts/preprocess_scrna_spatial.py` |
| [`spheroid/`](data/processed/spheroid/) | Keloid spheroid extractions (Choi et al.) and morphology sources | `keloidbench/scripts/data_extraction.py`, `build_morphology_table.py` |
| [`training/`](data/processed/training/) | Combined training corpus, 963 profiles, with splits and feature views ([composition below](#training-corpus)) | `keloidbench/scripts/build_training_corpus.py` |
| [`archive/`](data/processed/archive/) | Earlier versions removed from git on 2026-06-11 and 2026-07-09 | See its [README](data/processed/archive/README.md) |

Two things to know before using these cohorts:

- **GSE3189 is a melanoma study** (45 melanoma, 18 nevus, 7 normal) contributing 70 of the 154 microarray samples, so [`microarray/`](data/processed/microarray/) is not keloid-only.
- In `public_keloid/`, trust `public_keloid_sample_metadata.parquet` over the JSON sidecars. The parquet carries all 169 samples across 16 accessions, but `public_keloid_cohort_roles.json` lists only 15 (it omits GSE185309, 17 lockbox samples), and the `skipped` list in `public_keloid_summary.json` marks four accessions as failed that did contribute 44 samples between them.

#### Training corpus

The 963 profiles in [`training/`](data/processed/training/) come from six sources, only four of which have their own folder above:

| `source_dataset` | Profiles | Folder |
|---|---|---|
| `external_fibrosis` | 293 | assembled from sources outside `data/processed/` |
| `scrna_spatial` | 230 | [`scrna_spatial/`](data/processed/scrna_spatial/) |
| `public_keloid` | 169 | [`public_keloid/`](data/processed/public_keloid/) |
| `microarray` | 154 | [`microarray/`](data/processed/microarray/) |
| `arrayexpress` | 67 | assembled from sources outside `data/processed/` |
| `bulk_rnaseq` | 50 | [`bulk_rnaseq/`](data/processed/bulk_rnaseq/) |

> **Careful:** the `modalities` key in `training_corpus_summary.json` is grouped by **source**, not by assay. The manifest's own `modality` column is a genuinely different grouping — `external_fibrosis` 293, `scrna_pseudobulk` 246, `microarray` 232, `bulk_rnaseq` 188, `spatial_pseudobulk` 4. Both sum to 963. Group by whichever you actually mean.

### Tsao lab oxygen-tension RNA-seq

[`data/processed/tsaolab_o2/`](data/processed/tsaolab_o2/) holds 60 bulk RNA-seq samples (48 human, 12 mouse) from two Quintara orders, generated in the Tsao lab and released here with lab approval.

Five submitted cell lines (HL-1, NIH3T3, BJ, ATCC keloid, AC16) were each grown in four oxygen arms in triplicate — 5 × 4 × 3 = 60:

| Arm | `oxygen` | `treatment` | Human | Mouse |
|---|---|---|---|---|
| 20% O₂ throughout | `20%` | `none` | 12 | 3 |
| 1% O₂ throughout | `1%` | `hypoxia_1pct` | 12 | 3 |
| 20% → 1% O₂ | `20%->1%` | `normoxia_to_hypoxia_1pct` | 12 | 3 |
| 1% → 20% O₂ | `1%->20%` | `hypoxia_1pct_to_normoxia` | 12 | 3 |

| File | Contents |
|---|---|
| `tsaolab_o2_<species>_tpm_symbol.tsv.gz` | Gene symbol × sample TPM, not logged. **Start here.** |
| `tsaolab_o2_<species>_counts_symbol.tsv.gz` | Gene symbol × sample salmon estimated counts |
| `tsaolab_o2_<species>_tpm_ensembl.tsv.gz`, `_counts_ensembl.tsv.gz` | The same matrices keyed by Ensembl gene ID |
| `tsaolab_o2_<species>_genes.tsv.gz` | Ensembl gene ID → symbol and gene type |
| `tsaolab_o2_<species>_sample_metadata.csv` | One row per sample: cell line, oxygen arm, QC flags |
| `sample_sheets/` | One row per sample as submitted, from `make_sample_sheets.py` |
| `per_order/` | Per order and species: `*_harmonizable-USE.tsv.gz` log2(TPM+1) matrices, plus matching counts, gene annotation and metadata |
| `salmon_quant/` | Transcript-level salmon output (gzipped), readable by tximport |
| `qc/` | Read-level QC (`run_qc.tsv`), cell-identity QC (`qc_identity.csv`) and SNP genotyping (`genotype_match.tsv`) |

Human matrices carry 43,478 gene symbols across 48 samples; mouse, 34,195 across 12.

```python
import pandas as pd

tsao = "data/processed/tsaolab_o2"
tpm = pd.read_csv(f"{tsao}/tsaolab_o2_human_tpm_symbol.tsv.gz", sep="\t", index_col="gene_symbol")
meta = pd.read_csv(f"{tsao}/tsaolab_o2_human_sample_metadata.csv", index_col="sample_id")
```

#### Read this before analysing

Read the [dataset README](data/processed/tsaolab_o2/README.md) first. Three cell-line identities differ from the submission forms, each established by a different method:

- **Order 0013831748, samples 13–24** — all twelve submitted as NIH3T3 mouse fibroblasts — are HEK293-like **human** cells. Evidence: they map to the human transcriptome, carry adenovirus-5 E1 reads (967–2,401 per million), and are XIST-high with no real chrY signal. The HEK293 identity itself is an inference from those markers and still awaits STR profiling.
- **In the same order, 10 of the 12 BJ samples** carry HEK293-like admixture, estimated from their Ad5 E1 read fraction relative to the HEK-like samples (samples 27 and 35 are clean).
- **Order 1030390, samples 1–6** were submitted as the ATCC keloid line but are **BJ**. This one *is* from SNP genotyping against a 1,500-transcript mini-reference: they cluster with the BJ samples at 1.6–2.5% discordance and sit 72–77% from samples 7–12. Only samples 7–12 are genuinely keloid (ATCC CRL-1762, a single donor), and they were grown **only in the two arms that start at 1% O₂**.

The practical consequences: there are **6 real keloid samples, not 12**; keloid versus BJ is confounded with sequencing order; and the keloid samples do not span all four oxygen arms. 36 of the 48 human samples are flagged `ok` and 12 `suspect`; 14 are marked `eligible_for_keloid_pretraining`.

**The cardiomyocyte lines need the same scepticism.** AC16 is a fusion of adult ventricular cells with SV40-transformed fibroblasts, and in this dataset it carries essentially no cardiomyocyte identity: median TPM is 0.1 for TNNT2 and ~0 for MYH6, MYH7, TTN and NKX2-5, against COL1A1 ≈ 1,990 and VIM ≈ 6,660. Treat it as a fibroblast-like line, not a cardiomyocyte. HL-1 (mouse) does show the expected markers (Tnnt2 938, Myl7 3,610, Nkx2-5 311). Because AC16 is human and HL-1 mouse, cell identity is fully confounded with species, sequencing order and collection year — no cross-line comparison is identifiable here.

The raw FASTQs (about 258 GB) are too large for GitHub. They are on the lab Google Drive in `0013831748-2024 Raw_Data` and `1030390-2025 Raw_Data` — ask William for access. [`scripts/tsaolab_rnaseq/`](scripts/tsaolab_rnaseq/) rebuilds every file above from those FASTQs; see the [Quickstart](QUICKSTART.md#rebuild-the-tsao-lab-dataset).

## KeloidBench

[`keloidbench/`](keloidbench/) is the original keloid-only project, kept as published: an interpretable pipeline for keloid-versus-unaffected-skin transcriptomic classification, accepted as a poster at **NeurIPS GenAI4Health 2026**.

- [Paper on the NeurIPS site](https://neurips.cc/virtual/2026/loc/sydney/177621) · [PDF](keloidbench/KeloidBench_NeurIPS_GenAI4Health_2026_.pdf) · [Accepted poster](keloidbench/KeloidBench%20Poster_WW.pdf)
- [Archived interactive site](https://willibmwang.github.io/DermTrace/keloidbench/)
- [KeloidBench README](keloidbench/README.md) for results, caveats and reproduction

Run its code from that folder (`cd keloidbench`); it reads processed data from the shared [`data/`](data/) folder here. It is a research demonstration — its development results do not establish clinical performance.

## License

MIT — see [LICENSE](LICENSE). Copyright © 2026 willibmwang.

The processed public cohorts remain subject to the terms of their original GEO and ArrayExpress submissions; cite the original studies when you use them.

## Citation

If you use DermTrace or KeloidBench, please cite the KeloidBench paper:

```bibtex
@inproceedings{wang2026keloidbench,
  title     = {KeloidBench: Interpretable Machine Learning on Fibroblast Program
               Features for Keloid Detection},
  author    = {Wang, William and Tam, Joshua and Tsao, Hensin},
  booktitle = {NeurIPS Workshop on Generative AI for Health (GenAI4Health)},
  address   = {Sydney, Australia},
  year      = {2026},
  url       = {https://neurips.cc/virtual/2026/loc/sydney/177621}
}
```

When you use the public cohorts, cite the original GEO and ArrayExpress studies as well.

## Contact

William Wang (Yale University) · <william.b.wang@yale.edu>
Joshua Tam (Harvard-MIT HST) · Hensin Tsao (Harvard Medical School)

For data access, questions or collaboration, open an [issue](https://github.com/willibmwang/DermTrace/issues).
