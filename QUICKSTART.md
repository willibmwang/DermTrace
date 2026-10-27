# Quickstart

From `git clone` to a loaded expression matrix in about five minutes, then on to reproducing KeloidBench and rebuilding the RNA-seq pipeline.

- [1. Set up](#1-set-up) · 5 min
- [2. Load the data](#2-load-the-data) · 5 min
- [3. Reproduce KeloidBench](#3-reproduce-keloidbench) · 20 min
- [4. Rebuild the Tsao lab dataset](#4-rebuild-the-tsao-lab-dataset) · hours, needs FASTQs
- [5. Work on the website](#5-work-on-the-website) · 5 min

Not sure which you want? If you came for data, do 1 and 2 and stop. Everything after that is for rebuilding things that are already committed.

---

## 1. Set up

```bash
git clone https://github.com/willibmwang/DermTrace.git
cd DermTrace
```

The repository is about 400 MB; most of it is the committed Tsao lab matrices and the KeloidBench results. To skip the large blobs until you touch them, clone with `--filter=blob:none`.

For loading data you need **pandas** and **pyarrow** (the training feature views are Parquet):

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install pandas pyarrow
```

> There is no requirements file for the root-level scripts. [`scripts/tsaolab_rnaseq/`](scripts/tsaolab_rnaseq/) needs only `pandas`, `numpy`, `scipy` (for `genotype.py`) and `openpyxl` (for `make_sample_sheets.py`). Section 3 uses a different, fully pinned environment — keep them separate.

## 2. Load the data

Only two folders ship loadable matrices: `tsaolab_o2/` and `training/`. See the [data table](README.md#what-you-get-when-you-clone) for what the others contain.

### The Tsao lab oxygen-tension dataset

```python
import pandas as pd

tsao = "data/processed/tsaolab_o2"
tpm = pd.read_csv(f"{tsao}/tsaolab_o2_human_tpm_symbol.tsv.gz", sep="\t", index_col="gene_symbol")
meta = pd.read_csv(f"{tsao}/tsaolab_o2_human_sample_metadata.csv", index_col="sample_id")

print(tpm.shape)   # (43478, 48)
print(meta.shape)  # (48, 40)
```

`tpm.columns` and `meta.index` are the same 48 sample IDs in the same order, so you can slice directly:

```python
# Hypoxia response in the genuine keloid samples (order 1030390, samples 7-12)
keloid = meta[(meta.cell_line == "ATCC_keloid") & (meta.quality_flag == "ok")]
print(keloid.oxygen.value_counts())   # 1%: 3, 1%->20%: 3

hypoxic = tpm[keloid.index[keloid.oxygen == "1%"]]
reoxy   = tpm[keloid.index[keloid.oxygen == "1%->20%"]]
```

Swap `human` for `mouse` to get the 12 mouse samples (34,195 genes). Values are **TPM, not logged** — log them yourself if your method expects it.

#### Two sample-ID conventions

This trips people up. The dataset uses two ID spaces, and they are not interchangeable:

| Where | ID looks like | Column |
|---|---|---|
| Top-level matrices, `sample_metadata.csv` | `1030390-7` | `sample_id` |
| `per_order/`, `qc/`, `salmon_quant/` | `A7`, `Q13` | `sample_id` **there**, with `qb_id` holding the long form |

So the QC tables join to the metadata on a *different* column name:

```python
meta = pd.read_csv(f"{tsao}/tsaolab_o2_human_sample_metadata.csv")
qc   = pd.read_csv(f"{tsao}/qc/run_qc.tsv", sep="\t")
joined = meta.merge(qc, left_on="sample_id", right_on="qb_id")  # (48, 61)
```

The `per_order/*_harmonizable-USE.tsv.gz` files are already `log2(TPM+1)`, but they are keyed by Ensembl ID *and* by the short per-order IDs, so they need both a gene and a sample remap before they line up with the combined matrices.

**Before you interpret anything**, read [the three identity caveats](README.md#read-this-before-analysing). Six samples submitted as keloid are genetically BJ, twelve submitted as mouse NIH3T3 are human, and keloid versus BJ is confounded with sequencing order.

### The combined training corpus

963 profiles with 18 precomputed feature views, ready for modelling:

```python
import pandas as pd

tr = "data/processed/training"
manifest = pd.read_parquet(f"{tr}/profile_manifest.parquet")      # (963, 54)
X = pd.read_parquet(f"{tr}/features/rank_programs_only.parquet")  # (963, 7)

df = manifest.merge(X, on="sample_id")
print(df.source_dataset.value_counts())
print(df.keloid_vs_normal.value_counts(dropna=False))
```

There are **18** feature views in `data/processed/training/features/`, each with a matching `*_columns.txt` — from 6-column program scores (`rank_programs_only`) up to a 5,034-column shared-gene matrix (`shared_genes`). Predefined study-held-out splits are in `splits/splits.json` (a JSON **list**, not an object). `endpoint_eligibility_matrix.csv` tells you which profiles are valid for which task.

### Which folder should I use?

| If you want to... | Use |
|---|---|
| Train or evaluate a model | `training/` — this is what the frozen protocol consumes |
| Work with the in-house oxygen experiment | `tsaolab_o2/` |
| Look up public-cohort labels and provenance | `public_keloid/` (metadata only, no matrix) |
| Find a matrix that used to exist | `archive/2026-07-09/` — provenance only, superseded |
| — | `bulk_rnaseq/`, `microarray/`, `scrna_spatial/`, `spheroid/` are label-and-summary stubs; their data reaches you through `training/` |

## 3. Reproduce KeloidBench

KeloidBench has a pinned environment and a CI-verified command sequence. Use **Python 3.10** and run everything from the `keloidbench/` folder:

```bash
cd keloidbench
python3.10 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-revision.txt

python -m unittest discover -s tests -v              # 46 tests
python scripts/revision_report.py                    # primary all-study report
python scripts/run_marker_benchmark.py --report-only  # matched marker comparison
```

This is exactly what [`.github/workflows/revision-checks.yml`](.github/workflows/revision-checks.yml) runs, so it is the path most likely to work. It is fully offline and fast — the report takes seconds, the full re-fit under a minute on CPU. No GPU and no model endpoint are needed for scoring or report reproduction.

> **These commands rewrite tracked files in place.** `revision_report.py` and `run_marker_benchmark.py` overwrite their outputs under `keloidbench/results/reviewer_revision/`. In the pinned Python 3.10 environment they reproduce byte-identically and `git status` stays clean; under any other dependency set the numbers still agree to ~1e-13 but the bytes differ, which fails CI's diff gate. Check `git status` afterwards and `git checkout` the results folder if you only wanted to look.

To re-fit the matched benchmark rather than re-deriving it from saved predictions:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python scripts/run_marker_benchmark.py
```

Launch the interactive demo:

```bash
python demo/launch_demo.py --example confident_keloid
# or: python -m streamlit run demo/app.py
```

`--list-examples` shows all four (`confident_keloid`, `confident_unaffected`, `abstention`, `transfer_failure`) without starting a server.

> **Do not hand-edit** anything under `keloidbench/results/reviewer_revision/`, or `keloidbench/README.md` between its `<!-- REVISION_RESULTS_START/END -->` markers. Those are generated by `revision_report.py`, `run_marker_benchmark.py` and `generate_revision_docs.py`, and CI fails on any diff.

KeloidBench scripts find the shared `data/` folder from their own file path (`Path(__file__).resolve().parents[1].parent`), so data loading works from any working directory. The `cd keloidbench` matters for the *relative* paths in these commands.

**Read the results in context.** The headline is 8 studies / 101 profiles, mean macro F1 0.826 — but the worst study (GSE173900) scores 0.182, which fails the project's own 0.85/0.75 gates, and COL1A1 alone scores 0.914, beating every program view. The 0.918 seven-study figure was selected after observing that failure. See [keloidbench/README.md](keloidbench/README.md) for the full set of caveats.

## 4. Rebuild the Tsao lab dataset

Everything in `data/processed/tsaolab_o2/` is already committed. You only need this section to rebuild from raw FASTQs.

**You need, and the repository does not provide:**

- The **raw FASTQs** (~258 GB), on the lab Google Drive — ask William for access.
- A **reference folder** `$REF_DIR` that you build yourself. No script downloads or builds it. It must contain:
  - `gencode_human_v46/` with both `salmon_index/` and `gencode.v46.transcripts.fa.gz`
  - `gencode_mouse_vM35/` with both `salmon_index/` and `gencode.vM35.transcripts.fa.gz`
  - `adenovirus5_AC_000008.1/salmon_index` (only for the optional Ad5 screen)
- A conda env with **salmon 1.10.3**, **fastp** and **rclone** (rclone must be inside the env — it is called by absolute path).

The committed salmon runs used a transcriptome-only index with no decoys.

### The main path

```bash
DATA_DIR=/path/to/fastq/parent \
REF_DIR=/path/to/references \
RNASEQ_ENV=/opt/anaconda3/envs/rnaseq \
bash scripts/tsaolab_rnaseq/run_all.sh 0013831748
```

Then the same for order `1030390` — those are the only two accepted values. `run_all.sh` downloads the order, verifies MD5s, runs fastp and salmon, aggregates per-order tables, and — once *both* orders have finished — builds the combined dataset. It is safe to re-run; finished work is skipped. Progress goes to `$OUT_DIR/STATUS_<order>`.

`DATA_DIR` and `REF_DIR` are required. `OUT_DIR` defaults to `$DATA_DIR/processed`, `RNASEQ_ENV` to `/opt/anaconda3/envs/rnaseq`, `PYTHON` to `python3`.

### Two things `run_all.sh` will not do for you

**It does not produce the committed metadata.** `run_all.sh` calls `build_dermtrace.py` without `--qc` or `--runqc`, so seven QC columns are missing and every sample's `platform_id` silently defaults to `Element_AVITI` — wrong for the one NextSeq sample. To match what is committed, run the optional QC steps and then rebuild by hand:

```bash
# optional QC, in this order
bash   scripts/tsaolab_rnaseq/screen_ad5.sh   <sheet.csv> <raw_dir> <ad5.tsv>
python -I scripts/tsaolab_rnaseq/qc_identity.py   qc_identity.csv <sheet.csv>=<quant_dir>=<ad5.tsv>
python -I scripts/tsaolab_rnaseq/run_qc_summary.py run_qc.tsv     <sheet.csv>=<quant_dir>=<raw_dir>

# then rebuild with QC merged in (--qc/--runqc must come before the sheet=quant pairs)
REF_DIR=... python -I scripts/tsaolab_rnaseq/build_dermtrace.py out_dir \
  --qc qc_identity.csv --runqc run_qc.tsv \
  <sheet_1>=<quant_1> <sheet_2>=<quant_2>
```

Run the Python scripts with `python -I` as shown; they are not executable, and `build_dermtrace.py` and `qc_identity.py` import `aggregate.py` (which reads `REF_DIR` at import time, so an unset `REF_DIR` fails immediately).

**You cannot re-run the Python stages against the committed tree.** The scripts expect `quant_<order>/quant/<sample>/quant.sf` and write plain `.tsv`; the repository ships `salmon_quant/<order>/<sample>/` with `quant.sf.gz` and `.tsv.gz`. Decompress and rearrange to match, or start from FASTQs.

### Optional: SNP genotyping

This is what established the cell-line identities. It needs a *separate* conda env with minimap2, samtools and bcftools:

```bash
python -I scripts/tsaolab_rnaseq/genotype.py ref <transcripts.fa.gz> mini_ref.fa <quant_dir>...
GENO_ENV=/path/to/geno/env bash scripts/tsaolab_rnaseq/genotype_check.sh mini_ref.fa out_dir <sheet.csv>=<raw_dir>
python -I scripts/tsaolab_rnaseq/genotype.py compare out_dir/snps.tsv <sheet.csv> --out out_dir
```

`genotype.py ref` picks 1,500 transcripts expressed above 20 TPM everywhere; `compare` clusters donors at 35% discordance. `--out` is mandatory. You can add an external reference by passing `<name>=<R1 URL>` to `genotype_check.sh`, which is how the public ATCC CRL-1762 run was compared.

The sample sheets in `data/processed/tsaolab_o2/sample_sheets/` are already committed; `make_sample_sheets.py` only regenerates them from the original Quintara spreadsheets and is not needed.

## 5. Work on the website

The site source is in [`website/`](website/) and is published from the **`gh-pages` branch**, not from `main`. Nothing in CI builds or deploys it — the two are kept in sync by hand.

Preview locally:

```bash
python3 -m http.server 8765 --directory website
```

To publish, copy the changed files to the root of `gh-pages` and commit there.

To regenerate the KeloidBench demo data that the archived site uses:

```bash
cd keloidbench && python scripts/export_static_site_data.py   # writes website/assets/data/demo-data.json
```

> **Never regenerate `gh-pages:keloidbench/` from `main`.** It is a frozen snapshot of the KeloidBench site at commit `be60f64`, with its repository links deliberately pinned to that commit. `main/keloidbench/website/index.html` still points at the old `willibmwang/KeloidBench` repository, so re-syncing it would silently break the archive.

To preview the archived KeloidBench site instead:

```bash
python3 -m http.server 8765 --directory keloidbench/website
```

---

## Where to go next

| You want to... | Go to |
|---|---|
| Understand a dataset's caveats | [`data/processed/tsaolab_o2/README.md`](data/processed/tsaolab_o2/README.md) |
| Read the paper | [NeurIPS 2026 GenAI4Health](https://neurips.cc/virtual/2026/loc/sydney/177621) · [PDF](keloidbench/KeloidBench_NeurIPS_GenAI4Health_2026_.pdf) |
| See KeloidBench's results and limits | [`keloidbench/README.md`](keloidbench/README.md) |
| Reproduce the published numbers | [`keloidbench/REPRODUCIBILITY.md`](keloidbench/REPRODUCIBILITY.md) |
| Find a removed matrix | [`data/processed/archive/README.md`](data/processed/archive/README.md) |
| Ask a question | [open an issue](https://github.com/willibmwang/DermTrace/issues) |
