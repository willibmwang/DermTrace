# KeloidBench

KeloidBench is an interpretable research pipeline for keloid-versus-unaffected-skin transcriptomic classification. It combines fibroblast gene-program features, linear classifiers, study-held-out evaluation, adjustable abstention, and a separate narrative layer that cannot change predictions.

Accepted as a **poster at NeurIPS GenAI4Health**. [View the accepted poster](KeloidBench%20Poster_WW.pdf) · [Interactive website](https://willibmwang.github.io/KeloidBench/)

## Reviewer revision results

<!-- REVISION_RESULTS_START -->

The primary historical ensemble result includes **8 studies / 101 profiles**: **0.826 mean macro F1**, with **0.182 on the worst study**. The **0.918** seven-study result is an exploratory subset selected after observing GSE173900's failure.

Under the new matched protocol, COL1A1 alone scores **0.914** and the fixed rank-program view **0.896**. Inner-selected programs versus markers score **0.855 versus 0.829** (paired difference +0.025; descriptive study-bootstrap interval [-0.108, +0.148]). These results do **not** establish general superiority of program features.

The revised template audit passes automated screening on **101/101 profiles**, including **21 abstentions**. LLM narratives (Qwen3-8B writes the explanatory prose; code inserts every number, the decision, the QC note, and the disclaimer): 290/303 outputs passed all automated checks (96/101 profiles on every repeat; profile-cluster interval [0.917, 0.990]). All failures were contribution-direction errors. 215/303 were also English-only; the demo shows generated prose only when it passes every check and is English-only, and otherwise shows the template. The format and checks were developed on these same development profiles, so this is not independent validation. Independent human rating was not performed.

<!-- REVISION_RESULTS_END -->

The poster is preserved as originally uploaded and may contain submission-time claims. The revision reports below are the current evidence. This is a research demonstration; its development results do not establish clinical performance.

## Evidence and revision status

- [Full-cohort evaluation](results/reviewer_revision/evaluation.md), [per-study table](results/reviewer_revision/all_studies.csv), and [coverage by study/class](results/reviewer_revision/selective_by_study.csv).
- [Matched single-marker comparison](results/reviewer_revision/matched_comparison.md), [frozen protocol](protocols/reviewer_revision_v1.json), and [per-profile predictions](results/reviewer_revision/matched_predictions.csv).
- [Narrative audit status and human-review rubric](results/reviewer_revision/narrative_audit/README.md).
- [Reviewer response matrix](REVIEWER_RESPONSE.md) and [camera-ready results text](manuscript/reviewer_revision_results.md).
- [Reproduction guide](REPRODUCIBILITY.md) and [revision plan/status](REVIEW_REVISION_PLAN.md).

The live demonstration uses a frozen Stage-B scorer. The historical nested ensemble and the new matched marker benchmark are distinct evaluation artifacts. The browser website uses deterministic summaries. The optional model-backed demo labels generation and validation fallbacks explicitly.

## Run locally

Use Python 3.10 and the tested revision environment:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-revision.txt
python scripts/revision_report.py
python scripts/run_marker_benchmark.py --report-only
python -m streamlit run demo/app.py
```

See [demo/README.md](demo/README.md) for input formats and narrative endpoint configuration. No model endpoint is required for scoring, report reproduction, or deterministic summaries.

## Data

### Tsao lab oxygen-tension RNA-seq

The processed Tsao lab bulk RNA-seq is in [`data/processed/tsaolab_o2/`](data/processed/tsaolab_o2/). The dataset has 60 samples (48 human, 12 mouse) covering the ATCC keloid fibroblast line, BJ, AC16, HL-1 and HEK293-like cells. Each line was grown at 20% O₂, at 1% O₂, or switched between the two, in triplicate.

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

Read the [dataset README](data/processed/tsaolab_o2/README.md) before analysis. Samples 13–24 of order 0013831748 were submitted as NIH3T3 but are HEK293-like human cells. Ten of the 12 BJ samples carry HEK293-like admixture. All 12 keloid samples come from a single cell line.

The raw FASTQs (about 258 GB) are too large for GitHub. They are on the lab Google Drive in the folders `0013831748-2024 Raw_Data` and `1030390-2025 Raw_Data`; ask William for access. [`scripts/tsaolab_rnaseq/`](scripts/tsaolab_rnaseq/) rebuilds every file above from those FASTQs.

## Historical research

The repository retains prior experiments, processed data, and broader research directions. [PLAN.md](PLAN.md) is the historical encoder-decoder roadmap; [manuscript/draft.md](manuscript/draft.md) contains development history. Older result files are preserved for provenance and are superseded for current claims by the reviewer revision reports.
