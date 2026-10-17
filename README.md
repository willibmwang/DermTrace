# KeloidBench

KeloidBench is an interpretable research pipeline for keloid-versus-unaffected-skin transcriptomic classification. It combines fibroblast gene-program features, linear classifiers, study-held-out evaluation, adjustable abstention, and a separate narrative layer that cannot change predictions.

Accepted as a **poster at NeurIPS GenAI4Health**. [View the accepted poster](KeloidBench%20Poster_WW.pdf) · [Interactive website](https://willibmwang.github.io/KeloidBench/)

## Reviewer revision results

<!-- REVISION_RESULTS_START -->

The primary historical ensemble result includes **8 studies / 101 profiles**: **0.826 mean macro F1**, with **0.182 on the worst study**. The **0.918** seven-study result is an exploratory subset selected after observing GSE173900's failure.

Under the new matched protocol, COL1A1 alone scores **0.914** and the fixed rank-program view **0.896**. Inner-selected programs versus markers score **0.855 versus 0.829** (paired difference +0.025; descriptive study-bootstrap interval [-0.108, +0.148]). These results do **not** establish general superiority of program features.

The revised template audit passes automated screening on **101/101 profiles**, including **21 abstentions**. Frozen v1 (Qwen3-1.7B writes the whole narrative): 303/303 generations completed; 0/303 raw outputs passed all automated checks across 101 profiles (3 greedy requests per profile). In a post hoc diagnostic, removing echoed prompt labels raised this only to 6/303: 297/303 outputs restated numbers despite the instruction not to, and 609/933 numbers attached to a named program matched its contribution (255 reported its input value instead; 69 matched neither). Development-informed v2 (Qwen3-8B writes only the contribution prose; code adds the fixed opening, QC note, and disclaimer): 303/303 generations completed and 290/303 assembled outputs passed all automated checks, including a stricter clause-level contribution-role check (96/101 profiles on every repeat; profile-cluster interval [0.917, 0.990]). Post hoc, 75/303 outputs contained a non-English token (Qwen code-switching) that no check covered, so 215/303 were both check-passing and English-only; the demo now falls back to the template for such text. v2 was designed after inspecting v1 on the same packets and is not an independent validation. Independent human rating was not performed in this revision; automated screening is not human-validated faithfulness, and blank two-rater forms are released with each run.

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

## Historical research

The repository retains prior experiments, processed data, and broader research directions. [PLAN.md](PLAN.md) is the historical encoder-decoder roadmap; [manuscript/draft.md](manuscript/draft.md) contains development history. Older result files are preserved for provenance and are superseded for current claims by the reviewer revision reports.
