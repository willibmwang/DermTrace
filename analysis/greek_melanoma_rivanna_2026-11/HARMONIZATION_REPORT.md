# Harmonized quantification of all 42 eligible specimens

All Rivanna results below were run and **reported by the user**.

## Why harmonize

- In the mixed matrix, 33 columns come from the collaborator's pipeline and
  9 from ours, so pipeline origin cannot be told apart from specimen set.
- Re-quantifying all 42 specimens with one qualified pipeline removes that
  confound.
- It also gives a paired check of all 33 collaborator quantifications.

## Qualification chain behind the pipeline

| Level | Evidence (user-reported) |
|---|---|
| Real-tool qualification (r2, pilot) | 14/14 provisional-v1 gates. Transcript Pearson log1p TPM 0.998858; mapping-rate difference −0.0732 percentage points |
| Run-to-run (pilot) | Transcript Pearson 0.999890 and gene Pearson 0.999997 between the two runs. 0/18,529 expressed genes and 18/41,838 expressed transcripts differ by more than two-fold. Mapped fragments differ by 12 |
| Cross-specimen (four r2 specimens) | 56/56 gates (70/70 including pilot); lowest transcript Pearson 0.998858, lowest gene Pearson 0.999800. 28 distinct discordant genes, none recurring across specimens |
| Parser-specific qualification (r3) | r2 and r3 make identical eligibility decisions with identical resolved inputs for all 48 specimens; regression suites pass |

## Campaign

- **Reused outputs (14):** the five qualification specimens (r2) and the
  nine expansion specimens (r3). They were reused **only** after a fail-closed
  compatibility check:
  - structure valid;
  - provenance complete;
  - code identity r2 or r3;
  - every recorded configuration value identical to the pilot configuration;
  - collaborator index fingerprints present.
- **Quantified now (28):** one bounded Slurm array of 28 tasks,
  running the unchanged qualified `quantify_specimen.sbatch`. All 28 tasks
  completed with exit 0:0.
- **Recovery:** the first collection step stopped with exit 141 before
  writing anything. That was a shell SIGPIPE in the accounting display, not
  a data problem. It was replaced by a versioned recovery
  (`procedure-20261127r1`) without rerunning any quantification.

## Results (user-reported)

| Measure | Result |
|---|---|
| Compatible quantifications | **42/42** |
| Paired provisional-v1 concordance (33 collaborator specimens) | **33/33 PASS** |
| Acceptance gates | **462/462** (33 specimens × 14 gates) |
| Lowest transcript-level Pearson, log1p TPM (33 specimens) | **0.998545** |
| Lowest gene-level Pearson, log1p TPM (33 specimens) | **0.999699** |
| Harmonized matrix | **47,647 Ensembl genes × 42 specimens; 46,872 gene symbols** |
| Ensembl genes present in only one of the mixed and harmonized matrices | 108 |

**How the harmonized matrices were built.**
- They were built from scratch from all 42 outputs, with the collaborator's
  aggregation recipe (validated by rebuilding the committed 33 exactly).
- Columns are in the same order as the mixed matrix.
- They were written to a separate output. The collaborator matrices and the
  mixed matrices are kept unchanged.

**Provisional-v1, briefly.** Transcript-level Pearson of log1p TPM ≥ 0.99 and
absolute mapping-rate difference ≤ 1 percentage point, plus identity,
integrity and index-fingerprint gates. The thresholds were frozen before any
result and never changed.

**Primary representation.** harmonized42 is adopted as the primary research
representation. Earlier matrices are kept for comparison.

## What this establishes, and what it does not

**Established (computational reproducibility).** The qualified pipeline
reproduces the collaborator's quantification for every one of the 33 paired
specimens under provisional-v1. All 42 specimens now share one documented,
provenance-recorded pipeline.

**Not established:**
- That provisional-v1 is a validated equivalence margin. It is an
  engineering rule.
- The collaborator's original fastp version and Salmon k-mer length. Both
  remain unverified; we reproduce outputs with a candidate fastp (0.23.4)
  and k = 31.
- Anything about the biology or clinical meaning of the expression values.
  See `TECHNICAL_CONFOUNDING_REPORT.md` and `RESEARCH_ROADMAP.md`.
