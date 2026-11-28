# Qualification status of the Rivanna quantification scripts

This file gives the current status of the three levels of evidence. They
are kept separate on purpose: none of them implies the others.

All Rivanna results here were run and **reported by the user**. The
developer had no direct access to Rivanna. Local tests use stub tools only.

| Level | Version | Status | What it shows | What it does not show |
|---|---|---|---|---|
| 1. Real-tool qualification | r2 `rqc-src-20261124-r2-08dd537f` | **PASSED** on GM01 | Real fastp and Salmon on Rivanna reproduce the collaborator's committed GM01 output under provisional-v1 | That other specimens agree; the collaborator's fastp version or k; a validated equivalence margin |
| 2. Parser-specific qualification | r3 `rqc-src-20261125-r3-4432a801` | **PASSED** | The stricter manifest parser accepts the real `validated_pairs.tsv` and makes the same eligibility decision, with the same resolved inputs, as r2 for every specimen | Anything about expression values. It quantifies nothing |
| 3. Cross-specimen scientific validation | r2 code on GM02, GM08, GM11 and GM31 (plus GM01) | **PASSED**: 70 of 70 provisional-v1 gates across 5 specimens | The pipeline reproduces the collaborator's committed output across both sequencing runs and a range of QC profiles. The remaining differences are specimen-specific and do not move the frozen PCA | A validated equivalence margin; the collaborator's fastp version or k; anything about the nine specimens that have no collaborator output |

## 1. r2 real-tool qualification (GM01)

The full record is in `QUALIFICATION_GM01_r2.md`. Its file hashes are the
**deployment identity** of r2, the code that actually ran.

The repository holds r3, whose manifest SHA-256 is `4432a801…8f83`. r3
differs from r2 only in `rivanna_qc_tools.py` (`read_validated_pairs`) and
`tests/run_local_tests.sh`; see section 2.

**Production run, job 23402243:**
- Completed in 32 min 20 s with exit 0:0. `QUANT_OK`, and the output is
  structurally valid with complete provenance.
- All 14 provisional-v1 gates PASS.
- Transcript-level Pearson of log1p TPM: 0.998858 (threshold ≥ 0.99).
  Mapping-rate difference: −0.0732 percentage points (threshold ±1).
- Discordant expressed genes: 12 of 18,535. Discordant expressed
  transcripts: 129 of 41,896.

**Run-to-run check against the historical candidate run, job 23390924:**
- Transcript-level Pearson 0.999890; gene-level Pearson 0.999997.
- Discordant between the two runs: 0 of 18,529 genes and 18 of 41,838
  transcripts.
- Reads kept after trimming are identical, and mapped fragments differ
  by 12.

**Discordance against the collaborator is stable at gene level:**
- The same 12 genes are discordant in both runs (Jaccard 1.0).
- At transcript level the Jaccard is 0.8582: 115 shared, 5 only in the
  historical run, 14 only in the production run.
- So gene-level discordance comes from a fixed difference between the
  pipelines, not from run-to-run noise. Its cause has not been
  established. One candidate is the fastp version: 278,440 more reads are
  kept after trimming than in the collaborator's run.

## 2. r3 parser-specific qualification

Run with step 10 of the Rivanna procedure; nothing was quantified.

| Check | Result |
|---|---|
| Deployment manifests (r2, r3, procedure) | PASS |
| r2 and r3 identical except the manifest parser, tests, documentation and records | PASS. Every other script is byte-identical |
| r3 accepts the real `validated_pairs.tsv` | PASS. The manifest was accepted, and all 42 eligible rows carry core-ID labels. The report does not record the header line or the total row count |
| Eligibility for all 48 sample-sheet specimens, r2 vs r3 | Identical exit codes and byte-identical resolved output for all 48 |
| Eligible specimens | 42: 11-01 to 11-42 and 11-50, without 11-37 |
| Blocked specimens | Exactly GM37 (provenance adjudication pending) and GM43, GM44, GM46, GM47, GM48 (corrupt FASTQ) |
| Quantification regression suite (Rivanna's Python) | 188 passed, 0 failed |
| Legacy QC regression suite, including `ParserParity` | PASS |

The report is `eligibility_r2_vs_r3.json`, SHA-256
`98b9f440225c125597403a3bb4b01c84c48410fd14a318372326dee6d812360f`.
Its 42 resolved records were checked independently against the committed
sample sheet: core IDs, file names and MD5s all match.

**Why no new GM01 quantification was needed.**
- `quantify_specimen.sbatch` uses only what `resolve-specimen` prints.
- On the real manifest, that output is byte-identical under r2 and r3.
- So the input files, fastp and Salmon arguments, index and container are
  identical, and the r2 GM01 result carries over to r3 unchanged.
- r3 provenance records a different `tools_sha256`, so runs made with r2
  and with r3 can still be told apart.

**Behaviour change to be aware of.** r3 rejects the whole manifest on any
defective row. r2 rejected only the affected specimen. Any future edit to
`validated_pairs.tsv` that introduces a defect makes **no** specimen
eligible until the defect is fixed.

## 3. Cross-specimen scientific validation (passed)

The plan was fixed before any result: `PLAN_xspecimen_r2.md`, kept with the
Rivanna procedure outside this repository.

- **Specimens:** GM02, GM08, GM11 and GM31. They span sequencing runs
  LH00150:1058 and LH00150:1062, adapter content 45–68 %, duplication
  0.15–0.62 and insert peak 115–138.
- **Code:** frozen r2, with the same configuration as GM01. Only the output
  folders differed.
- **Jobs:** 23416488 to 23416491.

**Primary endpoint: met.**
- GM02, GM08, GM11 and GM31 pass 56 of 56 provisional-v1 gates. With GM01,
  70 of 70 gates pass across the five specimens.
- The thresholds were not changed.

**Secondary results (reported, not gated):**

| Measure | Result |
|---|---|
| Lowest transcript-level Pearson, log1p TPM, across the five | 0.998858 |
| Lowest gene-level Pearson, log1p TPM, across the five | 0.999800 |
| Discordant expressed genes (\|log2FC\| > 1) | 28 distinct genes across the five specimens; **none discordant in more than one specimen** |
| GM01's 12 discordant genes | None recurs in another specimen |
| Effect on the frozen exploratory PCA | Largest shift on PC1–PC3, over all five specimens: **0.0049 cohort SD** |
| Duplication rate vs transcript-level discordance | Spearman rho 0.90, exact permutation p 0.0833, n = 5. Exploratory only; it suggests a hypothesis and is not evidence |

**Interpretation.**
- No gene is discordant in more than one specimen. So these differences do
  not form a systematic, gene-specific bias that would need correcting.
- The effect on the frozen PCA is negligible.
- No feature was removed and no threshold was changed.

`QUALIFICATION_GM01_r2.md` keeps its original text as the record of its
own date, with an added supersession note. Its limitation 2 ("Only one specimen was tested") is answered by this
section.

**GM40** had no collaborator output to compare against. It was held back
from this validation, and it joins the production expansion described
below.

## 4. Production expansion (prepared; separate evidence)

Nine eligible specimens have no committed collaborator output: GM33, GM34,
GM35, GM36, GM38, GM39, GM40, GM41 and GM42.

- **How they are quantified:** with the parser-qualified r3 code and the
  r2-qualified configuration.
- **What each must show:** `QUANT_OK`, structural validation, and complete
  r3 provenance.
- **What cannot be shown:** concordance with the collaborator, because no
  reference output exists for them.
- **How the gene matrices are extended:**
  - the collaborator's 33 columns are kept byte-identical;
  - exactly nine columns are appended;
  - no gene is discarded.
- **How the result is recorded:** pipeline origin is a technical
  covariate in the expanded QC manifest and in the exploratory comparison
  of the 42-specimen and 33-specimen matrices.

Results of the expansion will be recorded separately once it has run.
