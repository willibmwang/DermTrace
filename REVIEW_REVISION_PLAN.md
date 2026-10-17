# KeloidBench reviewer revision plan

Prepared 2026-10-17 from the supplied NeurIPS GenAI4Health reviews and local repository inspection. Status: proposed work; no new experiments have been run for this plan.

## Execution update

The plan was pushed before implementation. The primary eight-study report has been reproduced, the matched marker/program experiment has run, the audit has been hardened, and the deterministic audit passes on 101 profiles. Current status and evidence are tracked in [REVIEWER_RESPONSE.md](REVIEWER_RESPONSE.md). The frozen local GPU LLM audit (Qwen3-1.7B, SLURM job 2865612) completed with 0/303 raw outputs passing all automated checks. This is reported as a negative primary result with a separately labeled post hoc failure diagnostic. A development-informed v2 protocol was frozen before generation: Qwen3-8B writes only the prose, and code inserts all facts. It passed 290/303 (215/303 English-only). Independent human rating was not performed and is stated as a limitation. The sections below preserve the original plan.

## Objective

Address the three concerns shared by both reviewers and the area chair: validate actual generated narratives, report cross-study performance without outcome-dependent exclusions, and compare program features against individual markers under a matched protocol. Keep the interpretable predictor and separate narrative layer as the core contribution.

## Findings that determine the plan

| Concern | Existing evidence | Implication |
| --- | --- | --- |
| Generalization | `results/accuracy_085_push/expanded_ensemble_v1.json` reports eight-study mean macro F1 **0.825951**, worst fold **0.181818**, and failed development gates. The seven-study report gives **0.917971**. | Make the full eight-study evaluation primary. Keep the seven-study analysis explicitly exploratory and selected after observing outcomes. |
| Exclusion independence | Appendix A3 explicitly says GSE173900 was excluded after outer-fold failure. Its diagnostic grid uses held-out labels. | Documentation cannot make this exclusion independent of test performance. A future label-independent eligibility rule requires separate validation. |
| Generative validation | The saved deterministic audit reports 101 profiles, 21 abstentions, and all eight checks passing. The runner supports `--arm llm`, but the generation function returns a deterministic fallback when endpoint/model settings are absent. | Require evidence of actual generation; an LLM-labeled output file alone is insufficient. |
| Marker comparison | Existing ablations compare marker unions and modules on the broader `keloid_binary` endpoint using weighted F1. The marker replication audit measures biological consistency. | Neither establishes superiority over individual-marker classifiers on the current primary endpoint and macro-F1 protocol. |
| Release consistency | Local `main` has diverged from GitHub and includes staged research artifacts. The README/older plan describe a broader encoder-decoder direction. | Integrate deliberately, preserve user work, and align the public entry points with the accepted contribution. |

## 1. Freeze the evaluation and correct the claims

Relevant code: `scripts/score_primary_ladder.py`, `scripts/train_accuracy_085_push.py`, `scripts/train_nested_loso.py`.

- Record a versioned protocol with endpoint, eligible cohorts, sample/donor identifiers, fixed outer folds, feature definitions, inner selection rules, calibration, thresholds, metric aggregation, seeds, and artifact hashes.
- Retain all eight eligible development studies in the primary report. Recompute from saved predictions and reconcile the result with the existing 0.825951 mean and 0.181818 worst fold before changing any model.
- Document the current mixed aggregation explicitly: donor-level scoring where appropriate, sample-level fallback for single-donor and mixed-class donor designs. Show sample-level results alongside it and use the same rules for every baseline.
- Report per-study scores, sample/donor counts, class counts, pooled results, worst-fold performance, uncertainty, and full-coverage versus selective results separately. For abstention, report coverage by study and class, including folds with no evaluable retained cases.
- Describe the seven-study result as a post hoc sensitivity analysis. Remove claims that its inclusion criteria were independent of test performance, that the failure proves a platform-wide limitation, or that an oracle diagnostic establishes a prospective performance ceiling.
- Keep GSE173900 as a documented failure case. Inspect provenance, labels, expression processing, and program coverage; label any resulting repair as development informed by observed outcomes.
- Preserve the reserved lockbox and its existing gate policy. Verify its access history before calling it untouched. Do not use it to choose revisions.
- If future cohort screening is pursued, use information available before outcome evaluation (for example, assay metadata and coverage), fit any learned screening only within training folds, and validate on independent data. A rule using held-out AUROC remains outcome-dependent even if written down in advance.

Completion evidence: one reproducible eight-study table, a clearly labeled seven-study sensitivity table, and consistent wording across manuscript, demo, website, and README.

## 2. Validate actual LLM narratives and the audit itself

Relevant code: `scripts/audit_copilot_faithfulness.py`, `demo/demo_utils.py`, `tests/test_demo_utils.py`.

- Make LLM audit mode fail immediately when generation is unconfigured. Record requested arm, actual mode, model/version, prompt hash, packet hash, generation settings, run identifier, and API failures per request. Never count template fallback as a successful generated narrative.
- Freeze the prompt and evaluation packets after development. Start with the existing 101-profile corpus, including all 21 current abstentions and the failed cohort; add separately labeled challenge cases for near-threshold decisions, poor coverage, missing evidence, and unseen cohorts.
- Separate prompt-development examples from the final audit. Preserve raw responses and report all final attempts, including failures and truncations. A practical initial design is three requests per packet to characterize repeatability, with repeated outputs treated as clustered observations rather than independent profiles.
- Strengthen semantic checks: numerical values must match their named fields, class statements must preserve the decision, contribution directions must be correct, and abstention/QC limitations must survive. Current numeric membership and class-word checks cannot establish those properties by themselves.
- Add targeted tests with swapped probabilities, contradictory class statements, suppressed abstentions, unsupported claims embedded in disclaimer sentences, empty responses, and fallback mislabeled as generation.
- Check the actual outgoing request for redaction and verify prediction immutability. Current redaction checking reconstructs an expected packet rather than inspecting the transmitted payload.
- Have two independent human raters review a preselected stratified sample using a written rubric. Report agreement and adjudication, plus all automatically flagged cases separately. Do not treat an LLM judge as the sole reference standard.
- Publish per-check and joint pass rates with denominators, generation-failure rates, subgroup results, uncertainty, and representative failures. Distinguish generation success, factual faithfulness, and usefulness.
- For usefulness, compare the frozen LLM arm against the deterministic template on a small blinded clarity/completeness rating task. Report parity or worse performance if observed; this is secondary to the reviewers' required faithfulness audit.
- After measuring raw generation, add a runtime validation/fallback policy for the demo. Report raw and displayed-output outcomes separately so fallback does not conceal generation errors.

Completion evidence: replayable generated-response artifacts, provenance, automated and human audit tables, and explicit separation from deterministic results. Execution needs a configured model endpoint and human raters; neither is assumed available by this plan.

## 3. Add matched single-marker baselines

Relevant code: `scripts/gene_modules.py`, `scripts/train_accuracy_085_push.py`, `scripts/train_nested_loso.py`.

- Freeze an initial marker panel before this new comparison: POSTN, COL1A1, COL3A1, ACTA2, FN1, and CTHRC1, with rationale drawn from the repository's existing marker definitions and literature records. Describe this as a revision-time specification, not a historical preregistration.
- Run each marker alone on `keloid_vs_unaffected_skin`, using identical outer folds, cohort membership, donor weighting, missingness policy, preprocessing boundaries, calibration, and threshold-selection protocol.
- Compare single markers with the fixed program feature views using the same linear classifier family and tuning budget. Separate this controlled feature comparison from the full program-ensemble comparison so model complexity is visible.
- Add an optional best-single-marker arm in which marker selection occurs only inside each outer training fold. Never choose the winning marker using the aggregate outer-test results.
- Record coverage and missing markers by cohort. Do not silently drop difficult cohorts or change eligible samples for a favorable comparison. Fit imputation/scaling on training data only.
- Save predictions and selected configurations for every arm. Report paired per-study differences, equal-study means, worst-fold scores, calibration, and coverage where abstention is used. Keep uncertainty at the study/donor level and acknowledge the small number of studies.
- Inspect whether data-derived gene panels or preprocessing were built using outer-held-out labels; if so, rebuild within folds or label that arm exploratory.
- Retain superiority claims only to the extent supported by the matched comparison. A similar-performing, interpretable program representation is a valid outcome.

Completion evidence: a reproducible all-eight-study comparison table and per-fold prediction files; the seven-study comparison is secondary.

## 4. Publish a reproducible research release

- Build revision work on a branch based on current GitHub `main`; explicitly bring over the local code and staged artifacts needed for the accepted paper. Preserve unrelated local work and check dependencies before integrating the diverged commits.
- Replace the README's lead description with the accepted pipeline, acceptance/poster link, actual evaluation scope, narrative mode, and reproducibility entry points. Preserve the earlier roadmap as historical context.
- Generate paper tables, demo metrics, website figures, and narrative evidence summaries from the same versioned results. Remove hard-coded headline values from `demo/demo_utils.py`, `demo/app.py`, and static exports.
- Keep the live model and nested LOSO evaluation artifacts separately identified; do not attach ensemble evaluation numbers to a different deployed classifier.
- Add portable commands and an environment specification for reproducing reports from saved outputs, rerunning the benchmark, and running deterministic or generative audits. Record data sources, preprocessing versions, and artifact checksums.
- Extend targeted tests for study/donor separation, train-only transformations, nested marker selection, metric aggregation, narrative audit failures, and artifact consistency. Keep network-dependent LLM experiments outside routine CI; use saved responses and mocked requests for CI.
- Add a reviewer-response matrix mapping each concern to the experiment, result artifact, code change, and manuscript section. Update poster materials if claims change, preserving the previously published version.

Completion evidence: a fresh checkout reproduces published tables, targeted checks pass, public claims match saved evidence, and each reviewer concern has a direct artifact link.

## Delivery order

1. **Protocol and reporting correction:** establish the primary eight-study result and freeze comparison rules.
2. **Audit correctness:** prevent false LLM runs and test the checker before collecting final responses.
3. **Matched marker benchmark and generative audit:** produce the two missing evidence packages under frozen protocols.
4. **Documentation and release:** integrate results into paper, demo, website, README, and reviewer-response matrix.

The camera-ready minimum is the honest full-cohort report, actual LLM audit, and matched marker comparison. Broader model searches, new modalities, platform adaptation, and a larger encoder-decoder redesign are follow-up research after these obligations are met.
