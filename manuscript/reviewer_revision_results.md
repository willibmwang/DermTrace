# Camera-ready revision results

The primary historical ensemble result includes **8 studies / 101 profiles**: **0.826 mean macro F1**, with **0.182 on the worst study**. The **0.918** seven-study result is an exploratory subset selected after observing GSE173900's failure.

Under the new matched protocol, COL1A1 alone scores **0.914** and the fixed rank-program view **0.896**. Inner-selected programs versus markers score **0.855 versus 0.829** (paired difference +0.025; descriptive study-bootstrap interval [-0.108, +0.148]). These results do **not** establish general superiority of program features.

The revised template audit passes automated screening on **101/101 profiles**, including **21 abstentions**. Qwen3-1.7B: 303/303 generations completed; 0/303 raw outputs passed all automated checks across 101 profiles (3 greedy requests per profile). In a post hoc diagnostic, removing echoed prompt labels raised this only to 6/303: 297/303 outputs restated numbers despite the instruction not to, and 609/933 numbers attached to a named program matched its contribution (255 reported its input value instead; 69 matched neither). Automated screening is not human-validated faithfulness. Independent human ratings remain pending.

## Evaluation protocol and scope

The primary endpoint is keloid versus unaffected skin. All eight eligible development studies remain in the historical ensemble evaluation. The primary metric gives equal weight to studies; within each study it averages probabilities by donor when multiple class-pure donors exist, and uses sample-level scoring for single-donor or mixed-class donor designs. Sample-level metrics and counts are reported alongside the primary metric. The existing mean-0.85/worst-0.75 development gates fail.

The GSE173900 exclusion was made after its outer-fold failure. It is not independent of test performance and cannot establish prospective platform eligibility. Oracle grid results are diagnostic and do not estimate achievable future performance. A future eligibility rule must use information available without held-out outcomes and undergo independent validation.

The matched revision benchmark uses the same outer test studies, fixed feature definitions, study/donor weighting, training-only median imputation and scaling, L2 logistic regression with C in {0.1, 1, 10}, inner-LOSO configuration selection, inner-OOF Platt calibration, and inner-selected thresholds for every arm. All other development samples are used for each outer training fit. This removes the historical train/validation reservation consistently across the new arms; the historical ensemble is a separate comparison and is not a matched arm.

Each fixed marker and each program view is evaluated separately. The inner-selected marker and program families choose their representation using outer-training data only. Selecting COL1A1 after inspecting the table is descriptive, not a new independent validation. Fixed source matrices include historical preprocessing; prospective raw-data generalization is not established.

## Matched results

| Arm | Mean macro F1 | Worst study |
|---|---:|---:|
| `inner_selected_marker` | 0.8295 | 0.1818 |
| `inner_selected_program` | 0.8548 | 0.4156 |
| `marker_ACTA2` | 0.3558 | 0.1429 |
| `marker_COL1A1` | 0.9144 | 0.4444 |
| `marker_COL3A1` | 0.8608 | 0.1818 |
| `marker_CTHRC1` | 0.8087 | 0.1818 |
| `marker_FN1` | 0.8181 | 0.3250 |
| `marker_POSTN` | 0.8260 | 0.3250 |
| `program_composition_only` | 0.7503 | 0.3250 |
| `program_fibrosis_only` | 0.8023 | 0.1818 |
| `program_fused_multiview` | 0.8356 | 0.3250 |
| `program_rank_programs_only` | 0.8965 | 0.4156 |
| `program_scar_discriminative` | 0.8419 | 0.4156 |

The inner-selected paired study-bootstrap interval is [-0.1081, +0.1480], based on only eight studies. Repeated development and overlapping training sets limit inference. No general superiority claim is supported.

## Narrative validation

The revised harness freezes the derived evidence packets, requires actual generation in LLM mode, and records every attempt, prompt/packet hash, returned model identifier, failure, and truncation flag. It preserves raw responses before runtime validation. The model is asked to retain an exact factual opening, the QC note, and a research-use disclaimer while generating the explanatory prose. Automated checks screen numeric field association, obvious class contradictions, abstention, contribution signs, prohibited claims, request redaction, length, and packet immutability. They are not a substitute for semantic human assessment.

The deterministic arm passed automated screening for all 101 profiles, including 21 abstentions. Qwen3-1.7B: 303/303 generations completed; 0/303 raw outputs passed all automated checks across 101 profiles (3 greedy requests per profile). In a post hoc diagnostic, removing echoed prompt labels raised this only to 6/303: 297/303 outputs restated numbers despite the instruction not to, and 609/933 numbers attached to a named program matched its contribution (255 reported its input value instead; 69 matched neither). Automated screening is not human-validated faithfulness. Independent human ratings remain pending. Blinded usefulness comparisons and human adjudication require independent raters.

The frozen local run is the primary generative result and is reported as a failure of the constrained-generation design with this model, not as a faithfulness rate. The post hoc diagnostic (`scripts/diagnose_narrative_audit.py`) was specified after inspecting outputs and does not change the frozen validator. It shows that echoed prompt labels explain the uniform anchor failure, but content failures remain: unrequested numbers, program input values described as contributions, and numbers matching no supplied field. The frozen `contribution_direction` check passed every output because it detects only explicit toward/against-keloid phrasing; it is therefore too weak to establish contribution faithfulness. After label normalization, the remaining decision-consistency flags are conservative restatements of the supplied class rather than observed class reversals. The demo shows a labeled template whenever raw generation fails validation, which would apply to every frozen-run output. Any revised prompt or stronger model requires a new frozen protocol version, and results on these same development packets would be development-informed.

## Provenance and limitations

GSE185309 was scored in an earlier experiment, as recorded in `results/breakthrough_sprint/breakthrough_lockbox_once.json` and the v2 frozen protocol. It was not rescored in this revision and must not be described as untouched. No additional reserved cohort was scored.

The live Stage-B scorer, historical nested ensemble, and matched revision classifiers are distinct evidence objects. The static website uses deterministic summaries. The preserved poster may contain older claims; current tables and this text supersede those claims.

Sources: `results/reviewer_revision/evaluation.json`, `matched_results.json`, `protocols/reviewer_revision_v1.json`, and `results/reviewer_revision/narrative_audit/`.
