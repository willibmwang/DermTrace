# Narrative audit status

The revised deterministic audit passed automated screening on **101/101 profiles**, including **21 abstentions**. This is template validation, not evidence that LLM-generated narratives are faithful.

## Frozen local LLM run (primary generative result)

Protocol [`narrative_local_v1`](../../../protocols/narrative_local_v1.json): pinned Qwen3-1.7B snapshot, greedy decoding, unchanged prompt, 101 frozen packets x 3 requests, SLURM job 2865612 on one A40 ([log](../../slurm/keloid-narrative-audit-2865612.out)). Outputs: [`qwen3_1_7b_v1/`](qwen3_1_7b_v1/).

- Generation succeeded for **303/303** requests; none truncated.
- **0/303** raw outputs passed all automated checks (profile-cluster bootstrap interval [0.0, 0.0]). No profile passed on any repeat.
- Every output failed `factual_anchor`, and therefore the dependent `numeric_fields_grounded` and `decision_consistent` checks. `research_use_stated` passed 156/303; `abstention_disclosed` passed 240/303 (all failures were the 63 abstention attempts). All other checks passed for every attempt.
- Greedy repeats were identical for 97/101 profiles, so repeats mostly measure determinism rather than independent variation.
- With this model, the demo's runtime validation would display the labeled template for every frozen-run output.

The frozen result is reported as is. The prompt was not tuned on these outputs.

## Post hoc failure diagnostic (not a pass rate)

[`scripts/diagnose_narrative_audit.py`](../../../scripts/diagnose_narrative_audit.py) was written after inspecting the outputs. It leaves the run and validator unchanged; outputs are in [`qwen3_1_7b_v1_diagnostics/`](qwen3_1_7b_v1_diagnostics/).

- **Formatting:** all 303 outputs echoed the prompt's `OPENING:` / `QC NOTE:` / `DISCLAIMER:` labels, so none began with the exact anchor sentence, although the anchor and QC note appeared verbatim in every output. 171 added an `END OF OUTPUT/RESPONSE` line, and 156 repeated the disclaimer.
- **Label-normalized rescoring:** removing only those labels and end markers yields **6/303** passes (cluster interval [0.000, 0.050]). Remaining failures are extra numbers (297), a disclaimer that is not final (147), and conservative decision flags (21). Manual inspection found that every decision flag restated the supplied class; none reversed it.
- **Numeric attribution (heuristic):** 297/303 outputs restated numbers despite the instruction not to. Of 933 numbers attached to a named program, 609 matched that program's contribution, 255 matched its input value while being described as a contribution, and 69 matched neither (for example, invented percentage shares). 120/303 outputs contained at least one such number.
- **Checker limitation:** the frozen `contribution_direction` check passed every output because it detects only explicit toward/against-keloid phrasing. It missed sign/value confusions such as "highest positive contribution (-0.156)". Contribution faithfulness therefore requires human rating or a stronger, separately validated checker.

**Conclusion (v1).** Constrained whole-narrative generation with this small local model is not faithful enough for display under the frozen criteria. This motivated the separately frozen v2 below.

## v2 local LLM run (current generative result; development-informed)

Protocol [`narrative_local_v2`](../../../protocols/narrative_local_v2.json) was frozen and pushed before generation (commit 2180d0d). Pinned Qwen3-8B, greedy decoding, 101 frozen packets x 3 requests, SLURM job 2866086 on one A40 ([log](../../slurm/keloid-narrative-v2-2866086.out)). Outputs: [`qwen3_8b_v2/`](qwen3_8b_v2/), including `raw_generation` for every attempt.

The model writes only the contribution prose. Code adds the fixed factual opening, transfer note, and disclaimer, so no number or class decision is model-generated. The live demo uses the same format.

- Generation succeeded for **303/303** requests, with no truncation.
- **290/303** assembled outputs passed all automated checks (profile-cluster interval [0.917, 0.990]). **96/101** profiles passed on every repeat.
- By decision: abstain 57/63, keloid 128/135, non-keloid 105/105.
- All 13 failures were `contribution_role` failures, in which a program was described in the wrong direction. Three also failed `contribution_direction`. Manual inspection confirmed that these were genuine errors.
- Greedy repeats were identical for 100/101 profiles.
- **Post hoc observation:** 75/303 outputs contained a Chinese token (`倾向`, code-switching), and all of them passed the frozen checks. Only **215/303** were both check-passing and English-only. The demo now also requires English-only prose before display.

**Caveat.** The format, model, and new check were chosen after inspecting v1 on the same development packets. v2 therefore shows that the constrained design is workable, not that it generalizes. The automated checks cover anchors, numbers, class statements, directions, prohibited claims, redaction, and immutability; they do not assess usefulness.

Independent human rating was not performed in this revision. Blank two-rater forms are released in each run directory.

The frozen packets identify the live Stage-B model's outputs; the historical LOSO ensemble is a separate evidence object. These public development profiles are not a new held-out clinical evaluation.

Reproduce the frozen local run on a GPU node (pinned snapshot from the protocol), then replay and diagnose it without a GPU:

```bash
sbatch scripts/run_narrative_audit_gpu.sbatch
python scripts/replay_narrative_audit.py \
  results/reviewer_revision/narrative_audit/qwen3_1_7b_v1/responses.jsonl
python scripts/diagnose_narrative_audit.py \
  results/reviewer_revision/narrative_audit/qwen3_1_7b_v1 --output /tmp/narrative_diagnostics
```

Alternatively, run with a configured OpenAI-compatible chat-completions endpoint (configure credentials outside the repository):

```bash
python scripts/audit_copilot_faithfulness.py --arm llm --repeats 3 \
  --output results/reviewer_revision/narrative_audit/llm_v1
python scripts/replay_narrative_audit.py \
  results/reviewer_revision/narrative_audit/llm_v1/responses.jsonl
```

Set `KELOIDBENCH_LLM_API_URL`, `KELOIDBENCH_LLM_MODEL`, and, when required, `KELOIDBENCH_LLM_API_KEY`. LLM mode refuses an unconfigured endpoint and never silently counts templates as generation. Output directories cannot be overwritten. Freeze a new protocol/run version before modifying prompts or evaluation packets. All raw attempts, API errors, truncation indicators, model identifiers, prompts, and packet hashes remain in the output. API credentials are not saved.

The prompt requires an exact labeled factual opening, the supplied QC note, and the research-use disclaimer. The LLM generates the explanatory prose between these anchors. Automated checks test these constraints, extra numbers, obvious contradictory class assertions, prohibited claims, contribution signs, packet immutability, and the outgoing payload. Passing these conservative checks is not a semantic guarantee. The demo separately screens raw prose and uses an explicitly labeled template fallback after generation/validation failure.

## Human review

Each run exports two blank rater forms and a separate key. Use independent raters who did not develop the prompt. Keep the key hidden until rating is complete. Rate the preselected study/decision-stratified sample; review every automatically flagged case separately, without folding that enriched sample into the random-sample pass rate.

For each item, compare the narrative with its redacted evidence packet:

- `factual_faithfulness`: pass only if every factual statement is supported; otherwise fail and quote the unsupported statement.
- `decision_preserved`, `abstention_preserved`, `qc_preserved`: pass/fail; use not-applicable for abstention on a non-abstained packet.
- `clarity_1_to_5`: 1 means confusing; 3 means understandable with effort; 5 means clear and concise.
- `completeness_1_to_5`: 1 omits key information; 3 conveys the decision but misses a limitation or explanation; 5 conveys the decision, uncertainty, QC limitation, and available contribution evidence.
- `error_notes`: classify numeric, class, abstention, contribution-direction, unsupported-mechanism, omission, or other error.

Report each rater's counts, raw agreement and chance-adjusted agreement, adjudicated outcomes, and denominators. Report request success separately from faithfulness among successful generations. Count API failures in the all-attempt success measure. Repeated outputs share a profile: use profile-cluster bootstrap intervals, and report study/decision breakdowns. Do not count repeated requests as independent profiles.

For a blinded template-versus-LLM usefulness comparison, mix both arms under randomly assigned item IDs and remove arm labels before rating. The per-run forms alone do not provide cross-arm blinding. This comparison and human adjudication were not performed in this revision.
