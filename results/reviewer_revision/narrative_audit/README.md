# Narrative audit status

The revised deterministic audit passed automated screening on **101/101 profiles**, including **21 abstentions**. This is template validation, not evidence that LLM-generated narratives are faithful.

**Actual LLM generation is pending:** no model endpoint was configured during this revision. Independent human review is also pending. No generative or human-reviewed pass rates are claimed.

The frozen packets identify the live Stage-B model's outputs; the historical LOSO ensemble is a separate evidence object. These public development profiles are not a new held-out clinical evaluation.

Run with a configured OpenAI-compatible chat-completions endpoint (configure credentials outside the repository):

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

For a blinded template-versus-LLM usefulness comparison, mix both arms under randomly assigned item IDs and remove arm labels before rating. The per-run forms alone do not provide cross-arm blinding. This comparison and human adjudication remain pending until actual LLM outputs and raters are available.
