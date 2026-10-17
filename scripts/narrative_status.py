"""Read narrative validation status from the published audit artifact."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / 'results/reviewer_revision/narrative_audit/qwen3_1_7b_v1/summary.json'
DIAGNOSTICS = ROOT / 'results/reviewer_revision/narrative_audit/qwen3_1_7b_v1_diagnostics/summary.json'
V2 = ROOT / 'results/reviewer_revision/narrative_audit/qwen3_8b_v2/summary.json'
HUMAN = ('Independent human rating was not performed in this revision; automated screening is not '
         'human-validated faithfulness, and blank two-rater forms are released with each run.')


def narrative_status():
    if not SUMMARY.exists():
        return ('Actual LLM audit is pending. Templates pass automated screening on 101/101 profiles; '
                'independent human ratings remain pending.')
    result = json.loads(SUMMARY.read_text())
    text = (f"Frozen v1 (Qwen3-1.7B writes the whole narrative): {result['generation_successes']}/{result['n_attempts']} generations completed; "
            f"{result['all_checks_pass']}/{result['n_attempts']} raw outputs passed all automated checks "
            f"across {result['n_profiles']} profiles ({result['repeats']} greedy requests per profile).")
    if DIAGNOSTICS.exists():
        d = json.loads(DIAGNOSTICS.read_text())
        norm, nums = d['label_normalized']['all_checks_pass'], d['numbers']
        attributed = nums['program_attributed_numbers']
        text += (f" In a post hoc diagnostic, removing echoed prompt labels raised this only to {norm['n']}/{norm['denominator']}: "
                 f"{nums['attempts_with_extra_numbers']['n']}/{d['n_attempts']} outputs restated numbers despite the instruction not to, and "
                 f"{attributed['matches_program_contribution']}/{attributed['n']} numbers attached to a named program matched its contribution "
                 f"({attributed['matches_program_input_value_only']} reported its input value instead; {attributed['matches_neither']} matched neither).")
    if V2.exists():
        v = json.loads(V2.read_text())
        lo, hi = v['profile_cluster_bootstrap_95_interval']
        rows = [json.loads(line) for line in (V2.parent / 'responses.jsonl').read_text().splitlines()]
        mixed = sum(not (row['raw_generation'] or '').isascii() for row in rows)
        clean = sum(row['all_checks_pass'] and (row['raw_generation'] or '').isascii() for row in rows)
        text += (f" Development-informed v2 (Qwen3-8B writes only the contribution prose; code adds the fixed opening, QC note, and disclaimer): "
                 f"{v['generation_successes']}/{v['n_attempts']} generations completed and {v['all_checks_pass']}/{v['n_attempts']} assembled outputs "
                 f"passed all automated checks, including a stricter clause-level contribution-role check "
                 f"({v['profiles_all_repeats_pass']}/{v['n_profiles']} profiles on every repeat; profile-cluster interval [{lo:.3f}, {hi:.3f}]). "
                 f"Post hoc, {mixed}/{v['n_attempts']} outputs contained a non-English token (Qwen code-switching) that no check covered, "
                 f"so {clean}/{v['n_attempts']} were both check-passing and English-only; the demo now falls back to the template for such text. "
                 "v2 was designed after inspecting v1 on the same packets and is not an independent validation.")
    return text + ' ' + HUMAN
