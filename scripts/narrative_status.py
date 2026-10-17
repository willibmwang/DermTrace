"""Read narrative validation status from the published audit artifact."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / 'results/reviewer_revision/narrative_audit/qwen3_1_7b_v1/summary.json'
DIAGNOSTICS = ROOT / 'results/reviewer_revision/narrative_audit/qwen3_1_7b_v1_diagnostics/summary.json'


def narrative_status():
    if not SUMMARY.exists():
        return ('Actual LLM audit is pending. Templates pass automated screening on 101/101 profiles; '
                'independent human ratings remain pending.')
    result = json.loads(SUMMARY.read_text())
    text = (f"Qwen3-1.7B: {result['generation_successes']}/{result['n_attempts']} generations completed; "
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
    return text + (' Automated screening is not human-validated faithfulness. Independent human ratings remain pending.')
