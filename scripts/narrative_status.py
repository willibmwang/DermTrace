"""Read narrative validation status from the published audit artifact."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / 'results/reviewer_revision/narrative_audit/qwen3_1_7b_v1/summary.json'


def narrative_status():
    if not SUMMARY.exists():
        return ('Actual LLM audit is pending. Templates pass automated screening on 101/101 profiles; '
                'independent human ratings remain pending.')
    result = json.loads(SUMMARY.read_text())
    return (f"Qwen3-1.7B: {result['generation_successes']}/{result['n_attempts']} generations completed; "
            f"{result['all_checks_pass']}/{result['n_attempts']} raw outputs passed all automated checks "
            f"across {result['n_profiles']} profiles ({result['repeats']} greedy requests per profile). "
            'Automated screening is not human-validated faithfulness. Independent human ratings remain pending.')
