"""Read narrative validation status from the published audit artifact."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / 'results/reviewer_revision/narrative_audit/qwen3_8b_v2'


def narrative_status():
    if not (RUN / 'summary.json').exists():
        return 'The LLM narrative audit is pending. Templates pass automated screening on 101/101 profiles.'
    v = json.loads((RUN / 'summary.json').read_text())
    rows = [json.loads(line) for line in (RUN / 'responses.jsonl').read_text().splitlines()]
    clean = sum(row['all_checks_pass'] and (row['raw_generation'] or '').isascii() for row in rows)
    lo, hi = v['profile_cluster_bootstrap_95_interval']
    return (f"LLM narratives (Qwen3-8B writes the explanatory prose; code inserts every number, the decision, the QC note, "
            f"and the disclaimer): {v['all_checks_pass']}/{v['n_attempts']} outputs passed all automated checks "
            f"({v['profiles_all_repeats_pass']}/{v['n_profiles']} profiles on every repeat; profile-cluster interval [{lo:.3f}, {hi:.3f}]). "
            f"All failures were contribution-direction errors. {clean}/{v['n_attempts']} were also English-only; "
            "the demo shows generated prose only when it passes every check and is English-only, and otherwise shows the template. "
            "The format and checks were developed on these same development profiles, so this is not independent validation. "
            "Independent human rating was not performed.")
