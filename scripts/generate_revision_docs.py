"""Generate current paper/README results from saved evaluation artifacts."""
import json
from pathlib import Path
import pandas as pd
from narrative_status import narrative_status
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/reviewer_revision'

def main():
    e=json.loads((OUT/'evaluation.json').read_text());p=e['primary'];s=e['post_hoc_seven_study_sensitivity']
    m=json.loads((OUT/'matched_results.json').read_text());arms=m['arms'];paired=m['paired_comparison'];lo,hi=paired['study_bootstrap_95_percentile_interval']
    rows=pd.read_csv(OUT/'matched_summary.csv')
    block=(f"The primary historical ensemble result includes **{p['n_studies']} studies / {p['n_profiles']} profiles**: "
        f"**{p['mean_macro_f1']:.3f} mean macro F1**, with **{p['worst_macro_f1']:.3f} on the worst study**. "
        f"The **{s['mean_macro_f1']:.3f}** seven-study result is an exploratory subset selected after observing GSE173900's failure.\n\n"
        f"Under the new matched protocol, COL1A1 alone scores **{arms['marker_COL1A1']['mean_macro_f1']:.3f}** and the fixed rank-program view **{arms['program_rank_programs_only']['mean_macro_f1']:.3f}**. "
        f"Inner-selected programs versus markers score **{arms['inner_selected_program']['mean_macro_f1']:.3f} versus {arms['inner_selected_marker']['mean_macro_f1']:.3f}** "
        f"(paired difference {paired['mean_delta']:+.3f}; descriptive study-bootstrap interval [{lo:+.3f}, {hi:+.3f}]). "
        "These results do **not** establish general superiority of program features.\n\n"
        "The revised template audit passes automated screening on **101/101 profiles**, including **21 abstentions**. "
        + narrative_status())
    readme=ROOT/'README.md';text=readme.read_text();a=text.index('<!-- REVISION_RESULTS_START -->')+len('<!-- REVISION_RESULTS_START -->');b=text.index('<!-- REVISION_RESULTS_END -->')
    readme.write_text(text[:a]+'\n\n'+block+'\n\n'+text[b:])
    prose=['# Camera-ready revision results','',block,'',
        '## Evaluation protocol and scope','',
        'The primary endpoint is keloid versus unaffected skin. All eight eligible development studies remain in the historical ensemble evaluation. The primary metric gives equal weight to studies; within each study it averages probabilities by donor when multiple class-pure donors exist, and uses sample-level scoring for single-donor or mixed-class donor designs. Sample-level metrics and counts are reported alongside the primary metric. The existing mean-0.85/worst-0.75 development gates fail.','',
        'The GSE173900 exclusion was made after its outer-fold failure. It is not independent of test performance and cannot establish prospective platform eligibility. Oracle grid results are diagnostic and do not estimate achievable future performance. A future eligibility rule must use information available without held-out outcomes and undergo independent validation.','',
        'The matched revision benchmark uses the same outer test studies, fixed feature definitions, study/donor weighting, training-only median imputation and scaling, L2 logistic regression with C in {0.1, 1, 10}, inner-LOSO configuration selection, inner-OOF Platt calibration, and inner-selected thresholds for every arm. All other development samples are used for each outer training fit. This removes the historical train/validation reservation consistently across the new arms; the historical ensemble is a separate comparison and is not a matched arm.','',
        'Each fixed marker and each program view is evaluated separately. The inner-selected marker and program families choose their representation using outer-training data only. Selecting COL1A1 after inspecting the table is descriptive, not a new independent validation. Fixed source matrices include historical preprocessing; prospective raw-data generalization is not established.','',
        '## Matched results','',
        '| Arm | Mean macro F1 | Worst study |','|---|---:|---:|']
    for r in rows.itertuples():prose.append(f'| `{r.arm}` | {r.mean_macro_f1:.4f} | {r.worst_macro_f1:.4f} |')
    prose += ['',f'The inner-selected paired study-bootstrap interval is [{lo:+.4f}, {hi:+.4f}], based on only eight studies. Repeated development and overlapping training sets limit inference. No general superiority claim is supported.','',
        '## Narrative validation','',
        'The revised harness freezes the derived evidence packets, requires actual generation in LLM mode, and records every attempt, prompt/packet hash, returned model identifier, failure, and truncation flag. It preserves raw responses before runtime validation. The model writes only the explanatory prose about which programs pushed the score toward keloid or non-keloid; code inserts the factual opening, QC note, and research-use disclaimer. Automated checks screen numeric field association, obvious class contradictions, abstention, contribution signs, prohibited claims, request redaction, length, and packet immutability. They are not a substitute for semantic human assessment.','',
        'The deterministic arm passed automated screening for all 101 profiles, including 21 abstentions. ' + narrative_status() + ' Blinded usefulness comparisons and human adjudication require independent raters.','',
        '## Provenance and limitations','',
        'GSE185309 was scored in an earlier experiment, as recorded in `results/breakthrough_sprint/breakthrough_lockbox_once.json` and the v2 frozen protocol. It was not rescored in this revision and must not be described as untouched. No additional reserved cohort was scored.','',
        'The live Stage-B scorer, historical nested ensemble, and matched revision classifiers are distinct evidence objects. The static website uses deterministic summaries. The preserved poster may contain older claims; current tables and this text supersede those claims.','',
        'Sources: `results/reviewer_revision/evaluation.json`, `matched_results.json`, `protocols/reviewer_revision_v1.json`, and `results/reviewer_revision/narrative_audit/`.']
    (ROOT/'manuscript/reviewer_revision_results.md').write_text('\n'.join(prose)+'\n')
    tex=['% Generated by scripts/generate_revision_docs.py; requires booktabs.',r'\subsection*{Matched single-marker comparison}',
         'All eight development studies are evaluated under a common revision-time nested protocol. The historical ensemble is a separate evidence object.',r'\begin{tabular}{lrr}',r'\toprule',r'Arm & Mean macro F1 & Worst study \\',r'\midrule']
    for r in rows.itertuples():tex.append(r.arm.replace('_',r'\_')+f' & {r.mean_macro_f1:.4f} & {r.worst_macro_f1:.4f} '+r'\\')
    tex += [r'\bottomrule',r'\end{tabular}',f'Inner-selected programs minus markers: {paired["mean_delta"]:+.4f}; descriptive study-bootstrap interval [{lo:+.4f}, {hi:+.4f}]. These results do not establish general program superiority.']
    (ROOT/'manuscript/reviewer_revision_results.tex').write_text('\n'.join(tex)+'\n')
if __name__=='__main__':main()
