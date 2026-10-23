"""Reproduce the primary all-study report with explicit keloid=0 encoding."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score, f1_score, brier_score_loss
from score_primary_ladder import score_from_predictions, bootstrap_ci

ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = ROOT.parent
OUT = ROOT / 'results/reviewer_revision'
RUN = ROOT / 'results/accuracy_085_push/expanded_ensemble_v1'
ENDPOINT = 'keloid_vs_unaffected_skin'
MARKERS = ['POSTN', 'COL1A1', 'COL3A1', 'ACTA2', 'FN1', 'CTHRC1']
PROGRAMS = ['fibrosis_only', 'rank_programs_only', 'composition_only', 'scar_discriminative', 'fused_multiview']

def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def summarize(predictions, manifest):
    folds = score_from_predictions(predictions, manifest)
    for row in folds:
        sub = predictions[predictions.accession.eq(row['accession'])]
        pred = np.where(sub.prob_keloid >= sub.threshold, 0, 1)
        row.update(n_keloid=int(sub.y_true.eq(0).sum()), n_unaffected=int(sub.y_true.eq(1).sum()),
                   sample_accuracy=float(accuracy_score(sub.y_true, pred)),
                   brier_score=float(brier_score_loss(sub.y_true.eq(0), sub.prob_keloid)))
    pred = np.where(predictions.prob_keloid >= predictions.threshold, 0, 1)
    return {'n_studies': len(folds), 'n_profiles': len(predictions),
            'n_donors': int(predictions[['accession','donor']].drop_duplicates().shape[0]),
            'mean_macro_f1': float(np.mean([r['macro_f1'] for r in folds])),
            'worst_macro_f1': float(min(r['macro_f1'] for r in folds)),
            'mean_sample_macro_f1': float(np.mean([r['sample_macro_f1'] for r in folds])),
            'pooled_accuracy': float(accuracy_score(predictions.y_true, pred)),
            'pooled_correct': int((predictions.y_true == pred).sum()),
            'pooled_sample_macro_f1': float(f1_score(predictions.y_true, pred, average='macro')),
            'study_bootstrap_ci': bootstrap_ci([r['macro_f1'] for r in folds]), 'folds': folds}

def selective_report(predictions, threshold=0.55):
    # Operational probability confidence; no tuning against held-out outcomes.
    records = []
    for acc, sub in predictions.groupby('accession'):
        keep = np.maximum(sub.prob_keloid, 1-sub.prob_keloid) >= threshold
        kept = sub[keep]
        evaluable = kept.y_true.nunique() == 2
        records.append({'accession': acc, 'n_total': len(sub), 'n_retained': len(kept),
                        'coverage': float(keep.mean()),
                        'keloid_coverage': float(keep[sub.y_true.eq(0)].mean()),
                        'unaffected_coverage': float(keep[sub.y_true.eq(1)].mean()),
                        'evaluable_both_classes': evaluable,
                        'sample_macro_f1': float(f1_score(kept.y_true, np.where(kept.prob_keloid >= kept.threshold,0,1), average='macro')) if evaluable else None})
    valid = [r['sample_macro_f1'] for r in records if r['sample_macro_f1'] is not None]
    return {'confidence_threshold': threshold, 'metric_grain': 'sample; both true classes required',
            'mean_coverage': float(np.mean([r['coverage'] for r in records])),
            'n_evaluable_folds': len(valid), 'mean_sample_macro_f1_evaluable': float(np.mean(valid)) if valid else None,
            'folds': records}

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    p = RUN / f'predictions_{ENDPOINT}.csv'
    predictions = pd.read_csv(p)
    manifest = pd.read_parquet(REPO_ROOT / 'data/processed/training/profile_manifest.parquet')
    labels = manifest.set_index('sample_id').loc[predictions.sample_id, ENDPOINT].map({'keloid':0,'non_keloid':1}).to_numpy()
    np.testing.assert_array_equal(labels, predictions.y_true)
    primary = summarize(predictions, manifest)
    historical = json.loads((ROOT / 'results/accuracy_085_push/expanded_ensemble_v1.json').read_text())
    np.testing.assert_allclose(primary['mean_macro_f1'], historical['mean_macro_f1'], atol=1e-12)
    secondary = summarize(predictions[~predictions.accession.eq('GSE173900')], manifest)
    result = {'status': 'recomputed from saved predictions; no new model training',
              'positive_label': 'keloid', 'positive_code': 0,
              'metric': 'equal-study macro F1; donor means for class-pure multi-donor studies; sample fallback otherwise',
              'primary': primary, 'post_hoc_seven_study_sensitivity': secondary,
              'primary_selective': selective_report(predictions),
              'source_sha256': {str(p.relative_to(ROOT)): sha256(p)},
              'limitations': ['Seven-study exclusion followed observation of test failure.',
                             'Study bootstrap has only eight clusters; estimates are descriptive.',
                             'Historical development data have been inspected repeatedly; this is not independent confirmation.',
                             'Live Stage-B model is distinct from this nested ensemble.']}
    (OUT/'evaluation.json').write_text(json.dumps(result,indent=2)+'\n')
    pd.DataFrame(primary['folds']).to_csv(OUT/'all_studies.csv',index=False)
    pd.DataFrame(result['primary_selective']['folds']).to_csv(OUT/'selective_by_study.csv',index=False)
    lines = ['# Reproduced eight-study evaluation', '', 'Keloid is code **0** in the saved predictions. Full coverage is primary.', '',
             '| Analysis | Studies | Profiles | Mean macro F1 | Worst fold |', '|---|---:|---:|---:|---:|']
    for name, item in [('Primary, all studies',primary),('Post hoc sensitivity',secondary)]:
        lines.append(f"| {name} | {item['n_studies']} | {item['n_profiles']} | {item['mean_macro_f1']:.4f} | {item['worst_macro_f1']:.4f} |")
    lines += ['', 'The primary result fails the existing 0.85 mean / 0.75 worst-fold gates. The subset is outcome-selected.', '',
              'The historical prediction cross-check reversed class encoding. That cross-check is now corrected; the main historical score, calculated from the training metrics, is unchanged.', '',
              'See `all_studies.csv` for per-study/sample metrics and `selective_by_study.csv` for coverage by class. Selective metrics use sample grain and are separate from the donor/sample primary metric.', '',
              'Uncertainty is descriptive: eight studies and repeated development do not establish prospective generalization.']
    (OUT/'evaluation.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({k:primary[k] for k in ('n_studies','n_profiles','mean_macro_f1','worst_macro_f1')}))
if __name__ == '__main__':
    main()
