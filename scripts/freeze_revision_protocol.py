"""Freeze revision-time inputs and matched comparison rules before training."""
import argparse
import json
from pathlib import Path
import pandas as pd
import pyarrow.parquet as pq
from revision_report import ROOT, OUT, RUN, ENDPOINT, MARKERS, PROGRAMS, sha256

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source-root', type=Path, default=ROOT)
args = parser.parse_args()
protocol_path = ROOT/'protocols/reviewer_revision_v1.json'
if protocol_path.exists():
    raise SystemExit('Protocol already frozen; create a new version rather than overwrite it.')
preds = pd.read_csv(RUN/f'predictions_{ENDPOINT}.csv')
manifest = pd.read_parquet(ROOT/'data/processed/training/profile_manifest.parquet').set_index('sample_id')
meta = manifest.loc[preds.sample_id].copy()
assert not meta.lockbox.any()
rows = []
sources = {}
for rel, group in meta.groupby('artifact_wide_path'):
    path = args.source_root / rel
    schema = pq.read_schema(path).names
    present = [g for g in MARKERS if g in schema]
    table = pd.read_parquet(path, columns=['sample_id', *present]).set_index('sample_id')
    part = table.reindex(index=group.index, columns=MARKERS)
    rows.append(part)
    sources[rel] = sha256(path)
markers = pd.concat(rows).loc[preds.sample_id]
markers.index.name = 'sample_id'
OUT.mkdir(parents=True, exist_ok=True)
markers.to_csv(OUT/'marker_inputs.csv')
meta[['accession','split_group',ENDPOINT,'platform_id','artifact_wide_path']].to_csv(OUT/'cohort_manifest.csv')
program = pd.read_parquet(ROOT/'data/processed/training/features/fused_multiview.parquet').set_index('sample_id')
assert set(meta.index) <= set(program.index)
# All non-held-out development profiles are used symmetrically in the new comparison.
# Historical folds reserved validation rows but the old inner selector used them;
# this explicit protocol removes that asymmetry for every matched arm.
splits = [{'accession':a, 'train_sample_ids':meta.index[~meta.accession.eq(a)].tolist(),
           'test_sample_ids':meta.index[meta.accession.eq(a)].tolist()} for a in sorted(meta.accession.unique())]
inputs = [OUT/'marker_inputs.csv', OUT/'cohort_manifest.csv'] + [ROOT/f'data/processed/training/features/{v}.parquet' for v in PROGRAMS]
protocol = {'version':'reviewer_revision_v1','specified_at':'2026-10-17',
 'status':'revision-time specification; development outcomes already known; not preregistration',
 'endpoint':ENDPOINT,'positive_code':0,'positive_label':'keloid','seed':13,
 'primary_studies':sorted(meta.accession.unique()), 'post_hoc_exclusion':['GSE173900'],
 'markers':MARKERS,'program_views':PROGRAMS,'splits':splits,
 'models':{'family':'L2 logistic regression','C':[0.1,1.0,10.0], 'max_iter':4000},
 'preprocessing':'fixed source expression inputs; median imputation and standard scaling fit on each training partition only',
 'source_preprocessing_limitation':'Source matrices may include cohort-level unsupervised harmonization; these are frozen preprocessed-input comparisons, not raw-data prospective validation.',
 'missingness':'preserve missing values; training median; all-missing training columns use zero; publish per-study coverage',
 'weighting':'equal study, then equal donor within study, then equal sample within donor; mean-one weights',
 'selection':'inner LOSO mean primary macro F1; lower-tail tie break; same C grid for each fixed arm',
 'calibration':'Platt fit on pooled inner-OOF probabilities; fixed for all arms; outer labels never used',
 'thresholds':[0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70],
 'threshold_selection':'inner-OOF equal-study primary macro F1; prefer 0.5 on ties',
 'metric':'donor probability means for multi-donor class-pure folds; sample fallback for single-donor or mixed-class designs',
 'comparison':'Each marker and each fixed program view uses identical new protocol. Best marker and best program view selected within outer training only. Historical ensemble reported separately, not as a matched arm.',
 'historical_protocol_difference':'New outer training uses all non-held-out development samples; historical train/validation reservation is removed consistently across arms.',
 'lockbox':'excluded; no new lockbox scoring; prior access history not certified',
 'uncertainty':'paired study bootstrap, descriptive with eight clusters; no superiority claim from point estimates alone',
 'source_sha256':sources,'input_sha256':{str(p.relative_to(ROOT)):sha256(p) for p in inputs}}
protocol_path.parent.mkdir(exist_ok=True)
protocol_path.write_text(json.dumps(protocol,indent=2)+'\n')
print('Frozen',len(meta),'profiles;',len(splits),'studies. Marker missing values:',int(markers.isna().sum().sum()))
