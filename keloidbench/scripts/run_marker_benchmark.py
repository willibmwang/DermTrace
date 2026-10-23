"""Run the frozen matched marker/program comparison without consulting outer labels."""
from __future__ import annotations
import json
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import f1_score
from revision_report import ROOT, REPO_ROOT, OUT, sha256, summarize, selective_report

PROTOCOL = ROOT/'protocols/reviewer_revision_v1.json'

def weights(meta):
    groups = meta.groupby(['accession','split_group']).sample_id.transform('size')
    donors = meta[['accession','split_group']].drop_duplicates().groupby('accession').size()
    w = 1.0 / (groups.to_numpy() * meta.accession.map(donors).to_numpy())
    return w / w.mean()

def fit_predict(x, y, meta, train, test, c):
    model = make_pipeline(SimpleImputer(strategy='median',keep_empty_features=True), StandardScaler(),
                          LogisticRegression(C=c,max_iter=4000,random_state=13))
    model.fit(x.loc[train],y.loc[train],logisticregression__sample_weight=weights(meta.loc[train]))
    return model.predict_proba(x.loc[test])[:,list(model.classes_).index(1)]

def primary_arrays(y, p, meta):
    # Internal y=1 denotes keloid; saved predictions explicitly convert to code 0.
    frame = pd.DataFrame({'y':np.asarray(y),'p':np.asarray(p),'donor':meta.split_group.to_numpy()})
    groups = frame.groupby('donor')
    if groups.ngroups > 1 and groups.y.nunique().max() == 1:
        donor = groups[['y','p']].mean()
        return donor.y.to_numpy(dtype=int), donor.p.to_numpy()
    return frame.y.to_numpy(dtype=int), frame.p.to_numpy()

def calibration_features(p):
    p = np.clip(p,1e-6,1-1e-6)
    return np.log(p/(1-p)).reshape(-1,1)

def select_config(x,y,meta,train,protocol):
    # This function receives only outer-training labels and metadata.
    assert set(x.index)==set(y.index)==set(meta.index)==set(train)
    candidates=[]
    for c in protocol['models']['C']:
        oof = pd.Series(index=train,dtype=float)
        for acc in sorted(meta.accession.unique()):
            valid=meta.index[meta.accession.eq(acc)].tolist()
            fit=meta.index[~meta.accession.eq(acc)].tolist()
            assert not set(meta.loc[fit,'split_group']) & set(meta.loc[valid,'split_group'])
            oof.loc[valid]=fit_predict(x,y,meta,fit,valid,c)
        assert oof.notna().all()
        calibrator=LogisticRegression(C=1.0,max_iter=4000,random_state=13)
        calibrator.fit(calibration_features(oof.to_numpy()),y.loc[train],sample_weight=weights(meta.loc[train]))
        calibrated=calibrator.predict_proba(calibration_features(oof.to_numpy()))[:,1]
        calibrated=pd.Series(calibrated,index=train)
        scored=[]
        for threshold in protocol['thresholds']:
            fold_scores=[]
            for _,group in meta.groupby('accession'):
                yy,pp=primary_arrays(y.loc[group.index],calibrated.loc[group.index],group)
                fold_scores.append(f1_score(yy,pp>=threshold,average='macro',zero_division=0))
            key=(float(np.mean(fold_scores)),float(min(fold_scores)),-abs(threshold-.5))
            scored.append((key,threshold))
        key,threshold=max(scored,key=lambda t:t[0])
        candidates.append((key,{'C':c,'threshold':threshold,'inner_mean_macro_f1':key[0],
                                'inner_worst_macro_f1':key[1],
                                'platt_coef':float(calibrator.coef_[0,0]),'platt_intercept':float(calibrator.intercept_[0])}))
    return max(candidates,key=lambda t:t[0])[1]

def apply_calibration(p,cfg):
    z=calibration_features(p).ravel()*cfg['platt_coef']+cfg['platt_intercept']
    return 1/(1+np.exp(-np.clip(z,-700,700)))

def main():
    protocol=json.loads(PROTOCOL.read_text())
    for rel,expected in protocol['input_sha256'].items():
        if sha256((REPO_ROOT if rel.startswith('data/') else ROOT)/rel)!=expected: raise RuntimeError(f'Frozen input changed: {rel}')
    meta=pd.read_csv(OUT/'cohort_manifest.csv').set_index('sample_id',drop=False)
    y=meta[protocol['endpoint']].eq('keloid').astype(int)
    markers=pd.read_csv(OUT/'marker_inputs.csv').set_index('sample_id')
    views={f'marker_{g}':markers[[g]] for g in protocol['markers']}
    views.update({f'program_{v}':pd.read_parquet(REPO_ROOT/f'data/processed/training/features/{v}.parquet').set_index('sample_id').loc[meta.index] for v in protocol['program_views']})
    coverage=[]
    for arm,x in views.items():
        for acc,group in meta.groupby('accession'):
            coverage.append({'arm':arm,'accession':acc,'n_profiles':len(group),'n_features':x.shape[1],
                             'missing_values':int(x.loc[group.index].isna().sum().sum())})
    pd.DataFrame(coverage).to_csv(OUT/'marker_program_coverage.csv',index=False)
    rows=[]; configs=[]
    for split in protocol['splits']:
        train,test=split['train_sample_ids'],split['test_sample_ids']
        assert not set(train)&set(test)
        assert not set(meta.loc[train,'accession'])&set(meta.loc[test,'accession'])
        assert not set(meta.loc[train,'split_group'])&set(meta.loc[test,'split_group'])
        for arm,x in views.items():
            cfg=select_config(x.loc[train],y.loc[train],meta.loc[train],train,protocol)
            p=apply_calibration(fit_predict(x,y,meta,train,test,cfg['C']),cfg)
            configs.append({'arm':arm,'accession':split['accession'],**cfg})
            for sid,prob in zip(test,p):
                rows.append({'arm':arm,'accession':split['accession'],'sample_id':sid,
                             'donor':meta.loc[sid,'split_group'],'y_true':1-int(y.loc[sid]),
                             'prob_keloid':float(prob),'threshold':cfg['threshold']})
        print('Completed',split['accession'],flush=True)
    preds=pd.DataFrame(rows); selections=pd.DataFrame(configs)
    # Select feature families from inner scores only, then reuse the corresponding outer predictions.
    for family in ['marker','program']:
        winners=[]
        for acc,sub in selections[selections.arm.str.startswith(family+'_')].groupby('accession'):
            winner=sub.sort_values(['inner_mean_macro_f1','inner_worst_macro_f1','arm'],ascending=[False,False,True]).iloc[0]
            selected=preds[preds.arm.eq(winner.arm)&preds.accession.eq(acc)].copy()
            selected['arm']='inner_selected_'+family
            winners.append(selected)
            configs.append({**winner.to_dict(),'arm':'inner_selected_'+family,'selected_feature_arm':winner.arm})
        preds=pd.concat([preds,*winners],ignore_index=True)
    preds.to_csv(OUT/'matched_predictions.csv',index=False)
    pd.DataFrame(configs).to_csv(OUT/'matched_selected_configs.csv',index=False)
    report_results(preds,protocol,meta)

def report_results(preds,protocol,meta):
    summaries=[]; fold_rows=[]; details={}
    for arm,sub in preds.groupby('arm'):
        result=summarize(sub,meta.reset_index(drop=True))
        sensitivity=summarize(sub[~sub.accession.eq('GSE173900')],meta.reset_index(drop=True))
        details[arm]={**result,'selective':selective_report(sub),'post_hoc_seven_study_sensitivity':sensitivity}
        summaries.append({'arm':arm,**{k:v for k,v in result.items() if not isinstance(v,(list,dict))}})
        fold_rows.extend({'arm':arm,**r} for r in result['folds'])
    table=pd.DataFrame(summaries).sort_values('arm')
    table.to_csv(OUT/'matched_summary.csv',index=False)
    pd.DataFrame([{'arm':arm, **{k:v for k,v in item['post_hoc_seven_study_sensitivity'].items() if not isinstance(v,(dict,list))}} for arm,item in details.items()]).to_csv(OUT/'matched_sensitivity_summary.csv',index=False)
    folds=pd.DataFrame(fold_rows);folds.to_csv(OUT/'matched_per_study.csv',index=False)
    pivot=folds.pivot(index='accession',columns='arm',values='macro_f1')
    delta=(pivot.inner_selected_program-pivot.inner_selected_marker).to_numpy()
    rng=np.random.default_rng(protocol['seed'])
    boot=np.mean(rng.choice(delta,size=(10000,len(delta)),replace=True),axis=1)
    paired={'contrast':'inner_selected_program minus inner_selected_marker',
            'mean_delta':float(delta.mean()),'study_bootstrap_95_percentile_interval':np.percentile(boot,[2.5,97.5]).tolist(),
            'per_study_delta':dict(zip(pivot.index,delta)),
            'limitation':'Descriptive paired resampling of eight studies; overlapping training sets and repeated development limit inference.'}
    (OUT/'matched_results.json').write_text(json.dumps({'protocol_sha256':sha256(PROTOCOL),'arms':details,'paired_comparison':paired},indent=2)+'\n')
    lines=['# Matched marker/program comparison','',
           'Revision-time nested LOSO on all eight development studies. Each fixed arm uses the same classifier family, C grid, training-only preprocessing, inner-OOF calibration, and threshold selection. Historical ensemble results are a separate evidence object.','',
           '| Arm | Mean primary macro F1 | Worst fold | Mean sample macro F1 |','|---|---:|---:|---:|']
    for r in table.itertuples(): lines.append(f'| {r.arm} | {r.mean_macro_f1:.4f} | {r.worst_macro_f1:.4f} | {r.mean_sample_macro_f1:.4f} |')
    lines += ['',f"Inner-selected program minus inner-selected marker: **{paired['mean_delta']:+.4f}**; descriptive paired study bootstrap interval **[{paired['study_bootstrap_95_percentile_interval'][0]:+.4f}, {paired['study_bootstrap_95_percentile_interval'][1]:+.4f}]**.",'',
              'No universal superiority claim is made. The development cohorts were already inspected; these results require independent confirmation. Missingness and per-study predictions are published alongside this report.']
    (OUT/'matched_comparison.md').write_text('\n'.join(lines)+'\n')
    print(table[['arm','mean_macro_f1','worst_macro_f1']].to_string(index=False))
if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report-only',action='store_true')
    if parser.parse_args().report_only:
        report_results(pd.read_csv(OUT/'matched_predictions.csv'),json.loads(PROTOCOL.read_text()),pd.read_csv(OUT/'cohort_manifest.csv').set_index('sample_id',drop=False))
    else:main()
