"""Freeze evidence packets, audit actual narratives, and export blinded human-review forms.

Automated passes are screening results, not human-validated faithfulness rates.
LLM mode fails closed when no endpoint/model is configured and records every attempt.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
import os
import random
import sys
from pathlib import Path
from urllib.parse import urlparse
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'demo'))
import demo_utils as du
from narrative_validation import validate_narrative

BASE=ROOT/'results/reviewer_revision/narrative_audit'
PACKETS=BASE/'frozen_packets.jsonl'

def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True).encode()).hexdigest()

def freeze_packets(path):
    if path.exists(): raise RuntimeError('Packets already frozen; choose a new output path for a new version.')
    manifest=du.load_profile_manifest(); _,bundle=du.load_product_bundle()
    feature=du.load_feature_views()[bundle['meta']['feature_set']]
    records=[]
    for sid in du.primary_example_ids(manifest):
        record=manifest[manifest.sample_id.eq(sid)].iloc[0]
        features=feature.loc[[sid]]
        score=du.score_sample(features,bundle,is_expression=False)
        qc=du.expression_qc(du.raw_expression_for_sample(sid,manifest),accession=str(record.accession))
        packet=du.build_evidence_packet(score,qc,du.program_score_table(features),accession=str(record.accession),sample_label=sid)
        records.append({'profile_id':sid,'accession':str(record.accession),'packet':packet,'packet_hash':digest(packet)})
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(''.join(json.dumps(r)+'\n' for r in records))
    print('Frozen',len(records),'packets;',sum(r['packet']['prediction']['decision']=='abstain' for r in records),'abstentions')

def run(args):
    if args.arm=='llm' and not all(os.environ.get(k,'').strip() for k in ['KELOIDBENCH_LLM_API_URL','KELOIDBENCH_LLM_MODEL']):
        raise RuntimeError('LLM audit requires a configured endpoint and model; deterministic fallback is forbidden.')
    if args.output.exists(): raise RuntimeError('Output directory exists; preserve it and choose a new run directory.')
    records=[json.loads(line) for line in args.packets.read_text().splitlines()]
    if args.limit: records=records[:args.limit]
    args.output.mkdir(parents=True)
    attempts=[]
    for record in records:
        packet=record['packet']; assert digest(packet)==record['packet_hash']
        for repeat in range(args.repeats):
            before=copy.deepcopy(packet);trace=[];narrative='';mode='generation_error';error=None
            try:
                if args.arm=='llm':
                    narrative,mode=du.generate_grounded_summary(packet,require_llm=True,request_recorder=trace.append,validate_output=False)
                    if mode!='grounded_llm': raise RuntimeError('Non-generative output in LLM audit')
                else: narrative,mode=du.deterministic_summary(packet),'deterministic'
            except Exception as exc:
                error=type(exc).__name__
            result=validate_narrative(narrative,packet)
            checks=result['checks']
            checks['packet_immutable']=packet==before
            payload=next((t for t in trace if 'messages' in t),None)
            external=du.external_evidence_packet(before)
            # Inspect the payload passed to requests.post, not a reconstructed stand-in.
            checks['request_redacted']=(payload is not None and record['profile_id'] not in json.dumps(payload)
                                        and json.loads(payload['messages'][0]['content'].split('EVIDENCE_JSON:\n',1)[1])==external) if args.arm=='llm' else True
            response=next((t['response_metadata'] for t in trace if 'response_metadata' in t),{})
            checks['not_truncated']=response.get('finish_reason')!='length' if args.arm=='llm' else True
            row={'profile_id':record['profile_id'],'accession':record['accession'],'repeat':repeat,
                 'decision':packet['prediction']['decision'],'requested_arm':args.arm,'actual_mode':mode,
                 'packet_hash':record['packet_hash'],'prompt_hash':digest(payload) if payload else None,
                 'generation_success':error is None,'error_type':error,'response_metadata':response,
                 'checks':checks,'all_checks_pass':error is None and all(checks.values()),'narrative':narrative,
                 'request_payload':payload}
            attempts.append(row)
            with (args.output/'responses.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
    names=attempts[0]['checks']
    by_profile={r['profile_id']:[] for r in attempts}
    for r in attempts:by_profile[r['profile_id']].append(r['all_checks_pass'])
    summary={'requested_arm':args.arm,'model':os.environ.get('KELOIDBENCH_LLM_MODEL') if args.arm=='llm' else None,
             'endpoint_host':urlparse(os.environ.get('KELOIDBENCH_LLM_API_URL','')).hostname if args.arm=='llm' else None,
             'n_profiles':len(records),'n_attempts':len(attempts),'repeats':args.repeats,
             'generation_successes':sum(r['generation_success'] for r in attempts),
             'all_checks_pass':sum(r['all_checks_pass'] for r in attempts),
             'profiles_all_repeats_pass':sum(all(v) for v in by_profile.values()),
             'checks':{k:{'n_pass':sum(r['checks'][k] for r in attempts),'denominator':len(attempts)} for k in names},
             'subgroups':{a:{'n_attempts':len([r for r in attempts if r['accession']==a]),
                             'all_checks_pass':sum(r['all_checks_pass'] for r in attempts if r['accession']==a)} for a in sorted({r['accession'] for r in attempts})},
             'human_review_status':'pending independent raters',
             'interpretation':'Automated screening only; repeats are clustered by profile. No human-validated generative pass rate is established.',
             'packets_sha256':hashlib.sha256(args.packets.read_bytes()).hexdigest(),
             'validator_sha256':hashlib.sha256((ROOT/'scripts/narrative_validation.py').read_bytes()).hexdigest()}
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    # Frozen stratification by study and decision, before any ratings are inspected.
    chosen=[];rng=random.Random(13)
    first=[r for r in attempts if r['repeat']==0]
    for acc in sorted({r['accession'] for r in first}):
        for decision in sorted({r['decision'] for r in first}):
            pool=[r for r in first if r['accession']==acc and r['decision']==decision]
            rng.shuffle(pool);chosen.extend(pool[:2])
    rng.shuffle(chosen)
    forms=[];key=[]
    for i,row in enumerate(chosen):
        blind=f'item-{i+1:03d}'
        forms.append({'item_id':blind,'narrative':row['narrative'],'packet':json.dumps(du.external_evidence_packet(next(r['packet'] for r in records if r['profile_id']==row['profile_id']))),
                      'factual_faithfulness':'','decision_preserved':'','abstention_preserved':'','qc_preserved':'','clarity_1_to_5':'','completeness_1_to_5':'','error_notes':''})
        key.append({'item_id':blind,'profile_id':row['profile_id'],'arm':args.arm})
    for rater in ['a','b']:pd.DataFrame(forms).to_csv(args.output/f'human_review_rater_{rater}.csv',index=False)
    pd.DataFrame(key).to_csv(args.output/'human_review_key.csv',index=False)
    print(json.dumps({k:summary[k] for k in ['requested_arm','n_profiles','n_attempts','all_checks_pass','human_review_status']}))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--arm',choices=['deterministic','llm'],default='deterministic')
    p.add_argument('--packets',type=Path,default=PACKETS)
    p.add_argument('--freeze-packets',action='store_true')
    p.add_argument('--repeats',type=int,default=1)
    p.add_argument('--limit',type=int)
    p.add_argument('--output',type=Path,default=BASE/'deterministic_v1')
    args=p.parse_args()
    if args.repeats<1: p.error('--repeats must be positive')
    if args.freeze_packets:freeze_packets(args.packets)
    else:run(args)
if __name__=='__main__':main()
