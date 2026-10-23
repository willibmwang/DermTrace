"""Post hoc failure diagnostics for a frozen narrative audit run.

This script never changes the frozen run, prompt, or validator. The primary
result remains the raw-output pass count in the run's summary.json. The
label-normalized rescoring below was specified after inspecting the outputs; it
separates formatting failures from content failures and is not a pass rate.
"""
from __future__ import annotations
import argparse
import collections
import json
import re
import sys
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'scripts'),str(ROOT/'demo')]
from audit_copilot_faithfulness import PACKETS,digest
from narrative_validation import DISCLAIMER,factual_anchor,validate_narrative

LABEL=re.compile(r'^\s*(?:OPENING|QC NOTE|DISCLAIMER)\s*:\s*',re.M)
END_MARKER=re.compile(r'^\s*END OF (?:OUTPUT|RESPONSE)\.?\s*$',re.M|re.I)
NUMBER=re.compile(r'(?<![A-Za-z])[-+]?\d+(?:\.\d+)?')

def normalize_labels(text):
    """Remove echoed prompt section labels and end markers only; prose is untouched."""
    return re.sub(r'\s+',' ',END_MARKER.sub('',LABEL.sub('',text))).strip()

def extra_numbers(narrative,packet):
    body=narrative.replace('\u2212','-').replace(factual_anchor(packet),'').replace(DISCLAIMER,'')
    body=body.replace(packet.get('quality_control',{}).get('transfer_note',''),'')
    packet_numbers={float(n) for n in NUMBER.findall(json.dumps(packet))}
    found=[float(n) for n in NUMBER.findall(LABEL.sub('',body))]
    # Contributions are stored at full precision; a number is "from packet" if it rounds from one.
    def grounded(x):
        decimals=len(str(x).split('.')[1]) if '.' in str(x) else 0
        return any(round(v,decimals)==x for v in packet_numbers)
    return found,[x for x in found if not grounded(x)]

def attributed_numbers(narrative,packet):
    """Heuristically pair each number with the nearest preceding named program in its sentence.

    Returns (n_attributed, n_contribution_match, n_value_only_match). A match equals
    the signed field within one unit of the last written digit (rounding or truncation).
    A value-only match reports the program's input value, not its contribution.
    """
    items=packet.get('top_supporting_contributions',[])+packet.get('top_opposing_contributions',[])
    keys={}
    for item in items:
        for alias in {item.get('display_name','').lower(),item.get('feature','').lower().replace('_',' '),
                      item.get('display_name','').split(' ')[0].lower()}:
            if alias:keys[alias]=item
    text=LABEL.sub('',narrative).replace('\u2212','-').replace('_',' ').lower()
    total=contribution=value_only=0
    for sentence in re.split(r'(?<=[.!?])\s+(?=[A-Za-z])',text):
        events=sorted([(m.start(),'f',keys[a]) for a in keys for m in re.finditer(re.escape(a),sentence)]
                      +[(m.start(),'n',m.group()) for m in NUMBER.finditer(sentence)],key=lambda e:e[0])
        current=None
        for _,kind,obj in events:
            if kind=='f':current=obj;continue
            if current is None:continue
            x=float(obj);d=len(obj.split('.')[1]) if '.' in obj else 0;tol=10**-d+1e-9
            hit=abs(x-float(current['contribution']))<=tol
            total+=1;contribution+=hit;value_only+=not hit and abs(x-float(current.get('value',np.nan)))<=tol
    return total,contribution,value_only

def bootstrap(values,seed=13,draws=10000):
    rng=np.random.default_rng(seed)
    return np.percentile(np.mean(rng.choice(values,size=(draws,len(values)),replace=True),axis=1),[2.5,97.5]).tolist()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('run',type=Path)
    p.add_argument('--packets',type=Path,default=PACKETS)
    p.add_argument('--output',type=Path)
    args=p.parse_args()
    out=args.output or args.run.with_name(args.run.name+'_diagnostics')
    if out.exists():raise RuntimeError('Diagnostics directory exists; choose a new path.')
    packets={r['profile_id']:r for r in map(json.loads,args.packets.read_text().splitlines())}
    rows=list(map(json.loads,(args.run/'responses.jsonl').read_text().splitlines()))
    records=[]
    for row in rows:
        record=packets[row['profile_id']];packet=record['packet']
        assert digest(packet)==record['packet_hash']==row['packet_hash']
        text=row['narrative'];anchor=factual_anchor(packet)
        qc=packet['quality_control'].get('transfer_note','')
        # Use only the validator checks the run recorded; later-added checks are not retrofitted.
        normalized={k:v for k,v in validate_narrative(normalize_labels(text),packet)['checks'].items() if k in row['checks']}
        found,ungrounded=extra_numbers(text,packet)
        attributed,contribution,value_only=attributed_numbers(text,packet)
        records.append({'profile_id':row['profile_id'],'accession':row['accession'],'repeat':row['repeat'],
            'decision':row['decision'],'raw_failed_checks':';'.join(k for k,v in row['checks'].items() if not v),
            'echoed_section_labels':bool(LABEL.search(text)),'end_marker':bool(END_MARKER.search(text)),
            'anchor_verbatim_anywhere':anchor in text,'qc_note_verbatim':bool(qc) and qc in text,
            'disclaimer_count':text.count(DISCLAIMER),'words':len(text.split()),
            'extra_numbers':len(found),'extra_numbers_not_in_packet':len(ungrounded),
            'attributed_numbers':attributed,'attributed_contribution_match':contribution,'attributed_value_only_match':value_only,
            'normalized_failed_checks':';'.join(k for k,v in normalized.items() if not v),
            'normalized_all_content_checks_pass':all(v for k,v in normalized.items())
                and row['checks']['packet_immutable'] and row['checks']['request_redacted'] and row['checks']['not_truncated']})
    df=pd.DataFrame(records)
    n=len(df)
    def count(mask):return {'n':int(mask.sum()),'denominator':n}
    norm_failed=collections.Counter(c for s in df.normalized_failed_checks for c in s.split(';') if c)
    by_profile=df.groupby('profile_id').normalized_all_content_checks_pass.mean().to_numpy()
    summary={'run':str(args.run.relative_to(ROOT)) if args.run.is_absolute() else str(args.run),
        'status':'post hoc diagnostic; not a primary or pre-specified pass rate',
        'raw_primary_all_checks_pass':int(sum(r['all_checks_pass'] for r in rows)),'n_attempts':n,
        'raw_failed_check_patterns':dict(collections.Counter(df.raw_failed_checks).most_common()),
        'format':{'echoed_section_labels':count(df.echoed_section_labels),'end_marker':count(df.end_marker),
            'anchor_verbatim_anywhere':count(df.anchor_verbatim_anywhere),'qc_note_verbatim':count(df.qc_note_verbatim),
            'disclaimer_repeated':count(df.disclaimer_count>1),'disclaimer_missing':count(df.disclaimer_count==0)},
        'numbers':{'attempts_with_extra_numbers':count(df.extra_numbers>0),
            'attempts_with_numbers_not_in_packet':count(df.extra_numbers_not_in_packet>0),
            'program_attributed_numbers':{'n':int(df.attributed_numbers.sum()),
                'matches_program_contribution':int(df.attributed_contribution_match.sum()),
                'matches_program_input_value_only':int(df.attributed_value_only_match.sum()),
                'matches_neither':int((df.attributed_numbers-df.attributed_contribution_match-df.attributed_value_only_match).sum()),
                'attempts_with_any_non_contribution_number':count(df.attributed_contribution_match<df.attributed_numbers),
                'method':'Heuristic: each number is paired with the nearest preceding named program in its sentence and compared with that program signed contribution, then its input value, within one unit of the last written digit.'}},
        'label_normalized':{'description':'Echoed OPENING/QC NOTE/DISCLAIMER labels and END OF OUTPUT/RESPONSE lines removed; validator unchanged.',
            'all_checks_pass':count(df.normalized_all_content_checks_pass),
            'profile_cluster_bootstrap_95_interval':bootstrap(by_profile),
            'failed_checks':{k:{'n':v,'denominator':n} for k,v in norm_failed.most_common()},
            'by_decision':{d:{'n_attempts':int(len(g)),'all_checks_pass':int(g.normalized_all_content_checks_pass.sum())} for d,g in df.groupby('decision')}},
        'repeat_identical_outputs':int(sum(len({r['narrative'] for r in rows if r['profile_id']==pid})==1 for pid in df.profile_id.unique())),
        'n_profiles':int(df.profile_id.nunique())}
    out.mkdir(parents=True)
    df.to_csv(out/'per_attempt.csv',index=False)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:summary[k] for k in ['raw_primary_all_checks_pass','n_attempts']}|{'normalized_pass':summary['label_normalized']['all_checks_pass']}))

if __name__=='__main__':main()
