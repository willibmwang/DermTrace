"""Verify saved audit responses without contacting a language-model endpoint."""
import argparse
import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'demo'))
from audit_copilot_faithfulness import PACKETS,digest
from narrative_validation import validate_narrative

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('responses',type=Path)
p.add_argument('--packets',type=Path,default=PACKETS)
args=p.parse_args()
packets={r['profile_id']:r for r in map(json.loads,args.packets.read_text().splitlines())}
rows=list(map(json.loads,args.responses.read_text().splitlines()))
for row in rows:
    record=packets[row['profile_id']]
    assert digest(record['packet'])==record['packet_hash']==row['packet_hash']
    checks=validate_narrative(row['narrative'],record['packet'])['checks']
    for key,value in checks.items():
        assert value==row['checks'][key],(row['profile_id'],key)
    if row['requested_arm']=='llm' and row['generation_success']:
        assert row['actual_mode']=='grounded_llm'
        assert row['request_payload'] and row['prompt_hash']==digest(row['request_payload'])
print('Replayed',len(rows),'saved responses; validation agrees.')
