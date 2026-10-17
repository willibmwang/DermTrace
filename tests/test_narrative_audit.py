import copy
import os
import sys
import unittest
from pathlib import Path
from unittest.mock import patch,Mock
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'scripts'),str(ROOT/'demo')]
from narrative_validation import factual_anchor, DISCLAIMER, validate_narrative
from demo_utils import generate_grounded_summary,deterministic_summary,build_evidence_packet

class NarrativeTests(unittest.TestCase):
    def setUp(self):
        self.packet={'sample':{'sample_id':'private-id'},
            'prediction':{'decision':'keloid','prob_keloid':.83,'confidence':.83,'confidence_threshold':.55},
            'quality_control':{'transfer_status':'unverified','transfer_note':'Unseen cohort; transfer is unverified.'},
            'top_supporting_contributions':[],'top_opposing_contributions':[]}
    def valid(self):
        return factual_anchor(self.packet)+' '+self.packet['quality_control']['transfer_note']+' '+DISCLAIMER
    def test_valid_and_empty(self):
        self.assertTrue(validate_narrative(self.valid(),self.packet)['automated_pass'])
        self.assertFalse(validate_narrative('',self.packet)['automated_pass'])
    def test_swapped_numbers_rejected(self):
        text=self.valid().replace('P(keloid)=0.830','P(keloid)=0.550')
        self.assertFalse(validate_narrative(text,self.packet)['checks']['numeric_fields_grounded'])
    def test_both_class_names_do_not_hide_contradiction(self):
        text=self.valid().replace(DISCLAIMER,'This is unaffected, although keloid was considered. '+DISCLAIMER)
        self.assertFalse(validate_narrative(text,self.packet)['checks']['decision_consistent'])
    def test_claim_in_research_sentence_not_stripped(self):
        text=self.valid().replace(DISCLAIMER,'For research use, this confirms a diagnosis. '+DISCLAIMER)
        self.assertFalse(validate_narrative(text,self.packet)['checks']['no_prohibited_claims'])
    def test_suppressed_abstention(self):
        text=self.valid();self.packet['prediction']['decision']='abstain'
        self.assertFalse(validate_narrative(text,self.packet)['checks']['abstention_disclosed'])
    def test_unconfigured_llm_audit_fails(self):
        with patch.dict(os.environ,{},clear=True):
            with self.assertRaises(RuntimeError):generate_grounded_summary(self.packet,require_llm=True)
    def test_raw_and_runtime_outputs_separate_and_request_captured(self):
        response=Mock();response.json.return_value={'model':'test-version','choices':[{'message':{'content':'Unsupported 99 diagnosis.'},'finish_reason':'stop'}]}
        before=copy.deepcopy(self.packet);trace=[]
        with patch.dict(os.environ,{'KELOIDBENCH_LLM_API_URL':'https://example.invalid','KELOIDBENCH_LLM_MODEL':'test'}),patch('demo_utils.requests.post',return_value=response) as post:
            raw,mode=generate_grounded_summary(self.packet,require_llm=True,validate_output=False,request_recorder=trace.append)
            self.assertEqual(mode,'grounded_llm');self.assertEqual(raw,'Unsupported 99 diagnosis.')
            self.assertEqual(trace[0],post.call_args.kwargs['json'])
            self.assertNotIn('private-id',str(trace[0]))
            displayed,mode=generate_grounded_summary(self.packet)
            self.assertEqual(mode,'validation_fallback')
            self.assertEqual(displayed,deterministic_summary(self.packet))
        self.assertEqual(before,self.packet)
    def test_generation_error_is_not_success(self):
        response=Mock();response.json.return_value={'choices':[{'message':{'content':''}}]}
        with patch.dict(os.environ,{'KELOIDBENCH_LLM_API_URL':'https://example.invalid','KELOIDBENCH_LLM_MODEL':'test'}),patch('demo_utils.requests.post',return_value=response):
            with self.assertRaises(ValueError):generate_grounded_summary(self.packet,require_llm=True,validate_output=False)
if __name__=='__main__':unittest.main()
