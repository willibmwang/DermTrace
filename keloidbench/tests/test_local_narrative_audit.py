import copy
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'scripts'), str(ROOT / 'demo')]
import audit_copilot_faithfulness as audit
import demo_utils as du
from narrative_validation import DISCLAIMER, factual_anchor


class LocalAuditTests(unittest.TestCase):
    def test_local_backend_records_raw_failure_without_templates(self):
        packet = {'sample': {'sample_id': 'private-id'},
                  'prediction': {'decision': 'keloid', 'prob_keloid': .83,
                                 'confidence': .83, 'confidence_threshold': .55},
                  'quality_control': {'transfer_note': 'Unseen cohort; transfer is unverified.'},
                  'top_supporting_contributions': [], 'top_opposing_contributions': []}
        payload = {'model': 'test-snapshot', 'messages': [
            {'role': 'user', 'content': du._copilot_prompt(du.external_evidence_packet(packet))}]}
        trace = [payload, {'response_metadata': {'model': 'test-snapshot', 'finish_reason': 'stop'}}]
        valid = factual_anchor(packet) + ' ' + packet['quality_control']['transfer_note'] + ' ' + DISCLAIMER
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            packets = root / 'packets.jsonl'
            packets.write_text(json.dumps({'profile_id': 'private-id', 'accession': 'study',
                                           'packet': packet, 'packet_hash': audit.digest(packet)}) + '\n')
            args = SimpleNamespace(arm='llm', local_model=Path('test-snapshot'), batch_size=2,
                                   packets=packets, limit=None, repeats=2, output=root / 'run',
                                   model_name='test', narrative_format='full', max_new_tokens=400)
            generations = iter([(valid, copy.deepcopy(trace), None),
                                ('Unsupported diagnosis 99.', copy.deepcopy(trace), None)])
            with patch('local_narrative_model.generate_local_attempts', return_value=generations):
                audit.run(args)
            summary = json.loads((args.output / 'summary.json').read_text())
            rows = list(map(json.loads, (args.output / 'responses.jsonl').read_text().splitlines()))
            self.assertEqual(summary['generation_successes'], 2)
            self.assertEqual(summary['all_checks_pass'], 1)
            self.assertEqual(summary['profiles_all_repeats_pass'], 0)
            self.assertEqual(rows[1]['narrative'], 'Unsupported diagnosis 99.')
            self.assertTrue(all(row['checks']['request_redacted'] for row in rows))
            self.assertEqual(summary['backend'], 'local_transformers')

    def test_prose_format_assembles_anchors_and_keeps_raw_generation(self):
        packet = {'sample': {'sample_id': 'private-id'},
                  'prediction': {'decision': 'non_keloid', 'prob_keloid': .3,
                                 'confidence': .7, 'confidence_threshold': .55},
                  'quality_control': {'transfer_status': 'represented', 'transfer_note': 'Included in development.'},
                  'top_supporting_contributions': [{'feature': 'fibroblast_composition_rank',
                      'display_name': 'Fibroblast composition', 'value': .1, 'contribution': .5}],
                  'top_opposing_contributions': [{'feature': 'keratinocyte_composition_rank',
                      'display_name': 'Keratinocyte composition', 'value': .2, 'contribution': -1.2}]}
        payload = {'model': 'test-snapshot', 'messages': [
            {'role': 'user', 'content': du._copilot_prose_prompt(du.external_evidence_packet(packet))}]}
        trace = [payload, {'response_metadata': {'model': 'test-snapshot', 'finish_reason': 'stop'}}]
        good = 'Fibroblast composition pushed the score toward keloid, while keratinocyte composition pushed it toward non-keloid.'
        swapped = 'Fibroblast composition pushed the score toward non-keloid.'
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            packets = root / 'packets.jsonl'
            packets.write_text(json.dumps({'profile_id': 'private-id', 'accession': 'study',
                                           'packet': packet, 'packet_hash': audit.digest(packet)}) + '\n')
            args = SimpleNamespace(arm='llm', local_model=Path('test-snapshot'), batch_size=2,
                                   packets=packets, limit=None, repeats=3, output=root / 'run',
                                   model_name='test', narrative_format='prose', max_new_tokens=160)
            generations = iter([(good, copy.deepcopy(trace), None), (swapped, copy.deepcopy(trace), None),
                                ('', copy.deepcopy(trace), None)])
            with patch('local_narrative_model.generate_local_attempts', return_value=generations):
                audit.run(args)
            rows = list(map(json.loads, (args.output / 'responses.jsonl').read_text().splitlines()))
            self.assertEqual([r['all_checks_pass'] for r in rows], [True, False, False])
            self.assertTrue(rows[0]['narrative'].startswith(factual_anchor(packet)))
            self.assertEqual(rows[0]['raw_generation'], good)
            self.assertFalse(rows[1]['checks']['contribution_role'])
            self.assertEqual(rows[2]['error_type'], 'EmptyGeneration')


if __name__ == '__main__':
    unittest.main()
