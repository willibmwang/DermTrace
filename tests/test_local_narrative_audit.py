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
                                   packets=packets, limit=None, repeats=2, output=root / 'run')
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


if __name__ == '__main__':
    unittest.main()
