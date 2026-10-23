import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'scripts'), str(ROOT / 'demo')]
from diagnose_narrative_audit import attributed_numbers, normalize_labels
from narrative_validation import DISCLAIMER, factual_anchor, validate_narrative

PACKET = {'prediction': {'decision': 'keloid', 'prob_keloid': .83, 'confidence': .83, 'confidence_threshold': .55},
          'quality_control': {'transfer_note': 'Unseen cohort; transfer is unverified.'},
          'top_supporting_contributions': [{'feature': 'fibroblast_composition_rank', 'display_name': 'Fibroblast composition',
                                            'value': -0.3985, 'contribution': 1.1829}],
          'top_opposing_contributions': [{'feature': 'keratinocyte_composition_rank', 'display_name': 'Keratinocyte composition',
                                          'value': 0.7377, 'contribution': -2.0602}]}


class NarrativeDiagnosticTests(unittest.TestCase):
    def test_label_normalization_only_removes_echoed_labels(self):
        raw = ('OPENING: ' + factual_anchor(PACKET) + '\nQC NOTE: Unseen cohort; transfer is unverified.\n'
               'Fibroblast composition supports the result.\nEND OF OUTPUT.\nDISCLAIMER: ' + DISCLAIMER)
        self.assertFalse(validate_narrative(raw, PACKET)['checks']['factual_anchor'])
        self.assertTrue(validate_narrative(normalize_labels(raw), PACKET)['automated_pass'])

    def test_attribution_separates_contributions_values_and_unsupported_numbers(self):
        text = ('Fibroblast composition contributed 1.18, while keratinocyte_composition_rank contributed 0.7377. '
                'Keratinocyte composition accounts for 57.5%.')
        self.assertEqual(attributed_numbers(text, PACKET), (3, 1, 1))
        self.assertEqual(attributed_numbers('Keratinocyte composition contributed −2.06.', PACKET), (1, 1, 0))


if __name__ == '__main__':
    unittest.main()
