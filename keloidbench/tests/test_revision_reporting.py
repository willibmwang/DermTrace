import sys
import unittest
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from score_primary_ladder import score_from_predictions

class ReportingTests(unittest.TestCase):
    def test_keloid_zero_encoding(self):
        frame = pd.DataFrame({'accession':['a']*4,'donor':['a','b','c','d'],
                              'y_true':[0,0,1,1],'prob_keloid':[.9,.8,.2,.1],'threshold':[.5]*4})
        self.assertEqual(score_from_predictions(frame,pd.DataFrame())[0]['macro_f1'],1.0)
    def test_mixed_donor_uses_sample_contrast(self):
        frame = pd.DataFrame({'accession':['a']*4,'donor':['a','a','b','b'],
                              'y_true':[0,1,0,1],'prob_keloid':[.9,.1,.8,.2],'threshold':[.5]*4})
        result = score_from_predictions(frame,pd.DataFrame())[0]
        self.assertEqual(result['grain'],'sample_mixed_donor')
        self.assertEqual(result['macro_f1'],1.0)
if __name__ == '__main__': unittest.main()
