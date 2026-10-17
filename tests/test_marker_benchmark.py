import sys
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from run_marker_benchmark import weights,primary_arrays,select_config,fit_predict

class MarkerTests(unittest.TestCase):
    def test_study_donor_weights(self):
        meta=pd.DataFrame({'sample_id':['a','b','c','d','e'],'accession':['A','A','A','B','B'],
                           'split_group':['a','a','b','c','c']})
        w=weights(meta)
        self.assertAlmostEqual(w[:3].sum(),w[3:].sum())
        self.assertAlmostEqual(w[:2].sum(),w[2])
    def test_mixed_labels_not_collapsed(self):
        meta=pd.DataFrame({'split_group':['a','a','b','b']})
        y,p=primary_arrays([1,0,1,0],[.9,.1,.8,.2],meta)
        self.assertEqual(len(y),4)
    def test_selection_rejects_labels_outside_training(self):
        x=pd.DataFrame({'x':[0,1]},index=['train','heldout'])
        y=pd.Series([0,1],index=x.index)
        with self.assertRaises(AssertionError): select_config(x,y,pd.DataFrame(index=x.index),['train'],{})
    def test_heldout_values_do_not_change_other_predictions(self):
        x=pd.DataFrame({'x':[0.,1.,2.,3.,np.nan,2.]},index=list('abcdef'))
        y=pd.Series([0,0,1,1,0,1],index=x.index)
        meta=pd.DataFrame({'sample_id':list('abcdef'),'accession':['A']*4+['B']*2,'split_group':list('abcdef')},index=x.index)
        before=fit_predict(x,y,meta,list('abcd'),list('ef'),1.)
        x.loc['f','x']=1e9
        after=fit_predict(x,y,meta,list('abcd'),list('ef'),1.)
        self.assertAlmostEqual(before[0],after[0])
if __name__=='__main__':unittest.main()
