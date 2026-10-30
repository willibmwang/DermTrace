import unittest
import hashlib
import json
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from experiments.expression_completion.loader import Profile, SplitAssignments, MaskPlan, align, fit_normalizer, make_episode, VisibleAnchorTotalNormalizer, select_training_anchors, make_mask_plan, AuthorAnnotatedCellPolicy
from experiments.expression_completion.anchor_normalization import FrozenAnchorNormalizer

G = tuple('ENSG'+str(i).zfill(11) for i in range(1, 7))

def profile(uid, v, partition='train', unit='counts', group='donor1'):
    return Profile(uid, 'source', uid, group, partition, unit, G, np.asarray(v, dtype=float), np.ones(6, bool))

class TestCompletion(unittest.TestCase):
    def setUp(self):
        self.p = profile('train1', [10,20,30,40,50,60])
        self.plan = MaskPlan(np.array([1,1,1,0,0,0],bool),np.array([0,0,0,1,0,0],bool),np.array([0,0,0,0,1,0],bool),np.array([1,1,0,0,0,0],bool))
        self.norm = fit_normalizer([self.p], G[:2], FrozenAnchorNormalizer)

    def test_group_and_sample_leakage_rejected(self):
        with self.assertRaises(ValueError):
            SplitAssignments([dict(source_id='a',sample_id='1',group_id='d',partition='train'),dict(source_id='a',sample_id='2',group_id='d',partition='test')])
        with self.assertRaises(ValueError):
            SplitAssignments([dict(source_id='a',sample_id='1',group_id='d',partition='train')]*2)

    def test_split_and_metadata_hashes_checked(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); sample=root/'metadata.csv'; sample.write_text('immutable')
            split=root/'fold.json'; split.write_text(json.dumps({'assignments':[dict(source_id='GSE144424',sample_id='1',group_id='d',partition='train')]}))
            digest=lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
            (root/'manifest.json').write_text(json.dumps({'input_metadata_sha256':{'cardiac':digest(sample)},'files':[{'file':split.name,'sha256':digest(split)}]}))
            with self.assertRaises(ValueError): SplitAssignments.read(split)
            self.assertEqual(SplitAssignments.read(split, metadata_paths={'cardiac':sample}).get('GSE144424','1'),('d','train'))
            sample.write_text('changed')
            with self.assertRaises(ValueError):SplitAssignments.read(split, metadata_paths={'cardiac':sample})
            sample.write_text('immutable');split.write_text('{}')
            with self.assertRaises(ValueError):SplitAssignments.read(split, metadata_paths={'cardiac':sample})

    def test_per_split_metadata_are_bound_and_consistent(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            for name in ('skin','cardiac','melanoma'):
                (root/name).write_text(name)
            sha=lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
            hashes={name:sha(root/name) for name in ('skin','cardiac','melanoma')}
            split=root/'fold.json'
            payload={'assignments':[dict(source_id='GSE163973',sample_id='KF1',group_id='d',partition='train')],
                     'input_metadata_sha256':{'skin':hashes['skin']}}
            split.write_text(json.dumps(payload))
            record={'file':split.name,'sha256':sha(split),'input_metadata_sha256':{'skin':hashes['skin']}}
            manifest=root/'manifest.json'
            manifest.write_text(json.dumps({'input_metadata_sha256':hashes,'files':[record]}))
            self.assertEqual(SplitAssignments.read(split,metadata_paths={'skin':root/'skin'}).get('GSE163973','KF1'),('d','train'))
            with self.assertRaises(ValueError):
                SplitAssignments.read(split,metadata_paths={'skin':root/'skin','cardiac':root/'cardiac'})
            record['input_metadata_sha256']['skin']='0'*64
            manifest.write_text(json.dumps({'input_metadata_sha256':hashes,'files':[record]}))
            with self.assertRaisesRegex(ValueError,'disagree'):
                SplitAssignments.read(split,metadata_paths={'skin':root/'skin'})

    def test_split_source_hashes_must_match_payload_and_cover_assignments(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            for name in ('skin','cardiac'):
                (root/name).write_text(name)
            digest=lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
            hashes={name:digest(root/name) for name in ('skin','cardiac')}
            rows=[dict(source_id='GSE163973',sample_id='KF1',group_id='skin1',partition='train'),
                  dict(source_id='GSE144424',sample_id='c1',group_id='donor1',partition='test')]
            split=root/'fold.json'
            manifest=root/'manifest.json'
            def write(payload, per_file):
                split.write_text(json.dumps(payload))
                manifest.write_text(json.dumps({'input_metadata_sha256':hashes,'files':[
                    {'file':split.name,'sha256':digest(split),'input_metadata_sha256':per_file}]}))
            write({'assignments':rows,'input_metadata_sha256':{'skin':hashes['skin']}}, {'skin':hashes['skin']})
            with self.assertRaisesRegex(ValueError,'cover'):
                SplitAssignments.read(split,metadata_paths={'skin':root/'skin'})
            write({'assignments':rows,'input_metadata_sha256':{'skin':hashes['skin']}},hashes)
            with self.assertRaisesRegex(ValueError,'cover'):
                SplitAssignments.read(split,metadata_paths={'skin':root/'skin','cardiac':root/'cardiac'})
            write({'assignments':rows,'input_metadata_sha256':hashes},hashes)
            self.assertEqual(SplitAssignments.read(split,metadata_paths={'skin':root/'skin','cardiac':root/'cardiac'}).get('GSE144424','c1'),('donor1','test'))

    def test_cell_qc_is_audited_but_does_not_auto_exclude_flags(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'qc.parquet'
            pd.DataFrame([dict(sample='KF1',barcode='bc',annotated=True,cell_type='Fibroblast',
                flag_low_counts=True,flag_low_features=False,flag_high_mito=False,flag_any=True)]).to_parquet(path)
            policy=AuthorAnnotatedCellPolicy(path,hashlib.sha256(path.read_bytes()).hexdigest())
            row=pd.Series(dict(source_metadata_row=2,cellType='Fibroblast'))
            self.assertTrue(policy('KF1','bc',row))
            self.assertEqual(policy.flags('KF1','bc'),('flag_low_counts',))
            with self.assertRaises(ValueError):AuthorAnnotatedCellPolicy(path,'0'*64)

    def test_train_only_and_measurement_family(self):
        with self.assertRaises(ValueError):
            fit_normalizer([self.p, profile('test',[1]*6,partition='test')], G[:2], FrozenAnchorNormalizer)
        with self.assertRaises(ValueError):
            fit_normalizer([self.p, profile('other',[1]*6,unit='other')], G[:2], FrozenAnchorNormalizer)

    def test_sparse_anchor_total_supports_zero_anchors(self):
        p=profile('sparse',[0,20,30,40,50,60])
        normalizer=fit_normalizer([p],G[:2],VisibleAnchorTotalNormalizer)
        obs,reward,ev=make_episode(p,self.plan,normalizer)
        self.assertTrue(np.isfinite(obs.context_values).all())
        self.assertEqual(normalizer.factor(pd.Series([0,20],index=G[:2])),20)
        with self.assertRaises(ValueError):
            normalizer.factor(pd.Series([0,0],index=G[:2]))

    def test_training_anchor_selection_and_seeded_masks(self):
        anchors=select_training_anchors(iter([self.p]),G,2)
        self.assertEqual(anchors,G[:2])
        first=make_mask_plan(self.p,anchors,seed=123)
        second=make_mask_plan(self.p,anchors,seed=123)
        np.testing.assert_array_equal(first.context,second.context)
        np.testing.assert_array_equal(first.reward,second.reward)
        np.testing.assert_array_equal(first.evaluation,second.evaluation)
        self.assertTrue(first.context[:2].all())
        with self.assertRaises(ValueError):
            select_training_anchors([profile('test',[1]*6,partition='test')],G,2)

    def test_measured_zero_and_missing_are_distinct(self):
        p=profile('p',[0,2,3,4,5,6])
        x=align(p,[G[0],'ENSG99999999999'])
        np.testing.assert_array_equal(x.values,[0,0]); np.testing.assert_array_equal(x.measured,[True,False])

    def test_disjoint_masks_and_anchors(self):
        bad=MaskPlan(self.plan.context,self.plan.context,self.plan.evaluation,self.plan.anchors)
        with self.assertRaises(ValueError): make_episode(self.p,bad,self.norm)
        bad=MaskPlan(self.plan.context,self.plan.reward,self.plan.evaluation,np.array([0,0,0,0,1,0],bool))
        with self.assertRaises(ValueError): make_episode(self.p,bad,self.norm)
        missing=self.p.measured.copy();missing[4]=False
        with self.assertRaises(ValueError): self.plan.validate(missing)

    def test_hidden_value_cannot_change_input_or_reward(self):
        a, r1, e1 = make_episode(self.p,self.plan,self.norm)
        altered=profile('train1',[10,20,30,40,5000,60])
        b, r2, e2 = make_episode(altered,self.plan,self.norm)
        np.testing.assert_array_equal(a.context_values,b.context_values)
        self.assertEqual(r1.reward_loss([0]),r2.reward_loss([0]))
        r1.stop(); r2.stop()
        self.assertNotEqual(e1.evaluate_once([0]),e2.evaluate_once([0]))
        self.assertEqual(a.context_values[3],0)
        self.assertEqual(a.context_values[4],0)
        self.assertFalse(hasattr(r1,'evaluate_once'))
        self.assertFalse(hasattr(e1,'reward_loss'))

    def test_reward_only_returns_values_after_stop(self):
        obs, reward, evaluation=make_episode(self.p,self.plan,self.norm)
        loss=reward.reward_loss([0]); self.assertGreater(loss,0)
        with self.assertRaises(ValueError):evaluation.evaluate_once([0])
        values=reward.stop();self.assertEqual(len(values),1)
        # Returned reward values use the same factor computed from C anchors.
        factor=self.norm.factor(pd.Series(self.p.values[:3],index=G[:3]))
        expected=self.norm.transform(pd.Series([self.p.values[3]],index=[G[3]]),factor)
        np.testing.assert_allclose(values,expected.to_numpy())
        with self.assertRaises(ValueError):reward.reward_loss([0])
        with self.assertRaises(ValueError):reward.stop()
        self.assertTrue(np.isfinite(evaluation.evaluate_once([0])))
        with self.assertRaises(ValueError):evaluation.evaluate_once([0])

    def test_anchor_order_does_not_change_membership(self):
        anchors = (G[1], G[0])
        norm = fit_normalizer([self.p], anchors, VisibleAnchorTotalNormalizer)
        plan = make_mask_plan(self.p, anchors, seed=4)
        obs, _, _ = make_episode(self.p, plan, norm)
        self.assertTrue(np.isfinite(obs.context_values).all())

    def test_reward_query_budget(self):
        _, reward, evaluator = make_episode(self.p, self.plan, self.norm, max_reward_queries=2)
        self.assertEqual(reward.queries_used, 0)
        with self.assertRaises(AttributeError):
            reward.max_queries = 100
        with self.assertRaises(AttributeError):
            reward.queries_used = 0
        reward.reward_loss([0])
        reward.reward_loss([1])
        self.assertEqual(reward.queries_used, 2)
        with self.assertRaises(ValueError):
            reward.reward_loss([2])
        reward.stop()
        self.assertTrue(np.isfinite(evaluator.evaluate_once([0])))
        with self.assertRaises(ValueError):
            make_episode(self.p, self.plan, self.norm, max_reward_queries=0)

if __name__ == '__main__': unittest.main()
