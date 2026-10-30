import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from experiments.expression_completion.splits import (
    collect_samples, create_splits, validate_assignments, write_splits,
)


class SplitTests(unittest.TestCase):
    def setUp(self):
        self.skin = [{'source_donor_key': x} for x in ('KF1','KF2','KF3','NF1','NF2','NF3')]
        self.cardiac = [dict(sample_id=f'c{d}_{j}', donor_id=f'd{d}',
                             published_primary_eligible='True', quality_flag='no_specific_paper_flag')
                        for d in range(4) for j in range(2)]
        self.melanoma = [dict(sample_id=f'm{i}', cell_line_family=fam, clean_hypoxia_eligible='True')
                         for i, fam in enumerate(('501MEL', 'IGR37_IGR39', 'IGR37_IGR39'))]

    def test_grouped_rotations_and_target_isolation(self):
        samples = collect_samples(self.skin, self.cardiac, self.melanoma)
        splits = create_splits(samples)
        self.assertEqual(len(splits), 31)
        for split in splits.values():
            rows = split['assignments']
            validate_assignments(rows)
            scenario = split['scenario']
            if '_to_' in scenario:
                source, target = scenario.split('_to_')
                names = {'skin': 'GSE163973', 'cardiac': 'GSE144424', 'melanoma': 'GSE95280_recount3'}
                self.assertEqual({r['partition'] for r in rows if r['source_id'] == names[target]}, {'test'})
                self.assertEqual({r['partition'] for r in rows if r['source_id'] == names[source]}, {'train','validation'})
            if scenario in ('within_keloid', 'within_normal_scar'):
                prefix = 'KF' if scenario == 'within_keloid' else 'NF'
                self.assertTrue(all(r['sample_id'].startswith(prefix) for r in rows))
        cardiac = [s for s in splits.values() if s['scenario'] == 'within_cardiac']
        for donor in range(4):
            self.assertEqual(sum(any(r['group_id'].endswith(f'd{donor}') and r['partition'] == 'test'
                                     for r in s['assignments']) for s in cardiac), 1)

    def test_bad_eligibility_or_identity_fails(self):
        self.cardiac[0]['donor_id'] = ''
        with self.assertRaises(ValueError):
            collect_samples(self.skin, self.cardiac, self.melanoma)
        self.cardiac[0]['donor_id'] = 'd0'
        self.cardiac[0]['quality_flag'] = 'published_expression_outlier'
        with self.assertRaises(ValueError):
            collect_samples(self.skin, self.cardiac, self.melanoma)
        self.cardiac[0]['quality_flag'] = 'no_specific_paper_flag'
        self.melanoma[0]['cell_line_family'] = 'unknown'
        with self.assertRaises(ValueError):
            collect_samples(self.skin, self.cardiac, self.melanoma)

    def test_duplicate_and_cross_partition_group_fail(self):
        rows = [dict(source_id='x', sample_id='a', group_id='g', partition='train'),
                dict(source_id='x', sample_id='b', group_id='g', partition='test')]
        with self.assertRaises(ValueError): validate_assignments(rows)
        rows[1]['partition'] = 'train'
        rows[1]['sample_id'] = 'a'
        with self.assertRaises(ValueError): validate_assignments(rows)

    def test_manifest_has_only_symbolic_input_names(self):
        with TemporaryDirectory() as temp:
            path = Path(temp)
            inputs = {x: path / (x + '.csv') for x in ('skin', 'cardiac', 'melanoma')}
            for p in inputs.values(): p.write_text('metadata')
            splits = create_splits(collect_samples(self.skin, self.cardiac, self.melanoma))
            manifest = write_splits(splits, inputs, path / 'out')
            self.assertEqual(set(manifest['input_metadata_sha256']), set(inputs))
            self.assertNotIn(temp, json.dumps(manifest))
            self.assertEqual(len(manifest['files']), len(splits))
            with self.assertRaises(ValueError):
                write_splits(splits, inputs, path / 'out')

    def test_single_source_folds_do_not_require_other_metadata(self):
        skin = create_splits(collect_samples(skin_rows=self.skin))
        self.assertEqual(len(skin), 9)
        for item in skin.values():
            self.assertEqual({r['source_id'] for r in item['assignments']}, {'GSE163973'})
            self.assertEqual({r['partition'] for r in item['assignments']}, {'train', 'validation', 'test'})
        cardiac = create_splits(collect_samples(cardiac_rows=self.cardiac))
        self.assertEqual(len(cardiac), 4)
        for item in cardiac.values():
            self.assertEqual({r['source_id'] for r in item['assignments']}, {'GSE144424'})
            self.assertEqual({r['partition'] for r in item['assignments']}, {'train', 'validation', 'test'})
        with self.assertRaisesRegex(ValueError, 'only two lineages'):
            create_splits(collect_samples(melanoma_rows=self.melanoma))
        with self.assertRaisesRegex(ValueError, 'source is empty'):
            collect_samples(skin_rows=self.skin, cardiac_rows=[])

    def test_each_split_pins_only_its_own_sources(self):
        with TemporaryDirectory() as temp:
            root = Path(temp)
            inputs = {x: root / (x + '.csv') for x in ('skin', 'cardiac', 'melanoma')}
            for name, path in inputs.items(): path.write_text(name)
            splits = create_splits(collect_samples(self.skin, self.cardiac, self.melanoma))
            manifest = write_splits(splits, inputs, root / 'out')
            for entry in manifest['files']:
                payload = json.loads((root / 'out' / entry['file']).read_text())
                source_map = {'GSE163973': 'skin', 'GSE144424': 'cardiac', 'GSE95280_recount3': 'melanoma'}
                expected = {source_map[row['source_id']] for row in payload['assignments']}
                self.assertEqual(set(entry['input_metadata_sha256']), expected)
                self.assertEqual(entry['input_metadata_sha256'], payload['input_metadata_sha256'])
                self.assertTrue(all(manifest['input_metadata_sha256'][key] == digest
                                    for key, digest in entry['input_metadata_sha256'].items()))


if __name__ == '__main__':
    unittest.main()
