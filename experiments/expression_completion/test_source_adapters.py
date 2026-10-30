"""Small processed-source fixtures exercise the actual adapter boundaries."""
import gzip
import hashlib
import io
import json
import tarfile
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmwrite
from scipy.sparse import csc_matrix

from experiments.expression_completion.loader import (
    SplitAssignments, cardiac_profiles, skin_pseudobulk_profiles,
    skin_cells, melanoma_profiles,
)

GENES = ('ENSG00000000001', 'ENSG00000000002', 'ENSG00000000003')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify(path, expected):
    if digest(path) != expected:
        raise ValueError('Source checksum mismatch')


def put_tar(tar, name, payload):
    item = tarfile.TarInfo(name)
    item.size = len(payload)
    tar.addfile(item, io.BytesIO(payload))


def mapping(source, assignments):
    return SplitAssignments([dict(source_id=source, sample_id=s, group_id=g, partition=p)
                             for s, g, p in assignments])


class SourceAdapterTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def test_cardiac_filters_flags_and_detects_split_donor(self):
        d = self.root
        frame = pd.DataFrame({'s1':[0.,2.,3.], 's2':[4.,5.,6.], 's3':[3.,3.,3.]}, index=GENES)
        frame.to_parquet(d/'counts.parquet')
        pd.DataFrame([
            dict(sample_id='s1', donor_id='D1', published_primary_eligible=True, quality_flag='no_specific_paper_flag'),
            dict(sample_id='s2', donor_id='D1', published_primary_eligible=False, quality_flag='no_specific_paper_flag'),
            dict(sample_id='s3', donor_id='D2', published_primary_eligible=False, quality_flag='published_expression_outlier'),
        ]).to_csv(d/'samples.csv', index=False)
        (d/'manifest.json').write_text(json.dumps({'outputs':{x:digest(d/x) for x in ('counts.parquet','samples.csv')}}))
        split = mapping('GSE144424', [('s1','d1','train'), ('s2','d1','train'), ('s3','d2','test')])
        primary = list(cardiac_profiles(d, split, verify))
        self.assertEqual([p.sample_id for p in primary], ['s1'])
        self.assertEqual(primary[0].values[0], 0.)
        self.assertTrue(primary[0].measured[0])
        self.assertEqual(len(list(cardiac_profiles(d, split, verify, primary_only=False))), 2)
        bad = mapping('GSE144424', [('s1','d1','train'), ('s2','other','test'), ('s3','d2','test')])
        with self.assertRaisesRegex(ValueError, 'donor crosses groups'):
            list(cardiac_profiles(d, bad, verify))
        with (d/'samples.csv').open('a') as stream:
            stream.write('\n')
        with self.assertRaisesRegex(ValueError, 'checksum'):
            list(cardiac_profiles(d, split, verify))

    def test_skin_pseudobulk_preserves_source_groups_and_omits_aggregates(self):
        d = self.root
        pd.DataFrame({'A':[0,8], 'B':[3,9], 'C':[7,10]}, index=['p1','p2']).to_parquet(d/'pseudobulk_sum_counts.parquet')
        pd.DataFrame([dict(sample_id='p1', source_donor_key='KF1', cell_type='Fibroblast'),
                      dict(sample_id='p2', source_donor_key='KF2', cell_type='Endothelial')]).to_parquet(d/'metadata.parquet')
        pd.DataFrame([dict(source_gene_id=GENES[0],gene_symbol='A',symbol_aggregates_multiple_features=False),
                      dict(source_gene_id=GENES[1],gene_symbol='B',symbol_aggregates_multiple_features=False),
                      dict(source_gene_id=GENES[2],gene_symbol='C',symbol_aggregates_multiple_features=True)]).to_parquet(d/'gene_annotations.parquet')
        names=('pseudobulk_sum_counts.parquet','metadata.parquet','gene_annotations.parquet')
        (d/'recovery_manifest.json').write_text(json.dumps({'output_sha256':{x:digest(d/x) for x in names}}))
        split=mapping('GSE163973',[('KF1','g1','train'),('KF2','g2','test')])
        profiles=list(skin_pseudobulk_profiles(d,split,verify))
        self.assertEqual([(p.sample_id,p.partition) for p in profiles],[('KF1','train'),('KF2','test')])
        self.assertEqual(profiles[0].genes, GENES[:2])
        np.testing.assert_array_equal(profiles[0].values,[0,3])
        self.assertEqual(len(list(skin_pseudobulk_profiles(d,split,verify,samples=['KF1']))),1)

    def test_skin_cells_uses_processed_archive_and_excludes_unannotated(self):
        d = self.root
        (d/'downloads').mkdir()
        pd.DataFrame({'source_gene_id':GENES}).to_parquet(d/'gene_annotations.parquet')
        pd.DataFrame([
            dict(source_donor_key='KF1', matrix_barcode='cell-1', source_metadata_row=0, cellType='Fibroblast'),
            dict(source_donor_key='KF1', matrix_barcode='empty-1', source_metadata_row=None, cellType=None),
        ]).to_parquet(d/'cell_barcode_crosswalk.parquet')
        mtx=io.BytesIO()
        mmwrite(mtx,csc_matrix(np.array([[0,7],[2,0],[4,0]],dtype=np.int64)))
        inner_data=io.BytesIO()
        with tarfile.open(fileobj=inner_data,mode='w:gz') as inner:
            put_tar(inner,'features.tsv.gz',gzip.compress(('\n'.join(g+'\tG\tGene Expression' for g in GENES)+'\n').encode()))
            put_tar(inner,'barcodes.tsv.gz',gzip.compress(b'cell-1\nempty-1\n'))
            put_tar(inner,'matrix.mtx.gz',gzip.compress(mtx.getvalue()))
        archive=d/'downloads/GSE163973_RAW.tar'
        with tarfile.open(archive,'w:') as outer:
            put_tar(outer,'GSM000_KL1_matrix.tar.gz',inner_data.getvalue())
        names=('gene_annotations.parquet','cell_barcode_crosswalk.parquet')
        (d/'recovery_manifest.json').write_text(json.dumps({'source_files':{archive.name:{'sha256':digest(archive)}},
            'output_sha256':{x:digest(d/x) for x in names}}))
        split=mapping('GSE163973',[('KF1','g1','test')])
        p=list(skin_cells(d,split,verify,lambda *_:True,samples=['KF1']))
        self.assertEqual(len(p),1)
        self.assertEqual((p[0].sample_id,p[0].group_id,p[0].partition),('KF1','g1','test'))
        self.assertEqual(p[0].cell_type,'Fibroblast')
        np.testing.assert_array_equal(p[0].values,[0,2,4])
        self.assertTrue(p[0].measured.all())
        self.assertEqual(list(skin_cells(d,split,verify,lambda *_:False,samples=['KF1'])),[])

    def test_melanoma_filters_perturbed_profiles_and_checks_lineage(self):
        d=self.root
        pd.DataFrame({'A':[0.,3.,5.], 'B':[2.,4.,6.], 'C':[8.,9.,10.]}, index=['m1','m2','m3']).to_parquet(d/'expression_native_symbol.parquet')
        pd.DataFrame([dict(sample_id='m1', clean_hypoxia_eligible=True, cell_line_family='family'),
                      dict(sample_id='m2', clean_hypoxia_eligible=True, cell_line_family='family'),
                      dict(sample_id='m3', clean_hypoxia_eligible=False, cell_line_family='family')]).to_parquet(d/'metadata.parquet')
        names=('expression_native_symbol.parquet','metadata.parquet')
        (d/'recount_source_manifest.json').write_text(json.dumps({'output_sha256':{x:digest(d/x) for x in names}}))
        def annotations(_):
            return {'gene_axis':'columns','species':'human','genes':[
                dict(mapping_basis='source_annotation', ensembl_gene_ids=[g], source_ensembl_ids=[g], native_label=s)
                for g,s in zip(GENES[:2],('A','B'))] + [dict(mapping_basis='source_annotation', ensembl_gene_ids=[GENES[2], 'ENSG00000000004'], source_ensembl_ids=[GENES[2]],native_label='C')]}
        split=mapping('GSE95280_recount3',[('m1','f','train'),('m2','f','train'),('m3','f','train')])
        p=list(melanoma_profiles(d,split,verify,annotations))
        self.assertEqual([x.sample_id for x in p],['m1','m2'])
        self.assertEqual(p[0].genes,GENES[:2])
        np.testing.assert_array_equal(p[0].values,[0,2])
        self.assertEqual(len(list(melanoma_profiles(d,split,verify,annotations,clean_only=False))),3)
        bad=mapping('GSE95280_recount3',[('m1','f1','train'),('m2','f2','test'),('m3','f1','train')])
        with self.assertRaisesRegex(ValueError,'lineage crosses groups'):
            list(melanoma_profiles(d,bad,verify,annotations))


if __name__ == '__main__':
    unittest.main()
