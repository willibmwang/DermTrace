"""Synthetic archive tests; no source expression or reads are required."""
import gzip
import io
from pathlib import Path
import tarfile
from tempfile import TemporaryDirectory
import unittest

from experiments.expression_completion.gse163973_reader import iter_matrices


def add(tar, name, data):
    info = tarfile.TarInfo(name)
    info.size = len(data)
    tar.addfile(info, io.BytesIO(data))


def archive(path, *, features=b'ENSG1\tA\tGene Expression\nENSG2\tB\tGene Expression\n',
            barcodes=b'cell1\ncell2\n', matrix=b'%%MatrixMarket matrix coordinate integer general\n2 2 2\n1 1 1\n2 2 2\n',
            missing=None, duplicate=None, duplicate_sample=False):
    data = io.BytesIO()
    with tarfile.open(fileobj=data, mode='w:gz') as nested:
        for name, content in [('features.tsv.gz', features), ('barcodes.tsv.gz', barcodes), ('matrix.mtx.gz', matrix)]:
            if name != missing:
                add(nested, 'matrix/' + name, gzip.compress(content))
            if name == duplicate:
                add(nested, 'other/' + name, gzip.compress(content))
    with tarfile.open(path, 'w') as outer:
        add(outer, 'GSM4994379_KL1_matrix.tar.gz', data.getvalue())
        if duplicate_sample:
            add(outer, 'GSM9999999_KL1_matrix.tar.gz', data.getvalue())


class ReaderTests(unittest.TestCase):
    def setUp(self):
        self.tmp = TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name) / 'processed.tar'

    def test_preserves_feature_id_symbol_and_counts(self):
        archive(self.path)
        sample, features, barcodes, matrix = next(iter_matrices(self.path, ['KF1']))
        self.assertEqual(sample, 'KF1')
        self.assertEqual(features.iloc[:, :2].values.tolist(), [['ENSG1', 'A'], ['ENSG2', 'B']])
        self.assertEqual(barcodes.tolist(), ['cell1', 'cell2'])
        self.assertEqual(matrix.toarray().tolist(), [[1, 0], [0, 2]])

    def test_missing_or_duplicate_component_and_sample_fail(self):
        for kwargs in ({'missing': 'features.tsv.gz'}, {'duplicate': 'matrix.mtx.gz'}, {'duplicate_sample': True}):
            archive(self.path, **kwargs)
            with self.assertRaises(ValueError):
                list(iter_matrices(self.path, ['KF1']))
        archive(self.path)
        with self.assertRaises(ValueError):
            list(iter_matrices(self.path, ['KF2']))

    def test_duplicate_ids_barcode_or_invalid_counts_fail(self):
        for kwargs in ({'features': b'ENSG1\tA\nENSG1\tB\n'},
                       {'barcodes': b'cell1\ncell1\n'},
                       {'matrix': b'%%MatrixMarket matrix coordinate integer general\n2 2 1\n1 1 -1\n'},
                       {'matrix': b'%%MatrixMarket matrix coordinate integer general\n2 2 2\n1 1 -1\n1 1 2\n'},
                       {'matrix': b'%%MatrixMarket matrix coordinate integer general\n2 2 2\n1 1 1\n1 1 2\n'}):
            archive(self.path, **kwargs)
            with self.assertRaises(ValueError):
                list(iter_matrices(self.path, ['KF1']))

    def test_wrong_dimensions_fail(self):
        archive(self.path, matrix=b'%%MatrixMarket matrix coordinate integer general\n3 2 1\n1 1 1\n')
        with self.assertRaises(ValueError):
            list(iter_matrices(self.path, ['KF1']))


if __name__ == '__main__':
    unittest.main()
