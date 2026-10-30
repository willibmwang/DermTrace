import hashlib
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest


class QcSafetyTests(unittest.TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.archive = self.root / 'archive.tar'
        self.crosswalk = self.root / 'crosswalk.parquet'
        self.metadata = self.root / 'metadata.parquet'
        for path in (self.archive, self.crosswalk, self.metadata):
            path.write_bytes(b'input')
        digest = hashlib.sha256(b'input').hexdigest()
        self.receipt = self.root / 'receipt.json'
        self.receipt.write_text(json.dumps({'source_files': {'GSE163973_RAW.tar': {'sha256': digest}},
                                           'output_sha256': {'cell_barcode_crosswalk.parquet': digest,
                                                             'source_cell_metadata.parquet': digest}}))
        self.output = self.root / 'output.json'

    def run_qc(self):
        command = [sys.executable, '-m', 'experiments.expression_completion.review_qc',
                   '--archive', str(self.archive), '--crosswalk', str(self.crosswalk),
                   '--cell-metadata', str(self.metadata), '--receipt', str(self.receipt),
                   '--output', str(self.output)]
        return subprocess.run(command, capture_output=True, text=True)

    def test_changed_input_fails_before_processing(self):
        self.archive.write_bytes(b'changed')
        result = self.run_qc()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('does not match', result.stderr)
        self.assertFalse(self.output.exists())

    def test_existing_output_is_not_overwritten(self):
        self.output.write_text('preserve me')
        result = self.run_qc()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('distinct new files', result.stderr)
        self.assertEqual(self.output.read_text(), 'preserve me')


if __name__ == '__main__':
    unittest.main()
