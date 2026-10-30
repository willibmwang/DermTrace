import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from experiments.expression_completion.source_integrity import AnnotationAwareVerifier


def digest(data):
    return hashlib.sha256(data).hexdigest()


class IntegrityTests(unittest.TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        root = Path(self.temp.name)
        self.original = b"original count data"
        self.annotated = self.original + b"annotation footer"
        self.path = root / "source.parquet"
        self.path.write_bytes(self.annotated)
        v = {"status": "annotated_and_verified", "original_sha256": digest(self.original),
             "original_size": len(self.original), "annotated_size": len(self.annotated),
             "annotated_sha256": digest(self.annotated), "values_axes_types_nulls_equal": True,
             "pandas_metadata_unchanged": True, "original_bytes_preserved": True,
             "annotation_roundtrip_exact": True}
        receipt = {"status": "complete", "complete": True, "table_count": 1,
                   "tables": [{"path": str(self.path), "original_sha256": digest(self.original),
                               "original_size": len(self.original), "verification": v}]}
        receipt_path = root / "receipt.json"
        receipt_path.write_text(json.dumps(receipt))
        self.verify = AnnotationAwareVerifier(receipt_path)

    def test_exact_and_verified_annotation(self):
        self.assertEqual(self.verify(self.path, digest(self.annotated))["status"], "exact")
        self.assertEqual(self.verify(self.path, digest(self.original))["status"], "verified_annotation_only")

    def test_corruption_or_different_annotation_fails(self):
        for changed in (b"X" + self.annotated[1:], self.original + b"different footer"):
            self.path.write_bytes(changed)
            with self.assertRaises(ValueError):
                self.verify(self.path, digest(self.original))

    def test_wrong_manifest_hash_and_missing_receipt_fail(self):
        with self.assertRaises(ValueError):
            self.verify(self.path, digest(b"different original"))
        unknown = self.path.parent / "unknown"
        unknown.write_bytes(self.annotated + b"unrecorded")
        with self.assertRaises(ValueError):
            self.verify(unknown, digest(self.original))

    def test_relocated_file_uses_both_hashes(self):
        moved = self.path.parent / "moved"
        moved.write_bytes(self.annotated)
        self.assertEqual(self.verify(moved, digest(self.original))["status"], "verified_annotation_only")
        moved.write_bytes(self.original + b"modified footer")
        with self.assertRaises(ValueError):
            self.verify(moved, digest(self.original))

    def test_duplicate_receipts_require_explicit_override(self):
        receipt_path = self.path.parent / "receipt.json"
        receipt = json.loads(receipt_path.read_text())
        duplicate = dict(receipt["tables"][0], path=str(self.path.parent / "other-original"))
        receipt["tables"].append(duplicate)
        receipt["table_count"] = 2
        receipt_path.write_text(json.dumps(receipt))
        moved = self.path.parent / "moved"
        moved.write_bytes(self.annotated)
        verify = AnnotationAwareVerifier(receipt_path)
        with self.assertRaises(ValueError):
            verify(moved, digest(self.original))
        verify = AnnotationAwareVerifier(receipt_path, path_overrides={moved:self.path})
        self.assertEqual(verify(moved, digest(self.original))["status"], "verified_annotation_only")

    def test_malformed_checksum_fails(self):
        with self.assertRaises(ValueError):
            self.verify(self.path, "no")


if __name__ == "__main__":
    unittest.main()
