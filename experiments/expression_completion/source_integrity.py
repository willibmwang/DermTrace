"""Verify immutable source files, including recorded annotation-only updates.

The historic expected hash must come from a trusted source manifest. A full-file
match is accepted directly. Otherwise the exact annotated bytes AND original
prefix must match a separately supplied annotation receipt. No files are changed.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path


def sha256_file(path: str | Path, limit: int | None = None) -> str:
    if limit is not None and (not isinstance(limit, int) or isinstance(limit, bool) or limit <= 0):
        raise ValueError("Prefix length must be a positive integer")
    digest = hashlib.sha256()
    remaining = limit
    with Path(path).open("rb") as stream:
        while remaining is None or remaining > 0:
            block = stream.read(1024 * 1024 if remaining is None else min(1024 * 1024, remaining))
            if not block:
                if remaining:
                    raise ValueError("File is shorter than the recorded original")
                break
            digest.update(block)
            if remaining is not None:
                remaining -= len(block)
    return digest.hexdigest()


class AnnotationAwareVerifier:
    """Callable compatible with ``verify(path, expected_sha256)``.

    The receipt records a specific annotation event; it is not permission to
    accept any other metadata footer. Both the receipt and source manifest are
    provenance inputs and must be protected by the surrounding experiment.
    """

    def __init__(self, receipt_path: str | Path, *, path_overrides=None):
        """Bind optional relocated input paths to exact paths in the receipt.

        ``path_overrides`` maps the actual local path to the path recorded in
        the trusted receipt. If absent, an unambiguous pair of historical and
        annotated hashes can identify a moved file without relying on its name.
        """
        self.receipt_path = Path(receipt_path).resolve(strict=True)
        self.receipt_sha256 = sha256_file(self.receipt_path)
        receipt = json.loads(self.receipt_path.read_text())
        if receipt.get("status") != "complete" or receipt.get("complete") is not True:
            raise ValueError("Annotation receipt is not complete")
        records = receipt.get("tables")
        if not isinstance(records, list) or receipt.get("table_count") != len(records):
            raise ValueError("Annotation receipt has an invalid table list")
        self._records = {}
        for record in records:
            path = Path(record["path"]).resolve()
            if path in self._records:
                raise ValueError("Annotation receipt repeats a source path")
            self._records[path] = record
        self._overrides = {Path(actual).resolve(): Path(recorded).resolve()
                           for actual, recorded in (path_overrides or {}).items()}
        if any(path not in self._records for path in self._overrides.values()):
            raise ValueError("Path override is not present in the annotation receipt")

    def __call__(self, path: str | Path, expected_sha256: str) -> dict:
        path = Path(path).resolve(strict=True)
        if not isinstance(expected_sha256, str) or len(expected_sha256) != 64 or any(c not in "0123456789abcdef" for c in expected_sha256):
            raise ValueError("Expected SHA-256 must be a lowercase hexadecimal digest")
        current = sha256_file(path)
        if current == expected_sha256:
            return {"status": "exact", "current_sha256": current}
        record = self._records.get(self._overrides.get(path, path))
        if record is None:
            # A relocated file must match both independent recorded hashes.
            # Identical duplicate records are ambiguous and require an override.
            matches = [r for r in self._records.values()
                       if r.get("original_sha256") == expected_sha256
                       and r.get("verification", {}).get("annotated_sha256") == current]
            if len(matches) != 1:
                raise ValueError(f"Source lacks a unique annotation receipt match: {path}")
            record = matches[0]
        v = record.get("verification", {})
        if (record.get("original_sha256") != expected_sha256
                or v.get("original_sha256") != expected_sha256
                or v.get("status") != "annotated_and_verified"
                or any(v.get(k) is not True for k in ("values_axes_types_nulls_equal", "pandas_metadata_unchanged", "original_bytes_preserved", "annotation_roundtrip_exact"))):
            raise ValueError(f"Annotation receipt does not verify this historical source: {path}")
        if (v.get("original_size") != record.get("original_size")
                or v.get("annotated_size") != path.stat().st_size
                or v.get("annotated_sha256") != current):
            raise ValueError(f"Annotated source no longer matches the verified receipt: {path}")
        if sha256_file(path, record["original_size"]) != expected_sha256:
            raise ValueError(f"Original source prefix does not match its historical checksum: {path}")
        return {"status": "verified_annotation_only", "current_sha256": current,
                "original_sha256": expected_sha256, "receipt_sha256": self.receipt_sha256}
