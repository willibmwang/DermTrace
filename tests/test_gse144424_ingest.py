"""Integrity checks for the independent GSE144424 curated package."""

from __future__ import annotations

import gzip
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from ingest_gse144424 import (  # noqa: E402
    ANNOTATION_COLUMNS, DEFAULT_PACKAGE, GENE_MAPPING, INPUT_URLS, curate, sample_table, sha256,
)


class GSE144424IngestionTests(unittest.TestCase):
    def test_native_counts_and_annotations_are_preserved_exactly(self):
        raw = pd.read_csv(
            DEFAULT_PACKAGE / "inputs/GSE144424_Counts_RNA_MCW_NEB.txt.gz", sep="\t"
        )
        counts = pd.read_parquet(DEFAULT_PACKAGE / "counts.parquet")
        expected = raw.set_index("Geneid").drop(columns=ANNOTATION_COLUMNS[1:])
        pd.testing.assert_frame_equal(counts, expected)
        pd.testing.assert_frame_equal(
            pd.read_csv(DEFAULT_PACKAGE / "genes.csv"), raw[ANNOTATION_COLUMNS]
        )

    def test_sample_pairing_and_source_exclusions_remain_explicit(self):
        samples = pd.read_csv(DEFAULT_PACKAGE / "samples.csv", dtype={"donor_id": str})
        self.assertEqual(samples.shape[0], 84)
        self.assertEqual(samples.geo_accession.nunique(), 84)
        self.assertEqual(samples.donor_id.nunique(), 15)
        self.assertEqual(samples.experiment_id.nunique(), 21)
        primary = samples.loc[samples.published_primary_eligible]
        self.assertEqual(len(primary), 60)
        self.assertTrue(primary.groupby("donor_id").experiment_id.nunique().eq(1).all())
        self.assertTrue(primary.groupby("donor_id").condition.nunique().eq(4).all())
        selected_replicates = primary.groupby("donor_id").replicate.first().to_dict()
        for donor in ("18511", "18852", "18855"):
            self.assertEqual(selected_replicates[donor], 2)
        flagged = samples.set_index("sample_id")
        self.assertEqual(flagged.loc["H18852_A", "quality_flag"], "published_expression_outlier")
        self.assertEqual(flagged.loc["H18852_D", "quality_flag"], "published_rin_undetermined")
        self.assertFalse(flagged.loc[["H18852_A", "H18852_D"], "published_primary_eligible"].any())

    def test_metadata_identity_disagreement_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            altered = Path(directory) / "altered.soft.gz"
            with gzip.open(DEFAULT_PACKAGE / "inputs/GSE144424_family.soft.gz", "rt") as handle:
                text = handle.read()
            text = text.replace("subject id: 18499", "subject id: 99999", 1)
            with gzip.open(altered, "wt") as handle:
                handle.write(text)
            with self.assertRaisesRegex(ValueError, "Title/subject mismatch"):
                sample_table(altered, DEFAULT_PACKAGE / "inputs/GSE144424_series_matrix.txt.gz")

    def test_rebuild_is_deterministic_and_changed_counts_fail(self):
        with tempfile.TemporaryDirectory() as directory:
            package = Path(directory)
            inputs = package / "inputs"
            inputs.mkdir()
            for filename in INPUT_URLS:
                shutil.copyfile(DEFAULT_PACKAGE / "inputs" / filename, inputs / filename)
            first = curate(package)
            first_manifest_hash = sha256(package / "manifest.json")
            second = curate(package)
            self.assertEqual(first, second)
            self.assertEqual(first_manifest_hash, sha256(package / "manifest.json"))
            for filename, digest in second["outputs"].items():
                self.assertEqual(sha256(package / filename), digest)
            raw_path = inputs / "GSE144424_Counts_RNA_MCW_NEB.txt.gz"
            raw = pd.read_csv(raw_path, sep="\t")
            raw.loc[0, "H18499_A"] = -1
            raw.to_csv(raw_path, sep="\t", index=False, compression="gzip")
            with self.assertRaisesRegex(ValueError, "Pinned source checksum mismatch"):
                curate(package)
            self.assertEqual(first_manifest_hash, sha256(package / "manifest.json"))

    def test_valid_count_change_fails_without_replacing_existing_package(self):
        with tempfile.TemporaryDirectory() as directory:
            package = Path(directory)
            inputs = package / "inputs"
            inputs.mkdir()
            for filename in INPUT_URLS:
                shutil.copyfile(DEFAULT_PACKAGE / "inputs" / filename, inputs / filename)
            before = curate(package)
            manifest_hash = sha256(package / "manifest.json")
            raw_path = inputs / "GSE144424_Counts_RNA_MCW_NEB.txt.gz"
            raw = pd.read_csv(raw_path, sep="\t")
            raw.loc[0, "H18499_A"] += 1  # Still a finite nonnegative integer.
            raw.to_csv(raw_path, sep="\t", index=False, compression="gzip")
            with self.assertRaisesRegex(ValueError, "Pinned source checksum mismatch"):
                curate(package)
            self.assertEqual(sha256(package / "manifest.json"), manifest_hash)
            for name, expected in before["outputs"].items():
                self.assertEqual(sha256(package / name), expected)

    def test_changed_mapping_or_existing_output_is_not_replaced(self):
        with tempfile.TemporaryDirectory() as directory:
            package = Path(directory) / "package"
            inputs = package / "inputs"
            inputs.mkdir(parents=True)
            for filename in INPUT_URLS:
                shutil.copyfile(DEFAULT_PACKAGE / "inputs" / filename, inputs / filename)
            mapping = Path(directory) / "mapping.tsv.gz"
            shutil.copyfile(GENE_MAPPING, mapping)
            curate(package, mapping_path=mapping)
            original_manifest = sha256(package / "manifest.json")
            with mapping.open("ab") as handle:
                handle.write(b"changed")
            with self.assertRaisesRegex(ValueError, "Pinned gene mapping checksum mismatch"):
                curate(package, mapping_path=mapping)
            self.assertEqual(sha256(package / "manifest.json"), original_manifest)
            shutil.copyfile(GENE_MAPPING, mapping)
            with (package / "genes.csv").open("a") as handle:
                handle.write("changed")
            changed_hash = sha256(package / "genes.csv")
            with self.assertRaisesRegex(ValueError, "Existing package output changed"):
                curate(package, mapping_path=mapping)
            self.assertEqual(sha256(package / "genes.csv"), changed_hash)
            self.assertEqual(sha256(package / "manifest.json"), original_manifest)

    def test_outputs_without_manifest_are_not_overwritten(self):
        with tempfile.TemporaryDirectory() as directory:
            package = Path(directory)
            inputs = package / "inputs"
            inputs.mkdir()
            for filename in INPUT_URLS:
                shutil.copyfile(DEFAULT_PACKAGE / "inputs" / filename, inputs / filename)
            (package / "genes.csv").write_text("previous output")
            with self.assertRaisesRegex(ValueError, "lack a manifest"):
                curate(package)
            self.assertEqual((package / "genes.csv").read_text(), "previous output")

    def test_manifest_exposes_chromosome_discrepancy_and_verified_markers(self):
        manifest = json.loads((DEFAULT_PACKAGE / "manifest.json").read_text())
        self.assertFalse(manifest["frozen_benchmark_member"])
        self.assertFalse(manifest["audit"]["geo_autosomal_only_claim_matches_file"])
        self.assertEqual(manifest["audit"]["autosomal_gene_count"], 18226)
        self.assertEqual(manifest["audit"]["nonautosomal_gene_count"], 855)
        marker_mapping = pd.read_csv(DEFAULT_PACKAGE / "marker_mapping.csv")
        source = pd.read_csv(ROOT / marker_mapping.mapping_source.iloc[0], sep="\t")
        for row in marker_mapping.itertuples():
            self.assertIn(
                row.ensembl_gene_id, source.loc[source.gene_symbol.eq(row.gene_symbol), "ensembl_gene_id"].tolist()
            )
        self.assertTrue(marker_mapping.present_in_counts.all())

    def test_relocated_checkouts_produce_identical_packages(self):
        with tempfile.TemporaryDirectory() as directory:
            manifests = []
            for name in ("first-checkout", "relocated-checkout"):
                checkout = Path(directory) / name
                script = checkout / "scripts/ingest_gse144424.py"
                script.parent.mkdir(parents=True)
                shutil.copyfile(ROOT / "scripts/ingest_gse144424.py", script)
                mapping = checkout / GENE_MAPPING.relative_to(ROOT)
                mapping.parent.mkdir(parents=True)
                shutil.copyfile(GENE_MAPPING, mapping)
                package = checkout / DEFAULT_PACKAGE.relative_to(ROOT)
                inputs = package / "inputs"
                inputs.mkdir(parents=True)
                for filename in INPUT_URLS:
                    shutil.copyfile(DEFAULT_PACKAGE / "inputs" / filename, inputs / filename)
                subprocess.run(
                    [sys.executable, str(script)], cwd=checkout, check=True,
                    capture_output=True, text=True,
                )
                manifests.append(json.loads((package / "manifest.json").read_text()))
                provenance = json.loads((package / "provenance.json").read_text())
                source = provenance["local_marker_mapping"]
                self.assertEqual(source["path"], GENE_MAPPING.relative_to(ROOT).as_posix())
                self.assertEqual(source["path_basis"], "repository_root")
                for filename in ("provenance.json", "marker_mapping.csv", "README.md"):
                    self.assertNotIn(str(checkout), (package / filename).read_text())
            self.assertEqual(manifests[0], manifests[1])

    def test_external_mapping_provenance_is_portable_and_hash_identified(self):
        with tempfile.TemporaryDirectory() as directory:
            manifests = []
            for name in ("first-location", "relocated"):
                location = Path(directory) / name
                location.mkdir()
                mapping = location / "custom-annotation.tsv.gz"
                shutil.copyfile(GENE_MAPPING, mapping)
                package = location / "package"
                inputs = package / "inputs"
                inputs.mkdir(parents=True)
                for filename in INPUT_URLS:
                    shutil.copyfile(DEFAULT_PACKAGE / "inputs" / filename, inputs / filename)
                manifests.append(curate(package, mapping_path=mapping))
                provenance = json.loads((package / "provenance.json").read_text())
                self.assertEqual(provenance["local_marker_mapping"], {
                    "path": mapping.name, "path_basis": "external_input_filename",
                    "sha256": sha256(mapping),
                })
                marker_mapping = pd.read_csv(package / "marker_mapping.csv")
                self.assertTrue(marker_mapping.mapping_source.eq(mapping.name).all())
                self.assertTrue(marker_mapping.mapping_source_basis.eq("external_input_filename").all())
                for filename in ("provenance.json", "marker_mapping.csv", "README.md"):
                    self.assertNotIn(str(location), (package / filename).read_text())
            self.assertEqual(manifests[0], manifests[1])


if __name__ == "__main__":
    unittest.main()
