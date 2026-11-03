"""Source, mapping and real-reader checks for the processed GSE145725 release."""
from __future__ import annotations

import gzip
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "keloidbench/scripts"))
import ingest_gse145725 as ingest
from expression_processing import METADATA_COLUMNS


class GSE145725Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.metadata, cls.native, cls.mapping = ingest.prepare(ingest.PACKAGE / "source")
        cls.wide = pd.read_parquet(ingest.PACKAGE / "microarray_expression_wide.parquet")
        cls.genes = (ingest.PACKAGE / "microarray_gene_vocab.txt").read_text().splitlines()

    def test_source_identity_labels_and_replication(self):
        self.assertEqual(self.metadata.sample_id.tolist(), ingest.SAMPLE_IDS)
        self.assertEqual(self.metadata.disease_label.value_counts().to_dict(), {"normal": 10, "keloid": 9})
        sizes = self.metadata.groupby("source_group_id").size().to_dict()
        self.assertEqual(len(sizes), 10)
        self.assertEqual(sizes.pop("keloid_fibroblast_3"), 1)
        self.assertEqual(set(sizes.values()), {2})
        self.assertNotIn("keloid_fibroblast_2", self.metadata.source_group_id.tolist())
        self.assertEqual(set(self.metadata.patient_id_basis), {"deposited_source_label_not_verified_donor"})

    def test_exposure_roles_and_boolean_types(self):
        self.assertEqual(set(self.metadata.oxygen), {"20%"})
        self.assertEqual(set(self.metadata.passage), {3})
        self.assertEqual(set(self.metadata.growth_medium_hydrocortisone_ug_ml), {0.4})
        self.assertEqual(set(self.metadata.cohort_role), {"fibroblast_sensitivity_only"})
        self.assertEqual(self.wide.lockbox.dtype, bool)
        self.assertFalse(self.wide.lockbox.any())
        self.assertFalse(self.wide.allow_accession_level_grouping.any())

    def test_ambiguous_and_invalid_symbols_are_not_manufactured(self):
        counts = self.mapping.mapping_status.value_counts().to_dict()
        self.assertEqual(counts, {"retained": 47239, "ambiguous": 1457, "unmapped": 688, "invalid_symbol": 111})
        self.assertTrue(self.mapping.loc[self.mapping.mapping_status.ne("retained"), "gene_symbol"].eq("").all())
        historical = self.mapping[self.mapping.original_gene_symbol.eq("C17orf78")]
        self.assertTrue(len(historical) > 0)
        self.assertEqual(set(historical.gene_symbol), {"C17ORF78"})
        self.assertNotIn("41154", self.genes)
        self.assertNotIn("IGL", self.genes)  # Do not turn IGL@ into a different symbol.

    def test_native_values_against_independent_probe_aggregation(self):
        # Read the deposited table independently of DermTrace's SeriesMatrix parser.
        text = gzip.decompress((ingest.PACKAGE / "source/GSE145725_series_matrix.txt.gz").read_bytes()).decode()
        text = text.split("!series_matrix_table_begin", 1)[1].split("!series_matrix_table_end", 1)[0].strip()
        probes = pd.read_csv(io.StringIO(text), sep="\t", index_col=0)
        mapping = self.mapping[self.mapping.mapping_status.eq("retained")].set_index("probe_id").gene_symbol
        rows = probes.loc[mapping.index].copy()
        rows.index = mapping.values
        # Independent numpy accumulator checks ALL 361,209 native gene/sample values.
        sums = np.zeros((len(self.genes), 19), dtype=float)
        counts = np.zeros(len(self.genes), dtype=int)
        lookup = {gene: i for i, gene in enumerate(self.genes)}
        for gene, row in zip(rows.index, rows.to_numpy()):
            idx = lookup[gene]
            sums[idx] += row
            counts[idx] += 1
        expected = (sums / counts[:, None]).T
        np.testing.assert_allclose(self.native.to_numpy(), expected, atol=1e-12, rtol=0)
        # Deposited-scale averages must remain inside the observed log2 RMA range.
        self.assertGreaterEqual(self.native.to_numpy().min(), 1.79995)
        self.assertLessEqual(self.native.to_numpy().max(), 14.1429)

    def test_standard_artifacts_alignment_and_scale(self):
        report = ingest.validate_release(ingest.PACKAGE, self.metadata, self.native)
        self.assertEqual(report["n_genes"], 19011)
        self.assertEqual(self.genes, sorted(set(self.genes)))
        for col in METADATA_COLUMNS:
            self.assertIn(col, self.wide)
        values = self.wide[self.genes].to_numpy()
        np.testing.assert_allclose(values.mean(axis=0), 0, atol=1e-12)
        np.testing.assert_allclose(values.std(axis=0, ddof=1), 1, atol=1e-12)

    def test_jsonl_is_aligned_subset_of_standardized_matrix(self):
        rows = [json.loads(line) for line in (ingest.PACKAGE / "microarray_encoder_decoder_samples.jsonl").read_text().splitlines()]
        self.assertEqual([row["sample_id"] for row in rows], ingest.SAMPLE_IDS)
        for idx, row in enumerate(rows):
            self.assertEqual(row["response"], self.metadata.iloc[idx].encoder_response)
            self.assertLessEqual(len(row["genes"]), 2048)
            for gene, value in row["genes"].items():
                self.assertEqual(value, self.wide.at[idx, gene])

    def test_real_corpus_reader_accepts_package_without_registry_edit(self):
        import build_training_corpus as corpus
        original = list(corpus.ARTIFACTS)
        artifact = corpus.ModalityArtifact("microarray", "microarray", ingest.PACKAGE)
        with patch.object(corpus, "ARTIFACTS", [artifact]):
            manifest, features, vocab, expression = corpus.load_artifacts(max_shared_genes=100)
        self.assertEqual(corpus.ARTIFACTS, original)
        self.assertEqual(manifest.sample_id.tolist(), ingest.SAMPLE_IDS)
        self.assertEqual(manifest.split_group.nunique(), 10)
        self.assertEqual(len(manifest), 19)
        self.assertEqual(set(manifest.keloid_vs_unaffected_skin), {"exclude"})
        self.assertEqual(vocab["modality_gene_counts"], {"microarray": 19011})
        np.testing.assert_array_equal(expression["microarray"].to_numpy(), self.wide[self.genes].to_numpy())
        self.assertTrue(features)

    def test_pinned_source_change_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            for name in ingest.INPUTS:
                (directory/name).symlink_to(ingest.PACKAGE / "source" / name)
            target = directory / "GSE145725_series_matrix.txt.gz"
            target.unlink()
            target.write_bytes(b"changed input")
            with self.assertRaisesRegex(ValueError, "changed pinned input"):
                ingest.prepare(directory)

    def test_release_refuses_overwrite_and_preserves_existing_bytes(self):
        manifest = ingest.PACKAGE / "manifest.json"
        before = ingest.sha256(manifest)
        with self.assertRaises(FileExistsError):
            ingest.build(ingest.PACKAGE / "source", ingest.PACKAGE)
        self.assertEqual(ingest.sha256(manifest), before)

    def test_existing_output_source_must_match_pinned_inputs(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp)
            (output / "source").mkdir()
            (output / "source/GSE145725_series_matrix.txt.gz").write_bytes(b"wrong source")
            with self.assertRaisesRegex(ValueError, "changed pinned input"):
                ingest.build(ingest.PACKAGE / "source", output)
            self.assertFalse((output / "manifest.json").exists())

    def test_manifest_hashes_match_every_saved_output(self):
        manifest = json.loads((ingest.PACKAGE / "manifest.json").read_text())
        for name, receipt in manifest["outputs"].items():
            path = ingest.PACKAGE / name
            self.assertEqual(ingest.sha256(path), receipt["sha256"], name)
            self.assertEqual(path.stat().st_size, receipt["bytes"], name)
        for name, digest in manifest["code"].items():
            self.assertEqual(ingest.sha256(REPO / "keloidbench/scripts" / name), digest, name)

    def test_rebuild_is_relocatable_and_numerically_identical(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "release"
            ingest.build(ingest.PACKAGE / "source", output)
            ingest.validate_release(output, self.metadata, self.native)
            summary = json.loads((output / "microarray_summary.json").read_text())
            self.assertTrue(all(not Path(value).is_absolute() for value in summary["outputs"].values()))
            np.testing.assert_array_equal(np.load(output / "microarray_X.npy"), np.load(ingest.PACKAGE / "microarray_X.npy"))


if __name__ == "__main__":
    unittest.main()
