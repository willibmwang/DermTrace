"""Scientific invariants for donor-level oxygen comparisons."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

from scripts.analyze_gse144424 import DATA, OUT, ROOT, autosomal_mask, bh_adjust, normalize_counts, paired_deltas, sha256


class DonorAnalysisTests(unittest.TestCase):
    def test_sequencing_depth_does_not_create_effect(self):
        counts = pd.DataFrame({"A": [10, 30, 100], "B": [20, 60, 200]}, index=["g1", "g2", "g3"])
        normalized, factors, _ = normalize_counts(counts, pd.Series(True, index=counts.index), "median_ratio")
        np.testing.assert_allclose(normalized["A"], normalized["B"])
        self.assertAlmostEqual(factors["B"] / factors["A"], 2)

    def test_pairs_stay_in_experiment_and_repeats_do_not_inflate_donors(self):
        columns, values, metadata = [], [], []
        for donor, repeat, baseline, change in [("d1", 1, 1, 2), ("d1", 2, 100, 4), ("d2", 1, 20, -1)]:
            for condition, value in zip("ABCD", [baseline, baseline + change, baseline + 1, baseline]):
                sample = f"{donor}_{repeat}_{condition}"
                columns.append(sample)
                values.append(value)
                metadata.append(dict(sample_id=sample, donor_id=donor, experiment_id=f"{donor}_{repeat}", condition=condition))
        expressions = pd.DataFrame([values], index=["gene"], columns=columns)
        _, donors = paired_deltas(expressions, pd.DataFrame(metadata))
        primary = donors.query("contrast == 'B-A'").set_index("donor_id")
        self.assertEqual(len(primary), 2)
        self.assertEqual(primary.loc["d1", "delta"], 3)
        self.assertEqual(primary.loc["d2", "delta"], -1)
        self.assertEqual(primary.loc["d1", "n_experiments"], 2)
        with self.assertRaisesRegex(ValueError, "Incomplete"):
            paired_deltas(expressions.iloc[:, :-1], pd.DataFrame(metadata[:-1]))

    def test_bh_preserves_original_marker_order(self):
        np.testing.assert_allclose(bh_adjust(np.array([.04, .001, .03, .2])), [.0533333333, .004, .0533333333, .2])

    def test_autosome_detection_uses_actual_featurecounts_annotation(self):
        observed = autosomal_mask(pd.Series(["chr1;chr1", "chrX", "chr2;chr3", "chr22", "chrM"]))
        self.assertEqual(observed.tolist(), [True, False, False, True, False])

    def test_ac16_reference_reconstructs_from_tracked_sources(self):
        provenance = json.loads((OUT / "reference/provenance.json").read_text())
        for source in [provenance["snapshot"], *provenance["source_files"]]:
            self.assertEqual(sha256(ROOT / source["path"]), source["sha256"])
        snapshot = pd.read_csv(ROOT / provenance["snapshot"]["path"])
        tpm = pd.read_csv(ROOT / provenance["source_files"][0]["path"], sep="\t", index_col=0)
        metadata = pd.read_csv(ROOT / provenance["source_files"][1]["path"])
        metadata = metadata.loc[metadata.cell_line.eq("AC16")]
        self.assertEqual(len(metadata), 12)
        for arm, sample_ids in provenance["selection"]["arms"].items():
            self.assertEqual(len(sample_ids), 3)
            self.assertEqual(set(sample_ids), set(metadata.loc[metadata.oxygen.eq(arm), "sample_id"]))
        self.assertEqual(len(snapshot), 32)
        self.assertFalse(snapshot.duplicated(["gene", "oxygen"]).any())
        for row in snapshot.itertuples():
            self.assertTrue(row.present)
            mean = tpm.loc[row.gene, provenance["selection"]["arms"][row.oxygen]].mean()
            baseline = tpm.loc[row.gene, provenance["selection"]["arms"]["20%"]].mean()
            np.testing.assert_allclose([row.mean_tpm, row.log2_ratio_mean_tpm_plus1_vs_20pct],
                                       [mean, np.log2((mean + 1) / (baseline + 1))], rtol=1e-12, atol=1e-12)

    def test_full_analysis_rebuilds_from_only_bundled_inputs_after_relocation(self):
        """Exercise hashing and the AC16 dependency, including on Python 3.10 CI."""
        with tempfile.TemporaryDirectory() as directory:
            checkout = Path(directory) / "relocated-checkout"
            data = checkout / DATA.relative_to(ROOT)
            output = checkout / OUT.relative_to(ROOT)
            data.mkdir(parents=True)
            output.mkdir(parents=True)
            script = checkout / "scripts/analyze_gse144424.py"
            script.parent.mkdir()
            shutil.copy2(ROOT / script.relative_to(checkout), script)
            for name in ("counts.parquet", "samples.csv", "genes.csv", "marker_mapping.csv", "manifest.json"):
                shutil.copy2(DATA / name, data / name)
            shutil.copy2(OUT / "analysis_plan.json", output / "analysis_plan.json")
            shutil.copytree(OUT / "reference", output / "reference")
            result = subprocess.run([sys.executable, str(script)], cwd=directory,
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            manifest = json.loads((output / "analysis_manifest.json").read_text())
            self.assertEqual(set(manifest["primary_direction_support_genes"]),
                             {"EGLN3", "SLC2A1", "VEGFA", "NDRG1"})
            for section in ("input_sha256", "output_sha256"):
                for relative, digest in manifest[section].items():
                    self.assertFalse(Path(relative).is_absolute())
                    self.assertEqual(sha256(checkout / relative), digest)
            for original in OUT.glob("*.csv"):
                pd.testing.assert_frame_equal(pd.read_csv(original), pd.read_csv(output / original.name),
                                              check_exact=False, rtol=1e-10, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
