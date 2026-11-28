"""Regression tests for the Step 17 recovery (r1): SIGPIPE-safe accounting display, deterministic header parsing,
strict concordance-report validation, upstream-failure propagation and exit-status preservation by the EXIT trap.

Run from the bundle directory:  python3 -I -m unittest discover -s tests -p 'test_*.py'
"""
import copy, gzip, json, os, shutil, stat, subprocess, sys, tempfile, unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
LIB, PYH = HERE / "collect_lib_r1.sh", HERE / "collect_r1.py"
BASH = shutil.which("bash")
COLLAB33 = [f"GM{i:02d}" for i in range(1, 33)] + ["GM50"]
NINE = ["GM33", "GM34", "GM35", "GM36", "GM38", "GM39", "GM40", "GM41", "GM42"]
GATES = (["test_output_complete_and_valid", "specimen_identity", "raw_input_reads_identical", "identical_transcript_ids", "target_counts_equal",
          "library_format_matches"] + [f"index_{c}_identical" for c in ("seq_hash", "name_hash", "seq_hash512", "name_hash512", "decoy_seq_hash", "decoy_name_hash")]
         + ["transcript_pearson_log1p_tpm_ge_0.99", "abs_mapping_diff_le_1pp"])


def sh(script, env=None, timeout=120):
    e = dict(os.environ); e.update(env or {})
    return subprocess.run([BASH, "-c", script], capture_output=True, text=True, env=e, timeout=timeout)


def py(*args):
    return subprocess.run([sys.executable, "-I", str(PYH), *map(str, args)], capture_output=True, text=True)


class Sacct(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.bin = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def stub(self, body):
        p = self.bin / "sacct"; p.write_text("#!/bin/bash\n" + body); p.chmod(p.stat().st_mode | stat.S_IEXEC)
        return {"PATH": f"{self.bin}:{os.environ['PATH']}"}

    def test_old_pattern_reproduces_sigpipe(self):
        env = self.stub("seq 1 300000 | sed 's/^/23433436_/'\n")
        r = sh("set -euo pipefail; sacct -j 1 -o JobID | head -40 >/dev/null; echo unreachable", env)
        self.assertEqual(r.returncode, 141, "the original step 17 pattern should die with SIGPIPE (141)")
        self.assertNotIn("unreachable", r.stdout)

    def test_large_output_no_sigpipe(self):
        env = self.stub("seq 1 300000 | sed 's/^/23433436_/'\n")
        r = sh(f"set -euo pipefail; . '{LIB}'; show_sacct 23433436 40 || exit 1; echo AFTER", env)
        self.assertEqual(r.returncode, 0, r.stderr)
        lines = r.stdout.splitlines()
        self.assertEqual(lines[0], "23433436_1"); self.assertEqual(lines[39], "23433436_40")
        self.assertIn("299960 further sacct lines not shown; all 300000 lines were read", lines[40])
        self.assertEqual(lines[-1], "AFTER"); self.assertEqual(len(lines), 42)

    def test_small_output_shown_in_full(self):
        env = self.stub("printf 'JobID\\nA\\nB\\n'\n")
        r = sh(f"set -euo pipefail; . '{LIB}'; show_sacct 1 40", env)
        self.assertEqual(r.returncode, 0); self.assertEqual(r.stdout.splitlines(), ["JobID", "A", "B"])

    def test_genuine_sacct_failure_detected(self):
        env = self.stub("echo 'sacct: error: Problem talking to the database' >&2; exit 1\n")
        r = sh(f"set -euo pipefail; . '{LIB}'; show_sacct 23433436 40 || exit 1; echo AFTER", env)
        self.assertEqual(r.returncode, 1); self.assertNotIn("AFTER", r.stdout)
        self.assertIn("STOP sacct failed for job 23433436 (exit 1)", r.stderr); self.assertIn("Problem talking to the database", r.stderr)

    def test_sacct_missing_detected(self):
        r = sh(f"set -euo pipefail; . '{LIB}'; PATH=/nonexistent; show_sacct 1 || exit 1; echo AFTER")
        self.assertEqual(r.returncode, 1); self.assertIn("STOP sacct failed", r.stderr)


class TrapStatus(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.p = Path(self.tmp.name) / "harmonized_r3.partial.123"; self.p.mkdir()

    def tearDown(self):
        self.tmp.cleanup()

    def run_with(self, tail):
        return sh(f"set -euo pipefail; . '{LIB}'; install_partial_trap '{self.p}'; {tail}")

    def test_explicit_exit_status_preserved(self):
        for code in (1, 3, 7, 141):
            r = self.run_with(f"exit {code}")
            self.assertEqual(r.returncode, code); self.assertIn(f"original exit status {code}", r.stderr)

    def test_set_e_failure_preserved(self):
        r = self.run_with("false; echo unreachable")
        self.assertEqual(r.returncode, 1); self.assertIn("incomplete build left at", r.stderr)

    def test_upstream_pipeline_failure_preserved(self):
        r = self.run_with("seq 1 300000 | head -1 > /dev/null; echo unreachable")
        self.assertEqual(r.returncode, 141); self.assertIn("original exit status 141", r.stderr)

    def test_success_after_adoption_is_quiet(self):
        r = self.run_with("PARTIAL_DIR=''; exit 0")
        self.assertEqual(r.returncode, 0); self.assertNotIn("NOTE", r.stderr)

    def test_scientific_exit_3_after_adoption(self):
        r = self.run_with("PARTIAL_DIR=''; exit 3")
        self.assertEqual(r.returncode, 3); self.assertNotIn("NOTE", r.stderr)


class Header(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.d = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def mat(self, name, header, body="ENSG1\t" + "\t".join(["1.0"] * 33) + "\n", raw=None):
        p = self.d / name
        if raw is not None:
            p.write_bytes(raw)
        else:
            with gzip.open(p, "wt") as fh:
                fh.write(header + body)
        return p

    def test_valid_33_and_42(self):
        p = self.mat("ok33.gz", "ensembl_gene_id\t" + "\t".join(COLLAB33) + "\n")
        r = py("header", "--matrix", p, "--expect", ",".join(COLLAB33)); self.assertEqual(r.returncode, 0, r.stderr); self.assertEqual(r.stdout.split(), COLLAB33)
        p = self.mat("ok42.gz", "gene_symbol\t" + "\t".join(COLLAB33 + NINE) + "\n")
        self.assertEqual(py("header", "--matrix", p, "--expect", ",".join(COLLAB33 + NINE)).returncode, 0)

    def test_large_matrix_header_only(self):
        p = self.d / "big.gz"
        with gzip.open(p, "wt") as fh:
            fh.write("ensembl_gene_id\t" + "\t".join(COLLAB33) + "\n")
            for i in range(60000):
                fh.write(f"ENSG{i:011d}\t" + "\t".join(["1.0000"] * 33) + "\n")
        r = py("header", "--matrix", p, "--expect", ",".join(COLLAB33)); self.assertEqual(r.returncode, 0, r.stderr)

    def test_malformed_and_mismatched(self):
        e33 = ",".join(COLLAB33)
        cases = {
            "not_gzip": (dict(raw=b"ensembl_gene_id\tGM01\n"), 2, "cannot read gzip header"),
            "truncated_gzip": (dict(raw=gzip.compress(b"ensembl_gene_id\tGM01\tGM02\n")[:12]), 2, "cannot read gzip header"),
            "empty": (dict(header="", body=""), 1, "missing or unterminated"),
            "no_newline": (dict(header="ensembl_gene_id\t" + "\t".join(COLLAB33), body=""), 1, "missing or unterminated"),
            "wrong_first_field": (dict(header="gene\t" + "\t".join(COLLAB33) + "\n"), 1, "malformed header"),
            "duplicate_id": (dict(header="ensembl_gene_id\t" + "\t".join(COLLAB33[:-1] + ["GM01"]) + "\n"), 1, "duplicate specimen IDs"),
            "missing_id": (dict(header="ensembl_gene_id\t" + "\t".join(COLLAB33[:-1]) + "\n"), 1, "missing ['GM50']"),
            "extra_id": (dict(header="ensembl_gene_id\t" + "\t".join(COLLAB33 + ["GM37"]) + "\n"), 1, "unexpected ['GM37']"),
            "order_differs": (dict(header="ensembl_gene_id\t" + "\t".join(COLLAB33[1:] + COLLAB33[:1]) + "\n"), 1, "order differs"),
            "empty_id": (dict(header="ensembl_gene_id\t" + "\t".join(COLLAB33 + [""]) + "\n"), 1, "empty specimen ID"),
        }
        for name, (kw, code, msg) in cases.items():
            p = self.mat(name + ".gz", **({"header": "x\n"} | kw) if "raw" in kw else kw)
            r = py("header", "--matrix", p, "--expect", e33)
            self.assertEqual(r.returncode, code, f"{name}: {r.stderr}"); self.assertIn(msg, r.stderr, name); self.assertEqual(r.stdout, "", name)
        r = py("header", "--matrix", self.d / "absent.gz", "--expect", e33); self.assertEqual(r.returncode, 2); self.assertIn("not found", r.stderr)

    def test_upstream_failure_propagates_through_command_substitution(self):
        r = sh(f"set -euo pipefail; ids=$('{sys.executable}' -I '{PYH}' header --matrix '{self.d}/absent.gz' --expect GM01); echo \"AFTER $ids\"")
        self.assertNotEqual(r.returncode, 0); self.assertNotIn("AFTER", r.stdout)


def report(sid, failed=(), rule=None, drop_gate=None, verdict=None):
    g = {k: k not in failed for k in GATES}
    if drop_gate:
        g.pop(drop_gate)
    return {"specimen": sid, "rule": rule or {"version": "provisional-v1", "min_pearson_log1p_tpm_transcript": 0.99, "max_abs_mapping_diff_pp": 1.0,
                                              "justification": "frozen provisional rule", "status": "PROVISIONAL"},
            "gates": g, "verdict": verdict or ("PASS" if all(g.values()) else "FAIL"),
            "transcript_level": {"pearson_log1p_tpm": 0.9989, "n_expressed_tpm_ge_1": 41890, "n_expressed_abs_log2fc_gt_1": 120},
            "gene_level_reported_not_gated": {"pearson_log1p_tpm": 0.9999, "n_expressed_tpm_ge_1": 18530, "n_expressed_abs_log2fc_gt_1": 12},
            "yield": {"percent_mapped_diff_pp": -0.07, "reads_after_trim_diff": 278440}}


class Reports(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.d = Path(self.tmp.name) / "concordance_33"; self.d.mkdir()

    def tearDown(self):
        self.tmp.cleanup()

    def build(self, mutate=None, codes=None):
        reps = {s: report(s) for s in COLLAB33}; rc = {s: 0 for s in COLLAB33}
        if mutate:
            mutate(reps, rc)
        for s, r in reps.items():
            if r is None:
                continue
            p = self.d / f"{s}_provisional-v1.json"
            p.write_text(r if isinstance(r, str) else json.dumps(r))
        (self.d / "codes.tsv").write_text(codes if codes is not None else "".join(f"{s}\t{c}\n" for s, c in rc.items()))
        return py("reports", "--dir", self.d, "--expect", ",".join(COLLAB33), "--expect-n", "33", "--exit-codes", self.d / "codes.tsv")

    def test_all_pass(self):
        r = self.build(); self.assertEqual(r.returncode, 0, r.stdout + r.stderr); self.assertIn("33/33 PASS", r.stdout)

    def test_scientific_failures_exit_3_and_classified(self):
        def m(reps, rc):
            reps["GM05"] = report("GM05", failed=("transcript_pearson_log1p_tpm_ge_0.99",)); rc["GM05"] = 1
            reps["GM07"] = report("GM07", failed=("raw_input_reads_identical", "abs_mapping_diff_le_1pp")); rc["GM07"] = 1
        r = self.build(m)
        self.assertEqual(r.returncode, 3, r.stdout)
        rows = {l.split("\t")[0]: l.split("\t") for l in r.stdout.splitlines() if l.startswith("GM")}
        self.assertEqual(rows["GM05"][1], "CONCORDANCE_FAIL"); self.assertEqual(rows["GM07"][1], "IDENTITY_OR_INTEGRITY_FAIL")
        self.assertIn("31/33 PASS", r.stdout); self.assertIn("not hidden", r.stdout)

    def test_evidence_problems_exit_1(self):
        bad = {
            "missing_report": lambda reps, rc: reps.__setitem__("GM09", None),
            "malformed_json": lambda reps, rc: reps.__setitem__("GM09", "{not json"),
            "wrong_specimen": lambda reps, rc: reps.__setitem__("GM09", report("GM10")),
            "threshold_changed": lambda reps, rc: reps.__setitem__("GM09", report("GM09", rule={"version": "provisional-v1", "min_pearson_log1p_tpm_transcript": 0.95, "max_abs_mapping_diff_pp": 1.0})),
            "rule_renamed": lambda reps, rc: reps.__setitem__("GM09", report("GM09", rule={"version": "alt-v1", "min_pearson_log1p_tpm_transcript": 0.99, "max_abs_mapping_diff_pp": 1.0})),
            "gate_missing": lambda reps, rc: reps.__setitem__("GM09", report("GM09", drop_gate="specimen_identity")),
            "verdict_inconsistent": lambda reps, rc: reps.__setitem__("GM09", report("GM09", failed=("abs_mapping_diff_le_1pp",), verdict="PASS")),
            "exit_code_inconsistent": lambda reps, rc: rc.__setitem__("GM09", 1),
            "compare_could_not_evaluate": lambda reps, rc: (reps.__setitem__("GM09", None), rc.__setitem__("GM09", 2)),
            "empty_section": lambda reps, rc: reps.__setitem__("GM09", {**report("GM09"), "yield": {}}),
        }
        for name, m in bad.items():
            for f in self.d.iterdir():
                f.unlink()
            r = self.build(m)
            self.assertEqual(r.returncode, 1, f"{name}: {r.stdout}")
            self.assertIn("EVIDENCE_INVALID", r.stdout, name); self.assertIn("not a scientific result", r.stdout, name)

    def test_extra_report_and_unexpected_codes(self):
        r = self.build(lambda reps, rc: reps.__setitem__("GM37", report("GM37")))
        self.assertEqual(r.returncode, 1); self.assertIn("unexpected report files", r.stdout)
        for f in self.d.iterdir():
            f.unlink()
        r = self.build(codes="".join(f"{s}\t0\n" for s in COLLAB33) + "GM40\t0\n")
        self.assertEqual(r.returncode, 1); self.assertIn("unexpected specimens", r.stdout)

    def test_missing_and_duplicate_exit_codes(self):
        r = self.build(codes="".join(f"{s}\t0\n" for s in COLLAB33[:-1]))
        self.assertEqual(r.returncode, 1); self.assertIn("no compare-quant exit code recorded", r.stdout)
        for f in self.d.iterdir():
            f.unlink()
        r = self.build(codes="".join(f"{s}\t0\n" for s in COLLAB33) + "GM01\t0\n")
        self.assertEqual(r.returncode, 2); self.assertIn("recorded twice", r.stderr)

    def test_expected_inventory_size_enforced(self):
        self.build()
        r = py("reports", "--dir", self.d, "--expect", ",".join(COLLAB33[:-1]), "--expect-n", "33", "--exit-codes", self.d / "codes.tsv")
        self.assertEqual(r.returncode, 2); self.assertIn("33 unique specimens", r.stderr)


if __name__ == "__main__":
    unittest.main()
