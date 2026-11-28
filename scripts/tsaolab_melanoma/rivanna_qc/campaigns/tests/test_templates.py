"""Tests for the campaign source integration: identity of committed files, template tokens and rendering."""
import csv, hashlib, os, re, subprocess, tempfile, unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
RENDER = HERE / "render_templates.sh"
PROCS = ["procedure-20261126", "procedure-20261127", "procedure-20261127r1"]
TOKENS = {"@MELANOMA_ROOT@", "@PYTHON_BIN@"}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def table():
    with open(HERE / "deployed_hashes.tsv", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class Identity(unittest.TestCase):
    def test_identical_files_match_deployed_hashes(self):
        rows = [r for r in table() if r["repo_form"] == "identical"]
        self.assertEqual(len(rows), 12)  # Published byte-identical campaign scripts and tests.
        for r in rows:
            self.assertEqual(sha(HERE / r["procedure"] / r["path"]), r["deployed_sha256"], r["path"])

    def test_every_repo_file_is_accounted_for(self):
        listed = {(r["procedure"], r["path"] + (".template" if r["repo_form"] == "template" else "")) for r in table() if r["repo_form"] != "excluded"}
        for p in PROCS:
            for f in (HERE / p).rglob("*"):
                if f.is_file() and "__pycache__" not in f.parts:
                    self.assertIn((p, str(f.relative_to(HERE / p))), listed, f)

    def test_excluded_files_are_absent(self):
        for r in table():
            if r["repo_form"] == "excluded":
                self.assertFalse((HERE / r["procedure"] / r["path"]).exists(), r["path"])


class Templates(unittest.TestCase):
    def test_only_the_two_site_tokens(self):
        for t in HERE.rglob("*.template"):
            found = set(re.findall(r"@[A-Z_]+@", t.read_text()))
            self.assertTrue(found and found <= TOKENS, (t, found))

    def test_no_private_paths_anywhere(self):
        for f in HERE.rglob("*"):
            if f.is_file() and f.suffix not in (".pyc",):
                pat = "|".join(["/scr" + "atch/[a-z]", "/apps/" + "software", "/Us" + "ers/"])   # split so this file cannot match itself
                self.assertIsNone(re.search(pat, f.read_text(errors="ignore")), f)

    def render(self, proc, out, verify=False, env=None):
        e = dict(os.environ, MELANOMA_ROOT="/site/project", PYTHON_BIN="/site/bin/python3"); e.update(env or {})
        return subprocess.run(["bash", str(RENDER), proc, str(out)] + (["--verify"] if verify else []), capture_output=True, text=True, env=e)

    def test_render_all_procedures(self):
        with tempfile.TemporaryDirectory() as t:
            for p in PROCS:
                r = self.render(p, Path(t) / p)
                self.assertEqual(r.returncode, 0, r.stderr)
                out = Path(t) / p / p
                for f in out.rglob("*"):
                    if f.is_file():
                        self.assertNotRegex(f.read_text(errors="ignore"), r"@[A-Z_]+@", f)
                        if f.suffix in (".sh", ".sbatch"):
                            self.assertEqual(subprocess.run(["bash", "-n", str(f)]).returncode, 0, f)
                            self.assertTrue(os.access(f, os.X_OK), f)
                self.assertFalse(any(out.rglob("*.template")))

    def test_render_deterministic_and_refuses_overwrite(self):
        with tempfile.TemporaryDirectory() as t:
            a, b = Path(t) / "a", Path(t) / "b"
            self.assertEqual(self.render("procedure-20261127", a).returncode, 0)
            self.assertEqual(self.render("procedure-20261127", b).returncode, 0)
            self.assertEqual((a / "procedure-20261127" / "RENDERED.sha256").read_text(), (b / "procedure-20261127" / "RENDERED.sha256").read_text())
            r = self.render("procedure-20261127", a); self.assertEqual(r.returncode, 1); self.assertIn("refusing to overwrite", r.stderr)

    def test_missing_or_unsafe_values_refused(self):
        with tempfile.TemporaryDirectory() as t:
            e = dict(os.environ); e.pop("MELANOMA_ROOT", None); e["PYTHON_BIN"] = "/x"
            r = subprocess.run(["bash", str(RENDER), "procedure-20261127", str(Path(t) / "o")], capture_output=True, text=True, env=e)
            self.assertNotEqual(r.returncode, 0)
            r = self.render("procedure-20261127", Path(t) / "o2", env={"MELANOMA_ROOT": "/a#b"}); self.assertEqual(r.returncode, 2)

    def test_verify_with_wrong_site_values_fails_only_on_templates(self):
        with tempfile.TemporaryDirectory() as t:
            r = self.render("procedure-20261127r1", Path(t) / "v", verify=True)
            self.assertEqual(r.returncode, 1)
            self.assertIn("DIFFERS   17r1_harmonize_collect.sh (template)", r.stdout)
            self.assertIn("MATCH     collect_r1.py (identical)", r.stdout)


if __name__ == "__main__":
    unittest.main()
