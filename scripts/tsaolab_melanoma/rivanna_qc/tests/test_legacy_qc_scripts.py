"""Unit and smoke tests for the five legacy Greek melanoma Rivanna QC Slurm scripts.

Standard library only. Every test builds a synthetic MELANOMA_ROOT (tiny gzip FASTQs plus stand-ins for the
historical outputs) and runs the real sbatch files with ``bash`` against stub ``module``, ``cutadapt`` and
``fastqc`` executables. No Slurm, Lmod, network, real data or real tools are used. The tests cover control
flow, input validation and output protection only. They say nothing about how the real Cutadapt or FastQC
behave on Rivanna.

Run:  python3 -m unittest discover -s scripts/tsaolab_melanoma/rivanna_qc/tests -p 'test_legacy_*.py' -v
"""

import gzip
import hashlib
import os
import re
import shutil
import stat
import subprocess
import tempfile
import time
import unittest
from pathlib import Path

QC_DIR = Path(__file__).resolve().parent.parent
BASH = shutil.which("bash") or "/bin/bash"

O5 = "run_cutadapt_pilot.sbatch"
O10 = "run_cutadapt_pilot_O10.sbatch"
FQ_ALL = "run_fastqc_all.sbatch"
FQ_TRIM = "run_fastqc_trimmed_pilot.sbatch"
GZIP = "validate_gzip.sbatch"
LEGACY = [O5, O10, FQ_ALL, FQ_TRIM, GZIP]
STEP = {O5: "cutadapt_pilot_O5", O10: "cutadapt_pilot_O10", FQ_ALL: "fastqc_all",
        FQ_TRIM: "fastqc_trimmed_pilot", GZIP: "gzip_validation"}
ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
RAW_11_01 = "20088XR-11-01_S0_L001_R{}_001.fastq.gz"
MANIFEST_SPECIMENS = [("20088XR-11-01", "S0", "L001"), ("20088XR-11-02", "S119", "L004"),
                      ("20088XR-11-03", "S215", "L006")]

STUB_MODULE = """#!/bin/bash
echo "$*" >> "$STUB_LOG_DIR/module.$$"
exit "${STUB_MODULE_RC:-0}"
"""
STUB_CUTADAPT = """#!/bin/bash
set -eu
if [ "${1:-}" = --version ]; then echo "${STUB_CUTADAPT_VERSION:-5.2}"; exit 0; fi
printf '%s\\n' "$@" > "$STUB_LOG_DIR/cutadapt.$$"
sleep "${STUB_SLEEP:-0}"
o=; p=; j=; args=("$@"); n=${#args[@]}
i=0; while [ $i -lt $n ]; do case "${args[$i]}" in
  -o) o=${args[$((i+1))]};; -p) p=${args[$((i+1))]};; --json) j=${args[$((i+1))]};; esac; i=$((i+1)); done
set -o noclobber
cat "${args[$((n-2))]}" > "$o"; cat "${args[$((n-1))]}" > "$p"; echo '{"stub": true}' > "$j"
"""
STUB_FASTQC = """#!/bin/bash
set -eu
if [ "${1:-}" = --version ]; then echo "${STUB_FASTQC_VERSION:-FastQC v0.12.1}"; exit 0; fi
printf '%s\\n' "$@" > "$STUB_LOG_DIR/fastqc.$$"
sleep "${STUB_SLEEP:-0}"
out=; files=()
while [ $# -gt 0 ]; do case "$1" in --outdir) out=$2; shift;; --threads) shift;; *) files+=("$1");; esac; shift; done
[ -d "$out" ] || { echo "stub fastqc: outdir missing" >&2; exit 2; }
set -o noclobber
for f in "${files[@]}"; do b=$(basename "$f"); b=${b%.gz}; b=${b%.fastq}; b=${b%.fq}
  echo stub > "$out/${b}_fastqc.html"; echo stub > "$out/${b}_fastqc.zip"; done
"""


def fastq_bytes(names, mate, seq="ACGTACGTAC"):
    out = []
    for n in names:
        out.append(f"@{n} {mate}:N:0:ACGT\n{seq}\n+\n{'I' * len(seq)}\n")
    return "".join(out).encode()


def write_fastq_gz(path, names, mate, corrupt_tail=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    data = gzip.compress(fastq_bytes(names, mate), mtime=0)
    if corrupt_tail:
        # Valid first gzip member followed by garbage: decompressing the first line works, gzip -t fails.
        big = gzip.compress(fastq_bytes([f"{names[0]}_{k}" for k in range(2000)], mate), mtime=0)
        data = data + big[: len(big) // 2] + b"\x00garbage"
    path.write_bytes(data)


def snapshot(root):
    """Map every path under root (directories included) to (kind, sha256, size, mtime_ns, mode)."""
    snap = {}
    for dirpath, dirnames, filenames in os.walk(root):
        for name in dirnames + filenames:
            p = Path(dirpath) / name
            st = p.lstat()
            if stat.S_ISLNK(st.st_mode):
                snap[str(p.relative_to(root))] = ("link", os.readlink(p), None, st.st_mtime_ns, st.st_mode)
            elif stat.S_ISDIR(st.st_mode):
                snap[str(p.relative_to(root))] = ("dir", None, None, None, st.st_mode)
            else:
                h = hashlib.sha256(p.read_bytes()).hexdigest()
                snap[str(p.relative_to(root))] = ("file", h, st.st_size, st.st_mtime_ns, st.st_mode)
    return snap


class Fixture:
    def __init__(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="legacy_qc_"))
        self.root = self.tmp / "melanoma_root"
        self.raw = self.root / "raw_fastq"
        self.stub = self.tmp / "bin"
        self.logs = self.tmp / "stub_logs"
        self.submit = self.tmp / "submit"
        for d in (self.raw, self.stub, self.logs, self.submit):
            d.mkdir(parents=True)
        for name, body in (("module", STUB_MODULE), ("cutadapt", STUB_CUTADAPT), ("fastqc", STUB_FASTQC)):
            p = self.stub / name
            p.write_text(body)
            p.chmod(0o755)
        for core, s, lane in MANIFEST_SPECIMENS:
            names = [f"LH00150:1058:22XX:{lane[-1]}:1101:{k}:1000" for k in range(3)]
            for mate in (1, 2):
                write_fastq_gz(self.raw / f"{core}_{s}_{lane}_R{mate}_001.fastq.gz", names, mate)
        self.write_manifest(self.default_manifest_rows())
        self.make_historical()

    def default_manifest_rows(self):
        rows = []
        for core, s, lane in MANIFEST_SPECIMENS:
            r1 = self.raw / f"{core}_{s}_{lane}_R1_001.fastq.gz"
            rows.append((core, str(r1), str(r1).replace("_R1_", "_R2_")))
        return rows

    def write_manifest(self, rows, header="sample\tr1_path\tr2_path", newline="\n", raw_text=None):
        """Default is the historical Rivanna format: header line (field 1 'sample'), then SAMPLE<TAB>R1<TAB>R2."""
        path = self.root / "validated_pairs.tsv"
        if raw_text is None:
            lines = ([header] if header is not None else []) + ["\t".join(r) for r in rows]
            raw_text = newline.join(lines) + newline
        path.write_bytes(raw_text.encode())

    def make_historical(self):
        """Stand-ins for every historical output the scripts must never touch."""
        r = self.root
        (r / "gzip_full_validation.tsv").write_text("FAIL\tx.fastq.gz\nPASS\ty.fastq.gz\n")
        (r / "gzip_validation_input.txt").write_text("x.fastq.gz\ny.fastq.gz\n")
        for d in ("trimmed_pilot_11_01", "trimmed_pilot_11_01_O10"):
            names = ["LH00150:1062:22YY:8:1101:1:1000"]
            for mate in (1, 2):
                write_fastq_gz(r / d / f"11-01_R{mate}.trimmed.fastq.gz", names, mate)
            (r / d / "cutadapt_report.json").write_text('{"historical": true}\n')
        (r / "fastqc_all").mkdir()
        (r / "fastqc_all" / "20088XR-11-01_S0_L001_R1_001_fastqc.html").write_text("historical\n")
        (r / "fastqc_trimmed_pilot").mkdir()
        (r / "fastqc_trimmed_pilot" / "11-01_R1.trimmed_fastqc.html").write_text("historical\n")

    def env(self, script, **over):
        env = {"PATH": f"{self.stub}:/usr/bin:/bin:/usr/sbin:/sbin", "HOME": str(self.tmp),
               "TMPDIR": str(self.tmp), "LC_ALL": "C", "MELANOMA_ROOT": str(self.root),
               "STUB_LOG_DIR": str(self.logs)}
        if script == FQ_ALL:
            env.update(SLURM_JOB_ID="900001", SLURM_ARRAY_JOB_ID="900000", SLURM_ARRAY_TASK_ID="0",
                       SLURM_ARRAY_TASK_MIN="0", SLURM_ARRAY_TASK_MAX=str(len(MANIFEST_SPECIMENS) - 1))
        else:
            env["SLURM_JOB_ID"] = "800001"
        for k, v in over.items():
            if v is None:
                env.pop(k, None)
            else:
                env[k] = v
        return env

    def run(self, script, **over):
        return subprocess.run([BASH, str(QC_DIR / script)], cwd=self.submit, env=self.env(script, **over),
                              capture_output=True, text=True, timeout=120)

    def popen(self, script, **over):
        return subprocess.Popen([BASH, str(QC_DIR / script)], cwd=self.submit, env=self.env(script, **over),
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    def run_dir(self, script, env):
        base = self.root / "qc_runs" / STEP[script]
        if script == FQ_ALL:
            return base / f"job_{env['SLURM_ARRAY_JOB_ID']}" / f"task_{int(env['SLURM_ARRAY_TASK_ID'])}"
        return base / f"job_{env['SLURM_JOB_ID']}"

    def tool_calls(self, tool):
        return sorted(self.logs.glob(f"{tool}.*"))

    def cleanup(self):
        shutil.rmtree(self.tmp, ignore_errors=True)


def historical_only(snap):
    return {k: v for k, v in snap.items() if not k.startswith("qc_runs")}


class Base(unittest.TestCase):
    def setUp(self):
        self.fx = Fixture()
        self.addCleanup(self.fx.cleanup)

    def assertSucceeded(self, res):
        self.assertEqual(res.returncode, 0, f"stdout:\n{res.stdout}\nstderr:\n{res.stderr}")

    def assertRefused(self, res, pattern):
        self.assertNotEqual(res.returncode, 0, f"expected refusal; stdout:\n{res.stdout}")
        self.assertRegex(res.stderr, pattern)

    def assertRefusedUntouched(self, script, pattern, **over):
        """Refusal before any output: the entire root must be byte-identical and no tool may run."""
        before = snapshot(self.fx.root)
        res = self.fx.run(script, **over)
        self.assertRefused(res, pattern)
        self.assertEqual(snapshot(self.fx.root), before, "files changed or were created on refusal")
        self.assertEqual(self.fx.tool_calls("cutadapt") + self.fx.tool_calls("fastqc"), [])
        return res


class StaticChecks(unittest.TestCase):
    def test_bash_syntax(self):
        for s in LEGACY:
            with self.subTest(script=s):
                r = subprocess.run([BASH, "-n", str(QC_DIR / s)], capture_output=True, text=True)
                self.assertEqual(r.returncode, 0, r.stderr)

    def test_guard_block_identical_in_all_scripts(self):
        blocks = {}
        for s in LEGACY:
            text = (QC_DIR / s).read_text()
            m = re.search(r"# >>> legacy-qc-guard v1 >>>\n.*?# <<< legacy-qc-guard v1 <<<\n", text, re.S)
            self.assertIsNotNone(m, f"guard block missing in {s}")
            blocks[s] = m.group(0)
        self.assertEqual(len(set(blocks.values())), 1, "guard blocks differ between scripts")

    def test_no_overwrite_bypass_and_no_historical_targets(self):
        for s in LEGACY:
            text = (QC_DIR / s).read_text()
            code = "\n".join(l for l in text.splitlines() if not l.lstrip().startswith("#"))
            with self.subTest(script=s):
                self.assertNotRegex(code, r"ALLOW_OVERWRITE|FORCE|--force|rm -rf|mkdir -p")
                self.assertNotRegex(code, r'\$ROOT/(trimmed_pilot_11_01[^"]*/[^"]*-o|fastqc_all|fastqc_trimmed_pilot|gzip_full_validation|gzip_validation_input)')
                self.assertNotIn(': >', code)
                self.assertIn("make_run_dir", code)

    def test_cutadapt_parameters_preserved_in_source(self):
        for s, o in ((O5, "5"), (O10, "10")):
            text = (QC_DIR / s).read_text()
            with self.subTest(script=s):
                self.assertIn("module load \"$1\"", text)
                self.assertIn("load_module cutadapt/5.2", text)
                for line in ("  -j 2 \\", f"  -O {o} \\", f"  -a {ADAPTER_R1} \\", f"  -A {ADAPTER_R2} \\"):
                    self.assertIn(line, text)

    def test_fastqc_pinned_everywhere(self):
        for s in (FQ_ALL, FQ_TRIM):
            text = (QC_DIR / s).read_text()
            with self.subTest(script=s):
                self.assertIn("load_module fastqc/0.12.1", text)
                self.assertNotRegex(text, r"load_module fastqc\s*$", )


class SuccessAndProtection(Base):
    """Every script: success writes only to a fresh qc_runs directory; historical outputs never change."""

    def test_success_writes_only_new_run_dir(self):
        for s in LEGACY:
            with self.subTest(script=s):
                before = snapshot(self.fx.root)
                env = self.fx.env(s)
                res = self.fx.run(s)
                self.assertSucceeded(res)
                rd = self.fx.run_dir(s, env)
                self.assertTrue((rd / "COMPLETE").is_file())
                self.assertFalse((rd / "INCOMPLETE").exists())
                info = (rd / "run_info.tsv").read_text()
                self.assertIn(f"script\t{s}\n", info)
                self.assertRegex(info, r"script_sha256\t[0-9a-f]{64}\n")
                after = snapshot(self.fx.root)
                self.assertEqual(historical_only(after), historical_only(before))
                new = set(after) - set(before)
                prefix = str(rd.relative_to(self.fx.root))
                self.assertTrue(all(p == prefix or p.startswith(prefix + "/") or prefix.startswith(p)
                                    for p in new), sorted(new))

    def test_refuses_existing_run_dir_byte_identical(self):
        for s in LEGACY:
            with self.subTest(script=s):
                env = self.fx.env(s)
                rd = self.fx.run_dir(s, env)
                rd.mkdir(parents=True)
                (rd / "earlier_result.txt").write_text("must survive\n")
                self.assertRefusedUntouched(s, r"already exists .*refusing to overwrite")

    def test_refuses_symlinked_output_parent(self):
        target = self.fx.root / "fastqc_all"
        (self.fx.root / "qc_runs").mkdir()
        (self.fx.root / "qc_runs" / "fastqc_trimmed_pilot").symlink_to(target)
        self.assertRefusedUntouched(FQ_TRIM, r"symlink")

    def test_concurrent_same_job_id_exactly_one_wins(self):
        for s in LEGACY:
            with self.subTest(script=s):
                before = historical_only(snapshot(self.fx.root))
                procs = [self.fx.popen(s, STUB_SLEEP="0.5") for _ in range(6)]
                outs = [p.communicate(timeout=120) for p in procs]
                results = [(p.returncode, o) for p, o in zip(procs, outs)]
                codes = [rc for rc, _ in results]
                self.assertEqual(codes.count(0), 1, results)
                for rc, (_, err) in results:
                    if rc != 0:
                        self.assertRegex(err, r"already exists .*refusing to overwrite")
                rd = self.fx.run_dir(s, self.fx.env(s))
                self.assertTrue((rd / "COMPLETE").is_file())
                self.assertFalse((rd / "INCOMPLETE").exists())
                self.assertEqual(historical_only(snapshot(self.fx.root)), before)

    def test_concurrent_distinct_job_ids_all_succeed(self):
        for s in (O5, FQ_TRIM, GZIP):
            with self.subTest(script=s):
                ids = [str(700000 + k) for k in range(4)]
                procs = [self.fx.popen(s, SLURM_JOB_ID=i, STUB_SLEEP="0.3") for i in ids]
                for p in procs:
                    p.communicate(timeout=120)
                self.assertEqual([p.returncode for p in procs], [0] * 4)
                for i in ids:
                    self.assertTrue((self.fx.root / "qc_runs" / STEP[s] / f"job_{i}" / "COMPLETE").is_file())

    def test_requires_slurm_job_id(self):
        for s in LEGACY:
            with self.subTest(script=s):
                over = {"SLURM_ARRAY_JOB_ID": None} if s == FQ_ALL else {"SLURM_JOB_ID": None}
                self.assertRefusedUntouched(s, r"is not set; submit this script with sbatch", **over)

    def test_rejects_non_numeric_job_id(self):
        self.assertRefusedUntouched(O5, r"not a non-negative integer", SLURM_JOB_ID="../../fastqc_all")

    def test_requires_absolute_existing_root(self):
        for s in LEGACY:
            with self.subTest(script=s):
                self.assertRefusedUntouched(s, r"set MELANOMA_ROOT", MELANOMA_ROOT=None)
                self.assertRefusedUntouched(s, r"must be an absolute path", MELANOMA_ROOT="melanoma_root")
                self.assertRefusedUntouched(s, r"not a directory", MELANOMA_ROOT=str(self.fx.tmp / "nope"))

    def test_failed_tool_leaves_incomplete_marker_and_no_complete(self):
        res = self.fx.run(O5, STUB_LOG_DIR=str(self.fx.tmp / "missing_dir"))  # stub cannot log -> exits 1
        self.assertNotEqual(res.returncode, 0)
        rd = self.fx.run_dir(O5, self.fx.env(O5))
        self.assertFalse((rd / "COMPLETE").exists())
        self.assertTrue((rd / "INCOMPLETE").is_file())
        # A later attempt with the same job ID must still refuse and leave the partial run alone.
        snap = snapshot(rd)
        self.assertRefused(self.fx.run(O5), r"already exists")
        self.assertEqual(snapshot(rd), snap)


class CutadaptPilots(Base):
    def test_exact_parameters_and_inputs(self):
        for s, o in ((O5, "5"), (O10, "10")):
            with self.subTest(script=s):
                for f in self.fx.logs.glob("*"):
                    f.unlink()
                env = self.fx.env(s)
                self.assertSucceeded(self.fx.run(s))
                rd = self.fx.run_dir(s, env)
                calls = self.fx.tool_calls("cutadapt")
                self.assertEqual(len(calls), 1)
                argv = calls[0].read_text().splitlines()
                raw = self.fx.raw
                self.assertEqual(argv, [
                    "-j", "2", "-O", o, "-a", ADAPTER_R1, "-A", ADAPTER_R2,
                    "--json", str(rd / "cutadapt_report.json"),
                    "-o", str(rd / "11-01_R1.trimmed.fastq.gz"),
                    "-p", str(rd / "11-01_R2.trimmed.fastq.gz"),
                    str(raw / RAW_11_01.format(1)), str(raw / RAW_11_01.format(2))])
                mods = "".join(p.read_text() for p in self.fx.logs.glob("module.*"))
                self.assertIn("load cutadapt/5.2", mods)

    def test_missing_r2_refused(self):
        (self.fx.raw / RAW_11_01.format(2)).unlink()
        for s in (O5, O10):
            with self.subTest(script=s):
                self.assertRefusedUntouched(s, r"missing, empty or unreadable input")

    def test_empty_input_refused(self):
        (self.fx.raw / RAW_11_01.format(1)).write_bytes(b"")
        self.assertRefusedUntouched(O10, r"missing, empty or unreadable input")

    def test_non_gzip_input_refused(self):
        (self.fx.raw / RAW_11_01.format(1)).write_text("@r1\nACGT\n+\nIIII\n")
        self.assertRefusedUntouched(O5, r"not gzip-compressed")

    def test_mate_header_mismatch_refused(self):
        write_fastq_gz(self.fx.raw / RAW_11_01.format(2), ["SOME:OTHER:READ"], 2)
        self.assertRefusedUntouched(O5, r"first read names differ between mates")

    def test_non_fastq_first_record_refused(self):
        (self.fx.raw / RAW_11_01.format(1)).write_bytes(gzip.compress(b"not a fastq\n"))
        self.assertRefusedUntouched(O10, r"not a FASTQ header")

    def test_wrong_cutadapt_version_refused(self):
        self.assertRefusedUntouched(O5, r"version mismatch: expected '5.2', got '5.1'", STUB_CUTADAPT_VERSION="5.1")

    def test_module_failure_refused(self):
        self.assertRefusedUntouched(O10, r"module load cutadapt/5.2 failed", STUB_MODULE_RC="1")

    def test_validation_reads_only_the_first_record(self):
        # Valid first record, corrupt gzip tail: a full decompression would fail, header validation must not.
        names = ["LH00150:1062:22YY:8:1101:1:1000"]
        write_fastq_gz(self.fx.raw / RAW_11_01.format(1), names, 1, corrupt_tail=True)
        write_fastq_gz(self.fx.raw / RAW_11_01.format(2), names, 2, corrupt_tail=True)
        self.assertNotEqual(subprocess.run(["gzip", "-t", str(self.fx.raw / RAW_11_01.format(1))],
                                           capture_output=True).returncode, 0)
        self.assertSucceeded(self.fx.run(O5))


class FastqcAll(Base):
    def run_task(self, task, **over):
        return self.fx.run(FQ_ALL, SLURM_ARRAY_TASK_ID=str(task), **over)

    def test_each_task_processes_its_row(self):
        for task, (core, s, lane) in enumerate(MANIFEST_SPECIMENS):
            with self.subTest(task=task):
                res = self.run_task(task)
                self.assertSucceeded(res)
                rd = self.fx.root / "qc_runs" / "fastqc_all" / "job_900000" / f"task_{task}"
                self.assertTrue((rd / f"{core}_{s}_{lane}_R1_001_fastqc.html").is_file())
                self.assertTrue((rd / f"{core}_{s}_{lane}_R2_001_fastqc.html").is_file())
                info = (rd / "run_info.tsv").read_text()
                self.assertIn(f"sample\t{core}\n", info)
                self.assertRegex(info, r"manifest_sha256\t[0-9a-f]{64}\n")
        mods = "".join(p.read_text() for p in self.fx.logs.glob("module.*"))
        self.assertIn("load fastqc/0.12.1", mods)

    def test_concurrent_tasks_of_one_array_job(self):
        procs = [self.fx.popen(FQ_ALL, SLURM_ARRAY_TASK_ID=str(t), STUB_SLEEP="0.3")
                 for t in range(len(MANIFEST_SPECIMENS))]
        for p in procs:
            p.communicate(timeout=120)
        self.assertEqual([p.returncode for p in procs], [0] * len(MANIFEST_SPECIMENS))

    def test_task_id_out_of_range(self):
        n = len(MANIFEST_SPECIMENS)
        self.assertRefusedUntouched(FQ_ALL, r"out of range: manifest has 3 data rows",
                                    SLURM_ARRAY_TASK_ID=str(n), SLURM_ARRAY_TASK_MAX=str(n))

    def test_array_max_beyond_manifest(self):
        self.assertRefusedUntouched(FQ_ALL, r"--array range ends at 42 but the manifest has only 3",
                                    SLURM_ARRAY_TASK_MAX="42")

    def test_missing_or_bad_task_id(self):
        self.assertRefusedUntouched(FQ_ALL, r"SLURM_ARRAY_TASK_ID is not set", SLURM_ARRAY_TASK_ID=None)
        self.assertRefusedUntouched(FQ_ALL, r"not a non-negative integer", SLURM_ARRAY_TASK_ID="-1")
        self.assertRefusedUntouched(FQ_ALL, r"not a non-negative integer", SLURM_ARRAY_TASK_ID="1a")

    def test_missing_manifest(self):
        (self.fx.root / "validated_pairs.tsv").unlink()
        self.assertRefusedUntouched(FQ_ALL, r"missing, empty or unreadable manifest")

    def test_manifest_format_errors(self):
        """Each defect is refused in both the historical (headered) and the headerless layout."""
        rows = self.fx.default_manifest_rows()
        cases = {
            "two_fields": ([rows[0][:2]] + rows[1:], r"line {L}: expected 3 tab-separated fields, found 2"),
            "four_fields": ([rows[0] + ("x",)] + rows[1:], r"line {L}: expected 3 tab-separated fields, found 4"),
            "trailing_tab": ([rows[0] + ("",)] + rows[1:], r"line {L}: expected 3 tab-separated fields, found 4"),
            "relative": ([(rows[0][0], "raw_fastq/" + Path(rows[0][1]).name, "raw_fastq/" + Path(rows[0][2]).name)]
                         + rows[1:], r"line {L}: R1 path is not absolute"),
            "not_mates": ([(rows[0][0], rows[0][1], rows[1][2])] + rows[1:], r"line {L}: R2 is not the mate of R1"),
            "swapped": ([(rows[0][0], rows[0][2], rows[0][1])] + rows[1:], r"line {L}: R1 name does not end in _R1_NNN"),
            "dup_label": ([rows[0], (rows[0][0],) + rows[1][1:], rows[2]], r"line {L1}: duplicate sample label"),
            "dup_file": ([rows[0], ("other",) + rows[0][1:], rows[2]], r"line {L1}: file listed twice"),
            "spaces": ([(rows[0][0] + " ",) + rows[0][1:]] + rows[1:], r"line {L}: empty field or surrounding whitespace"),
            "label_vs_file": ([("20088XR-11-09",) + rows[0][1:]] + rows[1:], r"line {L}: label 20088XR-11-09 is inconsistent"),
        }
        for layout, header in (("headered", "sample\tr1\tr2"), ("headerless", None)):
            first = 2 if header else 1
            for name, (rws, pattern) in cases.items():
                with self.subTest(layout=layout, case=name):
                    self.fx.write_manifest(rws, header=header)
                    self.assertRefusedUntouched(FQ_ALL, pattern.format(L=first, L1=first + 1))
        with self.subTest(case="header_only"):
            self.fx.write_manifest([], header="sample\tr1\tr2")
            self.assertRefusedUntouched(FQ_ALL, r"no data rows")
        with self.subTest(case="blank_line"):
            text = "sample\tr1\tr2\n" + "\t".join(rows[0]) + "\n\n" + "\t".join(rows[1]) + "\n"
            self.fx.write_manifest(None, raw_text=text)
            self.assertRefusedUntouched(FQ_ALL, r"line 3: blank line")
        for name, text in (("crlf_all", "\r\n".join(["sample\tr1\tr2"] + ["\t".join(r) for r in rows]) + "\r\n"),
                           ("crlf_one_row", "sample\tr1\tr2\n" + "\t".join(rows[0]) + "\r\n"
                            + "\n".join("\t".join(r) for r in rows[1:]) + "\n")):
            with self.subTest(case=name):
                self.fx.write_manifest(None, raw_text=text)
                self.assertRefusedUntouched(FQ_ALL, r"CR \(Windows\) line endings")

    def test_listed_file_missing_or_empty(self):
        rows = self.fx.default_manifest_rows()
        Path(rows[1][2]).unlink()
        self.assertRefusedUntouched(FQ_ALL, r"missing, empty or unreadable input", SLURM_ARRAY_TASK_ID="1")
        Path(rows[2][1]).write_bytes(b"")
        self.assertRefusedUntouched(FQ_ALL, r"missing, empty or unreadable input", SLURM_ARRAY_TASK_ID="2")

    def test_header_desync(self):
        rows = self.fx.default_manifest_rows()
        write_fastq_gz(Path(rows[0][2]), ["LH00150:1058:22XX:1:1101:999:1000"], 2)
        self.assertRefusedUntouched(FQ_ALL, r"first read names differ")

    def test_wrong_fastqc_version(self):
        self.assertRefusedUntouched(FQ_ALL, r"expected 'FastQC v0.12.1', got 'FastQC v0.11.9'",
                                    STUB_FASTQC_VERSION="FastQC v0.11.9")


class ManifestFormatRegression(Base):
    """Historical Rivanna validated_pairs.tsv (verified on Rivanna): header sample<TAB>r1_path<TAB>r2_path + 43 rows,
    3 tab-separated fields, absolute paths, LF. Headerless manifests remain supported."""

    N_REAL = 43
    REAL_HEADER = "sample\tr1_path\tr2_path"

    def make_real_manifest(self, header=REAL_HEADER):
        rows = []
        cores = [n for n in range(1, 51) if n not in (45, 49)][: self.N_REAL]
        for n in cores:
            core = f"20088XR-11-{n:02d}"
            names = [f"LH00150:1058:22XX:4:1101:{n}:1000"]
            r1 = self.fx.raw / f"{core}_S{100 + n}_L004_R1_001.fastq.gz"
            r2 = self.fx.raw / f"{core}_S{100 + n}_L004_R2_001.fastq.gz"
            if not r1.exists():
                write_fastq_gz(r1, names, 1)
                write_fastq_gz(r2, names, 2)
            rows.append((core, str(r1), str(r2)))
        self.fx.write_manifest(rows, header=header)
        return rows

    def run_task(self, task, **over):
        over.setdefault("SLURM_ARRAY_TASK_MAX", str(self.N_REAL - 1))
        return self.fx.run(FQ_ALL, SLURM_ARRAY_TASK_ID=str(task), **over)

    def info(self, task, job):
        text = (self.fx.root / "qc_runs" / "fastqc_all" / f"job_{job}" / f"task_{task}" / "run_info.tsv").read_text()
        return dict(l.split("\t", 1) for l in text.splitlines())

    def check_mapping(self, rows, header_flag, job):
        for task in (0, 1, self.N_REAL - 1):
            res = self.run_task(task, SLURM_ARRAY_JOB_ID=job)
            self.assertSucceeded(res)
            info = self.info(task, job)
            self.assertEqual(info["sample"], rows[task][0])
            self.assertEqual(info["input_r1"], rows[task][1])
            self.assertEqual(info["manifest_data_rows"], str(self.N_REAL))
            self.assertEqual(info["manifest_header_line"], header_flag)
            self.assertIn(f"Sample: {rows[task][0]}\n", res.stdout)

    def test_real_44_line_headered_format(self):
        rows = self.make_real_manifest()
        text = (self.fx.root / "validated_pairs.tsv").read_text()
        self.assertEqual(len(text.splitlines()), 44)
        self.check_mapping(rows, "1", "910000")

    def test_real_format_bounds(self):
        self.make_real_manifest()
        self.assertRefusedUntouched(FQ_ALL, r"out of range: manifest has 43 data rows \(valid 0-42\)",
                                    SLURM_ARRAY_TASK_ID="43", SLURM_ARRAY_TASK_MAX="43")
        self.assertRefusedUntouched(FQ_ALL, r"--array range ends at 43 but the manifest has only 43",
                                    SLURM_ARRAY_TASK_ID="0", SLURM_ARRAY_TASK_MAX="43")

    def test_headerless_variant_same_mapping(self):
        rows = self.make_real_manifest(header=None)
        self.check_mapping(rows, "0", "920000")

    def test_header_spellings(self):
        headers = ("sample\tr1\tr2", "SAMPLE\tR1\tR2", "sample_id\tr1_fastq\tr2_fastq", "label\tread1\tread2")
        for k, header in enumerate(headers):
            with self.subTest(header=header):
                rows = self.make_real_manifest(header=header)
                res = self.run_task(0, SLURM_ARRAY_JOB_ID=str(930000 + k))
                self.assertSucceeded(res)
                self.assertIn(f"Sample: {rows[0][0]}\n", res.stdout)

    def test_sbatch_array_header_matches_real_manifest(self):
        text = (QC_DIR / FQ_ALL).read_text()
        m = re.search(r"^#SBATCH --array=0-(\d+)", text, re.M)
        self.assertIsNotNone(m)
        self.assertEqual(int(m.group(1)), self.N_REAL - 1)


TOOLS_PY = QC_DIR / "rivanna_qc_tools.py"


@unittest.skipUnless(TOOLS_PY.exists(), "rivanna_qc_tools.py (quantification workstream) not present; parity not checked")
class ParserParity(Base):
    """C2/C3: run_fastqc_all.sbatch and rivanna_qc_tools.read_validated_pairs accept and reject the same manifests,
    and a malformed first row is never discarded as a header by either."""

    @classmethod
    def setUpClass(cls):
        import importlib.util
        spec = importlib.util.spec_from_file_location("rivanna_qc_tools", QC_DIR / "rivanna_qc_tools.py")
        cls.tools = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.tools)

    def tools_parse(self, path):
        import contextlib
        import io
        err = io.StringIO()
        try:
            with contextlib.redirect_stderr(err):
                rows = self.tools.read_validated_pairs(str(path))
            return True, rows, ""
        except SystemExit:
            return False, None, err.getvalue()

    def variants(self):
        rows = self.fx.default_manifest_rows()
        r0 = rows[0]
        dup_dir = self.fx.raw
        names = ["LH00150:1058:22XX:6:1101:7:1000"]
        for mate in (1, 2):
            write_fastq_gz(dup_dir / f"20088XR-11-01_S226_L006_R{mate}_001.fastq.gz", names, mate)
        dup_core = (str(dup_dir / "20088XR-11-01_S226_L006_R1_001.fastq.gz"),
                    str(dup_dir / "20088XR-11-01_S226_L006_R2_001.fastq.gz"))
        other = self.fx.tmp / "elsewhere"
        other.mkdir(exist_ok=True)
        H = "sample\tr1_path\tr2_path\n"
        body = lambda rs: "".join("\t".join(r) + "\n" for r in rs)
        return {
            # name: (manifest text, expected accepted?, line number that must be reported when rejected)
            "real_headered": (H + body(rows), True, None),
            "headerless": (body(rows), True, None),
            "header_uppercase": ("SAMPLE\tR1\tR2\n" + body(rows), True, None),
            "first_row_missing_gz": (H + body([(r0[0], r0[1][:-3], r0[2][:-3])] + rows[1:]), False, 2),
            "headerless_first_row_missing_gz": (body([(r0[0], r0[1][:-3], r0[2][:-3])] + rows[1:]), False, 1),
            "headerless_first_row_relative": (body([(r0[0], Path(r0[1]).name, Path(r0[2]).name)] + rows[1:]), False, 1),
            "headerless_first_row_two_fields": (body([r0[:2]] + rows[1:]), False, 1),
            "header_word_with_path_fields": ("sample\t" + r0[1] + "\t" + r0[2][:-3] + "\n" + body(rows[1:]), False, 1),
            "unknown_header_word": ("sample_name\tr1\tr2\n" + body(rows), False, 1),
            "dup_label": (H + body([r0, (r0[0],) + rows[1][1:], rows[2]]), False, 3),
            "dup_specimen_other_files": (H + body([r0, ("11-01-copy",) + dup_core] + rows[1:]), False, 3),
            "dup_specimen_short_label": (H + body([r0, ("11-01",) + dup_core] + rows[1:]), False, 3),
            "swapped_mates": (H + body([(r0[0], r0[2], r0[1])] + rows[1:]), False, 2),
            "mates_from_two_specimens": (H + body([(r0[0], r0[1], rows[1][2])] + rows[1:]), False, 2),
            "mates_in_two_dirs": (H + body([(r0[0], r0[1], str(other / Path(r0[2]).name))] + rows[1:]), False, 2),
            "label_other_specimen": (H + body([("11-09",) + r0[1:]] + rows[1:]), False, 2),
            "extra_column": (H + body([r0 + ("note",)] + rows[1:]), False, 2),
            "trailing_tab": (H + body([r0 + ("",)] + rows[1:]), False, 2),
            "padded_label": (H + body([(" " + r0[0],) + r0[1:]] + rows[1:]), False, 2),
            "crlf_all": ((H + body(rows)).replace("\n", "\r\n"), False, None),
            "crlf_one_row": (H + "\t".join(r0) + "\r\n" + body(rows[1:]), False, None),
            "blank_line": (H + body(rows[:1]) + "\n" + body(rows[1:]), False, 3),
            "header_only": (H, False, None),
        }

    def test_both_parsers_agree(self):
        for k, (name, (text, accept, line)) in enumerate(self.variants().items()):
            with self.subTest(case=name):
                self.fx.write_manifest(None, raw_text=text)
                ok, rows, terr = self.tools_parse(self.fx.root / "validated_pairs.tsv")
                res = self.fx.run(FQ_ALL, SLURM_ARRAY_JOB_ID=str(940000 + k), SLURM_ARRAY_TASK_ID="0",
                                  SLURM_ARRAY_TASK_MAX=None)
                self.assertEqual(ok, accept, f"rivanna_qc_tools: {terr}")
                self.assertEqual(res.returncode == 0, accept, f"run_fastqc_all.sbatch: {res.stderr}")
                if accept:
                    self.assertEqual(rows[0][0], self.fx.default_manifest_rows()[0][0])
                    self.assertIn(f"Sample: {rows[0][0]}\n", res.stdout)
                elif line is not None:
                    self.assertRegex(terr, rf"malformed line {line}\b")
                    self.assertRegex(res.stderr, rf"line {line}\b")
                if not accept:
                    self.assertNotRegex(res.stdout, r"FASTQC_COMPLETE")

    def test_rules_documented_in_both_parsers(self):
        script = (QC_DIR / FQ_ALL).read_text()
        tools = (QC_DIR / "rivanna_qc_tools.py").read_text()
        for word in ("sample|sample_id|specimen|label",):
            self.assertIn(word, script)
        self.assertIn('VP_HEADER_FIRST = {"sample", "sample_id", "specimen", "label"}', tools)
        self.assertIn("tests check parity", script)
        self.assertIn("test_legacy_qc_scripts.py checks parity", tools)


class FastqcTrimmedPilot(Base):
    def test_reads_historical_pilot_without_modifying_it(self):
        before = snapshot(self.fx.root / "trimmed_pilot_11_01")
        self.assertSucceeded(self.fx.run(FQ_TRIM))
        self.assertEqual(snapshot(self.fx.root / "trimmed_pilot_11_01"), before)
        mods = "".join(p.read_text() for p in self.fx.logs.glob("module.*"))
        self.assertIn("load fastqc/0.12.1", mods)

    def test_uses_completed_new_pilot_run(self):
        self.assertSucceeded(self.fx.run(O5, SLURM_JOB_ID="600001"))
        pilot = self.fx.root / "qc_runs" / "cutadapt_pilot_O5" / "job_600001"
        res = self.fx.run(FQ_TRIM, PILOT_DIR=str(pilot))
        self.assertSucceeded(res)

    def test_incomplete_new_pilot_run_refused(self):
        pilot = self.fx.root / "qc_runs" / "cutadapt_pilot_O5" / "job_600002"
        for mate in (1, 2):
            write_fastq_gz(pilot / f"11-01_R{mate}.trimmed.fastq.gz", ["r:1"], mate)
        self.assertRefusedUntouched(FQ_TRIM, r"without a COMPLETE marker", PILOT_DIR=str(pilot))

    def test_missing_inputs_and_relative_dir(self):
        self.assertRefusedUntouched(FQ_TRIM, r"must be an absolute path", PILOT_DIR="trimmed_pilot_11_01")
        (self.fx.root / "trimmed_pilot_11_01" / "11-01_R2.trimmed.fastq.gz").unlink()
        self.assertRefusedUntouched(FQ_TRIM, r"missing, empty or unreadable input")

    def test_wrong_fastqc_version(self):
        self.assertRefusedUntouched(FQ_TRIM, r"version mismatch", STUB_FASTQC_VERSION="FastQC v0.12.0")


class ValidateGzip(Base):
    def test_report_contents_and_historical_reports_untouched(self):
        bad = self.fx.raw / "20088XR-11-43_S148_L004_R1_001.fastq.gz"
        write_fastq_gz(bad, ["r:1"], 1, corrupt_tail=True)
        before = historical_only(snapshot(self.fx.root))
        res = self.fx.run(GZIP)
        self.assertSucceeded(res)
        rd = self.fx.run_dir(GZIP, self.fx.env(GZIP))
        names = sorted(p.name for p in self.fx.raw.glob("*.fastq.gz"))
        self.assertEqual((rd / "gzip_validation_input.txt").read_text().splitlines(), names)
        lines = (rd / "gzip_full_validation.tsv").read_text().splitlines()
        self.assertEqual(lines, sorted(lines))
        self.assertEqual(len(lines), len(names))
        self.assertIn(f"FAIL\t{bad.name}", lines)
        self.assertEqual(sum(l.startswith("PASS\t") for l in lines), len(names) - 1)
        self.assertFalse((rd / "gzip_full_validation.partial.tsv").exists())
        self.assertIn("FAIL: 1", res.stdout)
        self.assertEqual(historical_only(snapshot(self.fx.root)), before)

    def test_no_fastq_refused(self):
        for p in self.fx.raw.glob("*.fastq.gz"):
            p.unlink()
        self.assertRefusedUntouched(GZIP, r"no \*\.fastq\.gz files")

    def test_missing_raw_dir_refused(self):
        shutil.rmtree(self.fx.raw)
        self.assertRefusedUntouched(GZIP, r"raw FASTQ directory missing")


if __name__ == "__main__":
    unittest.main()
