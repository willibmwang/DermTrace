"""Hermetic tests for campaign_tools.py and confounding_analysis.py (synthetic fixtures; stdlib unittest).

Run from the procedure directory:  python3 -I -m unittest discover -s tests -p 'test_*.py'
Needs ../procedure-20261126/expand_dataset.py (or EXPAND_TOOL) and an r3 rivanna_qc_tools.py (RQC_TOOLS); the
compatibility test is skipped if RQC_TOOLS is not available.
"""
import csv, gzip, hashlib, json, math, os, random, shutil, subprocess, sys, tempfile, unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
CAMPAIGN, CONF = HERE / "campaign_tools.py", HERE / "confounding_analysis.py"
EXPAND = Path(os.environ.get("EXPAND_TOOL", HERE.parent / "procedure-20261126" / "expand_dataset.py"))
RQC = os.environ.get("RQC_TOOLS")
PY = sys.executable
R3_SHA = "d5cb81227ab3f8232f2d5c0ae300415feeeb9ca061d1a33cb4eb7dfda432fefb"
R2_SHA = "61888cb45011b4e3fd31c7a2d7c09fc654eb5362e2530dc278d64ed5b2c9aaee"
FPV = {"index_seq_hash": "a" * 64, "index_name_hash": "b" * 64, "index_seq_hash512": "c" * 128, "index_name_hash512": "d" * 128,
       "index_decoy_seq_hash": "e" * 64, "index_decoy_name_hash": "e" * 64}
CFG = {"salmon_version": "salmon 1.10.3", "salmon_quant_args": "-l A --validateMappings --gcBias --seqBias", "salmon_threads": "6",
       "fastp_args": "--detect_adapter_for_pe -l 25", "fastp_version": "fastp 0.23.4", "fastp_container": "f" * 64, "fastp_threads": "4",
       "index_kmer": "31", "index_fingerprints": "six SHA-256/SHA-512 hashes", "transcriptome": "GENCODE v46 transcripts"}
SYMS = ["MITF", "PMEL", "MLANA", "TYR", "KRT1", "KRT5", "KRT10", "KRT14", "PTPRC", "CD3E", "CD74", "COL1A1", "COL1A2", "DCN",
        "XIST", "RPS4Y1", "DDX3Y"] + [f"G{i}" for i in range(183)]


def run(*a):
    return subprocess.run([PY, "-I", *map(str, a)], capture_output=True, text=True)


def make_output(root, sid, tools_sha, rng, cfg=None, sif="f" * 64, fp=None, tx=None, scale=None):
    d = root / sid; (d / "quant" / "aux_info").mkdir(parents=True); (d / "fastp").mkdir()
    rows = [(t, 1000, 900.0, rng.uniform(0, 100) * (scale(t) if scale else 1.0)) for t in tx]
    s = sum(r[3] for r in rows)
    (d / "quant" / "quant.sf").write_text("Name\tLength\tEffectiveLength\tTPM\tNumReads\n" + "".join(f"{t}\t{l}\t{e}\t{v / s * 1e6:.6f}\t{v:.3f}\n" for t, l, e, v in rows))
    meta = {"num_valid_targets": len(rows), "num_processed": 490, "num_mapped": 170, "percent_mapped": 34.7, "frag_length_mean": 140.0, **(fp or FPV)}
    (d / "quant" / "aux_info" / "meta_info.json").write_text(json.dumps(meta))
    (d / "quant" / "lib_format_counts.json").write_text(json.dumps({"expected_format": "ISR", "compatible_fragment_ratio": 1.0}))
    (d / "fastp" / f"{sid}.json").write_text(json.dumps({"summary": {"before_filtering": {"total_reads": 1000}, "after_filtering": {"total_reads": 980, "q30_rate": 0.95}},
                                                       "adapter_cutting": {"adapter_trimmed_reads": 600}, "insert_size": {"peak": 128}, "duplication": {"rate": 0.4}}))
    r1 = root / f"{sid}_R1_001.fastq.gz"
    with gzip.open(r1, "wt") as fh:
        fh.write("@LH00150:1058:22XX:4:1101:1:1 1:N:0\nACGT\n+\nIIII\n")
    c = {k: {"value": v, "status": "x", "source": "x"} for k, v in (cfg or CFG).items()}
    (d / "provenance.json").write_text(json.dumps({"status": "complete", "specimen": sid, "core_id": sid, "slurm_job_id": "1", "script": "quantify_specimen.sbatch",
                                                   "tools_sha256": tools_sha, "lifted_exclusions": "none", "fastp_sif_sha256": sif, "r1": str(r1),
                                                   "salmon_step_modules": "gcc/14.2.0;salmon/1.10.3", "fastp_step_modules": "apptainer/1.4.5", "configuration": c}))
    return d


class Campaign(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory(); cls.T = Path(cls.tmp.name)
        cls.tx = []
        with gzip.open(cls.T / "fa.gz", "wt") as fh:
            for gi, sym in enumerate(SYMS):
                for k in range(2):
                    t = f"ENST{gi * 10 + k:011d}.1"; cls.tx.append(t)
                    fh.write(f">{t}|ENSG{gi:011d}.1|-|-|{sym}-20{k}|{sym}|1000|protein_coding|\nACGT\n")
        cls.common = ["--expect-seq-hash", "a" * 64, "--expect-name-hash", "b" * 64, "--expect-seq-hash512", "c" * 128,
                      "--expect-name-hash512", "d" * 128, "--expect-decoy-hash", "e" * 64, "--expect-sif-sha256", "f" * 64]

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    @unittest.skipUnless(RQC and Path(RQC).is_file(), "set RQC_TOOLS to an r2/r3 rivanna_qc_tools.py")
    def test_verify_compat(self):
        root = self.T / "compat"; root.mkdir(); rng = random.Random(1)
        ok = [make_output(root, s, sha, rng, tx=self.tx) for s, sha in (("GM01", R2_SHA), ("GM33", R3_SHA), ("GM34", R3_SHA))]
        def v(dirs, name):
            return run(CAMPAIGN, "verify-compat", "--rqc-tools", RQC, "--expect-targets", str(len(self.tx)), *self.common,
                       *sum((["--specimen", f"{d.name}={d}"] for d in dirs), []), "--json-out", self.T / f"{name}.json")
        r = v(ok, "c_ok"); self.assertEqual(r.returncode, 0, r.stdout + r.stderr); self.assertIn("COMPATIBILITY PASS", r.stdout)
        self.assertEqual(json.loads((self.T / "c_ok.json").read_text())["code_versions"], {"GM01": "r2", "GM33": "r3", "GM34": "r3"})
        cases = {"GM40": dict(cfg={**CFG, "salmon_quant_args": "-l A"}), "GM41": dict(sif="0" * 64),
                 "GM42": dict(fp={**FPV, "index_name_hash512": "0" * 128}), "GM35": dict(tools="9" * 64)}
        for sid, kw in cases.items():
            d = make_output(root, sid, kw.pop("tools", R3_SHA), rng, tx=self.tx, **kw)
            r = v(ok + [d], f"c_{sid}")
            self.assertEqual(r.returncode, 1, sid); self.assertIn(f"FAIL  {sid}", r.stdout)
        p = json.loads((ok[1] / "provenance.json").read_text()); p["status"] = "running"; (root / "GM36").mkdir()
        shutil.copytree(ok[1], root / "GM36", dirs_exist_ok=True); p["specimen"] = "GM36"; (root / "GM36" / "provenance.json").write_text(json.dumps(p))
        r = v(ok + [root / "GM36"], "c_running"); self.assertEqual(r.returncode, 1); self.assertIn("status is 'running'", r.stdout)

    def test_harmonized_and_qc(self):
        root = self.T / "harm"; root.mkdir(); rng = random.Random(2)
        ids = [f"GM{i:02d}" for i in range(1, 6)]
        dirs = [make_output(root, s, R3_SHA, rng, tx=self.tx) for s in ids]
        compat = self.T / "compat_h.json"
        compat.write_text(json.dumps({"code_versions": {s: "r3" for s in ids}}))
        out = self.T / "harm_out"
        r = run(CAMPAIGN, "harmonized", "--expand-tool", EXPAND, "--fasta", self.T / "fa.gz", *sum((["--specimen", f"{d.name}={d}"] for d in dirs), []),
                "--expect-order", ",".join(ids), "--compat-json", compat, "--out-dir", out)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        # cross-check against expand_dataset reproduce logic: rebuild one column independently
        sys.path.insert(0, str(EXPAND.parent)); import importlib; ed = importlib.import_module("expand_dataset")
        t2g, ginfo = ed.tx2gene(self.T / "fa.gz"); tpm, cnt = ed.aggregate(dirs[2] / "quant" / "quant.sf", t2g)
        with gzip.open(out / "tsaolab_melanoma_tpm_ensembl.tsv.gz", "rt") as fh:
            h = fh.readline().rstrip("\n").split("\t"); rows = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in fh}
        self.assertEqual(h[1:], ids)
        for g, v in tpm.items():
            self.assertEqual(rows[g][3], "%.4f" % v)
        self.assertEqual(run(CAMPAIGN, "harmonized", "--expand-tool", EXPAND, "--fasta", self.T / "fa.gz", "--specimen", f"GM01={dirs[0]}", "--out-dir", out).returncode, 2)
        qc = self.T / "harm_qc.tsv"
        r = run(CAMPAIGN, "qc-manifest", "--expand-tool", EXPAND, *sum((["--specimen", f"{d.name}={d}"] for d in dirs), []), "--compat-json", compat,
                "--matrices-dir", out, "--out-tsv", qc)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        q = list(csv.DictReader(open(qc), delimiter="\t"))
        self.assertEqual([x["sample_id"] for x in q], ids); self.assertEqual({x["code_version"] for x in q}, {"r3"})
        self.assertEqual(q[0]["sequencing_run"], "LH00150:1058"); self.assertEqual(q[0]["lane"], "L004")


class Confounding(unittest.TestCase):
    """Synthetic 42-specimen cohort with a planted lane effect and a planted marker-driven biological axis."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory(); T = cls.T = Path(cls.tmp.name); rng = random.Random(5)
        cols = [f"GM{i:02d}" for i in range(1, 33)] + ["GM50"] + ["GM33", "GM34", "GM35", "GM36", "GM38", "GM39", "GM40", "GM41", "GM42"]
        lanes = {c: ("L004" if i % 3 == 0 else "L006" if i % 3 == 1 else "L008") for i, c in enumerate(cols)}
        mel = {c: rng.gauss(0, 1) for c in cols}
        genes = [(f"ENSG{i:011d}", s) for i, s in enumerate(SYMS)]
        def value(c, gi, sym):
            base = 3 + (gi % 7) * 0.3 + rng.gauss(0, 0.15)
            if sym in ("MITF", "PMEL", "MLANA", "TYR"): base += 1.5 * mel[c]
            if sym in ("KRT1", "KRT5", "KRT10", "KRT14"): base -= 1.5 * mel[c]
            if gi % 5 == 0 and gi > 20: base += 1.2 * (lanes[c] == "L008")
            return math.expm1(max(0.1, base)) * 10
        vals = {c: {g: value(c, gi, s) for gi, (g, s) in enumerate(genes)} for c in cols}
        def write_rep(d, cs):
            d.mkdir(); hdr = "\t".join(cs)
            with gzip.open(d / "tsaolab_melanoma_tpm_ensembl.tsv.gz", "wt") as fh:
                fh.write(f"ensembl_gene_id\t{hdr}\n" + "".join(f"{g}\t" + "\t".join(f"{vals[c][g]:.4f}" for c in cs) + "\n" for g, _ in genes))
            with gzip.open(d / "tsaolab_melanoma_tpm_symbol.tsv.gz", "wt") as fh:
                fh.write(f"gene_symbol\t{hdr}\n" + "".join(f"{s}\t" + "\t".join(f"{vals[c][g]:.4f}" for c in cs) + "\n" for g, s in sorted(genes, key=lambda t: t[1])))
            with gzip.open(d / "tsaolab_melanoma_genes.tsv.gz", "wt") as fh:
                fh.write("ensembl_gene_id\tgene_symbol\tgene_type\n" + "".join(f"{g}\t{s}\tprotein_coding\n" for g, s in genes))
        write_rep(T / "orig", cols[:33]); write_rep(T / "mixed", cols); write_rep(T / "harm", cols)
        def qc(p, origin_new):
            with open(p, "w") as fh:
                fh.write("sample_id\tsequencing_run\tlane\tsex_inferred\tquality_flag\tpipeline_origin\tinsert_size_peak\tfastp_duplication_rate\tmapping_rate_pct\tpct_adapter_trimmed\treads_raw\n")
                for i, c in enumerate(cols):
                    fh.write("\t".join([c, "LH00150:1062" if lanes[c] == "L008" else "LH00150:1058", lanes[c], "male" if i % 2 else "female",
                                        "suspect" if c in ("GM10", "GM13", "GM15", "GM21", "GM40", "GM41") else "ok",
                                        origin_new if c in cols[33:] else ("collaborator_committed" if origin_new == "rivanna_qc_r3" else origin_new),
                                        str(120 + i % 20), f"{0.2 + 0.4 * (lanes[c] == 'L008') + rng.uniform(0, 0.1):.3f}", f"{33 + rng.uniform(-3, 3):.2f}",
                                        f"{55 + rng.uniform(0, 10):.2f}", str(60000000 + i * 100000)]) + "\n")
        qc(T / "qc_mixed.tsv", "rivanna_qc_r3"); qc(T / "qc_harm.tsv", "rivanna_qc_harmonized")
        cls.out = T / "conf"
        cls.r = run(CONF, "--original33", T / "orig", "--mixed42", T / "mixed", "--harmonized42", T / "harm", "--qc-mixed", T / "qc_mixed.tsv",
                    "--qc-harmonized", T / "qc_harm.tsv", "--permutations", "300", "--subsamples", "15", "--out-dir", cls.out)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_runs_and_outputs(self):
        self.assertEqual(self.r.returncode, 0, self.r.stdout + self.r.stderr)
        for f in ("pc_covariate_associations.tsv", "covariate_relations.tsv", "technical_vs_marker_r2.tsv", "clustering_stability.tsv",
                  "pca_stability.tsv", "cross_representation.tsv", "per_specimen_original_vs_harmonized.tsv", "suspect_identity_check.tsv", "summary.json"):
            self.assertTrue((self.out / f).is_file(), f)

    def test_planted_effects_recovered(self):
        a = list(csv.DictReader(open(self.out / "pc_covariate_associations.tsv"), delimiter="\t"))
        lane = [r for r in a if r["analysis"] == "harmonized42" and r["covariate"] == "lane"]
        self.assertTrue(any(float(r["stat"]) > 0.5 and float(r["p_perm"]) < 0.01 for r in lane), "planted lane effect not detected")
        tb = [r for r in csv.DictReader(open(self.out / "technical_vs_marker_r2.tsv"), delimiter="\t") if r["analysis"] == "harmonized42"]
        self.assertTrue(any(float(r["r2_markers"]) > 0.8 for r in tb), "planted marker axis not explained by marker scores")

    def test_identical_representations_agree_and_suspects_dropped(self):
        x = {r["comparison"]: r for r in csv.DictReader(open(self.out / "cross_representation.tsv"), delimiter="\t")}
        m = x["mixed42 vs harmonized42 (same 42 specimens; mixing effect)"]
        self.assertAlmostEqual(float(m["PC1_max_abs_r"]), 1.0, places=6); self.assertAlmostEqual(float(m["kmeans2_ari"]), 1.0, places=6)
        s = json.loads((self.out / "summary.json").read_text())
        self.assertEqual(s["representations"]["harmonized42_minus_suspect"]["n_specimens"], 36)
        self.assertEqual(s["representations"]["original33_minus_suspect"]["n_specimens"], 29)
        self.assertEqual(s["suspect_list"], ["GM10", "GM13", "GM15", "GM21", "GM40", "GM41"])

    def test_refuses_overwrite(self):
        r = run(CONF, "--original33", self.T / "orig", "--mixed42", self.T / "mixed", "--harmonized42", self.T / "harm", "--qc-mixed", self.T / "qc_mixed.tsv",
                "--qc-harmonized", self.T / "qc_harm.tsv", "--out-dir", self.out)
        self.assertEqual(r.returncode, 2)


if __name__ == "__main__":
    unittest.main()
