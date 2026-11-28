"""Hermetic tests for expand_dataset.py and mixed_matrix_analysis.py (synthetic fixtures; stdlib unittest).

Run: python3 -I -m unittest discover -s tests -p 'test_*.py'   (from the procedure directory)
"""
import csv, gzip, io, json, math, os, random, subprocess, sys, tempfile, unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
EXPAND = HERE / "expand_dataset.py"
MIXED = HERE / "mixed_matrix_analysis.py"
PY = sys.executable

# 12 genes, 30 transcripts; G3/G4 share symbol SHR; G11 is zero in all committed samples (added later);
# G12 has a brand-new symbol and is also zero in all committed samples.
GENES = [(f"ENSG{n:011d}", sym, typ) for n, sym, typ in [
    (1, "XIST", "lncRNA"), (2, "RPS4Y1", "protein_coding"), (3, "SHR", "protein_coding"), (4, "SHR", "lncRNA"),
    (5, "MITF", "protein_coding"), (6, "PMEL", "protein_coding"), (7, "KRT14", "protein_coding"), (8, "PTPRC", "protein_coding"),
    (9, "DDX3Y", "protein_coding"), (10, "ACTB", "protein_coding"), (11, "SHR", "protein_coding"), (12, "NEWSYM", "misc_RNA")]]
TX = []
for gi, (g, sym, typ) in enumerate(GENES):
    for k in range(2 + gi % 3):
        TX.append((f"ENST{gi * 10 + k:011d}.1", g, sym, typ if k == 0 else "retained_intron"))
COMMITTED = ["GM01", "GM02", "GM03"]
NEW = ["GM33", "GM34"]


def quant_rows(seed, zero_genes):
    rng = random.Random(seed)
    rows = []
    for t, g, *_ in TX:
        z = g in zero_genes
        rows.append((t, 1000, 900.0, 0.0 if z else round(rng.uniform(0, 5000), 6), 0.0 if z else round(rng.uniform(0, 300), 3)))
    s = sum(r[3] for r in rows)
    return [(t, l, e, tpm / s * 1e6, n) for t, l, e, tpm, n in rows]


def write_quant(p, rows, gz=False):
    text = "Name\tLength\tEffectiveLength\tTPM\tNumReads\n" + "".join(f"{t}\t{l}\t{e}\t{tpm:.6f}\t{n:.3f}\n" for t, l, e, tpm, n in rows)
    p.parent.mkdir(parents=True, exist_ok=True)
    if gz:
        with gzip.open(p, "wt") as fh:
            fh.write(text)
    else:
        p.write_text(text)


def naive_matrices(samples):
    """Independent reference implementation of build_dermtrace.py (plain dict sums of parsed quant text)."""
    tpm, cnt = {}, {}
    for sid, rows in samples.items():
        for t, l, e, x, n in rows:
            g = next(gg for tt, gg, *_ in TX if tt == t)
            tpm.setdefault(g, {}).setdefault(sid, 0.0); cnt.setdefault(g, {}).setdefault(sid, 0.0)
            tpm[g][sid] += float(f"{x:.6f}"); cnt[g][sid] += float(f"{n:.3f}")
    keep = sorted(g for g in cnt if sum(cnt[g].values()) > 0)
    sym = {g: s for g, s, _ in GENES}
    stpm, scnt = {}, {}
    for g in keep:
        for sid in samples:
            stpm.setdefault(sym[g], {}).setdefault(sid, 0.0); scnt.setdefault(sym[g], {}).setdefault(sid, 0.0)
            stpm[sym[g]][sid] += tpm[g][sid]; scnt[sym[g]][sid] += cnt[g][sid]
    return keep, tpm, cnt, stpm, scnt


def gzwrite(p, text):
    with gzip.open(p, "wt") as fh:
        fh.write(text)


class Expansion(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        T = cls.T = Path(cls.tmp.name)
        with gzip.open(T / "tx.fa.gz", "wt") as fh:
            for t, g, sym, typ in TX:
                fh.write(f">{t}|{g}.3|OTTHUMG|OTTHUMT|{sym}-201|{sym}|1000|{typ}|\nACGT\n")
        cls.q = {sid: quant_rows(i, {GENES[10][0], GENES[11][0]}) for i, sid in enumerate(COMMITTED)}
        cls.q.update({"GM33": quant_rows(33, set()), "GM34": quant_rows(34, {GENES[11][0]})})
        for sid, rows in cls.q.items():
            write_quant(T / "quant" / sid / "quant" / "quant.sf", rows)
            write_quant(T / "committed_quant" / f"{sid}.sf.gz", rows, gz=True)
        keep, tpm, cnt, stpm, scnt = naive_matrices({s: cls.q[s] for s in COMMITTED})
        C = cls.C = T / "committed"; C.mkdir()
        hdr = "\t".join(COMMITTED)
        gzwrite(C / "tsaolab_melanoma_tpm_ensembl.tsv.gz", f"ensembl_gene_id\t{hdr}\n" + "".join(f"{g}\t" + "\t".join(f"{tpm[g][s]:.4f}" for s in COMMITTED) + "\n" for g in keep))
        gzwrite(C / "tsaolab_melanoma_counts_ensembl.tsv.gz", f"ensembl_gene_id\t{hdr}\n" + "".join(f"{g}\t" + "\t".join(f"{cnt[g][s]:.3f}" for s in COMMITTED) + "\n" for g in keep))
        gzwrite(C / "tsaolab_melanoma_tpm_symbol.tsv.gz", f"gene_symbol\t{hdr}\n" + "".join(f"{k}\t" + "\t".join(f"{stpm[k][s]:.4f}" for s in COMMITTED) + "\n" for k in sorted(stpm)))
        gzwrite(C / "tsaolab_melanoma_counts_symbol.tsv.gz", f"gene_symbol\t{hdr}\n" + "".join(f"{k}\t" + "\t".join(f"{scnt[k][s]:.3f}" for s in COMMITTED) + "\n" for k in sorted(scnt)))
        info = {g: (s, t) for g, s, t in GENES}
        gzwrite(C / "tsaolab_melanoma_genes.tsv.gz", "ensembl_gene_id\tgene_symbol\tgene_type\n" + "".join(f"{g}\t{info[g][0]}\t{info[g][1]}\n" for g in keep))
        cls.keep = keep

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def run_expand(self, *args):
        return subprocess.run([PY, "-I", str(EXPAND), *args], capture_output=True, text=True)

    def build(self, name, new=NEW, extra=()):
        args = ["matrices", "--committed-dir", str(self.C), "--fasta", str(self.T / "tx.fa.gz"), "--out-dir", str(self.T / name), *extra]
        for s in new:
            args += ["--new", f"{s}={self.T / 'quant' / s}"]
        return self.run_expand(*args)

    def read(self, d, m):
        with gzip.open(self.T / d / f"tsaolab_melanoma_{m}.tsv.gz", "rt") as fh:
            return fh.read().rstrip("\n").split("\n")

    def test_reproduce_pass_and_fail(self):
        args = ["reproduce", "--committed-dir", str(self.C), "--fasta", str(self.T / "tx.fa.gz")]
        for s in COMMITTED:
            args += ["--sample", f"{s}={self.T / 'committed_quant' / (s + '.sf.gz')}"]
        r = self.run_expand(*args)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("REPRODUCTION PASS", r.stdout)
        r = self.run_expand("reproduce", "--committed-dir", str(self.C), "--fasta", str(self.T / "tx.fa.gz"), "--sample", f"GM01={self.T / 'committed_quant' / 'GM02.sf.gz'}")
        self.assertEqual(r.returncode, 1)

    def test_expansion_layout_and_integrity(self):
        r = self.build("exp1", extra=["--expect-new", "GM33,GM34"])
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        _, tpm, cnt, stpm, scnt = naive_matrices({s: self.q[s] for s in COMMITTED + NEW})
        for m, fmt, ref in (("tpm_ensembl", "%.4f", tpm), ("counts_ensembl", "%.3f", cnt), ("tpm_symbol", "%.4f", stpm), ("counts_symbol", "%.3f", scnt)):
            lines = self.read("exp1", m)
            with gzip.open(self.C / f"tsaolab_melanoma_{m}.tsv.gz", "rt") as fh:
                orig = fh.read().rstrip("\n").split("\n")
            self.assertEqual(lines[0], orig[0] + "\tGM33\tGM34", m)
            o = {l.split("\t", 1)[0]: l for l in orig[1:]}
            for l in lines[1:]:
                f = l.split("\t"); k = f[0]
                if k in o:
                    self.assertEqual("\t".join(f[:4]), o[k], f"{m} committed text changed for {k}")
                else:
                    self.assertEqual(f[1:4], [fmt % 0.0] * 3, f"{m} added row {k} not zero in committed columns")
                for j, s in enumerate(NEW, 4):
                    self.assertEqual(f[j], fmt % ref[k][s], f"{m} {k} {s}")
            keys = [l.split("\t", 1)[0] for l in lines[1:]]
            self.assertEqual(keys, sorted(keys))
            self.assertTrue(set(o) <= set(keys))
        ens = [l.split("\t", 1)[0] for l in self.read("exp1", "tpm_ensembl")[1:]]
        self.assertIn(GENES[10][0], ens)                       # G11 added (non-zero in GM33)
        self.assertIn(GENES[11][0], ens)                       # G12 added (non-zero in GM33, zero in GM34)
        syms = [l.split("\t", 1)[0] for l in self.read("exp1", "tpm_symbol")[1:]]
        self.assertIn("NEWSYM", syms)
        man = json.loads((self.T / "exp1" / "expansion_manifest.json").read_text())
        self.assertEqual(man["verdict"], "PASS")
        self.assertEqual(man["added_genes"], sorted([GENES[10][0], GENES[11][0]]))
        self.assertEqual(man["added_symbols"], ["NEWSYM"])  # SHR already existed: G11 joins the existing SHR row
        self.assertEqual(man["column_order"], COMMITTED + NEW)
        with gzip.open(self.T / "exp1" / "tsaolab_melanoma_genes.tsv.gz", "rt") as fh:
            g = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in list(fh)[1:]}
        self.assertEqual(g[GENES[11][0]][1:], ["NEWSYM", "misc_RNA"])  # first transcript in FASTA order

    def test_refusals(self):
        self.assertEqual(self.build("exp_ok").returncode, 0)
        r = self.build("exp_ok"); self.assertEqual(r.returncode, 2); self.assertIn("refusing to overwrite", r.stderr)
        r = self.build("exp_dup", new=["GM01"]); self.assertEqual(r.returncode, 2); self.assertIn("already in the committed matrix", r.stderr)
        r = self.build("exp_order", extra=["--expect-new", "GM34,GM33"]); self.assertEqual(r.returncode, 2); self.assertIn("!= expected", r.stderr)
        bad = self.T / "quant_bad" / "GM33" / "quant" / "quant.sf"
        write_quant(bad, self.q["GM33"][:-1] + [("ENST99999999999.1", 1, 1.0, 1.0, 1.0)])
        r = self.run_expand("matrices", "--committed-dir", str(self.C), "--fasta", str(self.T / "tx.fa.gz"), "--new", f"GM33={bad.parent.parent}", "--out-dir", str(self.T / "exp_badtx"))
        self.assertEqual(r.returncode, 2); self.assertIn("not in the GENCODE FASTA", r.stderr)
        self.assertFalse((self.T / "exp_badtx").exists())
        r = self.run_expand("matrices", "--committed-dir", str(self.C), "--fasta", str(self.T / "tx.fa.gz"), "--fasta-md5", "0" * 32, "--new", f"GM33={self.T / 'quant' / 'GM33'}", "--out-dir", str(self.T / "exp_md5"))
        self.assertEqual(r.returncode, 2); self.assertIn("MD5", r.stderr)

    def test_qc_manifest_and_mixed_analysis(self):
        self.assertEqual(self.build("exp_qc").returncode, 0)
        T = self.T; qroot = T / "qroot"
        for i, s in enumerate(NEW):
            d = qroot / s; (d / "fastp").mkdir(parents=True); (d / "quant" / "aux_info").mkdir(parents=True)
            write_quant(d / "quant" / "quant.sf", self.q[s])
            r1 = T / "raw" / f"{s}_R1_001.fastq.gz"; r1.parent.mkdir(exist_ok=True)
            with gzip.open(r1, "wt") as fh:
                fh.write(f"@LH00150:1070:22ABC:{6 + i}:1101:1:1 1:N:0\nACGT\n+\nIIII\n")
            (d / "provenance.json").write_text(json.dumps({"status": "complete", "specimen": s, "slurm_job_id": str(900 + i), "script": "quantify_specimen.sbatch",
                                                           "tools_sha256": "ab" * 32, "r1": str(r1), "configuration": {"fastp_version": {"value": "fastp 0.23.4"}}}))
            (d / "fastp" / f"{s}.json").write_text(json.dumps({"summary": {"before_filtering": {"total_reads": 1000}, "after_filtering": {"total_reads": 980, "read1_mean_length": 140, "read2_mean_length": 138, "q30_rate": 0.95}},
                                                               "adapter_cutting": {"adapter_trimmed_reads": 600}, "insert_size": {"peak": 130}, "duplication": {"rate": 0.4}}))
            (d / "quant" / "aux_info" / "meta_info.json").write_text(json.dumps({"num_processed": 490, "num_mapped": 170 + 40 * i, "percent_mapped": 25.0 + 15 * i, "frag_length_mean": 140.2}))
            (d / "quant" / "lib_format_counts.json").write_text(json.dumps({"expected_format": "ISR", "compatible_fragment_ratio": 1.0}))
        cols = ["sample_id", "core_id", "instrument", "sequencing_run", "flowcell", "lane", "reads_raw", "reads_after_trim", "pct_reads_kept", "pct_adapter_trimmed",
                "mean_read_length_after_trim", "insert_size_peak", "pct_q30_after_trim", "fastp_duplication_rate", "salmon_fragments", "salmon_mapped",
                "mapping_rate_pct", "library_type", "compatible_fragment_ratio", "frag_length_mean", "driver_ref_reads"]
        with open(T / "run_qc.tsv", "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for i, s in enumerate(COMMITTED):
                fh.write("\t".join([s, f"20088XR-11-{i+1:02d}", "Illumina_NovaSeqX", "LH00150:1062", "FC", "L008", "1000", "990", "99.0", "60.0", "140.0", "126", "97.0", str(0.3 + 0.1 * i), "495", "170", str(34.0 + i), "ISR", "1.0", "139.0", "0"]) + "\n")
        with open(T / "meta.csv", "w") as fh:
            fh.write("sample_id,quality_flag,qc_notes\n" + "".join(f"{s},ok,\n" for s in COMMITTED))
        with open(T / "sheet.csv", "w") as fh:
            fh.write("sample_id,core_id\n" + "".join(f"{s},20088XR-11-{s[2:]}\n" for s in COMMITTED + NEW))
        out = T / "qc42.tsv"
        r = self.run_expand("qc-manifest", "--expansion-manifest", str(T / "exp_qc" / "expansion_manifest.json"), "--expanded-dir", str(T / "exp_qc"),
                            "--run-qc", str(T / "run_qc.tsv"), "--metadata", str(T / "meta.csv"), "--sheet", str(T / "sheet.csv"), "--quant-root", str(qroot), "--out-tsv", str(out))
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        rows = list(csv.DictReader(open(out), delimiter="\t"))
        self.assertEqual([x["sample_id"] for x in rows], COMMITTED + NEW)
        new = {x["sample_id"]: x for x in rows if x["pipeline_origin"] == "rivanna_qc_r3"}
        self.assertEqual(set(new), set(NEW))
        self.assertEqual(new["GM33"]["sequencing_run"], "LH00150:1070"); self.assertEqual(new["GM33"]["lane"], "L006")
        self.assertEqual(new["GM33"]["quality_flag"], "suspect"); self.assertEqual(new["GM34"]["quality_flag"], "ok")   # <30 % mapping rule
        self.assertEqual(new["GM33"]["provenance_status"], "complete"); self.assertEqual(new["GM33"]["slurm_job_id"], "900")
        self.assertEqual(new["GM33"]["pct_adapter_trimmed"], "60.0")
        self.assertEqual({x["pipeline_origin"] for x in rows if x["sample_id"] in COMMITTED}, {"collaborator_committed"})
        # mixed analysis: frozen gene set = protein-coding committed genes; loadings arbitrary unit vectors; frozen scores
        # computed here from the committed columns, so the script's re-projection must reproduce them
        with gzip.open(T / "exp_qc" / "tsaolab_melanoma_tpm_ensembl.tsv.gz", "rt") as fh:
            h = fh.readline().rstrip("\n").split("\t")[1:]; m = {l.split("\t")[0]: [float(x) for x in l.rstrip("\n").split("\t")[1:]] for l in fh}
        fz = [g for g, s, t in GENES if t == "protein_coding" and g in self.keep][:6]
        rng = random.Random(7); load = {g: [rng.uniform(-1, 1) for _ in range(3)] for g in fz}
        mu = {g: sum(math.log1p(v) for v in m[g][:3]) / 3 for g in fz}
        scores = {s: {f"PC{k+1}": sum((math.log1p(m[g][i]) - mu[g]) * load[g][k] for g in fz) for k in range(3)} for i, s in enumerate(COMMITTED)}
        for s in scores: scores[s].update(PC4=0.0, PC5=0.0)
        (T / "frozen.json").write_text(json.dumps({"scores": scores, "score_sd": {f"PC{k}": 1.0 for k in range(1, 6)}}))
        with open(T / "ctx.tsv", "w") as fh:
            fh.write("# synthetic\ngene_id\tin_pca_gene_set\tpc1_loading\tpc2_loading\tpc3_loading\n")
            for g in self.keep:
                fh.write(f"{g}\t{int(g in fz)}\t" + "\t".join(f"{x:.6f}" for x in load.get(g, [0, 0, 0])) + "\n")
        r = subprocess.run([PY, "-I", str(MIXED), "--expanded-dir", str(T / "exp_qc"), "--committed-dir", str(self.C), "--qc-manifest", str(out),
                            "--cohort-context", str(T / "ctx.tsv"), "--frozen-scores", str(T / "frozen.json"), "--permutations", "200", "--out-dir", str(T / "mixed")],
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("INTEGRITY+REPRODUCTION PASS", r.stdout)
        s = json.loads((T / "mixed" / "summary.json").read_text())
        self.assertTrue(s["integrity"]["committed_values_identical"])
        self.assertLess(s["frozen_projection"]["max_abs_reproduction_error_33"], 1e-4)
        self.assertEqual(s["integrity"]["n_new_columns"], 2)


if __name__ == "__main__":
    unittest.main()
