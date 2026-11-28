#!/usr/bin/env python3
"""Harmonized-quantification campaign tools (Python >= 3.9, stdlib only).

Subcommands
  verify-compat   every specimen output is reusable: structurally valid (r3 validate-quant-dir rules), provenance
                  complete for the right specimen, produced by r2 or r3 quantify_specimen.sbatch (identical
                  quantification code), and EVERY recorded configuration value identical to the reference
                  configuration (Salmon/fastp versions and arguments, threads, container SHA-256, index
                  fingerprints and k, transcriptome, module sets); quant meta_info fingerprints equal the
                  expected collaborator values. Any difference fails closed (exit 1); nothing is reused silently.
  harmonized      42-specimen matrices built from scratch from r2/r3 outputs with the collaborator recipe (via
                  expand_dataset.py of procedure-20261126); compared row-by-row with the mixed matrix; never
                  overwrites, never intersects gene sets silently.
  qc-manifest     42-row technical QC/provenance manifest, all from r2/r3 outputs (same rules as the mixed one).
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.util
import io
import json
import os
import platform
import socket
import sys
from datetime import datetime, timezone
from pathlib import Path

R2_TOOLS_SHA = "61888cb45011b4e3fd31c7a2d7c09fc654eb5362e2530dc278d64ed5b2c9aaee"
R3_TOOLS_SHA = "d5cb81227ab3f8232f2d5c0ae300415feeeb9ca061d1a33cb4eb7dfda432fefb"
CODE = {R2_TOOLS_SHA: "r2", R3_TOOLS_SHA: "r3"}
CONFIG_KEYS = ["salmon_version", "salmon_quant_args", "salmon_threads", "fastp_args", "fastp_version", "fastp_container",
               "fastp_threads", "index_kmer", "index_fingerprints", "transcriptome"]
TOP_KEYS = ["fastp_sif_sha256", "salmon_step_modules", "fastp_step_modules", "script"]
FP = ["index_seq_hash", "index_name_hash", "index_seq_hash512", "index_name_hash512", "index_decoy_seq_hash", "index_decoy_name_hash"]
MIN_MAPPING_RATE = 30.0
SEX_GENES = dict(xist="XIST", chrY=["RPS4Y1", "DDX3Y", "KDM5D", "EIF1AY", "UTY", "USP9Y", "ZFY"])
INSTRUMENTS = {"LH": "Illumina_NovaSeqX", "A0": "Illumina_NovaSeq", "NB": "Illumina_NextSeq", "AV": "Element_AVITI"}


def fail(msg, code=2):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


def digest(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def write_new(path, data):
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    with os.fdopen(fd, "wb") as fh:
        fh.write(data if isinstance(data, bytes) else data.encode())


def parse_dirs(specs):
    out = []
    for s in specs:
        sid, _, p = s.partition("=")
        if not sid or not p:
            fail(f"expected ID=DIR, got {s!r}")
        out.append((sid, Path(p)))
    if len({s for s, _ in out}) != len(out):
        fail("duplicate specimen IDs")
    return out


def run_info():
    return {"utc": datetime.now(timezone.utc).isoformat(), "host": socket.gethostname(), "python": platform.python_version(),
            "script_sha256": digest(__file__)}


# ------------------------------------------------------------------------------------------ verify-compat
def provenance_signature(p):
    cfg = p.get("configuration", {})
    sig = {f"configuration.{k}": (cfg.get(k) or {}).get("value") for k in CONFIG_KEYS}
    sig.update({k: p.get(k) for k in TOP_KEYS})
    return sig


def cmd_verify(a):
    rq = load_module(a.rqc_tools, "rqc_tools")
    if digest(a.rqc_tools) not in CODE:
        fail("--rqc-tools is neither the r2 nor the r3 rivanna_qc_tools.py")
    expect_fp = dict(zip(FP, [a.expect_seq_hash, a.expect_name_hash, a.expect_seq_hash512, a.expect_name_hash512,
                              a.expect_decoy_hash, a.expect_decoy_hash]))
    specs = parse_dirs(a.specimen)
    rows, sigs = [], {}
    for sid, d in specs:
        probs = []
        p_struct, facts = rq.check_quant_dir(d, stage=False, specimen=sid, expect_targets=a.expect_targets, expect_libtype=a.expect_libtype)
        probs += [f"structure: {x}" for x in p_struct]
        prov = json.loads((d / "provenance.json").read_text()) if (d / "provenance.json").is_file() else {}
        code = CODE.get(prov.get("tools_sha256"), "unknown")
        if code == "unknown":
            probs.append(f"tools_sha256 {prov.get('tools_sha256')!r} is neither r2 nor r3")
        if prov.get("script") != "quantify_specimen.sbatch":
            probs.append(f"script {prov.get('script')!r}")
        if prov.get("lifted_exclusions") != "none":
            probs.append(f"lifted exclusions {prov.get('lifted_exclusions')!r}")
        if prov.get("fastp_sif_sha256") != a.expect_sif_sha256:
            probs.append(f"fastp SIF {prov.get('fastp_sif_sha256')!r} != expected")
        meta_p = d / "quant" / "aux_info" / "meta_info.json"
        meta = json.loads(meta_p.read_text()) if meta_p.is_file() else {}
        bad_fp = [k for k in FP if meta.get(k) != expect_fp[k]]
        if bad_fp:
            probs.append(f"meta_info fingerprints differ: {bad_fp}")
        sig = provenance_signature(prov)
        sigs[sid] = sig
        rows.append({"specimen": sid, "dir": str(d), "code": code, "slurm_job_id": prov.get("slurm_job_id"), "problems": probs,
                     "percent_mapped": facts.get("percent_mapped"), "num_processed": facts.get("num_processed")})
    ref_sid = a.reference_specimen or specs[0][0]
    if ref_sid not in sigs:
        fail(f"reference specimen {ref_sid} not among --specimen")
    ref = sigs[ref_sid]
    missing_ref = [k for k, v in ref.items() if v in (None, "")]
    for r in rows:
        diff = {k: (sigs[r["specimen"]].get(k), v) for k, v in ref.items() if sigs[r["specimen"]].get(k) != v}
        if diff:
            r["problems"].append("configuration differs from " + ref_sid + ": " + json.dumps(diff))
    verdict = "PASS" if not missing_ref and all(not r["problems"] for r in rows) else "FAIL"
    rep = {"kind": "harmonization-compatibility", "verdict": verdict, "reference_specimen": ref_sid, "reference_configuration": ref,
           "reference_configuration_missing_keys": missing_ref, "specimens": rows,
           "code_versions": {r["specimen"]: r["code"] for r in rows}, "run": run_info()}
    if a.json_out:
        write_new(a.json_out, json.dumps(rep, indent=1) + "\n")
    for r in rows:
        print(f"{'OK  ' if not r['problems'] else 'FAIL'}  {r['specimen']:<5} code={r['code']} job={r['slurm_job_id']} map={r['percent_mapped']}"
              + ("" if not r["problems"] else "  " + " | ".join(r["problems"])[:400]))
    if missing_ref:
        print(f"FAIL  reference configuration incomplete: {missing_ref}")
    print(f"COMPATIBILITY {verdict}: {len(rows)} specimen outputs; reference configuration from {ref_sid}")
    sys.exit(0 if verdict == "PASS" else 1)


# ------------------------------------------------------------------------------------------ harmonized
def cmd_harmonized(a):
    ed = load_module(a.expand_tool, "expand_dataset")
    out = Path(a.out_dir)
    if os.path.lexists(out):
        fail(f"output exists: {out}")
    if a.fasta_md5 and hashlib.md5(open(a.fasta, "rb").read()).hexdigest() != a.fasta_md5:
        fail("FASTA MD5 mismatch")
    t2g, ginfo = ed.tx2gene(a.fasta)
    specs = parse_dirs(a.specimen)
    ids = [s for s, _ in specs]
    if a.expect_order and ids != a.expect_order.split(","):
        fail("specimen order differs from --expect-order")
    agg = {sid: ed.aggregate(ed.quant_file(d), t2g) for sid, d in specs}
    genes = sorted({g for sid in ids for g, c in agg[sid][1].items()} )
    keep = [g for g in genes if sum(agg[sid][1].get(g, 0.0) for sid in ids) > 0]   # collaborator filter, over all 42
    syms = sorted({ginfo[g][0] for g in keep})
    sym_vals = {sid: (ed.symbol_sum(agg[sid][0], keep, ginfo), ed.symbol_sum(agg[sid][1], keep, ginfo)) for sid in ids}
    def mat(first, keys, f, fmt):
        return "\n".join(["\t".join([first] + ids)] + [k + "\t" + "\t".join(fmt % f(s, k) for s in ids) for k in keys]) + "\n"
    texts = {"tpm_ensembl": mat("ensembl_gene_id", keep, lambda s, g: agg[s][0].get(g, 0.0), "%.4f"),
             "counts_ensembl": mat("ensembl_gene_id", keep, lambda s, g: agg[s][1].get(g, 0.0), "%.3f"),
             "tpm_symbol": mat("gene_symbol", syms, lambda s, k: sym_vals[s][0].get(k, 0.0), "%.4f"),
             "counts_symbol": mat("gene_symbol", syms, lambda s, k: sym_vals[s][1].get(k, 0.0), "%.3f"),
             "genes": "ensembl_gene_id\tgene_symbol\tgene_type\n" + "".join(f"{g}\t{ginfo[g][0]}\t{ginfo[g][1]}\n" for g in keep)}
    cmp_ = {}
    if a.mixed_dir:
        for m in ("tpm_ensembl", "tpm_symbol"):
            h, rows, order = ed.read_matrix(Path(a.mixed_dir) / f"tsaolab_melanoma_{m}.tsv.gz")
            if h[1:] != ids:
                fail(f"mixed {m} column order differs from the harmonized order; refusing a misaligned comparison")
            mine = keep if m == "tpm_ensembl" else syms
            cmp_[m] = {"rows_harmonized": len(mine), "rows_mixed": len(order), "only_in_harmonized": len(set(mine) - set(order)),
                       "only_in_mixed": len(set(order) - set(mine)), "shared": len(set(mine) & set(order)),
                       "examples_only_in_mixed": sorted(set(order) - set(mine))[:10], "examples_only_in_harmonized": sorted(set(mine) - set(order))[:10]}
    sums = {sid: round(sum(agg[sid][0].values()), 2) for sid in ids}
    if any(abs(v - 1e6) > 10 for v in sums.values()):
        fail(f"TPM column sums not ~1e6: {sums}", 1)
    out.mkdir(parents=True)
    files = {}
    for m, t in texts.items():
        name = f"tsaolab_melanoma_{m}.tsv.gz"
        write_new(out / name, ed.gz_text(t))
        files[name] = digest(out / name)
    man = {"kind": "harmonized-matrices", "recipe": "collaborator build_dermtrace.py, applied to all 42 r2/r3 outputs",
           "column_order": ids, "code_versions": json.loads(Path(a.compat_json).read_text())["code_versions"] if a.compat_json else None,
           "rows": {"ensembl": len(keep), "symbol": len(syms)}, "tpm_column_sums": sums, "comparison_with_mixed": cmp_,
           "inputs": {sid: {"dir": str(d), "quant_sha256": digest(ed.quant_file(d))} for sid, d in specs},
           "fasta_sha256": digest(a.fasta), "outputs": files, "run": run_info()}
    write_new(out / "harmonized_manifest.json", json.dumps(man, indent=1) + "\n")
    print(f"harmonized: {len(keep)} genes, {len(syms)} symbols x {len(ids)} specimens -> {out}")
    for m, c in cmp_.items():
        print(f"  vs mixed {m}: shared {c['shared']}, only harmonized {c['only_in_harmonized']}, only mixed {c['only_in_mixed']}")


# ------------------------------------------------------------------------------------------ qc-manifest
def cmd_qc(a):
    ed = load_module(a.expand_tool, "expand_dataset")
    if os.path.lexists(a.out_tsv):
        fail(f"refusing to overwrite {a.out_tsv}")
    specs = parse_dirs(a.specimen)
    compat = json.loads(Path(a.compat_json).read_text())["code_versions"]
    _, sym_rows, _ = ed.read_matrix(Path(a.matrices_dir) / "tsaolab_melanoma_tpm_symbol.tsv.gz")
    hdr = gzip.open(Path(a.matrices_dir) / "tsaolab_melanoma_tpm_symbol.tsv.gz", "rt").readline().rstrip("\n").split("\t")
    eh, erows, _ = ed.read_matrix(Path(a.matrices_dir) / "tsaolab_melanoma_tpm_ensembl.tsv.gz")
    det = {s: 0 for s in eh[1:]}
    for l in erows.values():
        f = l.split("\t")
        for i, s in enumerate(eh[1:], 1):
            det[s] += float(f[i]) >= 1
    sv = lambda sym, sid: float(sym_rows[sym].split("\t")[hdr.index(sid)]) if sym in sym_rows else 0.0
    cols = ["sample_id", "pipeline_origin", "code_version", "slurm_job_id", "provenance_status", "instrument", "sequencing_run", "flowcell", "lane",
            "reads_raw", "reads_after_trim", "pct_adapter_trimmed", "insert_size_peak", "pct_q30_after_trim", "fastp_duplication_rate",
            "salmon_fragments", "salmon_mapped", "mapping_rate_pct", "library_type", "compatible_fragment_ratio", "frag_length_mean",
            "quality_flag", "qc_notes", "xist_tpm", "chrY_tpm", "sex_inferred", "genes_detected_tpm_ge_1", "quant_sf_sha256"]
    rows = []
    for sid, d in specs:
        prov = json.loads((d / "provenance.json").read_text())
        fp = json.loads(sorted((d / "fastp").glob("*.json"))[0].read_text()); s = fp["summary"]
        mi = json.loads((d / "quant" / "aux_info" / "meta_info.json").read_text())
        lib = json.loads((d / "quant" / "lib_format_counts.json").read_text())
        try:
            with gzip.open(prov["r1"], "rt") as fh:
                inst, run, flowcell, lane = fh.readline()[1:].split(":")[:4]
        except (OSError, KeyError, ValueError) as e:
            fail(f"{sid}: cannot read the first R1 read name from provenance r1={prov.get('r1')!r}: {e}")
        rate = mi["percent_mapped"]; xist = sv(SEX_GENES["xist"], sid); chry = sum(sv(g, sid) for g in SEX_GENES["chrY"])
        rows.append(dict(sample_id=sid, pipeline_origin="rivanna_qc_harmonized", code_version=compat.get(sid), slurm_job_id=prov.get("slurm_job_id"),
                         provenance_status=prov.get("status"), instrument=INSTRUMENTS.get(inst[:2], inst), sequencing_run=f"{inst}:{run}",
                         flowcell=flowcell, lane=f"L{int(lane):03d}", reads_raw=s["before_filtering"]["total_reads"],
                         reads_after_trim=s["after_filtering"]["total_reads"],
                         pct_adapter_trimmed=round(fp.get("adapter_cutting", {}).get("adapter_trimmed_reads", 0) / s["before_filtering"]["total_reads"] * 100, 2),
                         insert_size_peak=fp.get("insert_size", {}).get("peak"),
                         pct_q30_after_trim=round(s["after_filtering"]["q30_rate"] * 100, 2) if "q30_rate" in s["after_filtering"] else None,
                         fastp_duplication_rate=round(fp["duplication"]["rate"], 4) if "rate" in fp.get("duplication", {}) else None,
                         salmon_fragments=mi["num_processed"], salmon_mapped=mi["num_mapped"], mapping_rate_pct=round(rate, 2),
                         library_type=lib["expected_format"], compatible_fragment_ratio=round(lib["compatible_fragment_ratio"], 4),
                         frag_length_mean=round(mi.get("frag_length_mean", float("nan")), 1),
                         quality_flag="suspect" if rate < MIN_MAPPING_RATE else "ok",
                         qc_notes=f"low mapping rate ({rate:.1f}%)" if rate < MIN_MAPPING_RATE else "",
                         xist_tpm=round(xist, 2), chrY_tpm=round(chry, 2),
                         sex_inferred="male" if (chry >= 5 and xist < 5) else ("female" if (xist >= 5 and chry < 2) else "ambiguous"),
                         genes_detected_tpm_ge_1=det.get(sid), quant_sf_sha256=digest(d / "quant" / "quant.sf")))
    buf = io.StringIO(); w = csv.DictWriter(buf, fieldnames=cols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    write_new(a.out_tsv, buf.getvalue())
    print(f"harmonized QC manifest: {len(rows)} specimens -> {a.out_tsv}; suspect: {[r['sample_id'] for r in rows if r['quality_flag'] == 'suspect']}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    s = sp.add_parser("verify-compat")
    s.add_argument("--rqc-tools", required=True); s.add_argument("--specimen", action="append", required=True)
    s.add_argument("--reference-specimen"); s.add_argument("--expect-targets", type=int, default=253181); s.add_argument("--expect-libtype", default="ISR")
    for k in ("seq-hash", "name-hash", "seq-hash512", "name-hash512", "decoy-hash", "sif-sha256"):
        s.add_argument(f"--expect-{k}", required=True)
    s.add_argument("--json-out"); s.set_defaults(fn=cmd_verify)
    s = sp.add_parser("harmonized")
    s.add_argument("--expand-tool", required=True); s.add_argument("--fasta", required=True); s.add_argument("--fasta-md5")
    s.add_argument("--specimen", action="append", required=True); s.add_argument("--expect-order"); s.add_argument("--mixed-dir")
    s.add_argument("--compat-json"); s.add_argument("--out-dir", required=True); s.set_defaults(fn=cmd_harmonized)
    s = sp.add_parser("qc-manifest")
    s.add_argument("--expand-tool", required=True); s.add_argument("--specimen", action="append", required=True)
    s.add_argument("--compat-json", required=True); s.add_argument("--matrices-dir", required=True); s.add_argument("--out-tsv", required=True)
    s.set_defaults(fn=cmd_qc)
    a = ap.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()
