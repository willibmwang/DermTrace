#!/usr/bin/env python3
"""Expand the collaborator's DermTrace matrices with new Salmon quantifications (Python >= 3.9, stdlib only).

Mirrors scripts/tsaolab_melanoma/build_dermtrace.py (collaborator) exactly:
  tx2gene    GENCODE v46 FASTA header: transcript -> Ensembl gene (version stripped); gene_symbol = field 6,
             gene_type = field 8 of the FIRST transcript of the gene in FASTA order (pandas drop_duplicates)
  gene TPM / counts   sum of transcript TPM / NumReads per gene (pandas groupby-sum, Kahan-compensated)
  filter     keep genes with counts summed over ALL samples > 0
  symbol     sum of (unrounded, filtered) gene TPM / counts per gene_symbol
  layout     rows sorted by gene ID / symbol; TPM "%.4f", counts "%.3f"; tab-separated, gzip

Subcommands
  reproduce   rebuild columns from committed quant.sf files and compare them string-for-string with the committed
              matrices (proves the aggregation reproduces the collaborator's output)
  matrices    append new specimens: the committed 33 columns are copied verbatim (byte-identical text), exactly the
              new columns are appended in the given order, and genes that are zero in all committed samples but
              non-zero in a new one are ADDED as rows (zeros in the committed columns, exactly what the collaborator's
              builder would produce). No gene is discarded and no gene set is intersected.
  qc-manifest per-specimen technical QC + provenance manifest for committed (collaborator) and new (r3) specimens

Exit: 0 ok, 1 a check failed, 2 bad input. Outputs go to a NEW directory and are never overwritten.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
import platform
import socket
import sys
from datetime import datetime, timezone
from pathlib import Path

PREFIX = "tsaolab_melanoma"
MATRICES = {"tpm_ensembl": "%.4f", "counts_ensembl": "%.3f", "tpm_symbol": "%.4f", "counts_symbol": "%.3f"}
MIN_MAPPING_RATE = 30.0
SEX_GENES = dict(xist="XIST", chrY=["RPS4Y1", "DDX3Y", "KDM5D", "EIF1AY", "UTY", "USP9Y", "ZFY"])
INSTRUMENTS = {"LH": "Illumina_NovaSeqX", "A0": "Illumina_NovaSeq", "NB": "Illumina_NextSeq", "AV": "Element_AVITI"}


def fail(msg, code=2):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


def digest(p, algo="sha256"):
    h = hashlib.new(algo)
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def write_new(path: Path, data: bytes):
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    with os.fdopen(fd, "wb") as fh:
        fh.write(data)


def gz_text(text: str) -> bytes:
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb", mtime=0) as g:
        g.write(text.encode())
    return buf.getvalue()


class Kahan:
    """Compensated summation in input order (pandas groupby-sum uses Kahan summation)."""
    __slots__ = ("s", "c")

    def __init__(self):
        self.s = 0.0
        self.c = 0.0

    def add(self, x):
        y = x - self.c
        t = self.s + y
        self.c = (t - self.s) - y
        self.s = t


def tx2gene(fasta):
    t2g, ginfo = {}, {}
    with gzip.open(fasta, "rt") as fh:
        for line in fh:
            if line[0] == ">":
                f = line[1:].rstrip().split("|")
                g = f[1].split(".")[0]
                t2g[f[0]] = g
                ginfo.setdefault(g, (f[5], f[7]))  # first transcript in FASTA order
    return t2g, ginfo


def quant_file(p: Path) -> Path:
    if p.is_dir():
        for c in (p / "quant" / "quant.sf", p / "quant.sf", p / "quant.sf.gz"):
            if c.is_file():
                return c
        fail(f"no quant.sf under {p}")
    return p


def aggregate(qpath: Path, t2g):
    op = gzip.open if str(qpath).endswith(".gz") else open
    tpm, cnt = {}, {}
    with op(qpath, "rt", newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        if (next(rd, None) or [])[:5] != ["Name", "Length", "EffectiveLength", "TPM", "NumReads"]:
            fail(f"not a Salmon quant.sf: {qpath}")
        for n, r in enumerate(rd, 2):
            if len(r) != 5:
                fail(f"malformed row {n} in {qpath}")
            g = t2g.get(r[0])
            if g is None:
                fail(f"transcript {r[0]} of {qpath} not in the GENCODE FASTA")
            tpm.setdefault(g, Kahan()).add(float(r[3]))
            cnt.setdefault(g, Kahan()).add(float(r[4]))
    return {g: k.s for g, k in tpm.items()}, {g: k.s for g, k in cnt.items()}


def symbol_sum(values: dict, genes_in_order, ginfo):
    out = {}
    for g in genes_in_order:  # gene-ID-sorted order, as in pandas groupby on the filtered frame
        out.setdefault(ginfo[g][0], Kahan()).add(values.get(g, 0.0))
    return {s: k.s for s, k in out.items()}


def read_matrix(p):
    with gzip.open(p, "rt", newline="") as fh:
        lines = fh.read().split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    header = lines[0].split("\t")
    rows = {}
    order = []
    for l in lines[1:]:
        k, _, rest = l.partition("\t")
        if k in rows:
            fail(f"duplicate row {k} in {p}")
        rows[k] = l
        order.append(k)
    return header, rows, order


def committed_paths(d: Path):
    paths = {m: d / f"{PREFIX}_{m}.tsv.gz" for m in MATRICES}
    paths["genes"] = d / f"{PREFIX}_genes.tsv.gz"
    for k, p in paths.items():
        if not p.is_file():
            fail(f"committed file missing: {p}")
    return paths


def parse_specs(specs):
    out = []
    for s in specs:
        sid, _, p = s.partition("=")
        if not sid or not p:
            fail(f"expected ID=PATH, got {s!r}")
        out.append((sid, quant_file(Path(p))))
    if len({s for s, _ in out}) != len(out):
        fail("duplicate specimen IDs")
    return out


def check_fasta(a):
    if a.fasta_md5 and digest(a.fasta, "md5") != a.fasta_md5:
        fail(f"FASTA MD5 != {a.fasta_md5}")


# --------------------------------------------------------------------------------------------- reproduce
def cmd_reproduce(a):
    check_fasta(a)
    t2g, ginfo = tx2gene(a.fasta)
    paths = committed_paths(Path(a.committed_dir))
    mats = {m: read_matrix(paths[m]) for m in MATRICES}
    res = {}
    for sid, q in parse_specs(a.sample):
        tpm, cnt = aggregate(q, t2g)
        kept = mats["tpm_ensembl"][2]
        stpm, scnt = symbol_sum(tpm, kept, ginfo), symbol_sum(cnt, kept, ginfo)
        r = {}
        for m, vals in (("tpm_ensembl", tpm), ("counts_ensembl", cnt), ("tpm_symbol", stpm), ("counts_symbol", scnt)):
            header, rows, order = mats[m]
            if sid not in header:
                fail(f"{sid} is not a committed column of {m}")
            ci = header.index(sid)
            mism, maxdiff = [], 0.0
            for k in order:
                want = rows[k].split("\t")[ci]
                got = MATRICES[m] % vals.get(k, 0.0)
                if got != want:
                    mism.append((k, want, got))
                    maxdiff = max(maxdiff, abs(float(got) - float(want)))
            r[m] = {"rows": len(order), "exact_string_matches": len(order) - len(mism), "mismatches": len(mism),
                    "max_abs_diff": maxdiff, "examples": mism[:5]}
        dropped_nonzero = [g for g in set(cnt) - set(mats["counts_ensembl"][1]) if cnt[g] > 0 or tpm.get(g, 0) > 0]
        r["genes_absent_from_committed_but_nonzero_here"] = len(dropped_nonzero)
        res[sid] = r
        print(f"{sid}: " + "  ".join(f"{m} {v['exact_string_matches']}/{v['rows']}" for m, v in r.items() if isinstance(v, dict))
              + f"  absent-but-nonzero={len(dropped_nonzero)}")
    ok = all(v[m]["mismatches"] == 0 and v["genes_absent_from_committed_but_nonzero_here"] == 0 for v in res.values() for m in MATRICES)
    if a.json_out:
        write_new(Path(a.json_out), (json.dumps({"kind": "aggregation-reproduction", "verdict": "PASS" if ok else "FAIL",
                                                  "samples": res, "fasta_sha256": digest(a.fasta)}, indent=1) + "\n").encode())
    print(f"REPRODUCTION {'PASS' if ok else 'FAIL'}: every committed value rebuilt string-for-string" if ok else "REPRODUCTION FAIL")
    sys.exit(0 if ok else 1)


# ---------------------------------------------------------------------------------------------- matrices
def cmd_matrices(a):
    check_fasta(a)
    out = Path(a.out_dir)
    if os.path.lexists(out):
        fail(f"output directory exists (refusing to overwrite): {out}")
    t2g, ginfo = tx2gene(a.fasta)
    paths = committed_paths(Path(a.committed_dir))
    mats = {m: read_matrix(paths[m]) for m in MATRICES}
    base_cols = mats["tpm_ensembl"][0][1:]
    for m in MATRICES:
        if mats[m][0][1:] != base_cols:
            fail(f"committed {m} columns differ from tpm_ensembl")
    new = parse_specs(a.new)
    new_ids = [s for s, _ in new]
    if a.expect_new and new_ids != a.expect_new.split(","):
        fail(f"new specimens {new_ids} != expected {a.expect_new}")
    if set(new_ids) & set(base_cols):
        fail(f"new specimens already in the committed matrix: {sorted(set(new_ids) & set(base_cols))}")
    for k in mats["tpm_ensembl"][2]:
        if k not in ginfo:
            fail(f"committed gene {k} not in the GENCODE FASTA")
    if mats["tpm_ensembl"][2] != sorted(mats["tpm_ensembl"][2]) or mats["tpm_symbol"][2] != sorted(mats["tpm_symbol"][2]):
        fail("committed matrices are not sorted as expected; refusing to guess the row order")
    if set(mats["tpm_ensembl"][1]) != set(mats["counts_ensembl"][1]) or set(mats["tpm_symbol"][1]) != set(mats["counts_symbol"][1]):
        fail("committed TPM and counts matrices have different rows")

    agg = {sid: aggregate(q, t2g) for sid, q in new}
    committed_genes = set(mats["tpm_ensembl"][1])
    added = sorted({g for sid in new_ids for g, c in agg[sid][1].items() if c > 0} - committed_genes)
    all_genes = sorted(committed_genes | set(added))
    syms = {}
    for sid in new_ids:
        syms[sid] = (symbol_sum(agg[sid][0], all_genes, ginfo), symbol_sum(agg[sid][1], all_genes, ginfo))
    committed_syms = set(mats["tpm_symbol"][1])
    added_syms = sorted({ginfo[g][0] for g in added} - committed_syms)
    all_syms = sorted(committed_syms | set(added_syms))
    zeros = {m: "\t".join([f % 0.0] * len(base_cols)) for m, f in MATRICES.items()}

    def build(m, keys, values_for):
        header, rows, _ = mats[m]
        lines = ["\t".join(header + new_ids)]
        for k in keys:
            base = rows.get(k) or f"{k}\t{zeros[m]}"
            lines.append(base + "\t" + "\t".join(MATRICES[m] % values_for(sid, k) for sid in new_ids))
        return "\n".join(lines) + "\n"

    texts = {
        "tpm_ensembl": build("tpm_ensembl", all_genes, lambda s, g: agg[s][0].get(g, 0.0)),
        "counts_ensembl": build("counts_ensembl", all_genes, lambda s, g: agg[s][1].get(g, 0.0)),
        "tpm_symbol": build("tpm_symbol", all_syms, lambda s, k: syms[s][0].get(k, 0.0)),
        "counts_symbol": build("counts_symbol", all_syms, lambda s, k: syms[s][1].get(k, 0.0)),
    }
    # genes table: committed rows verbatim, added genes inserted in sorted position
    with gzip.open(paths["genes"], "rt", newline="") as fh:
        glines = fh.read().rstrip("\n").split("\n")
    grows = {l.split("\t", 1)[0]: l for l in glines[1:]}
    if set(grows) != committed_genes:
        fail("committed genes table rows != committed matrix rows")
    texts["genes"] = "\n".join([glines[0]] + [grows.get(g) or f"{g}\t{ginfo[g][0]}\t{ginfo[g][1]}" for g in all_genes]) + "\n"

    # ---------------- integrity checks (fail closed before anything is written)
    checks = {}
    for m in MATRICES:
        header, rows, order = mats[m]
        lines = texts[m].rstrip("\n").split("\n")
        h = lines[0].split("\t")
        n0 = len(base_cols)
        ok_cols = h[:n0 + 1] == header and h[n0 + 1:] == new_ids and len(h) == n0 + 1 + len(new_ids)
        orig_identical, added_zero, n_added = True, True, 0
        for l in lines[1:]:
            k = l.split("\t", 1)[0]
            prefix = "\t".join(l.split("\t")[:n0 + 1])
            if k in rows:
                orig_identical &= prefix == rows[k]
            else:
                n_added += 1
                added_zero &= prefix == f"{k}\t{zeros[m]}"
        keys = [l.split("\t", 1)[0] for l in lines[1:]]
        checks[m] = {"columns_committed_then_exactly_new_appended": ok_cols, "committed_rows_byte_identical": orig_identical,
                     "all_committed_rows_present": set(order) <= set(keys), "no_duplicate_rows": len(keys) == len(set(keys)),
                     "rows_sorted": keys == sorted(keys), "added_rows_zero_in_committed_columns": added_zero,
                     "rows_committed": len(order), "rows_added": n_added, "rows_total": len(keys)}
    checks["genes_table"] = {"rows": len(all_genes), "matches_ensembl_rows": [l.split("\t", 1)[0] for l in texts["genes"].split("\n")[1:-1]] == all_genes}
    colsums = {sid: {"tpm_sum": round(sum(agg[sid][0].values()), 2), "counts_sum": round(sum(agg[sid][1].values()), 1)} for sid in new_ids}
    checks["new_tpm_sums_near_1e6"] = all(abs(v["tpm_sum"] - 1e6) < 10 for v in colsums.values())
    bools = [v for c in checks.values() for v in (c.values() if isinstance(c, dict) else [c]) if isinstance(v, bool)]
    verdict = "PASS" if all(bools) else "FAIL"
    if verdict != "PASS":
        print(json.dumps(checks, indent=1), file=sys.stderr)
        fail("integrity checks failed; nothing written", 1)

    out.mkdir(parents=True)
    files = {}
    for m, t in texts.items():
        name = f"{PREFIX}_{m}.tsv.gz"
        write_new(out / name, gz_text(t))
        files[name] = digest(out / name)
    origin = {c: "collaborator_committed" for c in base_cols}
    origin.update({c: a.new_origin for c in new_ids})
    man = {"kind": "matrix-expansion", "verdict": verdict, "checks": checks, "added_genes": added, "added_symbols": added_syms,
           "new_specimens": new_ids, "column_order": base_cols + new_ids, "pipeline_origin": origin, "new_column_sums": colsums,
           "note_column_order": "new columns are appended after the committed columns (the collaborator's builder would order "
                                "columns by sample sheet); committed text is unchanged",
           "inputs": {"committed": {k: {"path": str(p), "sha256": digest(p)} for k, p in paths.items()},
                      "fasta": {"path": a.fasta, "sha256": digest(a.fasta)},
                      "new_quant": {sid: {"path": str(q), "sha256": digest(q)} for sid, q in new}},
           "outputs": files, "run": run_info()}
    write_new(out / "expansion_manifest.json", (json.dumps(man, indent=1) + "\n").encode())
    for m in MATRICES:
        c = checks[m]
        print(f"{m:<15} rows {c['rows_committed']} + {c['rows_added']} added = {c['rows_total']}; committed text identical={c['committed_rows_byte_identical']}; "
              f"columns {len(base_cols)} + {len(new_ids)} appended={c['columns_committed_then_exactly_new_appended']}")
    print(f"added genes: {len(added)}; added symbols: {len(added_syms)}; new column TPM sums: "
          + ", ".join(f"{s}={v['tpm_sum']}" for s, v in colsums.items()))
    print(f"EXPANSION {verdict} -> {out}")


def run_info():
    return {"utc": datetime.now(timezone.utc).isoformat(), "host": socket.gethostname(), "python": platform.python_version(),
            "script_sha256": digest(__file__)}


# ------------------------------------------------------------------------------------------- qc-manifest
def cmd_qc_manifest(a):
    out = Path(a.out_tsv)
    if os.path.lexists(out):
        fail(f"refusing to overwrite {out}")
    exp = json.loads(Path(a.expansion_manifest).read_text())
    _, sym_rows, _ = read_matrix(Path(a.expanded_dir) / f"{PREFIX}_tpm_symbol.tsv.gz")
    hdr = gzip.open(Path(a.expanded_dir) / f"{PREFIX}_tpm_symbol.tsv.gz", "rt").readline().rstrip("\n").split("\t")
    def symval(sym, sid):
        r = sym_rows.get(sym)
        return float(r.split("\t")[hdr.index(sid)]) if r else 0.0
    with open(a.run_qc, newline="") as fh:
        runqc = {r["sample_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    with open(a.metadata, newline="") as fh:
        meta = {r["sample_id"]: r for r in csv.DictReader(fh)}
    with open(a.sheet, newline="") as fh:
        sheet = {r["sample_id"]: r for r in csv.DictReader(fh)}
    cols = ["sample_id", "core_id", "pipeline_origin", "quantification_code", "slurm_job_id", "provenance_status",
            "instrument", "sequencing_run", "flowcell", "lane", "reads_raw", "reads_after_trim", "pct_reads_kept",
            "pct_adapter_trimmed", "mean_read_length_after_trim", "insert_size_peak", "pct_q30_after_trim",
            "fastp_duplication_rate", "fastp_version", "salmon_fragments", "salmon_mapped", "mapping_rate_pct",
            "library_type", "compatible_fragment_ratio", "frag_length_mean", "quality_flag", "qc_notes",
            "xist_tpm", "chrY_tpm", "sex_inferred", "genes_detected_tpm_ge_1", "quant_sf_sha256", "tools_sha256"]
    rows = []
    tpm_e = read_matrix(Path(a.expanded_dir) / f"{PREFIX}_tpm_ensembl.tsv.gz")
    eh = tpm_e[0]
    det = {sid: 0 for sid in eh[1:]}
    for l in tpm_e[1].values():
        f = l.split("\t")
        for i, sid in enumerate(eh[1:], 1):
            if float(f[i]) >= 1:
                det[sid] += 1
    def sex(x, y):
        return "male" if (y >= 5 and x < 5) else ("female" if (x >= 5 and y < 2) else "ambiguous")
    for sid in exp["column_order"]:
        xist = symval(SEX_GENES["xist"], sid)
        chry = sum(symval(g, sid) for g in SEX_GENES["chrY"])
        base = {"sample_id": sid, "core_id": sheet.get(sid, {}).get("core_id"), "pipeline_origin": exp["pipeline_origin"][sid],
                "xist_tpm": round(xist, 2), "chrY_tpm": round(chry, 2), "sex_inferred": sex(xist, chry), "genes_detected_tpm_ge_1": det[sid]}
        if exp["pipeline_origin"][sid] == "collaborator_committed":
            q, m = runqc[sid], meta[sid]
            base.update({k: q.get(k) for k in cols if k in q})
            base.update(quality_flag=m.get("quality_flag"), qc_notes=m.get("qc_notes"), quantification_code="collaborator run_quant.sh",
                        slurm_job_id="", provenance_status="committed outputs (no provenance schema)", fastp_version="unknown",
                        quant_sf_sha256="", tools_sha256="")
        else:
            d = Path(a.quant_root) / sid
            prov = json.loads((d / "provenance.json").read_text())
            fpj = sorted((d / "fastp").glob("*.json"))
            fp = json.loads(fpj[0].read_text()); s = fp["summary"]
            mi = json.loads((d / "quant" / "aux_info" / "meta_info.json").read_text())
            lib = json.loads((d / "quant" / "lib_format_counts.json").read_text())
            r1 = Path(prov["r1"])
            with gzip.open(r1, "rt") as fh:
                inst, run, flowcell, lane = fh.readline()[1:].split(":")[:4]
            rate = mi["percent_mapped"]
            notes = [f"low mapping rate ({rate:.1f}%)"] if rate < MIN_MAPPING_RATE else []
            cfg = prov.get("configuration", {})
            base.update(quantification_code=f"rivanna_qc r3 ({prov.get('script')})", slurm_job_id=prov.get("slurm_job_id"),
                        provenance_status=prov.get("status"), instrument=INSTRUMENTS.get(inst[:2], inst), sequencing_run=f"{inst}:{run}",
                        flowcell=flowcell, lane=f"L{int(lane):03d}", reads_raw=s["before_filtering"]["total_reads"],
                        reads_after_trim=s["after_filtering"]["total_reads"],
                        pct_reads_kept=round(s["after_filtering"]["total_reads"] / s["before_filtering"]["total_reads"] * 100, 2),
                        pct_adapter_trimmed=round(fp.get("adapter_cutting", {}).get("adapter_trimmed_reads", 0) / s["before_filtering"]["total_reads"] * 100, 2),
                        mean_read_length_after_trim=(round((s["after_filtering"]["read1_mean_length"] + s["after_filtering"]["read2_mean_length"]) / 2, 1)
                                                     if "read1_mean_length" in s["after_filtering"] and "read2_mean_length" in s["after_filtering"] else None),
                        insert_size_peak=fp.get("insert_size", {}).get("peak"), pct_q30_after_trim=round(s["after_filtering"]["q30_rate"] * 100, 2) if "q30_rate" in s["after_filtering"] else None,
                        fastp_duplication_rate=round(fp["duplication"]["rate"], 4) if "rate" in fp.get("duplication", {}) else None,
                        fastp_version=cfg.get("fastp_version", {}).get("value"), salmon_fragments=mi["num_processed"],
                        salmon_mapped=mi["num_mapped"], mapping_rate_pct=round(rate, 2), library_type=lib["expected_format"],
                        compatible_fragment_ratio=round(lib["compatible_fragment_ratio"], 4),
                        frag_length_mean=round(mi.get("frag_length_mean", float("nan")), 1),
                        quality_flag="suspect" if notes else "ok", qc_notes="; ".join(notes),
                        quant_sf_sha256=digest(d / "quant" / "quant.sf"), tools_sha256=prov.get("tools_sha256"))
        rows.append(base)
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
    w.writeheader(); w.writerows(rows)
    write_new(out, buf.getvalue().encode())
    n_new = sum(1 for r in rows if r["pipeline_origin"] != "collaborator_committed")
    print(f"QC manifest: {len(rows)} specimens ({len(rows) - n_new} collaborator, {n_new} new) -> {out}")
    for r in rows:
        if r["pipeline_origin"] != "collaborator_committed":
            print(f"  {r['sample_id']} job={r['slurm_job_id']} status={r['provenance_status']} run={r['sequencing_run']} lane={r['lane']} "
                  f"raw={r['reads_raw']} map={r['mapping_rate_pct']}% dup={r['fastp_duplication_rate']} flag={r['quality_flag']} sex={r['sex_inferred']}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    s = sp.add_parser("reproduce")
    s.add_argument("--committed-dir", required=True); s.add_argument("--fasta", required=True); s.add_argument("--fasta-md5")
    s.add_argument("--sample", action="append", required=True, help="ID=quant.sf(.gz) or quant dir (committed specimen)")
    s.add_argument("--json-out"); s.set_defaults(fn=cmd_reproduce)
    s = sp.add_parser("matrices")
    s.add_argument("--committed-dir", required=True); s.add_argument("--fasta", required=True); s.add_argument("--fasta-md5")
    s.add_argument("--new", action="append", required=True, help="ID=specimen output dir or quant.sf (append order)")
    s.add_argument("--expect-new", help="comma-separated expected new IDs in order")
    s.add_argument("--new-origin", default="rivanna_qc_r3")
    s.add_argument("--out-dir", required=True); s.set_defaults(fn=cmd_matrices)
    s = sp.add_parser("qc-manifest")
    for k in ("expansion-manifest", "expanded-dir", "run-qc", "metadata", "sheet", "quant-root", "out-tsv"):
        s.add_argument(f"--{k}", required=True)
    s.set_defaults(fn=cmd_qc_manifest)
    a = ap.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()
