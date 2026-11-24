#!/usr/bin/env python3
"""Standard-library helpers for the Greek melanoma Rivanna quantification scripts (Python >= 3.9).

Subcommands
  provenance          write a JSON provenance record (refuses to overwrite)
  resolve-specimen    resolve one specimen: sample sheet + validated_pairs.tsv eligibility + exclusion list
  verify-fingerprint  compare Salmon index fingerprints (index info.json or quant meta_info.json) to expected values
  validate-quant-dir  structural integrity check of one specimen output directory (staged or adopted)
  compare-quant       compare a reproduced specimen against the committed collaborator output (rule provisional-v1)
  audit-historical    read-only concordance audit of a quantification run that predates the provenance schema
                      (same provisional-v1 metrics; never adopts, never writes into the audited directory)

Exit codes: 0 success / PASS, 1 mismatch / FAIL, 2 required evidence unavailable or bad input.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import platform
import re
import socket
import sys
from datetime import datetime, timezone
from pathlib import Path

HEX32 = re.compile(r"^[0-9a-f]{32}$")
SAFE_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")   # no leading dot, no slash, no traversal
CORE_RE = re.compile(r"^20088XR-(\d{2}-\d{2})$")

# Field name in index info.json -> field name in quant aux_info/meta_info.json -> canonical name
FINGERPRINTS = [
    ("SeqHash", "index_seq_hash", "seq_hash"),
    ("NameHash", "index_name_hash", "name_hash"),
    ("SeqHash512", "index_seq_hash512", "seq_hash512"),
    ("NameHash512", "index_name_hash512", "name_hash512"),
    ("DecoySeqHash", "index_decoy_seq_hash", "decoy_seq_hash"),
    ("DecoyNameHash", "index_decoy_name_hash", "decoy_name_hash"),
]

# Built-in exclusions. They can be lifted only with an explicit, recorded adjudication reference.
BUILTIN_EXCLUSIONS = {
    "20088XR-11-37": "duplicate Drive objects; provenance adjudication pending",
    "20088XR-11-43": "corrupt FASTQ (gzip validation failed)",
    "20088XR-11-44": "corrupt FASTQ (gzip validation failed)",
    "20088XR-11-46": "corrupt FASTQ (gzip validation failed)",
    "20088XR-11-47": "corrupt FASTQ (gzip validation failed)",
    "20088XR-11-48": "corrupt FASTQ (gzip validation failed)",
}

# Provisional acceptance rule. Thresholds are frozen; any other values must carry a new rule version + justification.
RULE_V1 = {"version": "provisional-v1", "min_pearson_log1p_tpm_transcript": 0.99, "max_abs_mapping_diff_pp": 1.0}


def fail(msg: str, code: int = 2):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


def load_json(p) -> dict:
    try:
        return json.loads(Path(p).read_text())
    except (OSError, ValueError) as e:
        fail(f"cannot read JSON {p}: {e}")


def write_new(path: Path, text: str):
    """Write text to a path that must not exist (O_EXCL), so concurrent writers cannot clobber each other."""
    try:
        fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    except FileExistsError:
        fail(f"refusing to overwrite: {path}")
    with os.fdopen(fd, "w") as fh:
        fh.write(text)


# ------------------------------------------------------------------------------------------- provenance
def cmd_provenance(a):
    rec = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "hostname": socket.gethostname(),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "cwd": os.getcwd(),
    }
    for kv in a.set or []:
        if "=" not in kv:
            fail(f"--set expects key=value, got {kv!r}")
        k, v = kv.split("=", 1)
        rec[k] = v
    for kv in a.json or []:
        k, path = kv.split("=", 1)
        rec[k] = load_json(path)
    for kv in a.config or []:  # name=value::status::source
        try:
            name, rest = kv.split("=", 1)
            value, status, source = rest.split("::", 2)
        except ValueError:
            fail(f"--config expects name=value::status::source, got {kv!r}")
        if status not in ("verified-collaborator", "candidate", "inferred", "observed-this-run", "user-reported"):
            fail(f"unknown configuration status {status!r}")
        rec.setdefault("configuration", {})[name] = {"value": value, "status": status, "source": source}
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    write_new(out, json.dumps(rec, indent=2, sort_keys=True) + "\n")


# ------------------------------------------------------------------------------------- resolve-specimen
def specimen_ids(row):
    ids = {row["sample_id"], row["core_id"]}
    m = CORE_RE.match(row["core_id"])
    if m:
        ids.add(m.group(1))
    return ids


# validated_pairs.tsv rules. Keep identical to run_fastqc_all.sbatch (tests/test_legacy_qc_scripts.py checks parity).
VP_HEADER_FIRST = {"sample", "sample_id", "specimen", "label"}
VP_FASTQ_LIKE = re.compile(r"\.(fastq|fq)(\.gz)?$", re.I)
VP_R1_NAME = re.compile(r"^[^/]*_R1_\d{3}\.fastq\.gz$")
VP_FILE_CORE = re.compile(r"^(20088XR-\d{2}-\d{2})_")
VP_CORE_LABEL = re.compile(r"^(20088XR-)?\d{2}-\d{2}$")


def vp_is_header(fields):
    """Line 1 is a header only if it has 3 fields, field 1 is a known column name, and fields 2-3 are not paths."""
    return (len(fields) == 3 and fields[0].lower() in VP_HEADER_FIRST
            and not any("/" in x or VP_FASTQ_LIKE.search(x) for x in fields[1:]))


def read_validated_pairs(path):
    """Strictly parse validated_pairs.tsv; any defect fails the whole manifest (exit 2), nothing is skipped.

    Historical Rivanna file: a header line (field 1 'sample') + 43 rows. Headerless files are also accepted.
    Rows are SAMPLE<TAB>R1<TAB>R2: exactly 3 fields, no surrounding whitespace, LF line endings, no blank lines;
    R1/R2 are absolute paths in the same directory, R1 named *_R1_NNN.fastq.gz and R2 the same name with _R1_ ->
    _R2_; sample labels, files and filename specimens (20088XR-NN-NN) are unique; a label that is a core ID
    (20088XR-NN-NN or NN-NN) must name the same specimen as its files. Paths are returned verbatim; callers
    match on basenames and compare absolute paths with RAW_DIR.
    """
    data = Path(path).read_bytes()
    if b"\r" in data:
        fail(f"malformed line endings in {path}: CR (Windows) characters are not accepted")
    try:
        lines = data.decode("utf-8").split("\n")
    except UnicodeDecodeError:
        fail(f"malformed encoding in {path}: not UTF-8")
    if lines and lines[-1] == "":
        lines.pop()
    errs, rows, seen_label, seen_file, seen_core = [], [], {}, {}, {}
    for n, line in enumerate(lines, 1):
        f = line.split("\t")
        if n == 1 and vp_is_header(f):
            continue
        if line == "":
            errs.append(f"malformed line {n}: blank line")
            continue
        if len(f) != 3:
            errs.append(f"malformed line {n}: expected 3 tab-separated fields, found {len(f)}")
            continue
        if any(x == "" or x != x.strip() for x in f):
            errs.append(f"malformed line {n}: empty field or surrounding whitespace")
            continue
        label, r1, r2 = f
        for name, p in (("R1", r1), ("R2", r2)):
            if not os.path.isabs(p):
                errs.append(f"malformed line {n}: {name} path is not absolute: {p}")
        b1, b2 = os.path.basename(r1), os.path.basename(r2)
        if not VP_R1_NAME.match(b1):
            errs.append(f"malformed line {n}: R1 name does not end in _R1_NNN.fastq.gz: {b1}")
        if b1.count("_R1_") != 1 or b1.replace("_R1_", "_R2_") != b2 or os.path.dirname(r1) != os.path.dirname(r2):
            errs.append(f"malformed line {n}: R2 is not the mate of R1: {r1} | {r2}")
        for key, seen, what in ((label, seen_label, "duplicate sample label"), (r1, seen_file, "file listed twice"),
                                (r2, seen_file, "file listed twice")):
            if key in seen:
                errs.append(f"malformed line {n}: {what}: {key} (first on line {seen[key]})")
            else:
                seen[key] = n
        m = VP_FILE_CORE.match(b1)
        if m:
            core = m.group(1)
            if core in seen_core:
                errs.append(f"malformed line {n}: duplicate specimen {core} (first on line {seen_core[core]})")
            else:
                seen_core[core] = n
            if VP_CORE_LABEL.match(label) and label not in (core, core[len("20088XR-"):]):
                errs.append(f"malformed line {n}: label {label!r} is inconsistent with file specimen {core}")
        rows.append((label, r1, r2))
    if errs:
        fail(f"invalid validated_pairs {path}:\n  " + "\n  ".join(errs[:20]))
    if not rows:
        fail(f"no pairs in {path}")
    return rows


def cmd_resolve(a):
    with open(a.sheet, newline="") as fh:
        rows = list(csv.DictReader(fh))
    need = {"sample_id", "core_id", "r1", "r2", "md5_r1", "md5_r2"}
    if not rows or not need <= set(rows[0]):
        fail(f"sample sheet lacks columns {sorted(need)}: {a.sheet}")
    hits = [r for r in rows if a.specimen in specimen_ids(r)]
    if len(hits) != 1:
        fail(f"specimen {a.specimen!r} matched {len(hits)} rows in {a.sheet} (need exactly 1)")
    r = hits[0]
    for k in ("sample_id", "r1", "r2"):
        if not SAFE_NAME.match(r[k]):
            fail(f"unsafe {k} in sheet: {r[k]!r}")
    for k in ("md5_r1", "md5_r2"):
        if not HEX32.match(r[k].lower()):
            fail(f"malformed MD5 in sheet: {r[k]!r}")

    # exclusion list: built-in + caller-supplied; lifting a built-in exclusion needs a recorded reference
    lifted = {}
    for kv in a.lift_exclusion or []:
        k, _, ref = kv.partition("=")
        if not ref.strip():
            fail(f"--lift-exclusion needs CORE_ID=<adjudication reference>, got {kv!r}")
        lifted[k] = ref.strip()
    excluded = {k: v for k, v in BUILTIN_EXCLUSIONS.items() if k not in lifted}
    for x in (a.exclude or "").replace(",", " ").split():
        excluded.setdefault(x, "excluded by EXCLUDED_SPECIMENS")
    for x, why in excluded.items():
        if x in specimen_ids(r):
            fail(f"specimen {r['sample_id']} ({r['core_id']}) is excluded: {why}", 1)

    # eligibility: exactly one validated pair with the same R1/R2 filenames and a consistent sample label
    vp = read_validated_pairs(a.validated_pairs)
    match = [v for v in vp if os.path.basename(v[1]) == r["r1"] and os.path.basename(v[2]) == r["r2"]]
    if len(match) != 1:
        fail(f"{r['sample_id']} ({r['r1']}, {r['r2']}) has {len(match)} matching pairs in {a.validated_pairs}; "
             "only gzip-validated pairs are eligible", 1)
    label, v1, v2 = match[0]
    ids = specimen_ids(r)
    if not (label in ids or any(i in label for i in ids if len(i) >= 5)):
        fail(f"validated_pairs label {label!r} is inconsistent with {sorted(ids)}", 1)
    if a.raw_dir:
        for vpath, name in ((v1, r["r1"]), (v2, r["r2"])):
            raw = Path(a.raw_dir) / name
            if os.path.isabs(vpath) and os.path.exists(vpath) and raw.exists() and os.path.realpath(vpath) != os.path.realpath(raw):
                fail(f"validated path {vpath} and RAW_DIR file {raw} are different files", 1)
    print("\t".join([r["sample_id"], r["core_id"], r["r1"], r["r2"], r["md5_r1"].lower(), r["md5_r2"].lower(), label,
                     ";".join(f"{k}={v}" for k, v in lifted.items()) or "none"]))


# ----------------------------------------------------------------------------------- verify-fingerprint
def load_fingerprints(path: Path):
    """Return ({canonical: value}, source_kind, extra) from an index dir, info.json or meta_info.json."""
    if path.is_dir():
        info = path / "info.json"
        if not info.is_file():
            fail(f"index directory has no info.json: {path}")
        data = load_json(info)
        extra = {"k": data.get("k")}
        vi = path / "versionInfo.json"
        if vi.is_file():
            extra["versionInfo"] = load_json(vi)
        return {c: data.get(i) for i, _, c in FINGERPRINTS}, "index_info", extra
    data = load_json(path)
    if "index_seq_hash" in data:
        return ({c: data.get(m) for _, m, c in FINGERPRINTS}, "quant_meta_info",
                {"salmon_version": data.get("salmon_version"), "num_valid_targets": data.get("num_valid_targets")})
    return {c: data.get(i) for i, _, c in FINGERPRINTS}, "index_info", {"k": data.get("k")}


def cmd_fingerprint(a):
    observed, kind, extra = load_fingerprints(Path(a.observed))
    expected = {c: None for _, _, c in FINGERPRINTS}
    if a.reference:
        expected.update(load_fingerprints(Path(a.reference))[0])
    for _, _, c in FINGERPRINTS:
        v = getattr(a, f"expect_{c}")
        if v:
            if expected[c] and expected[c] != v:
                fail(f"--expect-{c} disagrees with --reference for {c}")
            expected[c] = v
    rows, status = [], {}
    for _, _, c in FINGERPRINTS:
        o, e = observed.get(c), expected.get(c)
        s = "UNAVAILABLE" if not o or not e else ("MATCH" if o.lower() == e.lower() else "MISMATCH")
        status[c] = s
        rows.append({"fingerprint": c, "status": s, "observed": o, "expected": e})
    required = [c.strip() for c in a.require.split(",") if c.strip()]
    unknown = set(required) - set(status)
    if unknown:
        fail(f"unknown fingerprint(s) in --require: {sorted(unknown)}")
    k_note = {"observed_k": extra.get("k"), "collaborator_k": "UNVERIFIED (not recorded in committed outputs)"}
    if a.expect_k is not None and extra.get("k") is not None:
        k_note["requested_k_status"] = ("MATCH" if int(extra["k"]) == a.expect_k else "MISMATCH") + f" (requested {a.expect_k})"
    k_bad = k_note.get("requested_k_status", "").startswith("MISMATCH")
    verdict = ("MISMATCH" if any(s == "MISMATCH" for s in status.values()) or k_bad
               else "UNAVAILABLE" if any(status[c] == "UNAVAILABLE" for c in required) else "MATCH")
    report = {"observed_source": str(a.observed), "observed_kind": kind, "reference": a.reference,
              "fingerprints": rows, "required": required, "kmer": k_note, "extra": extra, "verdict": verdict}
    for r in rows:
        print(f"{r['fingerprint']:<16} {r['status']}")
    print(f"{'k-mer':<16} observed={k_note['observed_k']} {k_note.get('requested_k_status', '')} (collaborator k UNVERIFIED)")
    print(f"VERDICT {verdict}")
    if a.json_out:
        write_new(Path(a.json_out), json.dumps(report, indent=2) + "\n")
    sys.exit({"MATCH": 0, "MISMATCH": 1, "UNAVAILABLE": 2}[verdict])


# ----------------------------------------------------------------------------------- validate-quant-dir
def read_quant(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(path, "rt", newline="") as fh:
            rd = csv.reader(fh, delimiter="\t")
            header = next(rd, None)
            if header is None or header[:5] != ["Name", "Length", "EffectiveLength", "TPM", "NumReads"]:
                fail(f"not a Salmon quant.sf: {path}")
            out = {}
            for n, row in enumerate(rd, 2):
                if len(row) != 5:
                    fail(f"truncated or malformed row {n} in {path}")
                if row[0] in out:
                    fail(f"duplicate transcript {row[0]} in {path}")
                out[row[0]] = (float(row[3]), float(row[4]))
    except (OSError, EOFError, ValueError) as e:
        fail(f"cannot read quant file {path}: {e}")
    if not out:
        fail(f"empty quant.sf: {path}")
    return out


def check_quant_dir(d: Path, stage: bool, specimen=None, expect_targets=None, expect_libtype=None):
    """Return (problems, facts). A non-staged directory must be an adopted, complete specimen output."""
    problems, facts = [], {}
    real = os.path.realpath(d)
    if not stage and ".partial" in os.path.basename(real):
        problems.append(f"path is a partial work directory: {real}")
    if not stage and d.is_symlink():
        problems.append("specimen directory is a symlink")
    q, meta, lib = d / "quant" / "quant.sf", d / "quant" / "aux_info" / "meta_info.json", d / "quant" / "lib_format_counts.json"
    for p in (q, meta, lib):
        if not p.is_file() or p.stat().st_size == 0:
            problems.append(f"missing or empty {p.relative_to(d)}")
    fastp = sorted((d / "fastp").glob("*.json")) if (d / "fastp").is_dir() else []
    if len(fastp) != 1:
        problems.append(f"expected exactly one fastp JSON, found {len(fastp)}")
    if problems:
        return problems, facts
    quant = read_quant(q)
    m, lf, fp = load_json(meta), load_json(lib), load_json(fastp[0])
    facts.update(n_rows=len(quant), num_valid_targets=m.get("num_valid_targets"), expected_format=lf.get("expected_format"),
                 percent_mapped=m.get("percent_mapped"), num_processed=m.get("num_processed"), num_mapped=m.get("num_mapped"),
                 compatible_fragment_ratio=lf.get("compatible_fragment_ratio"),
                 raw_reads=fp.get("summary", {}).get("before_filtering", {}).get("total_reads"),
                 reads_after_trim=fp.get("summary", {}).get("after_filtering", {}).get("total_reads"),
                 fastp_version_reported=fp.get("summary", {}).get("fastp_version"))
    if m.get("num_valid_targets") != len(quant):
        problems.append(f"quant.sf has {len(quant)} rows but meta_info num_valid_targets={m.get('num_valid_targets')}")
    if expect_targets is not None and len(quant) != expect_targets:
        problems.append(f"quant.sf has {len(quant)} rows, expected {expect_targets}")
    if expect_libtype and lf.get("expected_format") != expect_libtype:
        problems.append(f"library format {lf.get('expected_format')!r}, expected {expect_libtype!r}")
    if facts["raw_reads"] is None or facts["reads_after_trim"] is None:
        problems.append("fastp JSON lacks before/after filtering read counts")
    elif m.get("num_processed") is not None and 2 * m["num_processed"] != facts["reads_after_trim"]:
        problems.append(f"Salmon processed {m['num_processed']} fragments but fastp kept {facts['reads_after_trim']} reads")
    if not stage:
        prov = d / "provenance.json"
        if not prov.is_file():
            problems.append("missing provenance.json")
        else:
            pr = load_json(prov)
            facts["provenance_specimen"] = pr.get("specimen")
            if pr.get("status") != "complete":
                problems.append(f"provenance status is {pr.get('status')!r}, not 'complete'")
            if specimen and pr.get("specimen") != specimen and pr.get("core_id") != specimen:
                problems.append(f"provenance specimen {pr.get('specimen')!r}/{pr.get('core_id')!r} != {specimen!r}")
        if any(d.glob(".*.partial*")):
            problems.append("nested partial directory inside specimen output")
    return problems, facts


def cmd_validate(a):
    problems, facts = check_quant_dir(Path(a.dir), a.stage, a.specimen, a.expect_targets, a.expect_libtype)
    for k, v in facts.items():
        print(f"{k:<26} {v}")
    for p in problems:
        print(f"PROBLEM {p}")
    print("VALID" if not problems else "INVALID")
    sys.exit(0 if not problems else 1)


# ---------------------------------------------------------------------------------------- compare-quant
def pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    return sxy / math.sqrt(sxx * syy) if sxx > 0 and syy > 0 else float("nan")


def ranks(v):
    order = sorted(range(len(v)), key=v.__getitem__)
    r = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return r


def gene_map(fasta: Path):
    out = {}
    with (gzip.open(fasta, "rt") if fasta.suffix == ".gz" else open(fasta)) as fh:
        for line in fh:
            if line.startswith(">"):
                f = line[1:].rstrip().split("|")
                if len(f) < 2:
                    fail(f"not a GENCODE-style header: {line.strip()[:80]}")
                out[f[0]] = f[1].split(".")[0]
    return out


def level_stats(t: dict, r: dict, min_tpm: float):
    keys = sorted(set(t) | set(r))
    tv, rv = [t.get(k, 0.0) for k in keys], [r.get(k, 0.0) for k in keys]
    expressed = [i for i, k in enumerate(keys) if tv[i] >= min_tpm or rv[i] >= min_tpm]
    lfc = [math.log2((tv[i] + 1) / (rv[i] + 1)) for i in expressed]
    big = sum(1 for x in lfc if abs(x) > 1)
    return {"n_features": len(keys), "pearson_log1p_tpm": round(pearson([math.log1p(x) for x in tv], [math.log1p(x) for x in rv]), 6),
            "spearman_tpm": round(pearson(ranks(tv), ranks(rv)), 6), f"n_expressed_tpm_ge_{min_tpm:g}": len(expressed),
            "n_expressed_abs_log2fc_gt_1": big, "frac_expressed_abs_log2fc_gt_1": round(big / len(expressed), 6) if expressed else None,
            "max_abs_log2fc_expressed": round(max((abs(x) for x in lfc), default=0.0), 4)}


def concordance(test, tmeta, tlib, t_raw, t_after, ref, rmeta, rlib, rq, gm, rule, expect_targets, expect_libtype, min_tpm):
    """Identity/integrity gates, transcript- and gene-level statistics and mapping yield shared by compare-quant and
    audit-historical, so both apply exactly the same provisional-v1 computation. Returns (gates, tx, gene, yield)."""
    gates = {}
    gates["raw_input_reads_identical"] = t_raw is not None and str(t_raw) == str(rq.get("reads_raw"))
    gates["identical_transcript_ids"] = bool(test) and set(test) == set(ref)
    gates["target_counts_equal"] = (tmeta.get("num_valid_targets") == rmeta.get("num_valid_targets") == len(ref) == len(test)
                                    and (expect_targets is None or len(ref) == expect_targets))
    gates["library_format_matches"] = (tlib.get("expected_format") == rlib.get("expected_format") == expect_libtype)
    for _, key, c in FINGERPRINTS:
        gates[f"index_{c}_identical"] = bool(tmeta.get(key)) and tmeta.get(key) == rmeta.get(key)

    tx = level_stats({k: v[0] for k, v in test.items()}, {k: v[0] for k, v in ref.items()}, min_tpm) if test else {}
    missing = [t for t in set(test) | set(ref) if t not in gm]
    if missing:
        fail(f"{len(missing)} transcripts absent from gene map (e.g. {missing[0]})")

    def agg(q):
        out = {}
        for t, (tpm, _) in q.items():
            out[gm[t]] = out.get(gm[t], 0.0) + tpm
        return out
    gene = level_stats(agg(test), agg(ref), min_tpm) if test else {}

    yield_ = {
        "raw_reads_test": t_raw, "raw_reads_reference": int(rq["reads_raw"]),
        "reads_after_trim_test": t_after, "reads_after_trim_reference": int(rq["reads_after_trim"]),
        "num_processed_test": tmeta.get("num_processed"), "num_processed_reference": rmeta.get("num_processed"),
        "num_mapped_test": tmeta.get("num_mapped"), "num_mapped_reference": rmeta.get("num_mapped"),
        "percent_mapped_test": tmeta.get("percent_mapped"), "percent_mapped_reference": rmeta.get("percent_mapped"),
        "compatible_fragment_ratio_test": tlib.get("compatible_fragment_ratio"),
        "compatible_fragment_ratio_reference": rlib.get("compatible_fragment_ratio"),
    }
    for side in ("test", "reference"):
        raw, mapped = yield_[f"raw_reads_{side}"], yield_[f"num_mapped_{side}"]
        yield_[f"mapped_fragments_per_raw_pair_{side}"] = round(mapped / (raw / 2), 6) if raw and mapped is not None else None
    if tmeta.get("percent_mapped") is not None and rmeta.get("percent_mapped") is not None:
        yield_["percent_mapped_diff_pp"] = round(tmeta["percent_mapped"] - rmeta["percent_mapped"], 4)
    if t_after is not None and yield_["reads_after_trim_reference"]:
        yield_["reads_after_trim_diff"] = t_after - yield_["reads_after_trim_reference"]

    gates[f"transcript_pearson_log1p_tpm_ge_{rule['min_pearson_log1p_tpm_transcript']:g}"] = bool(tx) and tx["pearson_log1p_tpm"] >= rule["min_pearson_log1p_tpm_transcript"]
    gates[f"abs_mapping_diff_le_{rule['max_abs_mapping_diff_pp']:g}pp"] = ("percent_mapped_diff_pp" in yield_
                                                                          and abs(yield_["percent_mapped_diff_pp"]) <= rule["max_abs_mapping_diff_pp"])
    return gates, tx, gene, yield_


def rule_from_args(a):
    custom = a.min_pearson != RULE_V1["min_pearson_log1p_tpm_transcript"] or a.max_mapping_diff_pp != RULE_V1["max_abs_mapping_diff_pp"]
    if custom and (a.rule_version in (None, RULE_V1["version"]) or not a.rule_justification):
        fail("non-default thresholds require --rule-version (not provisional-v1) and --rule-justification")
    return {"version": a.rule_version or RULE_V1["version"], "min_pearson_log1p_tpm_transcript": a.min_pearson,
            "max_abs_mapping_diff_pp": a.max_mapping_diff_pp, "justification": a.rule_justification or "frozen provisional rule",
            "status": "PROVISIONAL engineering thresholds; not a validated equivalence margin"}


def read_run_qc_row(path, specimen):
    with open(path, newline="") as fh:
        rq = [r for r in csv.DictReader(fh, delimiter="\t") if specimen in (r.get("sample_id"), r.get("core_id"))]
    if len(rq) != 1:
        fail(f"specimen {specimen!r} matched {len(rq)} rows in {path}")
    return rq[0]


def cmd_compare(a):
    rule = rule_from_args(a)

    td = Path(a.test_dir)
    problems, tf = check_quant_dir(td, stage=False, specimen=a.specimen, expect_targets=a.expect_targets, expect_libtype=a.expect_libtype)
    test = read_quant(td / "quant" / "quant.sf") if (td / "quant" / "quant.sf").is_file() else {}
    tmeta = load_json(td / "quant" / "aux_info" / "meta_info.json") if (td / "quant" / "aux_info" / "meta_info.json").is_file() else {}
    tlib = load_json(td / "quant" / "lib_format_counts.json") if (td / "quant" / "lib_format_counts.json").is_file() else {}
    ref, rmeta, rlib = read_quant(Path(a.reference_quant)), load_json(a.reference_meta), load_json(a.reference_lib_format)
    rq = read_run_qc_row(a.reference_run_qc, a.specimen)

    gates = {"test_output_complete_and_valid": not problems}
    gates["specimen_identity"] = tf.get("provenance_specimen") in (rq.get("sample_id"), rq.get("core_id")) and a.specimen in (rq.get("sample_id"), rq.get("core_id"))
    shared, tx, gene, yield_ = concordance(test, tmeta, tlib, tf.get("raw_reads"), tf.get("reads_after_trim"), ref, rmeta, rlib, rq,
                                           gene_map(Path(a.gene_map_fasta)), rule, a.expect_targets, a.expect_libtype, a.min_tpm)
    gates.update(shared)
    verdict = "PASS" if all(gates.values()) else "FAIL"
    report = {"specimen": a.specimen, "test_dir": str(td), "reference_quant": a.reference_quant, "rule": rule,
              "test_problems": problems, "gates": gates, "transcript_level": tx, "gene_level_reported_not_gated": gene,
              "yield": yield_, "test_fastp_version_reported": tf.get("fastp_version_reported"), "verdict": verdict}
    print(f"rule {rule['version']} ({rule['status']})")
    for p in problems:
        print(f"PROBLEM {p}")
    for section, d in (("transcript", tx), ("gene", gene), ("yield", yield_)):
        for k, v in d.items():
            print(f"{section:<10} {k:<40} {v}")
    for k, v in gates.items():
        print(f"gate {k:<48} {'PASS' if v else 'FAIL'}")
    print(f"VERDICT {verdict} ({rule['version']})")
    if a.json_out:
        write_new(Path(a.json_out), json.dumps(report, indent=2) + "\n")
    sys.exit(0 if verdict == "PASS" else 1)


# ------------------------------------------------------------------------------------- audit-historical
# A run made before the provenance schema existed cannot pass validate-quant-dir (no provenance.json, ad hoc layout),
# and must not be made to. This audit reads it in place, never writes into it, labels every check PASS / FAIL /
# UNAVAILABLE, applies the same provisional-v1 computation as compare-quant, and never reports an adopted reproduction.
HIST_KINDS = {"quant": ("quant.sf", "quant.sf.gz"), "meta": ("meta_info.json",),
              "lib_format": ("lib_format_counts.json",), "cmd_info": ("cmd_info.json",)}
CORE_ANY = re.compile(r"20088XR-\d{2}-\d{2}")


def file_digest(p: Path, algo: str) -> str:
    h = hashlib.new(algo)
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def is_fastp_json(p: Path) -> bool:
    try:
        if p.stat().st_size > 200 * 1024 * 1024:
            return False
        d = json.loads(p.read_text())
    except (OSError, ValueError, UnicodeDecodeError):
        return False
    return isinstance(d, dict) and "before_filtering" in d.get("summary", {})


def discover_historical(root: Path):
    found = {k: [] for k in list(HIST_KINDS) + ["fastp", "provenance"]}
    known = {n for names in HIST_KINDS.values() for n in names}
    for dirpath, dirnames, filenames in os.walk(root, followlinks=False):
        dirnames.sort()
        for f in sorted(filenames):
            fp = Path(dirpath) / f
            for kind, names in HIST_KINDS.items():
                if f in names:
                    found[kind].append(fp)
            if f == "provenance.json":
                found["provenance"].append(fp)
            elif f.endswith(".json") and f not in known and is_fastp_json(fp):
                found["fastp"].append(fp)
    return found


def tri(ok):
    return "UNAVAILABLE" if ok is None else ("PASS" if ok else "FAIL")


def path_inside(child: Path, parent: Path) -> bool:
    c, p = os.path.realpath(child), os.path.realpath(parent)
    return c == p or c.startswith(p.rstrip(os.sep) + os.sep)


def cmd_audit(a):
    rule = rule_from_args(a)
    root = Path(a.test_dir)
    if not root.is_dir():
        fail(f"historical run directory not found: {root}")
    if a.json_out and (path_inside(Path(a.json_out).parent, root)):
        fail("--json-out must be outside the audited directory (the audit never writes into it)")

    found = discover_historical(root)
    overrides = {"quant": a.test_quant, "meta": a.test_meta, "lib_format": a.test_lib_format,
                 "fastp": a.test_fastp_json, "cmd_info": a.test_cmd_info}
    chosen, notes = {}, []
    for kind, ov in overrides.items():
        if ov:
            if not Path(ov).is_file():
                fail(f"--test-{kind.replace('_', '-')} is not a file: {ov}")
            chosen[kind] = Path(ov)
            notes.append(f"{kind}: explicit override {ov}")
        elif len(found[kind]) == 1:
            chosen[kind] = found[kind][0]
        elif len(found[kind]) > 1:
            fail(f"ambiguous: {len(found[kind])} candidate {kind} files under {root}: "
                 + ", ".join(str(x) for x in found[kind]) + f"; choose one with --test-{kind.replace('_', '-')}")
        elif kind in ("quant", "meta", "lib_format"):
            fail(f"required {kind} file not found under {root} (looked for {', '.join(HIST_KINDS[kind])}); "
                 f"give it with --test-{kind.replace('_', '-')}")
        else:
            chosen[kind] = None
            notes.append(f"{kind}: not found (checks depending on it are UNAVAILABLE)")
    if a.reference_dir:
        rd = Path(a.reference_dir)
        if a.json_out and path_inside(Path(a.json_out).parent, rd):
            fail("--json-out must be outside the reference directory")

    if a.gene_map_md5:
        got = file_digest(Path(a.gene_map_fasta), "md5")
        if got != a.gene_map_md5:
            fail(f"gene map FASTA MD5 {got} != expected {a.gene_map_md5}")

    # ---- inputs (all read-only)
    test, tmeta, tlib = read_quant(chosen["quant"]), load_json(chosen["meta"]), load_json(chosen["lib_format"])
    tfastp = load_json(chosen["fastp"]) if chosen["fastp"] else {}
    tcmd = load_json(chosen["cmd_info"]) if chosen["cmd_info"] else {}
    summ = tfastp.get("summary", {})
    t_raw = summ.get("before_filtering", {}).get("total_reads")
    t_after = summ.get("after_filtering", {}).get("total_reads")
    ref, rmeta, rlib = read_quant(Path(a.reference_quant)), load_json(a.reference_meta), load_json(a.reference_lib_format)
    rq = read_run_qc_row(a.reference_run_qc, a.specimen)

    # ---- reference bundle integrity (MANIFEST.sha256 next to the reference files, if present)
    ref_files = [Path(a.reference_quant), Path(a.reference_meta), Path(a.reference_lib_format), Path(a.reference_run_qc)]
    manifest_ok = None
    man = Path(a.reference_manifest) if a.reference_manifest else Path(a.reference_quant).parent / "MANIFEST.sha256"
    if man.is_file():
        listed = {}
        for line in man.read_text().splitlines():
            if line.strip():
                h, _, name = line.partition("  ")
                listed[os.path.basename(name.strip())] = h.strip()
        manifest_ok = all(listed.get(f.name) == file_digest(f, "sha256") for f in ref_files)

    # ---- historical integrity (provenance is NOT required: the run predates the schema)
    n = len(test)
    integrity = {
        "quant_rows_equal_num_valid_targets": tri(tmeta.get("num_valid_targets") == n),
        "quant_rows_equal_expected_targets": tri(a.expect_targets is None or n == a.expect_targets),
        "library_format_is_expected": tri(tlib.get("expected_format") == a.expect_libtype),
        "salmon_processed_equals_fastp_retained_pairs": tri(None if t_after is None or tmeta.get("num_processed") is None
                                                            else 2 * tmeta["num_processed"] == t_after),
        "no_partial_markers_in_artifact_paths": tri(not any(".partial" in str(p) for p in chosen.values() if p)),
        "reference_bundle_manifest_verified": tri(manifest_ok),
    }

    # ---- specimen identity from recorded inputs (no provenance): cmd_info mates and fastp command line
    texts = [json.dumps(tcmd.get(k, "")) for k in ("mates1", "mates2", "unmatedReads")] + [str(tfastp.get("command", ""))]
    blob = " ".join(texts)
    cores = set(CORE_ANY.findall(blob))
    sid = rq.get("sample_id") or ""
    sid_hit = bool(sid) and re.search(r"(?<![A-Za-z0-9])" + re.escape(sid) + r"(?![0-9])", blob) is not None
    other_gm = {m for m in re.findall(r"(?<![A-Za-z0-9])GM\d+(?![0-9])", blob) if m != sid}
    if cores:
        ident = cores == {rq.get("core_id")}
    elif sid_hit or other_gm:
        ident = sid_hit and not other_gm
    else:
        ident = None
    identity = {
        "specimen_in_reference_run_qc": tri(a.specimen in (rq.get("sample_id"), rq.get("core_id"))),
        "specimen_identity_from_recorded_inputs": tri(ident),
    }

    shared, tx, gene, yield_ = concordance(test, tmeta, tlib, t_raw, t_after, ref, rmeta, rlib, rq,
                                           gene_map(Path(a.gene_map_fasta)), rule, a.expect_targets, a.expect_libtype, a.min_tpm)
    gates = {k: tri(v) for k, v in shared.items()}
    if t_raw is None:
        gates["raw_input_reads_identical"] = "UNAVAILABLE"
    all_checks = {**{f"integrity.{k}": v for k, v in integrity.items()}, **{f"identity.{k}": v for k, v in identity.items()},
                  **{f"gate.{k}": v for k, v in gates.items()}}
    if any(v == "FAIL" for v in all_checks.values()):
        verdict, code = "HISTORICAL_DISCORDANT", 1
    elif any(v == "UNAVAILABLE" for v in all_checks.values()):
        verdict, code = "HISTORICAL_INCONCLUSIVE", 2
    else:
        verdict, code = "HISTORICAL_CONCORDANT", 0

    flags = lambda k: k in tcmd  # salmon records boolean flags as keys with [] values
    parameters_reported_not_gated = {
        "salmon_version_test": tmeta.get("salmon_version"), "salmon_version_reference": rmeta.get("salmon_version"),
        "salmon_libType": tcmd.get("libType"), "salmon_validateMappings": flags("validateMappings") if tcmd else None,
        "salmon_gcBias": flags("gcBias") if tcmd else None, "salmon_seqBias": flags("seqBias") if tcmd else None,
        "salmon_threads": tcmd.get("threads"), "salmon_index": tcmd.get("index"),
        "fastp_version_reported_by_run": summ.get("fastp_version"), "fastp_command": tfastp.get("command"),
        "expected_collaborator_salmon_args": "-l A --validateMappings --gcBias --seqBias",
        "expected_collaborator_fastp_args": "--detect_adapter_for_pe -l 25",
    }
    evidence_files = {}
    for kind, pth in list(chosen.items()) + [("reference_" + f.name, f) for f in ref_files] + [("gene_map_fasta", Path(a.gene_map_fasta))]:
        if pth:
            st = pth.stat()
            evidence_files[kind] = {"path": str(pth), "realpath": os.path.realpath(pth), "is_symlink": pth.is_symlink(),
                                    "bytes": st.st_size, "mtime_utc": datetime.fromtimestamp(st.st_mtime, timezone.utc).isoformat(),
                                    "sha256": file_digest(pth, "sha256")}
    prov = {"status": "absent: run predates the provenance schema; not fabricated",
            "provenance_json_found": [str(x) for x in found["provenance"]],
            "adoption": "none: historical audit only; validate-quant-dir not applied and not relaxed"}
    report = {"kind": "historical-audit", "specimen": a.specimen, "test_dir": str(root), "rule": rule, "verdict": verdict,
              "verdict_meaning": {"HISTORICAL_CONCORDANT": "every integrity, identity and provisional-v1 gate PASS",
                                  "HISTORICAL_DISCORDANT": "at least one check FAIL",
                                  "HISTORICAL_INCONCLUSIVE": "no FAIL, but at least one check UNAVAILABLE"}[verdict],
              "provenance": prov, "discovery_notes": notes, "integrity": integrity, "identity": identity, "gates": gates,
              "transcript_level": tx, "gene_level_reported_not_gated": gene, "yield": yield_,
              "parameters_reported_not_gated": parameters_reported_not_gated, "evidence_files": evidence_files,
              "audit_run": {"utc": datetime.now(timezone.utc).isoformat(), "host": socket.gethostname(),
                            "python": platform.python_version(), "tools_sha256": file_digest(Path(__file__), "sha256")}}
    print(f"HISTORICAL AUDIT (read-only; provenance absent; not an adopted reproduction)  rule {rule['version']} ({rule['status']})")
    for x in notes:
        print(f"note       {x}")
    for kind, ev in evidence_files.items():
        print(f"file       {kind:<40} {ev['sha256'][:16]}… {ev['path']}")
    for section, d in (("integrity", integrity), ("identity", identity)):
        for k, v in d.items():
            print(f"{section:<10} {k:<48} {v}")
    for section, d in (("transcript", tx), ("gene", gene), ("yield", yield_), ("param", parameters_reported_not_gated)):
        for k, v in d.items():
            print(f"{section:<10} {k:<48} {v}")
    for k, v in gates.items():
        print(f"gate       {k:<48} {v}")
    print(f"VERDICT {verdict} ({rule['version']}; historical, provenance absent)")
    if a.json_out:
        write_new(Path(a.json_out), json.dumps(report, indent=2) + "\n")
    sys.exit(code)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest="cmd", required=True)

    s = sp.add_parser("provenance")
    s.add_argument("--out", required=True)
    s.add_argument("--set", action="append", help="key=value (repeatable)")
    s.add_argument("--json", action="append", help="key=path.json: embed a JSON file (repeatable)")
    s.add_argument("--config", action="append", help="name=value::status::source (repeatable)")
    s.set_defaults(fn=cmd_provenance)

    s = sp.add_parser("resolve-specimen")
    s.add_argument("--sheet", required=True)
    s.add_argument("--validated-pairs", required=True, help="validated_pairs.tsv (SAMPLE<TAB>R1<TAB>R2, gzip-validated pairs)")
    s.add_argument("--specimen", required=True, help="sample_id (GM01), core_id (20088XR-11-01) or short id (11-01)")
    s.add_argument("--raw-dir", help="RAW_DIR; validated absolute paths must be the same files")
    s.add_argument("--exclude", default="", help="additional excluded ids (space/comma separated)")
    s.add_argument("--lift-exclusion", action="append", help="CORE_ID=<adjudication reference> for a built-in exclusion")
    s.set_defaults(fn=cmd_resolve)

    s = sp.add_parser("verify-fingerprint")
    s.add_argument("--observed", required=True, help="Salmon index dir, its info.json, or a quant aux_info/meta_info.json")
    s.add_argument("--reference", help="reference meta_info.json/info.json/index dir (e.g. committed collaborator run)")
    for _, _, c in FINGERPRINTS:
        s.add_argument(f"--expect-{c.replace('_', '-')}", dest=f"expect_{c}")
    s.add_argument("--expect-k", type=int, help="k requested for THIS build; the collaborator's k stays unverified")
    s.add_argument("--require", default=",".join(c for _, _, c in FINGERPRINTS))
    s.add_argument("--json-out")
    s.set_defaults(fn=cmd_fingerprint)

    s = sp.add_parser("validate-quant-dir")
    s.add_argument("--dir", required=True)
    s.add_argument("--stage", action="store_true", help="staged work dir (no provenance/status requirement)")
    s.add_argument("--specimen")
    s.add_argument("--expect-targets", type=int)
    s.add_argument("--expect-libtype")
    s.set_defaults(fn=cmd_validate)

    s = sp.add_parser("compare-quant")
    s.add_argument("--specimen", required=True, help="sample_id or core_id; must match test provenance and reference run_qc row")
    s.add_argument("--test-dir", required=True, help="adopted specimen output directory from quantify_specimen.sbatch")
    s.add_argument("--reference-quant", required=True)
    s.add_argument("--reference-meta", required=True)
    s.add_argument("--reference-lib-format", required=True)
    s.add_argument("--reference-run-qc", required=True, help="committed qc/run_qc.tsv (raw/trimmed read counts)")
    s.add_argument("--gene-map-fasta", required=True, help="GENCODE v46 transcripts FASTA(.gz) for gene-level aggregation")
    s.add_argument("--expect-targets", type=int, default=253181)
    s.add_argument("--expect-libtype", default="ISR")
    s.add_argument("--min-tpm", type=float, default=1.0)
    s.add_argument("--min-pearson", type=float, default=RULE_V1["min_pearson_log1p_tpm_transcript"])
    s.add_argument("--max-mapping-diff-pp", type=float, default=RULE_V1["max_abs_mapping_diff_pp"])
    s.add_argument("--rule-version")
    s.add_argument("--rule-justification")
    s.add_argument("--json-out")
    s.set_defaults(fn=cmd_compare)

    s = sp.add_parser("audit-historical", help="read-only audit of a pre-provenance-schema run (never adopts)")
    s.add_argument("--specimen", required=True)
    s.add_argument("--test-dir", required=True, help="historical run directory (read only; artifacts discovered by name)")
    s.add_argument("--test-quant"); s.add_argument("--test-meta"); s.add_argument("--test-lib-format")
    s.add_argument("--test-fastp-json"); s.add_argument("--test-cmd-info")
    s.add_argument("--reference-dir", help="directory holding the reference files (only used for path checks)")
    s.add_argument("--reference-quant", required=True)
    s.add_argument("--reference-meta", required=True)
    s.add_argument("--reference-lib-format", required=True)
    s.add_argument("--reference-run-qc", required=True)
    s.add_argument("--reference-manifest", help="MANIFEST.sha256 for the reference files (default: next to --reference-quant)")
    s.add_argument("--gene-map-fasta", required=True)
    s.add_argument("--gene-map-md5", help="expected MD5 of the gene map FASTA (GENCODE v46: f1549843717b0fecfa4336fe3d5cc85f)")
    s.add_argument("--expect-targets", type=int, default=253181)
    s.add_argument("--expect-libtype", default="ISR")
    s.add_argument("--min-tpm", type=float, default=1.0)
    s.add_argument("--min-pearson", type=float, default=RULE_V1["min_pearson_log1p_tpm_transcript"])
    s.add_argument("--max-mapping-diff-pp", type=float, default=RULE_V1["max_abs_mapping_diff_pp"])
    s.add_argument("--rule-version")
    s.add_argument("--rule-justification")
    s.add_argument("--json-out")
    s.set_defaults(fn=cmd_audit)

    a = p.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()
