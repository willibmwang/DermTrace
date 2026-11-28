#!/usr/bin/env python3
"""Step 17 recovery helpers (r1; Python >= 3.9, stdlib only). Replaces shell gzip|head pipelines and adds strict
evidence validation. Exit codes: 0 ok / all PASS, 1 failed check or invalid evidence, 3 valid evidence with
scientific concordance failures, 2 usage or unreadable input.

  header    read the header line of a gzip matrix deterministically (no pipes, no SIGPIPE) and require EXACTLY the
            expected specimen IDs, unique and in the expected order; prints them space-separated
  reports   validate the provisional-v1 concordance reports: exactly one valid report per expected specimen, no
            extra reports, frozen rule (provisional-v1; 0.99; 1.0 pp), the 14 expected gates, verdict consistent with
            the gates and with compare-quant's exit code. Each specimen is classified as
              PASS                        all 14 gates true
              CONCORDANCE_FAIL            only the quantitative gates failed (transcript Pearson, mapping difference)
              IDENTITY_OR_INTEGRITY_FAIL  an identity/integrity gate failed (specimen, raw reads, IDs, targets,
                                          library orientation, fingerprints, test-output validity)
              EVIDENCE_INVALID            report missing, unparsable, duplicated, for the wrong specimen or rule,
                                          incomplete, inconsistent, or compare-quant could not evaluate (exit 2)
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
import zlib
from pathlib import Path

RULE = {"version": "provisional-v1", "min_pearson_log1p_tpm_transcript": 0.99, "max_abs_mapping_diff_pp": 1.0}
QUANT_GATES = {"transcript_pearson_log1p_tpm_ge_0.99", "abs_mapping_diff_le_1pp"}
INTEGRITY_GATES = {"test_output_complete_and_valid", "specimen_identity", "raw_input_reads_identical", "identical_transcript_ids",
                   "target_counts_equal", "library_format_matches"} | {f"index_{c}_identical" for c in (
                       "seq_hash", "name_hash", "seq_hash512", "name_hash512", "decoy_seq_hash", "decoy_name_hash")}
ALL_GATES = QUANT_GATES | INTEGRITY_GATES          # 14


def fail(msg, code=1):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


def cmd_header(a):
    p = Path(a.matrix)
    if not p.is_file():
        fail(f"matrix not found: {p}", 2)
    try:
        with gzip.open(p, "rt", newline="") as fh:
            line = fh.readline()
    except (OSError, EOFError, zlib.error, UnicodeDecodeError) as e:
        fail(f"cannot read gzip header of {p}: {e}", 2)
    if not line.endswith("\n"):
        fail(f"{p}: header line missing or unterminated (empty or truncated file)")
    f = line.rstrip("\n").split("\t")
    if len(f) < 2 or f[0] not in ("ensembl_gene_id", "gene_symbol"):
        fail(f"{p}: malformed header (first field {f[0]!r}, {len(f) - 1} columns)")
    ids = f[1:]
    expect = a.expect.split(",")
    if len(expect) != len(set(expect)):
        fail("--expect contains duplicates", 2)
    dup = sorted({x for x in ids if ids.count(x) > 1})
    if dup:
        fail(f"{p}: duplicate specimen IDs in header: {dup}")
    if any(not x for x in ids):
        fail(f"{p}: empty specimen ID in header")
    if ids != expect:
        missing, extra = [x for x in expect if x not in ids], [x for x in ids if x not in expect]
        fail(f"{p}: header has {len(ids)} specimens; expected exactly {len(expect)} in order"
             + (f"; missing {missing}" if missing else "") + (f"; unexpected {extra}" if extra else "")
             + ("; order differs" if not missing and not extra else ""))
    print(" ".join(ids))


def classify(sid, path, rc):
    """Return (class, detail) for one specimen."""
    if rc is None:
        return "EVIDENCE_INVALID", "no compare-quant exit code recorded"
    if rc == 2:
        return "EVIDENCE_INVALID", "compare-quant could not evaluate (exit 2: missing or unreadable evidence)"
    if rc not in (0, 1):
        return "EVIDENCE_INVALID", f"unexpected compare-quant exit code {rc}"
    if not path.is_file():
        return "EVIDENCE_INVALID", "report missing"
    try:
        d = json.loads(path.read_text())
    except (OSError, ValueError) as e:
        return "EVIDENCE_INVALID", f"report unparsable: {e}"
    if not isinstance(d, dict):
        return "EVIDENCE_INVALID", "report is not a JSON object"
    if d.get("specimen") != sid:
        return "EVIDENCE_INVALID", f"report specimen {d.get('specimen')!r} != {sid!r}"
    rule = d.get("rule") or {}
    for k, v in RULE.items():
        if rule.get(k) != v:
            return "EVIDENCE_INVALID", f"rule {k}={rule.get(k)!r}, expected {v!r} (provisional-v1 is frozen)"
    g = d.get("gates")
    if not isinstance(g, dict) or set(g) != ALL_GATES or not all(isinstance(v, bool) for v in g.values()):
        got = sorted(g) if isinstance(g, dict) else g
        return "EVIDENCE_INVALID", f"gates are not exactly the 14 expected booleans: {got}"
    verdict = d.get("verdict")
    if verdict not in ("PASS", "FAIL") or (verdict == "PASS") != all(g.values()):
        return "EVIDENCE_INVALID", f"verdict {verdict!r} inconsistent with gates"
    if (rc == 0) != (verdict == "PASS"):
        return "EVIDENCE_INVALID", f"compare-quant exit {rc} inconsistent with verdict {verdict}"
    for sec in ("transcript_level", "gene_level_reported_not_gated", "yield"):
        if not isinstance(d.get(sec), dict) or not d[sec]:
            return "EVIDENCE_INVALID", f"section {sec} missing or empty"
    failed = sorted(k for k, v in g.items() if not v)
    if not failed:
        return "PASS", ""
    if set(failed) <= QUANT_GATES:
        return "CONCORDANCE_FAIL", ";".join(failed)
    return "IDENTITY_OR_INTEGRITY_FAIL", ";".join(failed)


def cmd_reports(a):
    d = Path(a.dir)
    expect = a.expect.split(",")
    if len(expect) != len(set(expect)) or (a.expect_n and len(expect) != a.expect_n):
        fail(f"--expect must list {a.expect_n} unique specimens", 2)
    codes = {}
    with open(a.exit_codes) as fh:
        for n, line in enumerate(fh, 1):
            f = line.rstrip("\n").split("\t")
            if len(f) != 2 or not f[1].lstrip("-").isdigit():
                fail(f"malformed exit-code line {n}: {line.rstrip()!r}", 2)
            if f[0] in codes:
                fail(f"exit code recorded twice for {f[0]}", 2)
            codes[f[0]] = int(f[1])
    found = sorted(x.name for x in d.glob("*_provisional-v1.json"))
    extra = [x for x in found if x[: -len("_provisional-v1.json")] not in expect]
    rows = []
    for sid in expect:
        cls, detail = classify(sid, d / f"{sid}_provisional-v1.json", codes.get(sid))
        rep = {}
        if cls != "EVIDENCE_INVALID":
            rep = json.loads((d / f"{sid}_provisional-v1.json").read_text())
        t, gl, y = rep.get("transcript_level", {}), rep.get("gene_level_reported_not_gated", {}), rep.get("yield", {})
        rows.append([sid, cls, detail or "-", codes.get(sid, ""), t.get("pearson_log1p_tpm", ""), gl.get("pearson_log1p_tpm", ""),
                     f"{t.get('n_expressed_abs_log2fc_gt_1', '')}/{t.get('n_expressed_tpm_ge_1', '')}" if t else "",
                     f"{gl.get('n_expressed_abs_log2fc_gt_1', '')}/{gl.get('n_expressed_tpm_ge_1', '')}" if gl else "",
                     y.get("percent_mapped_diff_pp", ""), y.get("reads_after_trim_diff", "")])
    unexpected_codes = sorted(set(codes) - set(expect))
    counts = {c: sum(1 for r in rows if r[1] == c) for c in ("PASS", "CONCORDANCE_FAIL", "IDENTITY_OR_INTEGRITY_FAIL", "EVIDENCE_INVALID")}
    hdr = ["specimen", "class", "detail", "compare_exit", "tx_pearson", "gene_pearson", "tx_discordant", "gene_discordant", "mapping_diff_pp", "trimmed_reads_diff"]
    text = "\t".join(hdr) + "\n" + "".join("\t".join(map(str, r)) + "\n" for r in rows)
    if a.out_tsv:
        fd = os.open(a.out_tsv, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
        with os.fdopen(fd, "w") as fh:
            fh.write(text)
    sys.stdout.write(text)
    valid = len(expect) - counts["EVIDENCE_INVALID"]
    print(f"REPORTS: {len(expect)} expected; {valid} valid and uniquely identified; " + ", ".join(f"{k}={v}" for k, v in counts.items()))
    if extra:
        print(f"ERROR: unexpected report files: {extra}")
    if unexpected_codes:
        print(f"ERROR: exit codes recorded for unexpected specimens: {unexpected_codes}")
    if counts["EVIDENCE_INVALID"] or extra or unexpected_codes:
        print("EVIDENCE INVALID: the concordance evidence is incomplete or malformed (not a scientific result)")
        sys.exit(1)
    if counts["PASS"] == len(expect):
        print(f"CONCORDANCE provisional-v1: {counts['PASS']}/{len(expect)} PASS")
        sys.exit(0)
    print(f"CONCORDANCE provisional-v1: {counts['PASS']}/{len(expect)} PASS; scientific failures are listed above, not hidden")
    sys.exit(3)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    s = sp.add_parser("header"); s.add_argument("--matrix", required=True); s.add_argument("--expect", required=True); s.set_defaults(fn=cmd_header)
    s = sp.add_parser("reports"); s.add_argument("--dir", required=True); s.add_argument("--expect", required=True)
    s.add_argument("--expect-n", type=int); s.add_argument("--exit-codes", required=True); s.add_argument("--out-tsv")
    s.set_defaults(fn=cmd_reports)
    a = ap.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()
