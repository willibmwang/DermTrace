#!/usr/bin/env python3
"""Mixed 42-specimen versus original 33-specimen matrix: exploratory technical comparison (stdlib only, Python >= 3.9).

Pipeline origin (collaborator_committed vs rivanna_qc_r3) is tracked as a TECHNICAL covariate next to sequencing run
and lane. Nothing is corrected, excluded, filtered away or modelled; no supervised model is trained.

  1. Integrity: the 33 committed columns of the expanded TPM matrix equal the committed matrix value for value.
  2. Frozen PCA (33 specimens; protein-coding, TPM>=1 in >=50 %, log1p, mean-centred, unscaled; 11,022 genes):
     the 33 committed columns are re-projected (must reproduce the frozen scores) and the nine new specimens are
     projected with the same means and loadings (a new specimen cannot move the frozen axes).
  3. Refitted PCA on 42 specimens, (a) on the frozen gene set and (b) with the same recipe re-applied to 42.
  4. Association of PCs and technical QC with pipeline origin, sequencing run and lane: eta-squared with Monte Carlo
     permutation p (exploratory). Origin is confounded with specimen identity and possibly with run/lane; a
     cross-tabulation is reported so that an origin association is not over-read.
  5. Context: the cross-specimen validation bounds the effect of the pipeline itself on frozen PC1-3 at
     --pipeline-bound-sd cohort SD; origin differences far above that bound cannot come from quantification alone.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import operator
import os
import platform
import random
import socket
import sys
from datetime import datetime, timezone
from pathlib import Path

PREFIX = "tsaolab_melanoma"


def fail(msg, code=2):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


def read_matrix(p):
    with gzip.open(p, "rt") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        rows = {}
        for l in fh:
            f = l.rstrip("\n").split("\t")
            rows[f[0]] = f[1:]
    return h[1:], rows


def jacobi(a, tol=1e-12, sweeps=100):
    n = len(a); a = [r[:] for r in a]; v = [[float(i == j) for j in range(n)] for i in range(n)]
    for _ in range(sweeps):
        off = sum(a[i][j] ** 2 for i in range(n) for j in range(i + 1, n))
        if off < tol * max(1e-300, sum(a[i][i] ** 2 for i in range(n))):
            break
        for p in range(n):
            for q in range(p + 1, n):
                if abs(a[p][q]) < 1e-300:
                    continue
                th = (a[q][q] - a[p][p]) / (2 * a[p][q])
                t = (1 if th >= 0 else -1) / (abs(th) + math.sqrt(th * th + 1))
                c = 1 / math.sqrt(t * t + 1); s = t * c
                for k in range(n):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p], a[k][q] = c * akp - s * akq, s * akp + c * akq
                for k in range(n):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k], a[q][k] = c * apk - s * aqk, s * apk + c * aqk
                for k in range(n):
                    vkp, vkq = v[k][p], v[k][q]
                    v[k][p], v[k][q] = c * vkp - s * vkq, s * vkp + c * vkq
    order = sorted(range(n), key=lambda i: -a[i][i])
    return [a[i][i] for i in order], [[v[k][i] for k in range(n)] for i in order]


def pca(columns, genes, data, k=5):
    """data[g] = list of log1p values per column; returns explained variance ratios and scores (n x k)."""
    n = len(columns)
    cols = [[0.0] * len(genes) for _ in range(n)]
    for gi, g in enumerate(genes):
        x = data[g]; m = sum(x) / n
        for i in range(n):
            cols[i][gi] = x[i] - m
    gram = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            gram[i][j] = gram[j][i] = sum(map(operator.mul, cols[i], cols[j]))
    lam, vec = jacobi(gram)
    tot = sum(max(0.0, l) for l in lam)
    scores = [[vec[c][i] * math.sqrt(max(0.0, lam[c])) for c in range(k)] for i in range(n)]
    return [max(0.0, l) / tot for l in lam[:k]], scores


def eta2(y, g):
    m = sum(y) / len(y); tot = sum((v - m) ** 2 for v in y)
    if tot == 0:
        return 0.0
    grp = {}
    for v, k in zip(y, g):
        grp.setdefault(k, []).append(v)
    return sum(len(v) * (sum(v) / len(v) - m) ** 2 for v in grp.values()) / tot


def perm_p(y, g, n=20000, seed=20261126):
    rng = random.Random(seed); obs = eta2(y, g); g2 = list(g); hit = 0
    for _ in range(n):
        rng.shuffle(g2)
        hit += eta2(y, g2) >= obs - 1e-12
    return (hit + 1) / (n + 1)


def write_new(path, text):
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    with os.fdopen(fd, "w") as fh:
        fh.write(text)


def tsv(path, rows):
    cols = list(rows[0]) if rows else []
    f = lambda v: "" if v is None else (f"{v:.4f}" if isinstance(v, float) else str(v))
    write_new(path, "\n".join(["\t".join(cols)] + ["\t".join(f(r.get(c)) for c in cols) for r in rows]) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--expanded-dir", required=True)
    ap.add_argument("--committed-dir", required=True)
    ap.add_argument("--qc-manifest", required=True)
    ap.add_argument("--cohort-context", required=True, help="gm01_cohort_context.tsv (frozen PCA gene set and loadings)")
    ap.add_argument("--frozen-scores", required=True, help="pca_frozen_scores.json")
    ap.add_argument("--pipeline-bound-sd", type=float, default=0.0049)
    ap.add_argument("--permutations", type=int, default=20000)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    out = Path(a.out_dir)
    if os.path.lexists(out):
        fail(f"output directory exists: {out}")

    man = json.loads((Path(a.expanded_dir) / "expansion_manifest.json").read_text())
    origin = man["pipeline_origin"]
    cols, rows = read_matrix(Path(a.expanded_dir) / f"{PREFIX}_tpm_ensembl.tsv.gz")
    ccols, crows = read_matrix(Path(a.committed_dir) / f"{PREFIX}_tpm_ensembl.tsv.gz")
    with gzip.open(Path(a.expanded_dir) / f"{PREFIX}_genes.tsv.gz", "rt") as fh:
        gtype = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2] for l in list(fh)[1:]}
    with open(a.qc_manifest, newline="") as fh:
        qc = {r["sample_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    ctx = {}
    with open(a.cohort_context, newline="") as fh:
        for r in csv.DictReader([l for l in fh if not l.startswith("#")], delimiter="\t"):
            if r["in_pca_gene_set"] == "1":
                ctx[r["gene_id"]] = [float(r[f"pc{i}_loading"]) for i in (1, 2, 3)]
    frozen = json.loads(Path(a.frozen_scores).read_text())
    new = [c for c in cols if origin[c] != "collaborator_committed"]

    # 1. integrity, independently of the expansion manifest
    n0 = len(ccols)
    integ = {"committed_columns_first_and_unchanged_order": cols[:n0] == ccols,
             "committed_values_identical": all(rows[g][:n0] == v for g, v in crows.items()),
             "committed_genes_all_present": set(crows) <= set(rows),
             "added_rows_zero_in_committed_columns": all(float(x) == 0.0 for g, v in rows.items() if g not in crows for x in v[:n0]),
             "n_columns": len(cols), "n_new_columns": len(new), "n_rows_committed": len(crows), "n_rows_expanded": len(rows)}

    lg = {g: [math.log1p(float(x)) for x in v] for g, v in rows.items()}
    idx33 = [cols.index(c) for c in ccols]
    # 2. frozen projection
    missing = [g for g in ctx if g not in lg]
    if missing:
        fail(f"{len(missing)} frozen PCA genes missing from the expanded matrix")
    mu = {g: sum(lg[g][i] for i in idx33) / n0 for g in ctx}
    proj = {}
    for ci, c in enumerate(cols):
        proj[c] = [sum((lg[g][ci] - mu[g]) * ctx[g][k] for g in ctx) for k in range(3)]
    repro = max(abs(abs(proj[c][k]) - abs(frozen["scores"][c][f"PC{k+1}"])) for c in ccols for k in range(3))
    sd = frozen["score_sd"]
    rng33 = {k: (min(proj[c][k] for c in ccols), max(proj[c][k] for c in ccols)) for k in range(3)}
    proj_rows = [{"specimen": c, "pipeline_origin": origin[c], "sequencing_run": qc.get(c, {}).get("sequencing_run"),
                  **{f"PC{k+1}": proj[c][k] for k in range(3)}, **{f"PC{k+1}_z_cohort_sd": proj[c][k] / sd[f"PC{k+1}"] for k in range(3)},
                  "inside_33_range_all_pc1_3": all(rng33[k][0] <= proj[c][k] <= rng33[k][1] for k in range(3))} for c in cols]
    mean_shift = {f"PC{k+1}": (sum(proj[c][k] for c in new) / len(new) - sum(proj[c][k] for c in ccols) / n0) / sd[f"PC{k+1}"] for k in range(3)} if new else {}

    # 3. refitted PCAs
    fz = sorted(ctx)
    ev_a, sc_a = pca(cols, fz, lg)
    pc_genes = sorted(g for g in rows if gtype.get(g) == "protein_coding" and sum(float(x) >= 1 for x in rows[g]) >= 0.5 * len(cols))
    ev_b, sc_b = pca(cols, pc_genes, lg)

    # 4. associations
    covs = {"pipeline_origin": [origin[c] for c in cols], "sequencing_run": [qc[c]["sequencing_run"] for c in cols],
            "lane": [qc[c]["lane"] for c in cols]}
    assoc = []
    for label, scores in (("frozen_projection", [proj[c] for c in cols]), ("refit_frozen_genes", sc_a), ("refit_recipe_42", sc_b)):
        for k in range(3):
            y = [s[k] for s in scores]
            for cv, g in covs.items():
                assoc.append({"pca": label, "pc": f"PC{k+1}", "covariate": cv, "eta2": eta2(y, g), "perm_p_exploratory": perm_p(y, g, a.permutations)})
    tech = []
    for col in ("mapping_rate_pct", "fastp_duplication_rate", "pct_adapter_trimmed", "insert_size_peak", "reads_raw", "genes_detected_tpm_ge_1", "pct_q30_after_trim"):
        vals = [(float(qc[c][col]), origin[c]) for c in cols if qc[c].get(col) not in (None, "", "nan")]
        if len(vals) == len(cols):
            y, g = [v for v, _ in vals], [o for _, o in vals]
            by = {o: sorted(v for v, oo in vals if oo == o) for o in set(g)}
            tech.append({"qc_measure": col, **{f"median_{o}": by[o][len(by[o]) // 2] for o in sorted(by)}, "eta2_origin": eta2(y, g),
                         "perm_p_exploratory": perm_p(y, g, a.permutations)})
    xtab = {}
    for c in cols:
        key = (origin[c], qc[c]["sequencing_run"], qc[c]["lane"])
        xtab[key] = xtab.get(key, 0) + 1
    xt_rows = [{"pipeline_origin": k[0], "sequencing_run": k[1], "lane": k[2], "n": v} for k, v in sorted(xtab.items())]
    only_new = sum(1 for g, v in rows.items() if all(float(x) < 1 for x in v[:n0]) and any(float(v[cols.index(c)]) >= 1 for c in new))

    summary = {"integrity": integ, "frozen_projection": {"max_abs_reproduction_error_33": repro, "new_specimens_mean_shift_in_cohort_sd": mean_shift,
               "new_inside_33_range": {r["specimen"]: r["inside_33_range_all_pc1_3"] for r in proj_rows if r["pipeline_origin"] != "collaborator_committed"},
               "pipeline_effect_bound_sd_from_cross_specimen_validation": a.pipeline_bound_sd},
               "refit": {"frozen_gene_set": {"n_genes": len(fz), "explained_variance": [round(x, 4) for x in ev_a]},
                         "recipe_on_42": {"n_genes": len(pc_genes), "explained_variance": [round(x, 4) for x in ev_b]}},
               "genes_detected_only_in_new_specimens_tpm_ge_1": only_new, "added_gene_rows": len(man.get("added_genes", [])),
               "caveats": ["exploratory; n_new = %d" % len(new), "pipeline origin is confounded with specimen identity and possibly run/lane",
                           "no batch correction, feature removal or supervised modelling performed"],
               "run": {"utc": datetime.now(timezone.utc).isoformat(), "host": socket.gethostname(), "python": platform.python_version()}}
    out.mkdir(parents=True)
    tsv(out / "frozen_projection.tsv", proj_rows); tsv(out / "associations_exploratory.tsv", assoc); tsv(out / "technical_qc_by_origin.tsv", tech)
    tsv(out / "origin_run_lane_crosstab.tsv", xt_rows)
    tsv(out / "refit_scores.tsv", [{"specimen": c, "pipeline_origin": origin[c], **{f"frozen_genes_PC{k+1}": sc_a[i][k] for k in range(3)},
                                    **{f"recipe42_PC{k+1}": sc_b[i][k] for k in range(3)}} for i, c in enumerate(cols)])
    write_new(out / "summary.json", json.dumps(summary, indent=2) + "\n")

    print("INTEGRITY", json.dumps(integ))
    print(f"FROZEN PROJECTION: 33 reproduce frozen scores within {repro:.4f}; new-specimen mean shift (cohort SD): "
          + ", ".join(f"{k}={v:+.3f}" for k, v in mean_shift.items()) + f"  [pipeline-only bound {a.pipeline_bound_sd} SD]")
    for r in proj_rows:
        if r["pipeline_origin"] != "collaborator_committed":
            print(f"  {r['specimen']} run={r['sequencing_run']} PC1={r['PC1_z_cohort_sd']:+.2f}SD PC2={r['PC2_z_cohort_sd']:+.2f}SD PC3={r['PC3_z_cohort_sd']:+.2f}SD inside33={r['inside_33_range_all_pc1_3']}")
    print(f"REFIT explained variance: frozen genes {summary['refit']['frozen_gene_set']['explained_variance'][:3]}; recipe-42 ({len(pc_genes)} genes) {summary['refit']['recipe_on_42']['explained_variance'][:3]}")
    print("ASSOCIATIONS (eta2, permutation p; exploratory)")
    for x in assoc:
        print(f"  {x['pca']:<20}{x['pc']:<5}{x['covariate']:<16} eta2={x['eta2']:.3f} p={x['perm_p_exploratory']:.4f}")
    print("TECHNICAL QC BY ORIGIN"); [print("  " + "  ".join(f"{k}={v:.4f}" if isinstance(v, float) else f"{k}={v}" for k, v in t.items())) for t in tech]
    print("ORIGIN x RUN x LANE"); [print(f"  {r['pipeline_origin']:<24}{r['sequencing_run']:<16}{r['lane']:<6}{r['n']}") for r in xt_rows]
    print(f"genes detected (TPM>=1) only in new specimens: {only_new}; added gene rows: {summary['added_gene_rows']}")
    print(f"outputs: {out}")
    ok = all(v for v in integ.values() if isinstance(v, bool)) and repro < 0.05
    print(f"INTEGRITY+REPRODUCTION {'PASS' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
