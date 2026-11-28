#!/usr/bin/env python3
"""Technical-confounding analysis across expression representations (exploratory; stdlib only, Python >= 3.9).

Representations: original33 (committed collaborator matrix), mixed42 (33 collaborator + 9 r3 columns) and
harmonized42 (all 42 from r2/r3). For each representation, and again without the pre-declared suspect specimens:

  PCA        protein-coding genes with TPM >= 1 in >= 50 % of the analysed specimens; log1p; per-gene mean centring;
             no scaling (the frozen exploratory recipe). PC1-5 scores and explained variance; top loadings.
  confounds  PC x covariate association: sequencing run, lane, insert size, duplication, mapping rate, adapter %,
             depth, inferred sex, quality flag, specimen set (original 33 vs expansion 9) and, in mixed42 only,
             pipeline origin. eta-squared (categorical) / Spearman (continuous); permutation p; BH q per analysis.
             Covariate-covariate relations are reported too (confounding among covariates themselves).
  tech/bio   per PC: R-squared from technical covariates (lane, insert, duplication, mapping) versus marker-based
             biological axes (melanocytic, keratinocyte, immune, stromal mean log1p TPM) and their unique shares.
             Markers can themselves carry technical effects; this separates covariance "as far as the data permit".
  stability  leave-one-out and repeated 80 % subsampling: |correlation| of matched PC scores with the full fit;
             k-means (k = 2..4) on PC1-5 with subsample stability (adjusted Rand index) and cluster-covariate ARI.
  cross-rep  same specimens, different quantification: original33 vs harmonized (33 shared) isolates the pipeline;
             mixed42 vs harmonized42 (42 shared) isolates the effect of mixing pipelines. A specimen-set association
             that persists in harmonized42 cannot be a pipeline effect.

No batch correction, feature removal or supervised modelling. Outputs go to a NEW directory.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import platform
import random
import socket
import sys
from datetime import datetime, timezone
from pathlib import Path

MARKERS = {
    "melanocytic": ["MITF", "PMEL", "MLANA", "TYR", "TYRP1", "DCT", "SOX10", "S100B"],
    "keratinocyte": ["KRT1", "KRT5", "KRT10", "KRT14", "KRT15", "LOR", "FLG", "DSG1"],
    "immune": ["PTPRC", "CD3E", "CD2", "CD8A", "CD74", "LYZ", "CD68"],
    "stromal": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"],
}
CATEGORICAL = ["sequencing_run", "lane", "sex_inferred", "quality_flag", "specimen_set", "pipeline_origin"]
CONTINUOUS = ["insert_size_peak", "fastp_duplication_rate", "mapping_rate_pct", "pct_adapter_trimmed", "reads_raw"]
TECH = ["insert_size_peak", "fastp_duplication_rate", "mapping_rate_pct"]
K = 5


def fail(msg, code=2):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


# ------------------------------------------------------------------------------------------ linear algebra
def jacobi(a, sweeps=100):
    n = len(a); a = [r[:] for r in a]; v = [[float(i == j) for j in range(n)] for i in range(n)]
    for _ in range(sweeps):
        off = sum(a[i][j] ** 2 for i in range(n) for j in range(i + 1, n))
        if off < 1e-22 * max(1e-300, sum(a[i][i] ** 2 for i in range(n))):
            break
        for p in range(n):
            for q in range(p + 1, n):
                if abs(a[p][q]) < 1e-300:
                    continue
                th = (a[q][q] - a[p][p]) / (2 * a[p][q]); t = (1 if th >= 0 else -1) / (abs(th) + math.sqrt(th * th + 1))
                c = 1 / math.sqrt(t * t + 1); s = t * c
                for k in range(n):
                    a[k][p], a[k][q] = c * a[k][p] - s * a[k][q], s * a[k][p] + c * a[k][q]
                for k in range(n):
                    a[p][k], a[q][k] = c * a[p][k] - s * a[q][k], s * a[p][k] + c * a[q][k]
                for k in range(n):
                    v[k][p], v[k][q] = c * v[k][p] - s * v[k][q], s * v[k][p] + c * v[k][q]
    order = sorted(range(n), key=lambda i: -a[i][i])
    return [a[i][i] for i in order], [[v[k][i] for k in range(n)] for i in order]


def pca_from_gram(G, idx, k=K):
    """PCA scores for samples idx from the uncentred inner-product matrix G (double centring = per-gene centring)."""
    n = len(idx); sub = [[G[i][j] for j in idx] for i in idx]
    rm = [sum(r) / n for r in sub]; tm = sum(rm) / n
    gc = [[sub[i][j] - rm[i] - rm[j] + tm for j in range(n)] for i in range(n)]
    lam, vec = jacobi(gc); tr = sum(gc[i][i] for i in range(n))
    if tr <= 0:
        fail("no variance among the analysed specimens (identical expression profiles?)")
    scores = [[vec[c][i] * math.sqrt(max(0.0, lam[c])) for c in range(k)] for i in range(n)]
    return [max(0.0, l) / tr for l in lam[:k]], scores, vec


def ols_r2(y, X):
    n, p = len(y), len(X[0])
    A = [[sum(X[r][i] * X[r][j] for r in range(n)) + (1e-9 if i == j else 0) for j in range(p)] for i in range(p)]
    b = [sum(X[r][i] * y[r] for r in range(n)) for i in range(p)]
    for c in range(p):  # Gauss-Jordan with partial pivoting
        piv = max(range(c, p), key=lambda r: abs(A[r][c])); A[c], A[piv] = A[piv], A[c]; b[c], b[piv] = b[piv], b[c]
        d = A[c][c] or 1e-12
        for r in range(p):
            if r != c:
                f = A[r][c] / d
                A[r] = [x - f * yv for x, yv in zip(A[r], A[c])]; b[r] -= f * b[c]
    beta = [b[i] / (A[i][i] or 1e-12) for i in range(p)]
    fit = [sum(X[r][i] * beta[i] for i in range(p)) for r in range(n)]
    m = sum(y) / n; sst = sum((v - m) ** 2 for v in y); sse = sum((v - f) ** 2 for v, f in zip(y, fit))
    r2 = 1 - sse / sst if sst > 0 else 0.0
    adj = 1 - (1 - r2) * (n - 1) / (n - p) if n > p else None
    return r2, adj


# ------------------------------------------------------------------------------------------ statistics
def ranks(v):
    o = sorted(range(len(v)), key=v.__getitem__); r = [0.0] * len(v); i = 0
    while i < len(o):
        j = i
        while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]:
            j += 1
        for t in range(i, j + 1):
            r[o[t]] = (i + j) / 2 + 1
        i = j + 1
    return r


def pearson(x, y):
    n = len(x); mx, my = sum(x) / n, sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y)); sxx = sum((a - mx) ** 2 for a in x); syy = sum((b - my) ** 2 for b in y)
    return sxy / math.sqrt(sxx * syy) if sxx > 0 and syy > 0 else float("nan")


def spearman(x, y):
    return pearson(ranks(x), ranks(y))


def eta2(y, g):
    m = sum(y) / len(y); tot = sum((v - m) ** 2 for v in y)
    if tot == 0:
        return 0.0
    grp = {}
    for v, k in zip(y, g):
        grp.setdefault(k, []).append(v)
    return sum(len(v) * (sum(v) / len(v) - m) ** 2 for v in grp.values()) / tot


def cramers_v(a, b):
    ca, cb = sorted(set(a)), sorted(set(b)); n = len(a)
    if len(ca) < 2 or len(cb) < 2:
        return 0.0
    t = {(x, y): 0 for x in ca for y in cb}
    for x, y in zip(a, b):
        t[(x, y)] += 1
    ra = {x: sum(t[(x, y)] for y in cb) for x in ca}; rb = {y: sum(t[(x, y)] for x in ca) for y in cb}
    chi = sum((t[(x, y)] - ra[x] * rb[y] / n) ** 2 / (ra[x] * rb[y] / n) for x in ca for y in cb if ra[x] * rb[y])
    return math.sqrt(chi / (n * (min(len(ca), len(cb)) - 1)))


def perm_p(stat_fn, y, g, n, rng):
    obs = stat_fn(y, g); g2 = list(g); hit = 0
    for _ in range(n):
        rng.shuffle(g2)
        hit += stat_fn(y, g2) >= obs - 1e-12
    return (hit + 1) / (n + 1)


def bh(ps):
    m = len(ps); order = sorted(range(m), key=lambda i: ps[i]); q = [0.0] * m; prev = 1.0
    for rank, i in reversed(list(enumerate(order, 1))):
        prev = min(prev, ps[i] * m / rank); q[i] = prev
    return q


def ari(a, b):
    n = len(a)
    if n < 2:
        return float("nan")
    c2 = lambda x: x * (x - 1) / 2
    t = {}
    for x, y in zip(a, b):
        t[(x, y)] = t.get((x, y), 0) + 1
    ra, rb = {}, {}
    for (x, y), v in t.items():
        ra[x] = ra.get(x, 0) + v; rb[y] = rb.get(y, 0) + v
    s = sum(c2(v) for v in t.values()); sa = sum(c2(v) for v in ra.values()); sb = sum(c2(v) for v in rb.values())
    exp = sa * sb / c2(n); mx = (sa + sb) / 2
    return (s - exp) / (mx - exp) if mx != exp else 1.0


def kmeans(X, k, rng, restarts=50, iters=100):
    best = None
    for _ in range(restarts):
        cent = [list(X[i]) for i in rng.sample(range(len(X)), k)]
        lab = [0] * len(X)
        for _ in range(iters):
            new = [min(range(k), key=lambda c: sum((a - b) ** 2 for a, b in zip(x, cent[c]))) for x in X]
            if new == lab and _:
                break
            lab = new
            for c in range(k):
                mem = [X[i] for i in range(len(X)) if lab[i] == c]
                if mem:
                    cent[c] = [sum(v) / len(mem) for v in zip(*mem)]
        inertia = sum(sum((a - b) ** 2 for a, b in zip(X[i], cent[lab[i]])) for i in range(len(X)))
        if best is None or inertia < best[0]:
            best = (inertia, lab)
    return best[1]


def silhouette(X, lab):
    n = len(X); d = [[math.dist(X[i], X[j]) for j in range(n)] for i in range(n)]; s = []
    for i in range(n):
        own = [d[i][j] for j in range(n) if lab[j] == lab[i] and j != i]
        if not own:
            s.append(0.0); continue
        a = sum(own) / len(own)
        b = min(sum(d[i][j] for j in range(n) if lab[j] == c) / sum(1 for j in range(n) if lab[j] == c) for c in set(lab) if c != lab[i])
        s.append((b - a) / max(a, b) if max(a, b) > 0 else 0.0)
    return sum(s) / n


# ------------------------------------------------------------------------------------------ data
def read_matrix(p):
    with gzip.open(p, "rt") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        rows = {l.split("\t", 1)[0]: l.rstrip("\n").split("\t")[1:] for l in fh}
    return h[1:], rows


def load_rep(d):
    d = Path(d)
    cols, tpm = read_matrix(d / "tsaolab_melanoma_tpm_ensembl.tsv.gz")
    _, sym = read_matrix(d / "tsaolab_melanoma_tpm_symbol.tsv.gz")
    with gzip.open(d / "tsaolab_melanoma_genes.tsv.gz", "rt") as fh:
        gtype = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2] for l in list(fh)[1:]}
    return {"cols": cols, "tpm": tpm, "sym": sym, "gtype": gtype}


def load_qc(p, origin_override=None):
    with open(p, newline="") as fh:
        rows = {r["sample_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    return rows


def covariates(qc, cols, expansion_set):
    cov = {}
    for c in cols:
        r = qc[c]
        cov[c] = {"sequencing_run": r.get("sequencing_run"), "lane": r.get("lane"), "sex_inferred": r.get("sex_inferred"),
                  "quality_flag": r.get("quality_flag"), "specimen_set": "expansion9" if c in expansion_set else "original33",
                  "pipeline_origin": r.get("pipeline_origin")}
        for k in CONTINUOUS:
            v = r.get(k)
            cov[c][k] = float(v) if v not in (None, "", "nan", "None") else None
    return cov


def markers(rep, cols):
    out = {}
    for name, genes in MARKERS.items():
        present = [g for g in genes if g in rep["sym"]]
        out[name] = {c: (sum(math.log1p(float(rep["sym"][g][rep["cols"].index(c)])) for g in present) / len(present) if present else None) for c in cols}
        out[name + "_n_genes"] = len(present)
    return out


# ------------------------------------------------------------------------------------------ analysis
def analyse(name, rep, cols, cov, a, rng):
    ci = [rep["cols"].index(c) for c in cols]
    genes = sorted(g for g, v in rep["tpm"].items() if rep["gtype"].get(g) == "protein_coding"
                   and sum(float(v[i]) >= 1 for i in ci) >= 0.5 * len(ci))
    X = [[math.log1p(float(rep["tpm"][g][i])) for g in genes] for i in ci]           # samples x genes
    n = len(cols)
    G = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            G[i][j] = G[j][i] = sum(p * q for p, q in zip(X[i], X[j]))
    ev, sc, vec = pca_from_gram(G, list(range(n)))
    # top loadings (centred data projected on eigenvectors)
    gm = [sum(X[i][g] for i in range(n)) / n for g in range(len(genes))]
    top = {}
    for k in range(3):
        load = [(sum((X[i][g] - gm[g]) * vec[k][i] for i in range(n)), genes[g]) for g in range(len(genes))]
        top[f"PC{k+1}"] = [gid for _, gid in sorted(load, key=lambda t: -abs(t[0]))[:15]]
    # associations
    assoc = []
    for k in range(K):
        y = [sc[i][k] for i in range(n)]
        for cv in CATEGORICAL:
            g = [cov[c][cv] for c in cols]
            if len(set(g)) < 2 or None in g:
                continue
            assoc.append({"analysis": name, "pc": f"PC{k+1}", "covariate": cv, "type": "eta2", "stat": eta2(y, g),
                          "p_perm": perm_p(eta2, y, g, a.permutations, rng)})
        for cv in CONTINUOUS:
            pts = [(y[i], cov[c][cv]) for i, c in enumerate(cols) if cov[c][cv] is not None]
            if len(pts) < 6:
                continue
            yy, xx = [p for p, _ in pts], [q for _, q in pts]
            f = lambda u, v: abs(spearman(u, v))
            assoc.append({"analysis": name, "pc": f"PC{k+1}", "covariate": cv, "type": "abs_spearman", "stat": f(yy, xx),
                          "p_perm": perm_p(f, yy, xx, a.permutations, rng)})
    for r, q in zip(assoc, bh([r["p_perm"] for r in assoc])):
        r["q_bh"] = q
    # covariate-covariate relations
    rel = []
    allc = [c for c in CATEGORICAL if len({cov[s][c] for s in cols}) > 1 and None not in {cov[s][c] for s in cols}] + CONTINUOUS
    for i, u in enumerate(allc):
        for v in allc[i + 1:]:
            pts = [(cov[s][u], cov[s][v]) for s in cols if cov[s][u] is not None and cov[s][v] is not None]
            if len(pts) < 6:
                continue
            uu, vv = [p for p, _ in pts], [q for _, q in pts]
            if u in CATEGORICAL and v in CATEGORICAL:
                st, ty = cramers_v(uu, vv), "cramers_v"
            elif u in CATEGORICAL:
                st, ty = eta2(vv, uu), "eta2"
            elif v in CATEGORICAL:
                st, ty = eta2(uu, vv), "eta2"
            else:
                st, ty = spearman(uu, vv), "spearman"
            rel.append({"analysis": name, "a": u, "b": v, "type": ty, "stat": st, "n": len(pts)})
    # technical versus marker-based biological axes
    mk = markers(rep, cols)
    lanes = sorted({cov[c]["lane"] for c in cols})
    techbio = []
    for k in range(K):
        idx = [i for i, c in enumerate(cols) if all(cov[c][t] is not None for t in TECH) and all(mk[m][c] is not None for m in MARKERS)]
        y = [sc[i][k] for i in idx]
        tech = [[1.0] + [float(cov[cols[i]]["lane"] == l) for l in lanes[1:]] + [cov[cols[i]][t] for t in TECH] for i in idx]
        bio = [[1.0] + [mk[m][cols[i]] for m in MARKERS] for i in idx]
        both = [t + b[1:] for t, b in zip(tech, bio)]
        rt, at = ols_r2(y, tech); rb, ab = ols_r2(y, bio); rf, af = ols_r2(y, both)
        techbio.append({"analysis": name, "pc": f"PC{k+1}", "n": len(idx), "r2_technical": rt, "adj_r2_technical": at,
                        "r2_markers": rb, "adj_r2_markers": ab, "r2_both": rf, "unique_technical": rf - rb, "unique_markers": rf - rt,
                        "shared": rt + rb - rf})
    # stability
    def matched(sub_scores, sub_idx):
        res = []
        for k in range(3):
            full_k = [sc[i][k] for i in sub_idx]
            res.append(max(abs(pearson(full_k, [s[j] for s in sub_scores])) for j in range(K)))
        return res
    loo = []
    for drop in range(n):
        idx = [i for i in range(n) if i != drop]
        _, s2, _ = pca_from_gram(G, idx)
        loo.append(matched(s2, idx))
    subs, clus = [], {}
    labels_full = {kk: kmeans([s[:K] for s in sc], kk, random.Random(a.seed + kk)) for kk in (2, 3, 4)}
    for kk in (2, 3, 4):
        clus[kk] = {"silhouette": silhouette([s[:K] for s in sc], labels_full[kk]), "aris": []}
    srng = random.Random(a.seed)
    for _ in range(a.subsamples):
        idx = sorted(srng.sample(range(n), max(4, int(round(0.8 * n)))))
        _, s2, _ = pca_from_gram(G, idx)
        subs.append(matched(s2, idx))
        for kk in (2, 3, 4):
            lab = kmeans([s[:K] for s in s2], kk, random.Random(a.seed + 100 + kk), restarts=20)
            clus[kk]["aris"].append(ari(lab, [labels_full[kk][i] for i in idx]))
    q = lambda v, p: sorted(v)[max(0, min(len(v) - 1, int(p * (len(v) - 1))))]
    stab = {f"PC{k+1}": {"loo_min_abs_r": min(r[k] for r in loo), "loo_median_abs_r": q([r[k] for r in loo], 0.5),
                         "subsample_median_abs_r": q([r[k] for r in subs], 0.5), "subsample_p05_abs_r": q([r[k] for r in subs], 0.05)}
            for k in range(3)}
    cl = []
    for kk in (2, 3, 4):
        row = {"analysis": name, "k": kk, "silhouette": clus[kk]["silhouette"], "subsample_median_ari": q(clus[kk]["aris"], 0.5),
               "subsample_p05_ari": q(clus[kk]["aris"], 0.05), "sizes": json.dumps(sorted([labels_full[kk].count(c) for c in set(labels_full[kk])], reverse=True))}
        for cv in CATEGORICAL:
            g = [cov[c][cv] for c in cols]
            if len(set(g)) > 1 and None not in g:
                row[f"ari_vs_{cv}"] = ari(labels_full[kk], g)
        cl.append(row)
    return {"name": name, "cols": cols, "n_genes": len(genes), "explained_variance": ev, "scores": sc, "top_loadings": top,
            "assoc": assoc, "relations": rel, "techbio": techbio, "stability": stab, "clusters": cl, "labels": labels_full,
            "markers_n_genes": {m: mk[m + "_n_genes"] for m in MARKERS}}


def cross(A, B, label):
    shared = [c for c in A["cols"] if c in B["cols"]]
    ia, ib = [A["cols"].index(c) for c in shared], [B["cols"].index(c) for c in shared]
    r = {"comparison": label, "n_shared": len(shared)}
    for k in range(3):
        r[f"PC{k+1}_max_abs_r"] = max(abs(pearson([A["scores"][i][k] for i in ia], [B["scores"][j][m] for j in ib])) for m in range(K))
    for kk in (2, 3, 4):
        r[f"kmeans{kk}_ari"] = ari([A["labels"][kk][i] for i in ia], [B["labels"][kk][j] for j in ib])
    return r


def per_specimen_corr(r1, r2, cols):
    genes = sorted(set(r1["tpm"]) & set(r2["tpm"]))
    out = []
    for c in cols:
        i, j = r1["cols"].index(c), r2["cols"].index(c)
        x = [math.log1p(float(r1["tpm"][g][i])) for g in genes]; y = [math.log1p(float(r2["tpm"][g][j])) for g in genes]
        out.append({"specimen": c, "pearson_log1p_gene_tpm": pearson(x, y), "n_genes": len(genes)})
    return out


def write_new(path, text):
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    with os.fdopen(fd, "w") as fh:
        fh.write(text)


def tsv(path, rows):
    if not rows:
        write_new(path, ""); return
    cols = []
    for r in rows:
        cols += [k for k in r if k not in cols]
    f = lambda v: "" if v is None else (f"{v:.4g}" if isinstance(v, float) else str(v))
    write_new(path, "\n".join(["\t".join(cols)] + ["\t".join(f(r.get(c)) for c in cols) for r in rows]) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--original33", required=True); ap.add_argument("--mixed42", required=True); ap.add_argument("--harmonized42", required=True)
    ap.add_argument("--qc-mixed", required=True, help="expansion_qc_provenance_42.tsv (collaborator QC for 33, r3 for 9)")
    ap.add_argument("--qc-harmonized", required=True, help="harmonized_qc_provenance_42.tsv (all r2/r3)")
    ap.add_argument("--expansion", default="GM33,GM34,GM35,GM36,GM38,GM39,GM40,GM41,GM42")
    ap.add_argument("--suspect", default="GM10,GM13,GM15,GM21,GM40,GM41", help="pre-declared suspect specimens (fixed list)")
    ap.add_argument("--permutations", type=int, default=2000); ap.add_argument("--subsamples", type=int, default=100)
    ap.add_argument("--seed", type=int, default=20261127); ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    out = Path(a.out_dir)
    if os.path.lexists(out):
        fail(f"output exists: {out}")
    reps = {"original33": load_rep(a.original33), "mixed42": load_rep(a.mixed42), "harmonized42": load_rep(a.harmonized42)}
    qcm, qch = load_qc(a.qc_mixed), load_qc(a.qc_harmonized)
    exp = set(a.expansion.split(",")); sus = a.suspect.split(",")
    if reps["mixed42"]["cols"] != reps["harmonized42"]["cols"] or reps["mixed42"]["cols"][:33] != reps["original33"]["cols"]:
        fail("column orders differ between representations")
    all42 = reps["mixed42"]["cols"]
    missing = [s for s in sus if s not in all42]
    if missing:
        fail(f"suspect specimens not in the matrices: {missing}")
    sus_check = [{"specimen": s, "flag_mixed_qc": qcm[s]["quality_flag"], "mapping_mixed_qc": qcm[s]["mapping_rate_pct"],
                  "flag_harmonized_qc": qch[s]["quality_flag"], "mapping_harmonized_qc": qch[s]["mapping_rate_pct"]} for s in sus]
    flagged_m = sorted(s for s in all42 if qcm[s]["quality_flag"] == "suspect"); flagged_h = sorted(s for s in all42 if qch[s]["quality_flag"] == "suspect")
    rng = random.Random(a.seed)
    plan = [("original33", "original33", reps["original33"]["cols"], qcm), ("mixed42", "mixed42", all42, qcm), ("harmonized42", "harmonized42", all42, qch),
            ("harmonized33", "harmonized42", reps["original33"]["cols"], qch)]
    plan += [(n + "_minus_suspect", r, [c for c in cols if c not in sus], q) for n, r, cols, q in plan]
    res = {}
    for name, rname, cols, qc in plan:
        print(f"... {name}: {len(cols)} specimens", file=sys.stderr)
        res[name] = analyse(name, reps[rname], cols, covariates(qc, cols, exp), a, rng)
    crosses = [cross(res["original33"], res["harmonized33"], "original33 vs harmonized33 (same 33 specimens; pipeline only)"),
               cross(res["mixed42"], res["harmonized42"], "mixed42 vs harmonized42 (same 42 specimens; mixing effect)"),
               cross(res["original33_minus_suspect"], res["harmonized33_minus_suspect"], "original33 vs harmonized33, minus suspect"),
               cross(res["mixed42_minus_suspect"], res["harmonized42_minus_suspect"], "mixed42 vs harmonized42, minus suspect")]
    psc = per_specimen_corr(reps["original33"], reps["harmonized42"], reps["original33"]["cols"])
    out.mkdir(parents=True)
    tsv(out / "pc_covariate_associations.tsv", [r for v in res.values() for r in v["assoc"]])
    tsv(out / "covariate_relations.tsv", [r for v in res.values() for r in v["relations"]])
    tsv(out / "technical_vs_marker_r2.tsv", [r for v in res.values() for r in v["techbio"]])
    tsv(out / "clustering_stability.tsv", [r for v in res.values() for r in v["clusters"]])
    tsv(out / "pca_stability.tsv", [{"analysis": k, "pc": pc, **s} for k, v in res.items() for pc, s in v["stability"].items()])
    tsv(out / "cross_representation.tsv", crosses); tsv(out / "per_specimen_original_vs_harmonized.tsv", psc)
    tsv(out / "suspect_identity_check.tsv", sus_check)
    tsv(out / "pc_scores.tsv", [{"analysis": k, "specimen": c, **{f"PC{j+1}": v["scores"][i][j] for j in range(K)}} for k, v in res.items() for i, c in enumerate(v["cols"])])
    summary = {"representations": {k: {"n_specimens": len(v["cols"]), "n_genes": v["n_genes"], "explained_variance": [round(x, 4) for x in v["explained_variance"]],
                                       "top_loadings": v["top_loadings"], "markers_n_genes": v["markers_n_genes"]} for k, v in res.items()},
               "suspect_list": sus, "flagged_suspect_mixed_qc": flagged_m, "flagged_suspect_harmonized_qc": flagged_h, "cross_representation": crosses,
               "per_specimen_original_vs_harmonized_min": min(r["pearson_log1p_gene_tpm"] for r in psc),
               "caveats": ["exploratory; small n; associations generate hypotheses only", "markers may carry technical effects",
                           "no batch correction, feature removal or supervised modelling"],
               "run": {"utc": datetime.now(timezone.utc).isoformat(), "host": socket.gethostname(), "python": platform.python_version(),
                       "permutations": a.permutations, "subsamples": a.subsamples, "seed": a.seed}}
    write_new(out / "summary.json", json.dumps(summary, indent=2) + "\n")
    print("SUSPECT IDENTITIES (fixed list):", sus, "| flagged by mixed QC:", flagged_m, "| flagged by harmonized QC:", flagged_h)
    for s in sus_check:
        print(f"  {s['specimen']}: mixed {s['flag_mixed_qc']} ({s['mapping_mixed_qc']}%)  harmonized {s['flag_harmonized_qc']} ({s['mapping_harmonized_qc']}%)")
    print("EXPLAINED VARIANCE PC1-3:"); [print(f"  {k:<30} n={len(v['cols'])} genes={v['n_genes']} {[round(x, 3) for x in v['explained_variance'][:3]]}") for k, v in res.items()]
    print("STRONGEST PC-COVARIATE ASSOCIATIONS (q_bh < 0.10; exploratory):")
    for k, v in res.items():
        for r in sorted(v["assoc"], key=lambda r: r["q_bh"]):
            if r["q_bh"] < 0.10:
                print(f"  {k:<30} {r['pc']:<4} {r['covariate']:<24} {r['type']:<12} {r['stat']:.3f} p={r['p_perm']:.4f} q={r['q_bh']:.3f}")
    print("SPECIMEN SET (original33 vs expansion9) ASSOCIATION, mixed vs harmonized (pipeline is identical in harmonized):")
    for k in ("mixed42", "harmonized42", "mixed42_minus_suspect", "harmonized42_minus_suspect"):
        for r in res[k]["assoc"]:
            if r["covariate"] in ("specimen_set", "pipeline_origin") and r["pc"] in ("PC1", "PC2", "PC3"):
                print(f"  {k:<30} {r['pc']:<4} {r['covariate']:<16} eta2={r['stat']:.3f} p={r['p_perm']:.4f}")
    print("TECHNICAL vs MARKER R2 (PC1-3):")
    for k in ("original33", "harmonized42", "mixed42"):
        for r in res[k]["techbio"][:3]:
            print(f"  {k:<14} {r['pc']}: tech {r['r2_technical']:.2f} (unique {r['unique_technical']:.2f})  markers {r['r2_markers']:.2f} (unique {r['unique_markers']:.2f})")
    print("CROSS-REPRESENTATION:"); [print("  " + "  ".join(f"{kk}={vv:.3f}" if isinstance(vv, float) else f"{kk}={vv}" for kk, vv in c.items())) for c in crosses]
    print(f"per-specimen original vs harmonized gene log1p Pearson: min {summary['per_specimen_original_vs_harmonized_min']:.5f}")
    print("PCA STABILITY (PC1-3 subsample median |r|):"); [print(f"  {k:<30} " + " ".join(f"{pc}={s['subsample_median_abs_r']:.3f}" for pc, s in v["stability"].items())) for k, v in res.items()]
    print(f"outputs: {out}")


if __name__ == "__main__":
    main()
