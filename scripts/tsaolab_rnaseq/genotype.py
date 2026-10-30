"""Cell-line genotype check from RNA-seq SNPs (with genotype_check.sh).

  ref      pick one transcript per gene that is well expressed in every given salmon run and write a small
           FASTA to align to (no paralogs or rRNA, so every read has one place to go)
  compare  turn the joint SNP calls into pairwise genotype discordance and group samples into donors

Sites whose allele fraction is ambiguous (allelic imbalance, admixture) are left uncalled. Replicates of one culture then
disagree at ~0-2% of sites, the same line sequenced by another lab at ~5-8%, a sample with ~10-20% admixture from another
donor or species at ~1-17%, and unrelated donors at ~60-80%.

Usage: python -I genotype.py ref <transcripts.fa.gz> <out.fa> <salmon_quant_dir> [...]
       python -I genotype.py compare <snps.tsv> <sample_sheet.csv> [...] --out <dir>
"""
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

N_TRANSCRIPTS = 1500
MIN_TPM = 20.0
MIN_DEPTH = 20  # reads at a site, per sample, to call its genotype
HOM_MAX_AF = 0.05  # alternate-allele fraction up to this is 0/0 (and from 1 - this, 1/1)
HET_AF = (0.25, 0.75)  # 0/1; fractions between these bands are left uncalled
SAME_DONOR_PCT = 35.0  # discordance cut-off for grouping samples into one donor


def make_ref(fasta, out, quant_dirs):
    tpm = pd.DataFrame({d: pd.read_csv(Path(d) / "quant.sf", sep="\t", index_col=0).TPM for d in quant_dirs})
    seqs, gene, name = {}, {}, None
    with gzip.open(fasta, "rt") as fh:
        for line in fh:
            if line[0] == ">":
                f = line[1:].split("|")
                name, gene[f[0]], seqs[f[0]] = f[0], f[5], []
            else:
                seqs[name].append(line.strip())
    low = tpm.min(axis=1)
    low = low[low > MIN_TPM].to_frame("tpm")
    low["gene"] = low.index.map(gene)
    low = low[~low.gene.str.match(r"^(MT-|RPL|RPS|HLA-)")]
    best = low.sort_values("tpm", ascending=False).groupby("gene").head(1).index[:N_TRANSCRIPTS]
    with open(out, "w") as fh:
        for t in best:
            fh.write(f">{t}_{gene[t]}\n{''.join(seqs[t])}\n")
    print(f"{len(best)} transcripts -> {out}")


def genotypes(snps):
    """0/1/2 alternate-allele calls from allele depths (AD); NaN below MIN_DEPTH or at ambiguous fractions."""
    df = pd.read_csv(snps, sep="\t")
    df.columns = [c.split("]")[-1].split(":")[0] for c in df.columns]
    calls = {}
    for s in df.columns[2:]:
        ad = df[s].str.split(":").str[1].str.split(",", expand=True).replace(".", np.nan).astype(float)
        dp = ad.sum(axis=1)
        alt = (1 - ad[0] / dp).values
        het = (alt >= HET_AF[0]) & (alt <= HET_AF[1])
        calls[s] = np.select([alt <= HOM_MAX_AF, het, alt >= 1 - HOM_MAX_AF], [0, 1, 2], np.nan)
        calls[s][(dp < MIN_DEPTH).values] = np.nan
    return pd.DataFrame(calls)


def compare(snps, sheets, out):
    gt = genotypes(snps)
    names = list(gt.columns)
    disc = pd.DataFrame(0.0, index=names, columns=names)
    sites = pd.DataFrame(0, index=names, columns=names)
    for a in names:
        for b in names:
            ok = gt[a].notna() & gt[b].notna() & ((gt[a] > 0) | (gt[b] > 0))
            sites.loc[a, b] = ok.sum()
            disc.loc[a, b] = (gt.loc[ok, a] != gt.loc[ok, b]).mean() * 100
    donor = fcluster(linkage(squareform(disc.values, checks=False), "average"), SAME_DONOR_PCT, "distance")

    labels = pd.concat([pd.read_csv(s, dtype=str) for s in sheets]).set_index("sample_id")
    match = pd.DataFrame({"sample_id": names, "donor_group": donor})
    for col in ("cell_line", "submitted_as"):  # tissue sample sheets (tsaolab_melanoma) have neither column
        match[col] = match.sample_id.map(labels[col]).fillna("reference") if col in labels else "n/a"
    members = match.groupby("donor_group").sample_id.apply(lambda s: " ".join(s))
    match["donor_group_members"] = match.donor_group.map(members)
    for col, fn in (("within_group_discordance_pct", np.median), ("nearest_other_group_pct", np.min)):
        vals = []
        for r in match.itertuples():
            same = match.sample_id[(match.donor_group == r.donor_group) & (match.sample_id != r.sample_id)]
            other = match.sample_id[match.donor_group != r.donor_group]
            pool = same if col.startswith("within") else other
            vals.append(round(fn(disc.loc[r.sample_id, pool]), 1) if len(pool) else np.nan)
        match[col] = vals
    match["sites_called"] = gt.notna().sum().values

    out.mkdir(parents=True, exist_ok=True)
    disc.round(1).to_csv(out / "genotype_discordance.tsv", sep="\t")
    match.to_csv(out / "genotype_match.tsv", sep="\t", index=False)
    print(f"{len(gt)} SNPs, {match.donor_group.nunique()} donor groups")
    print(match.drop(columns="donor_group_members").to_string(index=False))
    for g, m in members.items():
        print(f"group {g}: {m}")


if __name__ == "__main__":
    cmd, args = sys.argv[1], sys.argv[2:]
    if cmd == "ref":
        make_ref(args[0], args[1], args[2:])
    elif cmd == "compare":
        i = args.index("--out")
        compare(args[0], args[1:i], Path(args[i + 1]))
    else:
        sys.exit(__doc__)
