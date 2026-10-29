"""Somatic hotspot and driver-gene variants from RNA-seq of FFPE melanoma, without a genome alignment.

  ref   write a small transcript reference of melanoma driver genes (MANE Select transcript of each), their
        pseudogenes and closest paralogs (so reads from those have somewhere better to go), plus a CDS table.
        run_quant.sh aligns every trimmed read pair to it with minimap2.
  call  joint bcftools mpileup over the per-sample BAMs, then per-site allele counts -> annotated variants:
        transcript position -> CDS codon -> amino-acid change; hotspot codon membership; cohort-wide flags
        (likely germline polymorphism; possible FFPE/RNA artefact classes).

RNA caveats: only expressed alleles are seen (a nonsense allele may be under-represented by NMD), promoter
mutations (TERT) are invisible, A>G can be ADAR editing, and low-fraction C>T / G>A in FFPE can be cytosine
deamination, so those are flagged rather than dropped.

Usage: python -I variants.py ref  <MANE summary.txt.gz> <gencode pc_transcripts.fa.gz> <gencode transcripts.fa.gz> <out_dir>
       python -I variants.py call <ref_dir> <sites.tsv from bcftools query> <sample_sheet.csv> <out_dir>
"""
import gzip
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Gene -> hotspot codons (MANE Select protein numbering). Empty tuple: tumour suppressor, any coding change is reported.
DRIVERS = {
    "BRAF": (466, 469, 581, 594, 597, 600, 601), "NRAS": (12, 13, 61), "KRAS": (12, 13, 61, 146), "HRAS": (12, 13, 61),
    "KIT": (502, 503, 557, 558, 559, 560, 576, 642, 816, 822), "NF1": (), "PTEN": (), "TP53": (), "CDKN2A": (),
    "CDK4": (24,), "MAP2K1": (56, 57, 121, 124), "MAP2K2": (), "RAC1": (29,), "IDH1": (132,), "IDH2": (140, 172),
    "CTNNB1": (33, 37, 41, 45), "SF3B1": (625, 666), "PPP6C": (264,), "ARID2": (), "ARID1A": (), "RB1": (),
    "DDX3X": (), "KNSTRN": (24,), "STK19": (89,), "MITF": (318,), "GNAQ": (183, 209), "GNA11": (183, 209),
    "CYSLTR2": (129,), "PLCB4": (630,), "BAP1": (), "EIF1AX": tuple(range(1, 16)), "PIK3CA": (542, 545, 1047),
    "AKT1": (17,), "PREX2": (), "RASA2": (), "PTPN11": (), "FBXW7": (), "SMARCA4": (), "SETD2": (), "KMT2C": (),
    "KMT2D": (), "ERBB2": (), "EGFR": (), "MTOR": (), "TERT": (), "RAF1": (), "MAP3K1": (), "IDH1": (132,),
}
# Paralogs / family members close enough in sequence to steal or donate reads. Listed only to absorb reads.
DECOYS = ["ARAF", "RRAS", "RRAS2", "MRAS", "RALA", "RALB", "RAP1A", "RAP1B", "RAP2A", "CDK6", "AKT2", "AKT3",
          "RAC2", "RAC3", "RHOG", "GNA14", "GNA15", "GNAI1", "GNAI2", "GNAS", "ERBB3", "ERBB4", "PIK3CB",
          "PIK3CD", "ARID1B", "DDX3Y", "SMARCA2", "CSF1R", "PDGFRA", "FLT3", "RAF1", "MAP2K2", "IDH3A"]
PSEUDO = re.compile(r"^(" + "|".join(DRIVERS) + r")P\d*$")  # e.g. PTENP1, BRAFP1, KRASP1
MIN_DEPTH = 10
MIN_ALT_READS = 3
MIN_VAF_HOTSPOT = 0.02
MIN_VAF = 0.05
GERMLINE_FRAC = 0.25  # present at >= 30% VAF in this share of samples -> likely polymorphism

CODON = {}
_B = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for _i, _a in enumerate(_AA):
    CODON[_B[_i // 16] + _B[_i // 4 % 4] + _B[_i % 4]] = _a


def read_fasta(path):
    name, seqs = None, {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line[0] == ">":
                name = line[1:].rstrip()
                seqs[name] = []
            else:
                seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def make_ref(mane_path, pc_path, tx_path, out):
    out.mkdir(parents=True, exist_ok=True)
    mane = pd.read_csv(mane_path, sep="\t")
    mane = mane[mane.MANE_status == "MANE Select"]
    mane_enst = dict(zip(mane.Ensembl_nuc.str.split(".").str[0], mane.symbol))
    want = set(DRIVERS) | set(DECOYS)
    rows, fa = [], []
    for hdr, seq in read_fasta(pc_path).items():
        f = hdr.split("|")
        enst, gene = f[0].split(".")[0], f[5]
        if gene in want and mane_enst.get(enst) == gene:
            cds = re.search(r"CDS:(\d+)-(\d+)", hdr)
            if not cds:
                continue
            rows.append(dict(name=f"{gene}|{f[0]}", gene=gene, transcript=f[0], role="driver" if gene in DRIVERS else "decoy",
                             cds_start=int(cds[1]), cds_end=int(cds[2]), length=len(seq)))
            fa.append((f"{gene}|{f[0]}", seq))
    found = {r["gene"] for r in rows}
    missing = want - found
    if missing:
        print("no MANE Select protein-coding transcript for:", " ".join(sorted(missing)))
    for hdr, seq in read_fasta(tx_path).items():  # pseudogenes of the driver genes, all transcripts
        f = hdr.split("|")
        if PSEUDO.match(f[5]):
            rows.append(dict(name=f"{f[5]}|{f[0]}", gene=f[5], transcript=f[0], role="pseudogene",
                             cds_start=0, cds_end=0, length=len(seq)))
            fa.append((f"{f[5]}|{f[0]}", seq))
    with open(out / "driver_ref.fa", "w") as fh:
        for n, s in fa:
            fh.write(f">{n}\n{s}\n")
    pd.DataFrame(rows).to_csv(out / "driver_ref.tsv", sep="\t", index=False)
    print(f"{len(fa)} sequences ({sum(len(s) for _, s in fa) / 1e6:.2f} Mb) -> {out / 'driver_ref.fa'}")


def annotate(seq, cds_start, cds_end, pos, ref, alt):
    """Protein-level consequence of ref->alt at 1-based transcript position pos. Returns (codon_no, hgvsp, consequence)."""
    if cds_start == 0:
        return np.nan, "", "pseudogene"
    if len(ref) == len(alt) == 1:
        if pos < cds_start:
            return np.nan, "", "5'UTR"
        if pos > cds_end:
            return np.nan, "", "3'UTR"
        off = pos - cds_start
        codon_no, frame = off // 3 + 1, off % 3
        start = cds_start - 1 + (codon_no - 1) * 3
        codon = seq[start:start + 3]
        if len(codon) < 3:
            return codon_no, "", "incomplete_codon"
        assert codon[frame] == ref, (pos, codon, ref)
        mut = codon[:frame] + alt + codon[frame + 1:]
        a, b = CODON.get(codon, "X"), CODON.get(mut, "X")
        cons = "synonymous" if a == b else "stop_gained" if b == "*" else "stop_lost" if a == "*" else "missense"
        return codon_no, f"p.{a}{codon_no}{b}", cons
    # indel: bcftools gives REF/ALT with a shared first base
    first = pos + 1
    if first < cds_start or pos > cds_end:
        return np.nan, "", "UTR_indel"
    d = len(alt) - len(ref)
    codon_no = (first - cds_start) // 3 + 1
    kind = "inframe" if d % 3 == 0 else "frameshift"
    return codon_no, f"p.{codon_no}{'ins' if d > 0 else 'del'}{abs(d)}{'fs' if kind == 'frameshift' else ''}", f"{kind}_{'insertion' if d > 0 else 'deletion'}"


def call(ref_dir, sites_path, sheet, out):
    out.mkdir(parents=True, exist_ok=True)
    refs = pd.read_csv(ref_dir / "driver_ref.tsv", sep="\t").set_index("name")
    seqs = {}
    for n, s in read_fasta_plain(ref_dir / "driver_ref.fa").items():
        seqs[n] = s
    samples = pd.read_csv(sheet, dtype=str)
    core2sid = dict(zip(samples.core_id, samples.sample_id))

    # sites.tsv: CHROM POS REF ALT INDEL then one AD column per sample (header from bcftools query -H)
    df = pd.read_csv(sites_path, sep="\t", dtype=str)
    df.columns = [re.sub(r"^\[\d+\]", "", c).replace("# ", "").replace(":AD", "") for c in df.columns]
    snames = [c for c in df.columns if c not in ("CHROM", "POS", "REF", "ALT", "INDEL")]
    snames_out = [core2sid.get(s, s) for s in snames]
    df["POS"] = df.POS.astype(int)
    recs = []
    for r in df.itertuples(index=False):
        alts = [a for a in r.ALT.split(",") if a != "<*>"]
        if not alts:
            continue
        ad = np.array([[float(x) if x != "." else 0.0 for x in getattr(r, s).split(",")] for s in snames])  # samples x alleles
        depth = ad.sum(axis=1)
        for j, alt in enumerate(alts, start=1):
            altc = ad[:, j]
            vaf = np.divide(altc, depth, out=np.zeros_like(altc), where=depth > 0)
            if (altc >= MIN_ALT_READS).sum() == 0:
                continue
            recs.append((r.CHROM, r.POS, r.REF, alt, r.INDEL == "1", depth, altc, vaf))
    rows = []
    for chrom, pos, ref, alt, is_indel, depth, altc, vaf in recs:
        gene, tx = chrom.split("|")
        info = refs.loc[chrom]
        codon_no, hgvsp, cons = annotate(seqs[chrom], info.cds_start, info.cds_end, pos, ref, alt)
        hot = DRIVERS.get(gene)
        is_hot = bool(hot) and codon_no in hot and cons not in ("synonymous",)
        ts = f"{ref}>{alt}" if not is_indel else "indel"
        artefact = ("possible_RNA_editing" if ts == "A>G" else
                    "possible_deamination" if ts in ("C>T", "G>A") else "")
        n_called = int(((depth >= MIN_DEPTH) & (altc >= MIN_ALT_READS) & (vaf >= (MIN_VAF_HOTSPOT if is_hot else MIN_VAF))).sum())
        n_high = int(((depth >= MIN_DEPTH) & (vaf >= 0.3)).sum())
        n_cov = int((depth >= MIN_DEPTH).sum())
        germline = n_cov > 0 and n_high / max(n_cov, 1) >= GERMLINE_FRAC and n_high >= 3
        if n_called == 0:
            continue
        rows.append(dict(
            gene=gene, transcript=tx, role=info.role, tx_pos=pos, ref=ref, alt=alt, change=ts,
            codon=codon_no, hgvsp=hgvsp, consequence=cons, hotspot=is_hot,
            n_samples_called=n_called, n_samples_covered=n_cov,
            likely_polymorphism=germline, artefact_class=artefact,
            **{f"{s}_vaf": round(v, 3) for s, v in zip(snames_out, vaf)},
            **{f"{s}_alt": int(a) for s, a in zip(snames_out, altc)},
            **{f"{s}_depth": int(d) for s, d in zip(snames_out, depth)},
        ))
    allv = pd.DataFrame(rows)
    if allv.empty:
        print("no variants")
        return
    allv.to_csv(out / "driver_variants_all.tsv", sep="\t", index=False)

    # long table: one row per sample x variant actually called in that sample
    long = []
    for r in allv.itertuples(index=False):
        d = r._asdict()
        for s in snames_out:
            vaf, alt_n, dep = d[f"{s}_vaf"], d[f"{s}_alt"], d[f"{s}_depth"]
            thr = MIN_VAF_HOTSPOT if r.hotspot else MIN_VAF
            if dep >= MIN_DEPTH and alt_n >= MIN_ALT_READS and vaf >= thr:
                long.append(dict(sample_id=s, gene=r.gene, hgvsp=r.hgvsp, consequence=r.consequence, change=r.change,
                                 tx_pos=r.tx_pos, ref=r.ref, alt=r.alt, vaf=vaf, alt_reads=alt_n, depth=dep,
                                 hotspot=r.hotspot, likely_polymorphism=r.likely_polymorphism, artefact_class=r.artefact_class))
    long = pd.DataFrame(long)
    long.to_csv(out / "driver_variants_by_sample.tsv", sep="\t", index=False)

    # somatic candidates: coding, non-synonymous, driver gene, not a cohort-wide polymorphism
    som = long[(long.consequence.isin(["missense", "stop_gained", "stop_lost"]) | long.consequence.str.contains("frameshift|inframe"))
               & ~long.likely_polymorphism & long.gene.isin(DRIVERS)]
    som.to_csv(out / "driver_variants_somatic_candidates.tsv", sep="\t", index=False)

    # per-sample summary of hotspot status + coverage at the key hotspots
    key = [("BRAF", 600), ("NRAS", 61), ("NRAS", 12), ("NRAS", 13), ("KIT", 576), ("KIT", 642), ("KIT", 816),
           ("GNAQ", 209), ("GNA11", 209), ("MAP2K1", 124), ("RAC1", 29), ("CTNNB1", 45), ("IDH1", 132), ("SF3B1", 625)]
    summ = []
    cov = depth_at_codons(ref_dir, sites_path, snames, key)
    for i, s in enumerate(snames_out):
        mine = som[som.sample_id == s]
        hs = mine[mine.hotspot]
        summ.append(dict(
            sample_id=s,
            hotspot_calls="; ".join(f"{r.gene} {r.hgvsp} ({r.vaf:.0%}, {r.alt_reads}/{r.depth})" for r in hs.itertuples()) or "none",
            other_driver_calls="; ".join(f"{r.gene} {r.hgvsp} ({r.vaf:.0%})" for r in mine[~mine.hotspot].itertuples()) or "none",
            n_somatic_candidates=len(mine),
            **{f"depth_{g}_{c}": cov[(g, c)][i] for g, c in key},
        ))
    summ = pd.DataFrame(summ)
    for g, c in (("BRAF", 600), ("NRAS", 61)):
        col = f"depth_{g}_{c}"
        summ[f"{g}_{c}_status"] = np.where(summ.hotspot_calls.str.contains(f"{g} p.[A-Z*]{c}"), "mutant",
                                 np.where(summ[col] >= 30, "wild-type", np.where(summ[col] >= MIN_DEPTH, "wild-type (low depth)", "not covered")))
    summ.to_csv(out / "hotspot_summary.tsv", sep="\t", index=False)
    print(f"{len(allv)} variant sites, {len(som)} somatic-candidate calls in {som.sample_id.nunique()} samples")
    print(summ[["sample_id", "BRAF_600_status", "NRAS_61_status", "hotspot_calls"]].to_string(index=False))


def read_fasta_plain(path):
    name, seqs = None, {}
    for line in Path(path).read_text().splitlines():
        if line[:1] == ">":
            name = line[1:]
            seqs[name] = []
        else:
            seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def depth_at_codons(ref_dir, sites_path, snames, key):
    """Max depth per sample over the three bases of each hotspot codon (0 if the site has no pileup line)."""
    refs = pd.read_csv(ref_dir / "driver_ref.tsv", sep="\t")
    refs = refs[refs.role == "driver"].set_index("gene")
    want = {}
    for g, c in key:
        info = refs.loc[g]
        start = info.cds_start + (c - 1) * 3
        for p in (start, start + 1, start + 2):
            want[(info.name, p)] = (g, c)
    out = {k: np.zeros(len(snames), dtype=int) for k in key}
    with open(sites_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            k = (f[0], int(f[1]))
            if k in want:
                dep = np.array([sum(float(x) for x in c.split(",") if x != ".") for c in f[5:]])
                out[want[k]] = np.maximum(out[want[k]], dep.astype(int))
    return out


if __name__ == "__main__":
    cmd, a = sys.argv[1], [Path(x) for x in sys.argv[2:]]
    if cmd == "ref":
        make_ref(*a)
    elif cmd == "call":
        call(*a)
    else:
        sys.exit(__doc__)
