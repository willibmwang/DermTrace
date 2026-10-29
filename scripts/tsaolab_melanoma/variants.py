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
       python -I variants.py call <ref_dir> <sites_dedup.tsv> <sites_raw.tsv> <sample_sheet.csv> <out_dir>   (via call_variants.sh)
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
MIN_DEPTH = 8  # unique molecules at a site, per sample
MIN_ALT_READS = 5  # duplicate-containing reads
MIN_ALT_MOLECULES = 2  # after duplicate removal
WT_MIN_MOLECULES = 20  # unique molecules over a hotspot codon to call it wild-type with some confidence
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


def load_sites(path, core2sid):
    """sites.tsv from bcftools query -H -> dict (chrom, pos, ref, alt) -> (depth, alt_reads) arrays, plus sample names."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    df.columns = [re.sub(r"^\[\d+\]", "", c.lstrip("# ")).replace(":AD", "") for c in df.columns]
    snames = [core2sid.get(c, c) for c in df.columns if c not in ("CHROM", "POS", "REF", "ALT", "INDEL")]
    sites, cover = {}, {}
    for r in df.itertuples(index=False):
        ad = np.array([[float(x) if x != "." else 0.0 for x in v.split(",")] for v in r[5:]])
        depth = ad.sum(axis=1)
        cover[(r.CHROM, int(r.POS))] = np.maximum(cover.get((r.CHROM, int(r.POS)), 0), depth)
        for j, alt in enumerate(a for a in r.ALT.split(",") if a != "<*>"):
            sites[(r.CHROM, int(r.POS), r.REF, alt)] = (depth, ad[:, j + 1])
    return snames, sites, cover


def call(ref_dir, dedup_sites, raw_sites, sheet, out):
    out.mkdir(parents=True, exist_ok=True)
    refs = pd.read_csv(ref_dir / "driver_ref.tsv", sep="\t").set_index("name")
    seqs = read_fasta_plain(ref_dir / "driver_ref.fa")
    samples = pd.read_csv(sheet, dtype=str)
    core2sid = dict(zip(samples.core_id, samples.sample_id))
    snames, uniq, ucov = load_sites(dedup_sites, core2sid)
    snames_raw, raw, rcov = load_sites(raw_sites, core2sid)
    assert snames == snames_raw, (snames, snames_raw)
    n = len(snames)

    rows = []
    for (chrom, pos, ref, alt), (udep, ualt) in uniq.items():
        if (ualt >= MIN_ALT_MOLECULES).sum() == 0:
            continue
        rdep, ralt = raw.get((chrom, pos, ref, alt), (np.zeros(n), np.zeros(n)))
        gene, tx = chrom.split("|")
        info = refs.loc[chrom]
        is_indel = len(ref) != len(alt)
        codon_no, hgvsp, cons = annotate(seqs[chrom], info.cds_start, info.cds_end, pos, ref, alt)
        hot = DRIVERS.get(gene)
        is_hot = bool(hot) and codon_no in hot and cons != "synonymous"
        ts = "indel" if is_indel else f"{ref}>{alt}"
        artefact = ("possible_RNA_editing" if ts == "A>G" else "possible_deamination" if ts in ("C>T", "G>A") else "")
        if is_indel and in_repeat(seqs[chrom], pos, ref, alt):
            artefact = "repeat_context_indel"
        uvaf = np.divide(ualt, udep, out=np.zeros(n), where=udep > 0)
        rvaf = np.divide(ralt, rdep, out=np.zeros(n), where=rdep > 0)
        called = (ualt >= MIN_ALT_MOLECULES) & (ralt >= MIN_ALT_READS) & (rvaf >= (MIN_VAF_HOTSPOT if is_hot else MIN_VAF)) & (udep >= MIN_DEPTH)
        n_cov = int((udep >= MIN_DEPTH).sum())
        n_high = int(((udep >= MIN_DEPTH) & (uvaf >= 0.3)).sum())
        germline = n_high >= 3 and n_high / max(n_cov, 1) >= GERMLINE_FRAC
        if called.sum() == 0:
            continue
        rows.append(dict(
            gene=gene, transcript=tx, role=info.role, tx_pos=pos, ref=ref, alt=alt, change=ts, codon=codon_no,
            hgvsp=hgvsp, consequence=cons, hotspot=is_hot, n_samples_called=int(called.sum()), n_samples_covered=n_cov,
            likely_polymorphism=germline, artefact_class=artefact,
            **{f"{s}_vaf": round(v, 3) for s, v in zip(snames, uvaf)},
            **{f"{s}_alt": int(a) for s, a in zip(snames, ualt)},
            **{f"{s}_depth": int(d) for s, d in zip(snames, udep)},
            **{f"{s}_raw_alt": int(a) for s, a in zip(snames, ralt)},
            **{f"{s}_raw_depth": int(d) for s, d in zip(snames, rdep)},
            called_mask=called,
        ))
    if not rows:
        print("no variants")
        return
    allv = pd.DataFrame(rows)
    allv.drop(columns="called_mask").to_csv(out / "driver_variants_all.tsv", sep="\t", index=False)

    long = []
    for r in allv.itertuples(index=False):
        d = r._asdict()
        for i, s in enumerate(snames):
            if d["called_mask"][i]:
                long.append(dict(sample_id=s, gene=r.gene, hgvsp=r.hgvsp, consequence=r.consequence, change=r.change,
                                 tx_pos=r.tx_pos, ref=r.ref, alt=r.alt, vaf=d[f"{s}_vaf"], alt_molecules=d[f"{s}_alt"],
                                 depth_molecules=d[f"{s}_depth"], raw_alt_reads=d[f"{s}_raw_alt"], raw_depth=d[f"{s}_raw_depth"],
                                 hotspot=r.hotspot, likely_polymorphism=r.likely_polymorphism, artefact_class=r.artefact_class))
    long = pd.DataFrame(long)
    long.to_csv(out / "driver_variants_by_sample.tsv", sep="\t", index=False)
    coding = long.consequence.isin(["missense", "stop_gained", "stop_lost"]) | long.consequence.str.contains("frameshift|inframe")
    som = long[coding & ~long.likely_polymorphism & long.gene.isin(DRIVERS)]
    som.to_csv(out / "driver_variants_somatic_candidates.tsv", sep="\t", index=False)

    # per-sample hotspot status with unique-molecule and raw depth over each key codon
    key = [("BRAF", 600), ("BRAF", 601), ("NRAS", 61), ("NRAS", 12), ("NRAS", 13), ("KRAS", 12), ("HRAS", 61), ("KIT", 576),
           ("KIT", 642), ("KIT", 816), ("GNAQ", 209), ("GNA11", 209), ("MAP2K1", 124), ("RAC1", 29), ("CTNNB1", 45),
           ("IDH1", 132), ("SF3B1", 625), ("PPP6C", 264), ("TP53", 273), ("CDKN2A", 80)]
    drv = refs[refs.role == "driver"].reset_index().set_index("gene")
    summ = []
    for i, s in enumerate(snames):
        mine = som[som.sample_id == s]
        hs = mine[mine.hotspot]
        row = dict(sample_id=s,
                   hotspot_calls="; ".join(f"{r.gene} {r.hgvsp} ({r.vaf:.0%}, {r.alt_molecules}/{r.depth_molecules} molecules, {r.raw_alt_reads}/{r.raw_depth} reads)" for r in hs.itertuples()) or "none",
                   other_driver_calls="; ".join(f"{r.gene} {r.hgvsp} ({r.vaf:.0%}, {r.alt_molecules}/{r.depth_molecules})" for r in mine[~mine.hotspot].itertuples()) or "none",
                   n_somatic_candidates=len(mine))
        for g, c in key:
            info = drv.loc[g]
            ps = [(info["name"], info.cds_start + (c - 1) * 3 + k) for k in range(3)]
            row[f"{g}_{c}_molecules"] = int(max(ucov.get(p, np.zeros(n))[i] if np.ndim(ucov.get(p, 0)) else 0 for p in ps))
            row[f"{g}_{c}_reads"] = int(max(rcov.get(p, np.zeros(n))[i] if np.ndim(rcov.get(p, 0)) else 0 for p in ps))
        summ.append(row)
    summ = pd.DataFrame(summ)
    for g, c in (("BRAF", 600), ("NRAS", 61), ("KIT", 576), ("GNAQ", 209), ("GNA11", 209)):
        m, rd = summ[f"{g}_{c}_molecules"], summ[f"{g}_{c}_reads"]
        mutant = summ.hotspot_calls.str.contains(rf"{g} p\.[A-Z*]{c}[A-Z*]")
        summ[f"{g}_{c}_status"] = np.select([mutant, m >= WT_MIN_MOLECULES, m >= MIN_DEPTH], ["mutant", "wild-type", "wild-type (low depth)"], "not covered")
    summ.to_csv(out / "hotspot_summary.tsv", sep="\t", index=False)
    print(f"{len(allv)} variant sites, {len(som)} somatic-candidate calls in {som.sample_id.nunique()} samples")
    print(summ[["sample_id", "BRAF_600_status", "BRAF_600_molecules", "NRAS_61_status", "NRAS_61_molecules", "hotspot_calls"]].to_string(index=False))


def in_repeat(seq, pos, ref, alt):
    """True when an indel's inserted/deleted unit is one copy of a tandem repeat (>= 3 units at the site), where
    alignment of short, damaged reads slips."""
    unit = (alt[1:] if len(alt) > len(ref) else ref[1:])
    if not unit or len(unit) > 12:
        return False
    k = 0
    i = pos  # 0-based index right after the shared anchor base
    while seq[i + k * len(unit): i + (k + 1) * len(unit)] == unit:
        k += 1
    return k >= 3


def read_fasta_plain(path):
    name, seqs = None, {}
    for line in Path(path).read_text().splitlines():
        if line[:1] == ">":
            name = line[1:]
            seqs[name] = []
        else:
            seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


if __name__ == "__main__":
    cmd, a = sys.argv[1], [Path(x) for x in sys.argv[2:]]
    if cmd == "ref":
        make_ref(*a)
    elif cmd == "call":
        call(*a)
    else:
        sys.exit(__doc__)
