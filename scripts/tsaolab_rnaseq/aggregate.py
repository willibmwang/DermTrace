"""Summarize salmon transcript quantifications to genes and write, per species:

  <project>_<species>_harmonizable-USE.tsv   ensembl_gene_id x sample, log2(TPM + 1)
                                            (same layout as Korean_Keloid_NHF_harmonizable-USE.tsv)
  <project>_<species>_counts.tsv            ensembl_gene_id x sample, salmon estimated counts
  <project>_<species>.meta-USE.csv          one row per sample, incl. QC
  <project>_<species>_genes.tsv             ensembl_gene_id -> gene_symbol, gene_type

Usage: python -I aggregate.py <sample_sheet.csv> <quant_out_dir> <dest_dir>
"""
import gzip
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REF = Path(os.environ["REF_DIR"])  # folder holding gencode_human_v46/ and gencode_mouse_vM35/
FASTA = {
    "human": REF / "gencode_human_v46" / "gencode.v46.transcripts.fa.gz",
    "mouse": REF / "gencode_mouse_vM35" / "gencode.vM35.transcripts.fa.gz",
}
MIN_MAPPING_RATE = 50.0  # percent; below this a sample is flagged "suspect"


def tx2gene(fasta):
    rows = []
    with gzip.open(fasta, "rt") as fh:
        for line in fh:
            if line[0] == ">":
                f = line[1:].rstrip().split("|")
                rows.append((f[0], f[1], f[5], f[7]))
    df = pd.DataFrame(rows, columns=["tx", "gene_version", "gene_symbol", "gene_type"])
    # Strip versions; PAR_Y copies collapse onto the X-linked gene ID.
    df["ensembl_gene_id"] = df.gene_version.str.split(".").str[0]
    return df


def main(sheet, qdir, dest):
    samples = pd.read_csv(sheet, dtype=str)
    status = pd.read_csv(qdir / "sample_status.tsv", sep="\t", dtype=str).set_index("sample_id")
    project = samples.project.iloc[0]
    dest.mkdir(parents=True, exist_ok=True)

    for species, sdf in samples.groupby("species"):
        t2g = tx2gene(FASTA[species])
        tpm, counts, meta = {}, {}, []
        for r in sdf.itertuples():
            qf = qdir / "quant" / r.sample_id / "quant.sf"
            if not qf.exists():
                print(f"[{r.sample_id}] not quantified: {status.detail.get(r.sample_id, '')}")
                continue
            q = pd.read_csv(qf, sep="\t", index_col=0)
            gene = t2g.set_index("tx").loc[q.index, "ensembl_gene_id"].values
            tpm[r.sample_id] = q.TPM.groupby(gene).sum()
            counts[r.sample_id] = q.NumReads.groupby(gene).sum()
            info = json.loads((qdir / "quant" / r.sample_id / "aux_info" / "meta_info.json").read_text())
            lib = json.loads((qdir / "quant" / r.sample_id / "lib_format_counts.json").read_text())
            rate = info["percent_mapped"]
            meta.append(dict(
                **r._asdict(),
                num_fragments=info["num_processed"],
                num_mapped=info["num_mapped"],
                mapping_rate_pct=round(rate, 2),
                library_type=lib["expected_format"],
                quality_flag="ok" if rate >= MIN_MAPPING_RATE else "suspect",
                notes="" if rate >= MIN_MAPPING_RATE else f"low mapping rate ({rate:.1f}%)",
            ))
        if not meta:
            print(f"{project}/{species}: no quantified samples yet")
            continue

        tpm = pd.DataFrame(tpm).fillna(0)
        counts = pd.DataFrame(counts).fillna(0)
        keep = counts.sum(axis=1) > 0  # drop genes with no reads in any sample
        prefix = dest / f"{project}_{species}"
        logtpm = np.log2(tpm[keep] + 1).rename_axis("ensembl_gene_id")
        logtpm.to_csv(f"{prefix}_harmonizable-USE.tsv", sep="\t")
        counts[keep].round(3).rename_axis("ensembl_gene_id").to_csv(f"{prefix}_counts.tsv", sep="\t")
        meta = pd.DataFrame(meta).drop(columns="Index")
        meta.to_csv(f"{prefix}.meta-USE.csv", index=False)
        (t2g.drop_duplicates("ensembl_gene_id")[["ensembl_gene_id", "gene_symbol", "gene_type"]]
            .set_index("ensembl_gene_id").loc[logtpm.index].to_csv(f"{prefix}_genes.tsv", sep="\t"))
        print(f"{project}/{species}: {len(meta)}/{len(sdf)} samples, {keep.sum()} genes -> {prefix}_*")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]))
