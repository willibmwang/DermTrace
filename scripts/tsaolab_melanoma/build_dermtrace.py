"""DermTrace-format dataset from the 20088-11 salmon runs.

Outputs in <dest>/ (genes x samples, NON-log; DermTrace loaders apply log1p):
  tsaolab_melanoma_tpm_symbol.tsv.gz / _counts_symbol.tsv.gz / _tpm_ensembl.tsv.gz / _counts_ensembl.tsv.gz
  tsaolab_melanoma_genes.tsv.gz                 ensembl_gene_id, gene_symbol, gene_type
  tsaolab_melanoma_sample_metadata.csv          METADATA_COLUMNS + specimen/QC/sex/hotspot columns

Usage: python -I build_dermtrace.py <dest_dir> <sample_sheet.csv> <quant_dir> [--runqc run_qc.tsv] [--hotspots hotspot_summary.tsv]
"""
import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PREFIX = "tsaolab_melanoma"
FASTA = Path(__file__).resolve().parents[3] / "references" / "gencode_human_v46" / "gencode.v46.transcripts.fa.gz"
MIN_MAPPING_RATE = 30.0  # FFPE total RNA: much of the library is intronic, so the transcriptome mapping rate is lower than for fresh cells
KB_COLUMNS = [
    "sample_id", "accession", "sample_title", "platform_id", "modality", "source_dataset",
    "disease_domain", "disease_label", "keloid_vs_normal", "lesional_status", "scar_type",
    "cell_type", "treatment", "patient_id", "contrast_type", "eligible_for_keloid_pretraining",
    "encoder_task", "encoder_prompt", "encoder_response", "eval_grain", "cohort_role",
    "allow_accession_level_grouping", "lockbox",
]
SEX_GENES = dict(xist="XIST", chrY=["RPS4Y1", "DDX3Y", "KDM5D", "EIF1AY", "UTY", "USP9Y", "ZFY"])


def tx2gene(fasta):
    rows = []
    with gzip.open(fasta, "rt") as fh:
        for line in fh:
            if line[0] == ">":
                f = line[1:].rstrip().split("|")
                rows.append((f[0], f[1].split(".")[0], f[5], f[7]))
    return pd.DataFrame(rows, columns=["tx", "ensembl_gene_id", "gene_symbol", "gene_type"]).set_index("tx")


def infer_sex(xist, chry):
    if chry >= 5 and xist < 5:
        return "male"
    if xist >= 5 and chry < 2:
        return "female"
    return "ambiguous"


def main(dest, sheet, qdir, runqc=None, hotspots=None):
    dest.mkdir(parents=True, exist_ok=True)
    samples = pd.read_csv(sheet, dtype=str)
    t2g = tx2gene(FASTA)
    runqc = pd.read_csv(runqc, sep="\t").set_index("sample_id") if runqc else None
    hot = pd.read_csv(hotspots, sep="\t").set_index("sample_id") if hotspots else None
    tpm, counts, meta, skipped = {}, {}, [], []
    for r in samples.itertuples():
        qd = qdir / "quant" / r.sample_id
        if not (qd / "quant.sf").exists():
            skipped.append(r.sample_id)
            continue
        q = pd.read_csv(qd / "quant.sf", sep="\t", index_col=0)
        gene = t2g.loc[q.index, "ensembl_gene_id"].values
        tpm[r.sample_id] = q.TPM.groupby(gene).sum()
        counts[r.sample_id] = q.NumReads.groupby(gene).sum()
        info = json.loads((qd / "aux_info" / "meta_info.json").read_text())
        rate = info["percent_mapped"]
        meta.append((r, info["num_processed"], rate))
    if skipped:
        print("not quantified yet, left out:", " ".join(skipped))
    tpm = pd.DataFrame(tpm).fillna(0)
    counts = pd.DataFrame(counts).fillna(0)
    keep = counts.sum(axis=1) > 0
    tpm, counts = tpm[keep].rename_axis("ensembl_gene_id"), counts[keep].rename_axis("ensembl_gene_id")
    genes = t2g.drop_duplicates("ensembl_gene_id").set_index("ensembl_gene_id").loc[tpm.index, ["gene_symbol", "gene_type"]]
    sym_tpm = tpm.groupby(genes.gene_symbol.values).sum()

    rows = []
    for r, nfrag, rate in meta:
        xist = float(sym_tpm.loc[SEX_GENES["xist"], r.sample_id]) if SEX_GENES["xist"] in sym_tpm.index else 0.0
        chry = float(sym_tpm.reindex(SEX_GENES["chrY"]).fillna(0)[r.sample_id].sum())
        notes = []
        if rate < MIN_MAPPING_RATE:
            notes.append(f"low mapping rate ({rate:.1f}%)")
        extra = {}
        if runqc is not None and r.sample_id in runqc.index:
            q = runqc.loc[r.sample_id]
            extra = dict(instrument=q.instrument, sequencing_run=q.sequencing_run, lane=q.lane, reads_raw=q.reads_raw,
                         pct_adapter_trimmed=q.pct_adapter_trimmed, insert_size_peak=q.insert_size_peak,
                         fastp_duplication_rate=q.fastp_duplication_rate, library_type=q.library_type)
        if hot is not None and r.sample_id in hot.index:
            h = hot.loc[r.sample_id]
            extra.update(braf_v600_status=h.BRAF_600_status, nras_q61_status=h.NRAS_61_status,
                         hotspot_calls=h.hotspot_calls, other_driver_calls=h.other_driver_calls)
        rows.append(dict(
            sample_id=r.sample_id, accession="TSAO_20088-11", sample_title=r.core_id,
            platform_id=extra.get("instrument", "Illumina_NovaSeqX"), modality="bulk_rnaseq",
            source_dataset="Tsao_lab_Greek_melanoma_20088-11", disease_domain="melanoma", disease_label="melanoma",
            keloid_vs_normal="unknown", lesional_status="lesional", scar_type="unknown", cell_type="tissue",
            treatment="unknown", patient_id="unknown", contrast_type="unknown",
            eligible_for_keloid_pretraining=False, encoder_task="unknown", encoder_prompt="unknown",
            encoder_response="unknown", eval_grain="sample", cohort_role="unassigned",
            allow_accession_level_grouping=False, lockbox=False,
            species="human", tissue=r.tissue, specimen=r.specimen, core_id=r.core_id, project=r.project,
            num_fragments=nfrag, mapping_rate_pct=round(rate, 2),
            quality_flag="suspect" if notes else "ok", qc_notes="; ".join(notes),
            xist_tpm=round(xist, 2), chrY_tpm=round(chry, 2), sex_inferred=infer_sex(xist, chry), **extra,
        ))
    p = dest / PREFIX
    tpm.to_csv(f"{p}_tpm_ensembl.tsv.gz", sep="\t", float_format="%.4f")
    counts.to_csv(f"{p}_counts_ensembl.tsv.gz", sep="\t", float_format="%.3f")
    sym_tpm.rename_axis("gene_symbol").to_csv(f"{p}_tpm_symbol.tsv.gz", sep="\t", float_format="%.4f")
    counts.groupby(genes.gene_symbol.values).sum().rename_axis("gene_symbol").to_csv(f"{p}_counts_symbol.tsv.gz", sep="\t", float_format="%.3f")
    genes.to_csv(f"{p}_genes.tsv.gz", sep="\t")
    md = pd.DataFrame(rows)
    md = md[KB_COLUMNS + [c for c in md.columns if c not in KB_COLUMNS]]
    md.to_csv(f"{p}_sample_metadata.csv", index=False)
    print(f"{len(md)} samples, {len(tpm)} genes ({sym_tpm.shape[0]} symbols) -> {p}_*")
    print(md[["sample_id", "mapping_rate_pct", "sex_inferred", "quality_flag"]].to_string(index=False))


if __name__ == "__main__":
    args, opts = sys.argv[1:], {}
    while "--runqc" in args or "--hotspots" in args:
        for k in ("--runqc", "--hotspots"):
            if k in args:
                i = args.index(k)
                opts[k[2:]] = Path(args[i + 1])
                del args[i:i + 2]
    main(Path(args[0]), Path(args[1]), Path(args[2]), opts.get("runqc"), opts.get("hotspots"))
