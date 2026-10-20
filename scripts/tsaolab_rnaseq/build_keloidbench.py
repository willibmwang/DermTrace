"""Combine the salmon results of every Quintara order into one KeloidBench/DermTrace-ready dataset.

KeloidBench loaders expect genes x samples matrices of NON-log values (they apply log1p themselves)
and sample metadata with the METADATA_COLUMNS schema (scripts/expression_processing.py).
Outputs, per species, in <dest>/:

  tsaolab_o2_<species>_tpm_symbol.tsv        gene_symbol x sample, TPM (duplicate symbols summed)
  tsaolab_o2_<species>_counts_symbol.tsv     gene_symbol x sample, salmon estimated counts
  tsaolab_o2_<species>_tpm_ensembl.tsv       ensembl_gene_id x sample, TPM
  tsaolab_o2_<species>_counts_ensembl.tsv    ensembl_gene_id x sample, counts
  tsaolab_o2_<species>_sample_metadata.csv   KeloidBench METADATA_COLUMNS + experiment/QC columns
  tsaolab_o2_<species>_genes.tsv             ensembl_gene_id, gene_symbol, gene_type

Usage: python -I build_keloidbench.py <dest_dir> [--qc qc_identity.csv] <sample_sheet.csv>=<quant_dir> [...]
"""
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))  # python -I does not add the script dir
from aggregate import FASTA, MIN_MAPPING_RATE, tx2gene  # noqa: E402

PREFIX = "tsaolab_o2"
TREATMENT = {  # 20% O2 is the untreated baseline arm
    "20%": "none",
    "1%": "hypoxia_1pct",
    "20%->1%": "normoxia_to_hypoxia_1pct",
    "1%->20%": "hypoxia_1pct_to_normoxia",
}
KB_COLUMNS = [
    "sample_id", "accession", "sample_title", "platform_id", "modality", "source_dataset",
    "disease_domain", "disease_label", "keloid_vs_normal", "lesional_status", "scar_type",
    "cell_type", "treatment", "patient_id", "contrast_type", "eligible_for_keloid_pretraining",
    "encoder_task", "encoder_prompt", "encoder_response", "eval_grain", "cohort_role",
    "allow_accession_level_grouping", "lockbox",
]


def kb_metadata(r, qc):
    keloid = r.disease_label == "keloid"
    # skin fibroblasts that passed identity QC are the only samples offered to keloid pretraining
    human_fib = r.species == "human" and r.cell_type == "fibroblast" and qc["quality_flag"] == "ok"
    return dict(
        sample_id=r.qb_id,
        accession=f"TSAO_QB{r.project}",
        sample_title=r.sample_name,
        platform_id="Element_AVITI",
        modality="bulk_rnaseq",
        source_dataset=f"Tsao_lab_Quintara_{r.project}",
        disease_domain="keloid" if keloid else "unknown",
        disease_label=r.disease_label if r.cell_type == "fibroblast" else "unknown",
        # only skin fibroblasts carry a keloid-vs-normal label; cardiomyocytes stay out of keloid tasks
        keloid_vs_normal=("keloid" if keloid else "normal") if r.cell_type == "fibroblast" else "unknown",
        lesional_status="unknown",
        scar_type="keloid" if keloid else "unknown",
        cell_type=r.cell_type,
        treatment=TREATMENT[r.oxygen],
        patient_id=r.cell_line,  # one donor/cell line per culture; replicates share it
        contrast_type="oxygen_tension",
        eligible_for_keloid_pretraining=bool(human_fib),
        encoder_task="unknown", encoder_prompt="unknown", encoder_response="unknown",
        eval_grain="sample",
        cohort_role="unassigned",
        allow_accession_level_grouping=False,
        lockbox=False,
        # experiment + QC (extra columns; KeloidBench ignores them)
        species=r.species, submitted_as=r.submitted_as, cell_line=r.cell_line,
        group=f"{r.project}_{r.group}", oxygen=r.oxygen,
        collection_date=r.collection_date, order_id=r.project, **qc,
    )


def main(dest, runs, qc_path=None):
    dest.mkdir(parents=True, exist_ok=True)
    qc = pd.read_csv(qc_path, dtype={"sample_id": str}).set_index("qb_id") if qc_path else None
    qc_cols = ["xist_tpm", "chrY_tpm", "ad5_e1_per_million", "hek_admixture_est_pct", "sex_inferred"]
    samples = []
    for sheet, qdir in runs:
        s = pd.read_csv(sheet, dtype=str)
        s["qdir"] = str(qdir)
        samples.append(s)
    samples = pd.concat(samples, ignore_index=True)

    for species, sdf in samples.groupby("species"):
        t2g = tx2gene(FASTA[species]).set_index("tx")
        tpm, counts, meta, skipped = {}, {}, [], []
        for r in sdf.itertuples():
            qd = Path(r.qdir) / "quant" / r.sample_id
            if not (qd / "quant.sf").exists():
                skipped.append(r.qb_id)
                continue
            q = pd.read_csv(qd / "quant.sf", sep="\t", index_col=0)
            gene = t2g.loc[q.index, "ensembl_gene_id"].values
            tpm[r.qb_id] = q.TPM.groupby(gene).sum()
            counts[r.qb_id] = q.NumReads.groupby(gene).sum()
            info = json.loads((qd / "aux_info" / "meta_info.json").read_text())
            rate = info["percent_mapped"]
            flags, notes = ["ok" if rate >= MIN_MAPPING_RATE else "suspect"], []
            if rate < MIN_MAPPING_RATE:
                notes.append(f"low mapping rate ({rate:.1f}%)")
            extra = {}
            if qc is not None:
                q = qc.loc[r.qb_id]
                flags.append(q.qc_flag)
                if isinstance(q.qc_note, str) and q.qc_note:
                    notes.append(q.qc_note)
                extra = {c: q.get(c) for c in qc_cols}
            meta.append(kb_metadata(r, dict(
                num_fragments=info["num_processed"], mapping_rate_pct=round(rate, 2),
                quality_flag="suspect" if "suspect" in flags else "ok", qc_notes="; ".join(notes), **extra,
            )))
        if skipped:
            print(f"{species}: not quantified yet, left out: {', '.join(skipped)}")
        if not meta:
            continue

        tpm = pd.DataFrame(tpm).fillna(0)
        counts = pd.DataFrame(counts).fillna(0)
        keep = counts.sum(axis=1) > 0
        tpm, counts = tpm[keep].rename_axis("ensembl_gene_id"), counts[keep].rename_axis("ensembl_gene_id")
        genes = (t2g.drop_duplicates("ensembl_gene_id").set_index("ensembl_gene_id")
                 .loc[tpm.index, ["gene_symbol", "gene_type"]])
        sym = genes.gene_symbol

        p = dest / f"{PREFIX}_{species}"
        tpm.to_csv(f"{p}_tpm_ensembl.tsv", sep="\t", float_format="%.4f")
        counts.to_csv(f"{p}_counts_ensembl.tsv", sep="\t", float_format="%.3f")
        tpm.groupby(sym.values).sum().rename_axis("gene_symbol").to_csv(
            f"{p}_tpm_symbol.tsv", sep="\t", float_format="%.4f")
        counts.groupby(sym.values).sum().rename_axis("gene_symbol").to_csv(
            f"{p}_counts_symbol.tsv", sep="\t", float_format="%.3f")
        genes.to_csv(f"{p}_genes.tsv", sep="\t")
        meta = pd.DataFrame(meta)
        meta = meta[KB_COLUMNS + [c for c in meta.columns if c not in KB_COLUMNS]]
        meta.to_csv(f"{p}_sample_metadata.csv", index=False)
        print(f"{species}: {len(meta)} samples, {len(tpm)} genes ({sym.nunique()} symbols) -> {p}_*")


if __name__ == "__main__":
    args = sys.argv[2:]
    qc_path = None
    if args[:1] == ["--qc"]:
        qc_path, args = Path(args[1]), args[2:]
    main(Path(sys.argv[1]), [tuple(map(Path, a.split("=", 1))) for a in args], qc_path)
