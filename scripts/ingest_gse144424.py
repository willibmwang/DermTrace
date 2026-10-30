#!/usr/bin/env python3
"""Curate the deposited GSE144424 counts and metadata without fitting a model."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PACKAGE = ROOT / "data/curated/gse144424_cardiomyocyte_v1"
GENE_MAPPING = ROOT / "data/processed/tsaolab_o2/tsaolab_o2_human_genes.tsv.gz"
ANNOTATION_COLUMNS = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
MARKERS = ["CA9", "EGLN3", "SLC2A1", "LDHA", "PDK1", "VEGFA", "BNIP3", "NDRG1"]
REPLICATED_DONORS = {"18511", "18852", "18855"}
PAPER_URL = "https://elifesciences.org/articles/57345"
INPUT_URLS = {
    "GSE144424_Counts_RNA_MCW_NEB.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144424/suppl/GSE144424_Counts_RNA_MCW_NEB.txt.gz",
    "GSE144424_family.soft.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144424/soft/GSE144424_family.soft.gz",
    "GSE144424_series_matrix.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144424/matrix/GSE144424_series_matrix.txt.gz",
}
# Immutable inputs for this version. A changed deposit or annotation requires a
# new package version; recomputing a hash from a changed file is not validation.
INPUT_SHA256 = {
    "GSE144424_Counts_RNA_MCW_NEB.txt.gz": "cad9ac4c6514550ea9bfb2b491cc2934f6952894d7cbd17338d5054d03da6f7c",
    "GSE144424_family.soft.gz": "156b46d9e11fc4a59219530808c617887ac622e711d137964e1ccb88b0766907",
    "GSE144424_series_matrix.txt.gz": "40dda761c345cb757eaa4d2c0eefcebad948c134137ffdcdc0af855384c90ec4",
}
GENE_MAPPING_SHA256 = "b463ad00c0f6d6b3de0ec78fc8de3c3df44cb7451a56dd6e8d2a1b3f2dd43ccb"
OUTPUT_NAMES = ("counts.parquet", "genes.csv", "samples.csv", "marker_mapping.csv", "provenance.json", "README.md")
CONDITIONS = {
    "A": ("normoxia", 10, 0, 0, 6),
    "B": ("hypoxia", 1, 6, 0, 6),
    "C": ("re-oxygenation1", 10, 6, 6, 12),
    "D": ("re-oxygenation2", 10, 6, 24, 30),
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_soft_samples(path: Path) -> list[dict[str, list[str]]]:
    """Retain repeated SOFT fields rather than silently overwriting them."""
    samples = []
    current = None
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("^SAMPLE = "):
                current = {"entity_accession": [line.strip().split(" = ", 1)[1]]}
                samples.append(current)
            elif line.startswith("^"):
                current = None
            elif current is not None and line.startswith("!Sample_"):
                key, value = line.rstrip("\n").split(" = ", 1)
                current.setdefault(key.removeprefix("!Sample_"), []).append(value)
    return samples


def single(record: dict[str, list[str]], field: str) -> str:
    values = record[field]
    if len(values) != 1:
        raise ValueError(f"Expected one {field}, found {values}")
    return values[0]


def sample_table(soft_path: Path, series_path: Path) -> pd.DataFrame:
    rows = []
    for record in read_soft_samples(soft_path):
        title = single(record, "title")
        if single(record, "geo_accession") != single(record, "entity_accession"):
            raise ValueError(f"SOFT entity/accession mismatch: {title}")
        match = re.fullmatch(r"(\d+)(?:_([23]))?_([ABCD])_RNA-seq", title)
        if match is None:
            raise ValueError(f"Unexpected RNA sample title: {title}")
        donor, replicate_text, condition = match.groups()
        replicate = int(replicate_text or "1")
        characteristics = {}
        for entry in record["characteristics_ch1"]:
            key, value = entry.split(": ", 1)
            if key.lower() in characteristics:
                raise ValueError(f"Duplicate characteristic {key} in {title}")
            characteristics[key.lower()] = value
        if donor != characteristics["subject id"]:
            raise ValueError(f"Title/subject mismatch: {title}")
        label, final_oxygen, hypoxia_hours, reoxygenation_hours, elapsed_hours = CONDITIONS[condition]
        if characteristics["condition"] != label:
            raise ValueError(f"Title/condition mismatch: {title}")
        if single(record, "organism_ch1") != "Homo sapiens":
            raise ValueError(f"Unexpected organism: {title}")
        matrix_id = "H" + title.removesuffix("_RNA-seq")
        flag, notes = "no_specific_paper_flag", "No sample-specific quality conclusion is inferred."
        if matrix_id == "H18852_A":
            flag = "published_expression_outlier"
            notes = "Paper excluded this expression outlier; recorded RIN=6. Retained here with a flag."
        elif matrix_id == "H18852_D":
            flag = "published_rin_undetermined"
            notes = "Paper reports undetermined RIN; this experiment was not selected for published DE/eQTL."
        rows.append({
            "sample_id": matrix_id,
            "matrix_id": matrix_id,
            "geo_accession": single(record, "geo_accession"),
            "source_title": title,
            "source_subject_id": characteristics["subject id"],
            "donor_id": donor,
            "replicate": replicate,
            "experiment_id": f"{donor}_rep{replicate}",
            "condition": condition,
            "condition_label": label,
            "organism": single(record, "organism_ch1"),
            "cell_type": characteristics["cell type"],
            "population": characteristics["population"],
            "sex": characteristics["sex"],
            "platform_id": single(record, "platform_id"),
            "oxygen_reference_pct": 10,
            "oxygen_final_pct": final_oxygen,
            "hypoxia_hours": hypoxia_hours,
            "reoxygenation_hours": reoxygenation_hours,
            "elapsed_experiment_hours": elapsed_hours,
            "published_primary_eligible": replicate == (2 if donor in REPLICATED_DONORS else 1),
            "quality_flag": flag,
            "qc_notes": notes,
            "treatment_protocol": single(record, "treatment_protocol_ch1"),
        })
    table = pd.DataFrame(rows)
    for column in ("sample_id", "matrix_id", "geo_accession"):
        if table[column].duplicated().any():
            raise ValueError(f"Duplicate {column} in SOFT metadata")
    if len(table) != 84 or table.donor_id.nunique() != 15 or table.experiment_id.nunique() != 21:
        raise ValueError("Expected 84 profiles from 15 donors and 21 experiments")
    if not table.groupby("experiment_id").condition.apply(lambda values: set(values) == set(CONDITIONS)).all():
        raise ValueError("Every differentiation must contain all four conditions")
    if not table.groupby("experiment_id").size().eq(4).all():
        raise ValueError("Every differentiation must contain exactly four profiles")
    replicated = set(table.groupby("donor_id").replicate.nunique().loc[lambda n: n > 1].index)
    if replicated != REPLICATED_DONORS:
        raise ValueError("Unexpected replicated donor identities")
    primary = table.loc[table.published_primary_eligible]
    if len(primary) != 60 or not primary.groupby("donor_id").size().eq(4).all():
        raise ValueError("Published selection must have four profiles for each of 15 donors")
    with gzip.open(series_path, "rt", encoding="utf-8") as handle:
        series_rows = {
            row[0]: row[1:] for row in csv.reader(handle, delimiter="\t")
            if row and row[0] in {"!Sample_title", "!Sample_geo_accession"}
        }
    series_pairs = list(zip(series_rows["!Sample_geo_accession"], series_rows["!Sample_title"], strict=True))
    if len(series_pairs) != 84 or set(series_pairs) != set(zip(table.geo_accession, table.source_title)):
        raise ValueError("Series matrix and SOFT do not identify the same 84 samples")
    return table


def curate(package: Path = DEFAULT_PACKAGE, mapping_path: Path = GENE_MAPPING) -> dict:
    """Validate source identities and write only the new curated package."""
    mapping_path = mapping_path.resolve()
    try:
        mapping_reference = mapping_path.relative_to(ROOT).as_posix()
        mapping_basis = "repository_root"
    except ValueError:
        mapping_reference = mapping_path.name
        mapping_basis = "external_input_filename"
    inputs = package / "inputs"
    for name, expected in INPUT_SHA256.items():
        if sha256(inputs / name) != expected:
            raise ValueError(f"Pinned source checksum mismatch for {name}; use a new package version for changed inputs")
    mapping_hash = sha256(mapping_path)
    if mapping_hash != GENE_MAPPING_SHA256:
        raise ValueError("Pinned gene mapping checksum mismatch; use a new package version for changed annotation")
    manifest_path = package / "manifest.json"
    if manifest_path.exists():
        existing = json.loads(manifest_path.read_text())
        recorded = {Path(item["path"]).name: item["sha256"] for item in existing.get("inputs", [])}
        if (existing.get("dataset_id") != "gse144424_cardiomyocyte_v1"
                or recorded != INPUT_SHA256
                or existing.get("local_marker_mapping_sha256") != mapping_hash
                or set(existing.get("outputs", {})) != set(OUTPUT_NAMES)):
            raise ValueError("Existing package manifest disagrees with pinned version; refuse replacement")
        for name, expected in existing["outputs"].items():
            if sha256(package / name) != expected:
                raise ValueError(f"Existing package output changed: {name}; refuse replacement")
        # Preserve the existing immutable release and its recorded builder hash.
        return existing
    if any((package / name).exists() for name in OUTPUT_NAMES):
        raise ValueError("Existing outputs lack a manifest; refuse same-version replacement")
    sources = [
        {"path": f"inputs/{name}", "url": url, "sha256": sha256(inputs / name)}
        for name, url in INPUT_URLS.items()
    ]
    raw = pd.read_csv(inputs / "GSE144424_Counts_RNA_MCW_NEB.txt.gz", sep="\t")
    if raw.columns[:6].tolist() != ANNOTATION_COLUMNS or raw.Geneid.duplicated().any():
        raise ValueError("Invalid featureCounts annotation or duplicate gene identifiers")
    if not raw.Geneid.str.fullmatch(r"ENSG\d{11}").all():
        raise ValueError("Expected unversioned Ensembl gene identifiers")
    genes = raw[ANNOTATION_COLUMNS].copy()
    counts = raw.drop(columns=ANNOTATION_COLUMNS).copy()
    values = counts.to_numpy()
    if not np.issubdtype(values.dtype, np.number):
        raise ValueError("Counts must be numeric")
    if not np.isfinite(values).all() or (values < 0).any() or not np.equal(values, np.floor(values)).all():
        raise ValueError("Counts must be finite, nonnegative integers")
    counts = counts.astype("int64")
    counts.index = pd.Index(raw.Geneid, name="Geneid")
    samples = sample_table(inputs / "GSE144424_family.soft.gz", inputs / "GSE144424_series_matrix.txt.gz")
    if len(counts.columns) != 84 or set(counts.columns) != set(samples.matrix_id):
        raise ValueError("Count columns and sample metadata are not a bijection")
    samples = samples.set_index("matrix_id", drop=False).loc[counts.columns].reset_index(drop=True)
    samples["deposited_gene_count_sum"] = counts.sum(axis=0).to_numpy()
    mapping = pd.read_csv(mapping_path, sep="\t", dtype=str)
    marker_rows = []
    for marker in MARKERS:
        selected = mapping.loc[mapping.gene_symbol.eq(marker), "ensembl_gene_id"].drop_duplicates()
        if len(selected) != 1:
            raise ValueError(f"Expected one source-mapped Ensembl ID for {marker}, found {selected.tolist()}")
        gene_id = selected.iloc[0]
        if gene_id not in counts.index:
            raise ValueError(f"Source-mapped marker absent from counts: {marker} ({gene_id})")
        marker_rows.append({
            "gene_symbol": marker, "ensembl_gene_id": gene_id,
            "mapping_source": mapping_reference, "mapping_source_basis": mapping_basis,
            "mapping_source_sha256": mapping_hash, "present_in_counts": True,
        })
    marker_mapping = pd.DataFrame(marker_rows)
    chromosomes = {key: int(value) for key, value in genes.Chr.value_counts().sort_index().items()}
    autosomal = genes.Chr.isin([f"chr{i}" for i in range(1, 23)])
    audit = {
        "n_genes": len(genes), "n_profiles": len(samples), "n_donors": samples.donor_id.nunique(),
        "n_differentiation_experiments": samples.experiment_id.nunique(),
        "n_published_primary_profiles": int(samples.published_primary_eligible.sum()),
        "matrix_metadata_bijection": True, "series_soft_sample_identity_match": True,
        "finite_nonnegative_integer_counts": True, "chromosome_gene_counts": chromosomes,
        "autosomal_gene_count": int(autosomal.sum()), "nonautosomal_gene_count": int((~autosomal).sum()),
        "geo_autosomal_only_claim_matches_file": bool(autosomal.all()),
        "flagged_profiles_retained": samples.loc[samples.quality_flag.ne("no_specific_paper_flag"), "sample_id"].tolist(),
    }
    provenance = {
        "study": "GSE144424", "paper_url": PAPER_URL, "input_files": sources,
        "local_marker_mapping": {
            "path": mapping_reference, "path_basis": mapping_basis, "sha256": mapping_hash,
        },
        "deposited_units": "Unnormalized, WASP-filtered featureCounts gene counts; not TPM or RUV-corrected values.",
        "genome_build": "hg19",
        "metadata_discrepancy": "GEO says autosomal-only, but the deposited matrix contains chrX, chrY and chrM. All source rows are preserved.",
        "published_selection": "Replicate 2 for donors 18511, 18852 and 18855; sole experiment for the other 12 donors. Selection is documented by the paper, not fitted here.",
        "qc_source": "Paper RNA/DNA extraction and RNA-seq analysis sections: 18852A is an expression outlier with RIN=6; 18852D has undetermined RIN.",
        "pairing": "Condition subcultures were split from the same differentiation; repeated differentiations share a donor and are not independent donors.",
        "timing_limit": "A/B were harvested at 6 h, C at 12 h and D at 30 h after experiment start; later time-matched continued-normoxia controls are not described.",
        "scope": "Standalone exploratory package, excluded from every existing frozen release, model and benchmark. No normalization or predictive training is performed by ingestion.",
        "audit": audit,
    }
    package.mkdir(parents=True, exist_ok=True)
    counts.to_parquet(package / "counts.parquet", compression="zstd")
    genes.to_csv(package / "genes.csv", index=False)
    samples.to_csv(package / "samples.csv", index=False)
    marker_mapping.to_csv(package / "marker_mapping.csv", index=False)
    (package / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    input_lines = "\n".join(f"- `{item['path']}` — SHA-256 `{item['sha256']}` — {item['url']}" for item in sources)
    readme = f"""# GSE144424 cardiomyocyte expression package v1

Standalone exploratory data addition to DermTrace. Existing frozen datasets, models and benchmarks are unchanged. This package does not establish predictive improvement or a clinical claim.

## Material and units

The deposited matrix contains {len(genes):,} Ensembl genes × 84 profiles from 15 human donors and 21 differentiation experiments. Donors 18511, 18852 and 18855 each have three experiments. Each experiment has four condition subcultures. The matrix contains unnormalized, WASP-filtered featureCounts counts aligned to hg19, with native source annotations preserved. It is not TPM, TMM-normalized or RUV-corrected expression.

| Condition | Exposure | Harvest after experiment start |
|---|---|---:|
| A | 10% oxygen reference | 6 h |
| B | 1% oxygen for 6 h | 6 h |
| C | 1% for 6 h, then 10% for 6 h | 12 h |
| D | 1% for 6 h, then 10% for 24 h | 30 h |

The source describes iPSC-derived cardiomyocyte cultures (reported purity 40–97%, median 81%). Condition subcultures came from the same differentiation. Replicate differentiations are nested within donor. Later time-matched normoxia controls are not described, limiting a causal interpretation of C/D versus A. This experiment has a 10% reference, unlike the Tsao 20% reference; do not pool absolute expression values across studies.

## QC and provenance decisions

All 84 deposited profiles remain present. `H18852_A` is flagged as the paper's expression outlier (RIN=6); `H18852_D` is flagged for undetermined RIN. An unflagged row is not a new claim that it passed an independent QC assessment.

`published_primary_eligible` selects experiment 2 for each of the three repeated donors and the sole experiment for each other donor: 60 profiles, 15 donors. The paper used this selection for its DE/eQTL analyses. Its differential model also used TMM/voom/limma and four RUV factors; those analyses are not reproduced by this ingestion script.

GEO metadata says the counts contain only autosomal genes, but the actual file contains {int(autosomal.sum()):,} autosomal genes and {int((~autosomal).sum()):,} nonautosomal genes (chrX, chrY, chrM). All rows are retained and their chromosomes are explicit in `genes.csv`.

Eight marker Ensembl IDs are verified from the existing Tsao human annotation, using exact gene-symbol matches and requiring one unique ID per marker and presence in the deposited counts. The local mapping is used only for gene identity, not expression values or outcomes; its gene-type labels are not used for filtering.

## Files

- `counts.parquet`: genes × 84 sample IDs; index `Geneid`; native integer counts.
- `genes.csv`: the six original annotation columns (`Geneid`, `Chr`, `Start`, `End`, `Strand`, `Length`).
- `samples.csv`: GEO accession, native matrix ID, donor, differentiation, condition, timing, published eligibility and explicit QC flags. Count-library totals refer only to the deposited gene matrix.
- `marker_mapping.csv`: eight exact symbol–Ensembl mappings and their annotation source hash.
- `provenance.json`: source URLs/hashes, decisions and audit results.
- `manifest.json`: input/output hashes and validation summary.

## Source hashes

{input_lines}
- Local gene mapping: `{mapping_reference}` (path basis: `{mapping_basis}`) — SHA-256 `{mapping_hash}`.
- Design, pairing and published exclusions: {PAPER_URL}.

Mapping paths with basis `repository_root` are relative to the checkout root. Custom mappings outside the checkout are recorded by filename with basis `external_input_filename`; supply that file through `curate(..., mapping_path=...)` and verify its SHA-256 to reproduce the package.

From the repository root, use Python 3.10 with `keloidbench/requirements-revision.txt`, then run `python scripts/ingest_gse144424.py`. The script is local-only, performs no downloads, and writes only this package. Raw reads were not downloaded.
"""
    (package / "README.md").write_text(readme)
    manifest = {
        "dataset_id": "gse144424_cardiomyocyte_v1", "status": "exploratory_only",
        "frozen_benchmark_member": False, "inputs": sources,
        "local_marker_mapping_sha256": mapping_hash, "audit": audit,
        "outputs": {name: sha256(package / name) for name in OUTPUT_NAMES},
        "ingestion_script_sha256": sha256(Path(__file__)),
    }
    (package / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package", type=Path, default=DEFAULT_PACKAGE)
    args = parser.parse_args()
    manifest = curate(args.package)
    print(json.dumps(manifest["audit"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
