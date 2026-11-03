#!/usr/bin/env python3
"""Build the versioned GSE145725 package from pinned, processed GEO inputs.

No downloads or raw-array processing occur here. The deposited log2 RMA scale
is retained in the native table. The compatibility exports use DermTrace's
existing writer and its accession-wise standardization.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from expression_processing import METADATA_COLUMNS, is_hgnc_like_symbol, module_coverage_report, write_expression_artifacts
from preprocess_microarray import expression_to_genes, infer_labels, parse_series_matrix

REPO_ROOT = Path(__file__).resolve().parents[2]
PACKAGE = REPO_ROOT / "data/processed/microarray/GSE145725/v1"
PREFIX = "microarray"
INPUTS = {
    "GSE145725_series_matrix.txt.gz": "62fcdb50492ed769bdd986e69c5a7558fae308dc664d36bb641b084bc3c25b9c",
    "GPL16043_full.txt.gz": "43eb01e4a5a1b11ee7bc87790a379692c13ba755f51620fb4bd32040b2e58425",
    "GSM4331585_processed.rma.chp.gz": "ebcf2dae0210ac65adfa0ee44db427ec38776838959c424aba1f26ffd9274631",
}
PLATFORM_TEXT_SHA256 = "6cc8e8fcdfc057411b6c69befd6533fdf9b5a4a8afc61089f7a26665ef884c54"
SAMPLE_IDS = [f"GSM{i}" for i in range(4331585, 4331604)]
SOURCE_GROUPS = ["keloid_fibroblast_1", "keloid_fibroblast_3", "keloid_fibroblast_4", "keloid_fibroblast_5", "keloid_fibroblast_6"]
SOURCE_GROUPS += [f"normal_fibroblast_{i}" for i in range(1, 6)]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_sources(source_dir: Path) -> None:
    for filename, expected in INPUTS.items():
        path = source_dir / filename
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f"Missing or changed pinned input: {filename}")


def platform_mapping(path: Path) -> pd.DataFrame:
    """Preserve every probe; do not choose a gene from multi-symbol annotations."""
    raw = gzip.decompress(path.read_bytes())
    if hashlib.sha256(raw).hexdigest() != PLATFORM_TEXT_SHA256:
        raise ValueError("Changed decompressed GPL16043 annotation")
    text = raw.decode("utf-8")
    table = text.split("!platform_table_begin", 1)[1].split("!platform_table_end", 1)[0].strip()
    frame = pd.read_csv(io.StringIO(table), sep="\t", dtype=str, keep_default_na=False)
    required = ["ID", "Gene Symbol", "Ensembl", "Entrez Gene"]
    if not set(required).issubset(frame):
        raise ValueError("GPL16043 annotation schema changed")
    audit = frame[required].rename(columns={"ID": "probe_id", "Gene Symbol": "original_gene_symbol", "Ensembl": "original_ensembl", "Entrez Gene": "original_entrez_gene"}).copy()
    if len(audit) != 49495 or audit.probe_id.duplicated().any():
        raise ValueError("Unexpected platform probe count or duplicate probe ID")
    def status(value: str) -> str:
        value = value.strip()
        if value in {"", "---"}:
            return "unmapped"
        if "///" in value:
            return "ambiguous"
        return "retained" if is_hgnc_like_symbol(value.upper()) else "invalid_symbol"
    audit["mapping_status"] = audit.original_gene_symbol.map(status)
    # Uppercase matches the shared vocabulary convention, without alias updates
    # or stripping punctuation to manufacture a different gene symbol.
    audit["gene_symbol"] = audit.original_gene_symbol.str.strip().str.upper().where(audit.mapping_status.eq("retained"), "")
    if audit.mapping_status.value_counts().to_dict() != {"retained": 47239, "ambiguous": 1457, "unmapped": 688, "invalid_symbol": 111}:
        raise ValueError("Unexpected GPL16043 mapping coverage")
    return audit


def prepare(source_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    verify_sources(source_dir)
    series = parse_series_matrix(source_dir / "GSE145725_series_matrix.txt.gz")
    if series.platform_id != "GPL16043" or series.expression.shape != (49495, 19):
        raise ValueError("Unexpected series platform or dimensions")
    if series.expression.columns.tolist() != SAMPLE_IDS or series.metadata.sample_id.tolist() != SAMPLE_IDS:
        raise ValueError("Unexpected sample IDs/order or metadata alignment")
    if series.expression.index.duplicated().any() or not np.isfinite(series.expression.to_numpy()).all():
        raise ValueError("Invalid deposited expression")
    metadata = infer_labels(series.metadata)
    if metadata.disease_label.value_counts().to_dict() != {"normal": 10, "keloid": 9}:
        raise ValueError("Unexpected disease labels")
    if sorted(metadata.patient_id.unique()) != sorted(SOURCE_GROUPS):
        raise ValueError("Unexpected recorded source groups")
    if not metadata.treatment_protocol_ch1.str.contains("normoxia (20% O2, 5% CO2)", regex=False).all():
        raise ValueError("Unexpected oxygen treatment")
    metadata["replicate_id"] = metadata["title"].str.extract(r"_rep([12])$", expand=False).astype(int)
    if metadata[["patient_id", "replicate_id"]].duplicated().any():
        raise ValueError("Duplicate source/replicate")
    metadata["source_group_id"] = metadata["patient_id"]
    metadata["patient_id_basis"] = "deposited_source_label_not_verified_donor"
    metadata["oxygen"] = "20%"
    metadata["passage"] = 3
    metadata["growth_medium_hydrocortisone_ug_ml"] = 0.4
    metadata["eval_grain"] = "profile"
    metadata["cohort_role"] = "fibroblast_sensitivity_only"
    metadata["lockbox"] = False
    metadata["allow_accession_level_grouping"] = False
    mapping = platform_mapping(source_dir / "GPL16043_full.txt.gz")
    if set(series.expression.index) != set(mapping.probe_id):
        raise ValueError("Series and platform probe IDs differ")
    keep = mapping[mapping.mapping_status.eq("retained")]
    expr = expression_to_genes(series.expression, dict(zip(keep.probe_id, keep.gene_symbol)), series.platform_id)
    expr = expr.reindex(columns=sorted(expr.columns)).reset_index(drop=True)
    if expr.shape != (19, 19011) or not np.isfinite(expr.to_numpy()).all():
        raise ValueError("Unexpected mapped gene expression")
    return metadata.reset_index(drop=True), expr, mapping


def write_native(path: Path, metadata: pd.DataFrame, expr: pd.DataFrame) -> None:
    frame = pd.concat([metadata[METADATA_COLUMNS], expr], axis=1)
    table = pa.Table.from_pandas(frame, preserve_index=False)
    attrs = dict(table.schema.metadata or {})
    attrs[b"dermtrace.expression_scale"] = b"deposited_log2_RMA; arithmetic_mean_of_unambiguous_probes"
    attrs[b"dermtrace.gene_annotation"] = b"GPL16043_full.txt.gz; probe_mapping.tsv.gz; deposited_historical_symbols"
    pq.write_table(table.replace_schema_metadata(attrs), path)


def validate_release(directory: Path, metadata: pd.DataFrame, native: pd.DataFrame) -> dict:
    genes = (directory / "microarray_gene_vocab.txt").read_text().splitlines()
    wide = pd.read_parquet(directory / "microarray_expression_wide.parquet")
    saved_meta = pd.read_parquet(directory / "microarray_sample_metadata.parquet")
    saved_native = pd.read_parquet(directory / "microarray_expression_native.parquet")
    labels = pd.read_csv(directory / "microarray_encoder_labels.csv")
    matrix = np.load(directory / "microarray_X.npy", allow_pickle=False)
    if genes != native.columns.tolist() or genes != sorted(set(genes)):
        raise ValueError("Gene vocabulary/order mismatch")
    for frame in (wide, saved_meta, saved_native, labels):
        if frame.sample_id.tolist() != SAMPLE_IDS or frame.sample_id.duplicated().any():
            raise ValueError("Output sample alignment mismatch")
    for frame in (wide, saved_meta, saved_native):
        pd.testing.assert_frame_equal(frame[METADATA_COLUMNS], metadata[METADATA_COLUMNS], check_dtype=False)
    if wide.columns.tolist()[-len(genes):] != genes:
        raise ValueError("Wide gene order mismatch")
    np.testing.assert_array_equal(saved_native[genes].to_numpy(), native.to_numpy())
    expected = ((native - native.mean()) / native.std(ddof=1).replace(0, np.nan)).fillna(0)
    np.testing.assert_allclose(wide[genes].to_numpy(), expected.to_numpy(), atol=1e-12, rtol=0)
    np.testing.assert_array_equal(matrix, expected.to_numpy(dtype=np.float32))
    if not np.isfinite(matrix).all() or matrix.dtype != np.float32:
        raise ValueError("Invalid compatibility matrix")
    if labels.encoder_response.tolist() != metadata.encoder_response.tolist():
        raise ValueError("Label mismatch")
    coverage = module_coverage_report(native, accession="GSE145725")
    if coverage["all_zero_programs"]:
        raise ValueError("No mapped program signal")
    return {"n_profiles": len(wide), "n_genes": len(genes), "n_recorded_source_groups": int(metadata.patient_id.nunique()), "all_values_finite": True, "native_values_preserved": True, "standard_writer_reproduced": True, "n_modules_with_signal": coverage["n_modules_with_signal"]}


def build(source_dir: Path, output: Path) -> dict:
    """Validate in a temporary directory before publishing a new release."""
    if output.exists() and any(p.name not in {"source", "README.md"} for p in output.iterdir()):
        raise FileExistsError(f"Refusing to replace an existing release: {output}")
    if (output / "source").exists():
        verify_sources(output / "source")
    metadata, native, mapping = prepare(source_dir)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".gse145725-", dir=output.parent) as temp:
        stage = Path(temp)
        if not (output / "source").exists():
            (stage / "source").mkdir()
            for name in INPUTS:
                shutil.copyfile(source_dir / name, stage / "source" / name)
            for name in ("processed_chp_matrix_comparison.json", "sources.json"):
                if (source_dir / name).exists():
                    shutil.copyfile(source_dir / name, stage / "source" / name)
            verify_sources(stage / "source")
        summary = write_expression_artifacts(out_dir=stage, prefix=PREFIX, metadata=metadata, expr=native, dataset_summaries=[{"accession": "GSE145725", "platform_id": "GPL16043", "n_samples": 19, "n_probe_rows": 49495, "n_gene_features": 19011, "labels": {"normal": 10, "keloid": 9}}])
        write_native(stage / "microarray_expression_native.parquet", metadata, native)
        # Gzip timestamps are fixed for reproducible annotation/source records.
        with (stage / "probe_mapping.tsv.gz").open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
                zipped.write(mapping.to_csv(sep="\t", index=False).encode())
        metadata.to_csv(stage / "sample_metadata.tsv", sep="\t", index=False)
        summary["outputs"] = {key: Path(value).name for key, value in summary["outputs"].items()}
        summary.update({"native_expression": "microarray_expression_native.parquet", "expression_scale": "accession_zscore_ddof1", "native_scale": "deposited_log2_RMA", "probe_mapping_policy": "exclude_ambiguous_unmapped_invalid; uppercase_only; average_retained_probes_on_deposited_scale", "cohort_role": "fibroblast_sensitivity_only", "n_recorded_source_groups": 10, "automatically_registered_for_training": False})
        (stage / "microarray_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        report = validate_release(stage, metadata, native)
        (stage / "validation.json").write_text(json.dumps(report, indent=2) + "\n")
        receipt = {"release": "GSE145725/v1", "inputs": {name: {"sha256": value, "bytes": (source_dir/name).stat().st_size} for name, value in INPUTS.items()}, "code": {name: sha256(Path(__file__).parent / name) for name in ["ingest_gse145725.py", "preprocess_microarray.py", "expression_processing.py", "gene_modules.py"]}, "outputs": {path.name: {"sha256": sha256(path), "bytes": path.stat().st_size} for path in sorted(stage.iterdir()) if path.is_file()}, "software": {"numpy": np.__version__, "pandas": pd.__version__, "pyarrow": pa.__version__}}
        (stage / "manifest.json").write_text(json.dumps(receipt, indent=2) + "\n")
        # Recheck sources before publishing. Original/archived artifacts are never written.
        verify_sources(source_dir)
        output.mkdir(parents=True, exist_ok=True)
        for path in stage.iterdir():
            shutil.move(str(path), output / path.name)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=PACKAGE / "source")
    parser.add_argument("--out-dir", type=Path, default=PACKAGE)
    args = parser.parse_args()
    print(json.dumps(build(args.source_dir, args.out_dir), indent=2))


if __name__ == "__main__":
    main()
