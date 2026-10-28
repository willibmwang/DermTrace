"""Donor-paired, exploratory validation of the fixed AC16 marker panel.

Reads the separate curated GSE144424 package and bundled AC16 reference;
never loads a classifier, benchmark, or split. Run after ingestion.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/curated/gse144424_cardiomyocyte_v1"
OUT = ROOT / "analysis/gse144424_oxygen_validation_2026-10-27"
CONTRASTS = {"B-A": ("B", "A"), "C-A": ("C", "A"), "D-A": ("D", "A"),
             "C-B": ("C", "B"), "D-B": ("D", "B")}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def autosomal_mask(chromosomes: pd.Series) -> pd.Series:
    """featureCounts can repeat a chromosome for each exon; reject mixed loci."""
    autosomes = {f"chr{i}" for i in range(1, 23)}
    return chromosomes.map(lambda value: len(set(value.split(";"))) == 1
                           and value.split(";")[0] in autosomes)


def normalize_counts(counts: pd.DataFrame, autosomal: pd.Series,
                     method: str) -> tuple[pd.DataFrame, pd.Series, int]:
    """Normalize within this study and sample selection, without gene selection by effect."""
    reference = counts.loc[autosomal].to_numpy(dtype=float)
    if not np.isfinite(reference).all() or (reference < 0).any():
        raise ValueError("Counts must be finite and nonnegative")
    if method == "median_ratio":
        positive = (reference > 0).all(axis=1)
        if not positive.any():
            raise ValueError("No all-positive autosomal normalization genes")
        logs = np.log(reference[positive])
        ratios = np.exp(logs - logs.mean(axis=1, keepdims=True))
        factors = np.median(ratios, axis=0)
        factors /= np.exp(np.log(factors).mean())
        n_reference = int(positive.sum())
    elif method == "cpm":
        factors = reference.sum(axis=0) / 1e6
        n_reference = len(reference)
    else:
        raise ValueError(f"Unknown normalization: {method}")
    if not np.isfinite(factors).all() or (factors <= 0).any():
        raise ValueError("Invalid normalization factors")
    factors = pd.Series(factors, index=counts.columns, name="size_factor")
    return counts.div(factors, axis=1), factors, n_reference


def paired_deltas(log_expression: pd.DataFrame, samples: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Pair within differentiation first; then give every donor exactly one value."""
    records = []
    if samples["sample_id"].duplicated().any():
        raise ValueError("Duplicate sample identifiers")
    if set(log_expression.columns) != set(samples["sample_id"]):
        raise ValueError("Expression/sample axes disagree")
    for (donor, experiment), group in samples.groupby(["donor_id", "experiment_id"], sort=True):
        if len(group) != 4 or set(group["condition"]) != set("ABCD"):
            raise ValueError(f"Incomplete or duplicated condition quartet: {experiment}")
        ids = group.set_index("condition")["sample_id"]
        for contrast, (after, before) in CONTRASTS.items():
            effects = log_expression[ids[after]] - log_expression[ids[before]]
            for gene, effect in effects.items():
                records.append({"donor_id": str(donor), "experiment_id": experiment,
                                "gene": gene, "contrast": contrast, "delta": float(effect)})
    experiments = pd.DataFrame(records)
    donors = experiments.groupby(["gene", "contrast", "donor_id"], sort=True)["delta"].agg(
        delta="mean", n_experiments="size").reset_index()
    if donors.groupby(["gene", "contrast"]).size().ne(samples["donor_id"].nunique()).any():
        raise ValueError("Donor aggregation changed the independent sample count")
    return experiments, donors


def bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    order = np.argsort(pvalues)
    ranked = pvalues[order] * len(pvalues) / np.arange(1, len(pvalues) + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1].clip(0, 1)
    result = np.empty(len(pvalues))
    result[order] = adjusted
    return result


def summarize_donors(donors: pd.DataFrame) -> pd.DataFrame:
    records = []
    for (gene, contrast), group in donors.groupby(["gene", "contrast"], sort=True):
        effects = group.sort_values("donor_id")["delta"].to_numpy()
        n = len(effects)
        if n < 2 or not np.isfinite(effects).all():
            raise ValueError("Insufficient finite donor contrasts")
        rng = np.random.default_rng(144424)
        boot = np.median(effects[rng.integers(0, n, size=(20_000, n))], axis=1)
        lower, upper = np.quantile(boot, [0.025, 0.975])
        positive, negative = int((effects > 0).sum()), int((effects < 0).sum())
        pvalue = (binomtest(positive, positive + negative, p=0.5, alternative="greater").pvalue
                  if positive + negative else 1.0)
        records.append({"gene": gene, "contrast": contrast, "n_donors": n,
                        "n_positive": positive, "n_zero": int((effects == 0).sum()),
                        "n_negative": negative, "mean_delta": effects.mean(),
                        "median_delta": np.median(effects), "ci_low": lower, "ci_high": upper,
                        "p_direction": pvalue if contrast == "B-A" else np.nan, "q_bh": np.nan})
    result = pd.DataFrame(records)
    primary = result["contrast"].eq("B-A")
    result.loc[primary, "q_bh"] = bh_adjust(result.loc[primary, "p_direction"].to_numpy())
    return result


def reoxygenation_proximity(donors: pd.DataFrame) -> pd.DataFrame:
    wide = donors.pivot(index=["gene", "donor_id"], columns="contrast", values="delta")
    records = []
    for gene, group in wide.groupby(level="gene", sort=True):
        for condition, hours in (("C", 6), ("D", 24)):
            hypoxia, reoxygenation = group["B-A"], group[f"{condition}-A"]
            closer = reoxygenation.abs() < hypoxia.abs()
            induced = hypoxia > 0
            records.append({"gene": gene, "reoxygenation_hours": hours, "n_donors": len(group),
                            "n_closer_to_reference": int(closer.sum()),
                            "n_hypoxia_positive": int(induced.sum()),
                            "n_closer_among_hypoxia_positive": int((closer & induced).sum()),
                            "n_still_above_reference": int((reoxygenation > 0).sum()),
                            "n_overshoot_from_hypoxia_positive": int((induced & (reoxygenation < 0)).sum())})
    return pd.DataFrame(records)


def main() -> None:
    plan = json.loads((OUT / "analysis_plan.json").read_text())
    markers = plan["markers"]
    counts = pd.read_parquet(DATA / "counts.parquet")
    samples = pd.read_csv(DATA / "samples.csv", dtype={"donor_id": str, "source_subject_id": str})
    genes = pd.read_csv(DATA / "genes.csv").set_index("Geneid").loc[counts.index]
    mapping = pd.read_csv(DATA / "marker_mapping.csv").set_index("gene_symbol").loc[markers]
    if mapping.index.has_duplicates or mapping["ensembl_gene_id"].duplicated().any():
        raise ValueError("Marker mapping is not one-to-one")
    marker_ids = mapping["ensembl_gene_id"].tolist()
    if not set(marker_ids).issubset(counts.index):
        raise ValueError("A prespecified marker is missing; report rather than silently drop it")
    autosomal = autosomal_mask(genes["Chr"])
    primary = samples.loc[samples["published_primary_eligible"]].copy()
    eligible_repeats = samples.loc[~(samples["donor_id"].eq("18852") & samples["replicate"].eq(1))].copy()
    selections = [("primary", "median_ratio", primary, 60),
                  ("repeat_averaged_qc", "median_ratio", eligible_repeats, 80),
                  ("repeat_averaged_all", "median_ratio", samples, 84),
                  ("primary_cpm", "cpm", primary, 60)]
    all_experiments, all_donors, all_summaries, all_factors, all_proximity = [], [], [], [], []
    audit = []
    for label, method, selected, n_profiles in selections:
        if len(selected) != n_profiles or selected["donor_id"].nunique() != 15:
            raise ValueError(f"Unexpected selection size for {label}")
        normalized, factors, n_reference = normalize_counts(
            counts.loc[:, selected["sample_id"]], autosomal, method)
        marker_values = normalized.loc[marker_ids].copy()
        marker_values.index = pd.Index(markers, name="gene")
        experiments, donors = paired_deltas(np.log2(marker_values + 1), selected)
        summaries = summarize_donors(donors)
        proximity = reoxygenation_proximity(donors)
        for condition in "AB":
            ids = selected.loc[selected["condition"].eq(condition), "sample_id"]
            medians = marker_values.loc[:, ids].median(axis=1)
            summaries[f"median_normalized_{condition}"] = summaries["gene"].map(medians)
        # The fixed count threshold has no equivalent meaning in CPM units.
        summaries["low_abundance_both_A_B"] = (
            ((summaries["median_normalized_A"] < 10) & (summaries["median_normalized_B"] < 10))
            if method == "median_ratio" else np.nan
        )
        factors_frame = factors.rename_axis("sample_id").reset_index()
        for frame in (experiments, donors, summaries, factors_frame, proximity):
            frame.insert(0, "normalization", method)
            frame.insert(0, "analysis", label)
        all_experiments.append(experiments)
        all_donors.append(donors)
        all_summaries.append(summaries)
        all_factors.append(factors_frame)
        all_proximity.append(proximity)
        audit.append({"analysis": label, "normalization": method, "n_profiles": len(selected),
                      "n_donors": selected["donor_id"].nunique(),
                      "n_experiments": selected["experiment_id"].nunique(),
                      "n_autosomal_reference_genes": n_reference,
                      "min_factor": factors.min(), "max_factor": factors.max()})
        if label == "primary":
            expression = marker_values.T.rename_axis("sample_id").reset_index().melt(
                id_vars="sample_id", var_name="gene", value_name="normalized_count")
            expression["log2_normalized_count_plus1"] = np.log2(expression["normalized_count"] + 1)
            expression.merge(selected, on="sample_id", validate="many_to_one").to_csv(
                OUT / "primary_donor_expression.csv", index=False)
    for name, frames in (("experiment_deltas", all_experiments), ("donor_deltas", all_donors),
                         ("marker_summary", all_summaries), ("normalization_factors", all_factors),
                         ("reoxygenation_summary", all_proximity)):
        pd.concat(frames, ignore_index=True).to_csv(OUT / f"{name}.csv", index=False)
    main_summary = all_summaries[0].query("contrast == 'B-A'").set_index("gene").loc[markers]
    tsao_path = OUT / "reference/human_AC16_marker_contrasts.csv"
    tsao = pd.read_csv(tsao_path).query("oxygen == '20%->1%'").set_index("gene").loc[markers]
    comparison = main_summary.copy()
    comparison["ac16_log2_tpm_mean_plus1_contrast"] = tsao["log2_ratio_mean_tpm_plus1_vs_20pct"]
    comparison["same_median_direction_as_ac16"] = (
        np.sign(comparison["median_delta"]) == np.sign(comparison["ac16_log2_tpm_mean_plus1_contrast"]))
    comparison.to_csv(OUT / "ac16_direction_comparison.csv")
    inputs = [DATA / name for name in ("counts.parquet", "samples.csv", "genes.csv", "marker_mapping.csv", "manifest.json")]
    inputs.extend([OUT / "analysis_plan.json", Path(__file__).resolve(), tsao_path,
                   OUT / "reference/provenance.json"])
    outputs = sorted(OUT.glob("*.csv"))
    result = {"analysis": plan["version"], "scope": "exploratory; no predictive fit or benchmark change",
              "normalization_audit": audit,
              "input_sha256": {str(p.relative_to(ROOT)): sha256(p) for p in inputs},
              "output_sha256": {str(p.relative_to(ROOT)): sha256(p) for p in outputs},
              "primary_direction_support_genes": main_summary.index[(main_summary["q_bh"] <= .05)
                  & (main_summary["median_delta"] > 0)].tolist(),
              "primary_same_direction_genes": comparison.index[comparison["same_median_direction_as_ac16"]].tolist()}
    (OUT / "analysis_manifest.json").write_text(json.dumps(result, indent=2) + "\n")
    print(main_summary[["n_positive", "n_donors", "median_delta", "ci_low", "ci_high", "q_bh"]].to_string())


if __name__ == "__main__":
    main()
