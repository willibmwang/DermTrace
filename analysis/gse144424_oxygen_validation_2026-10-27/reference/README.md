# AC16 marker-direction reference

`human_AC16_marker_contrasts.csv` is the compact, unchanged snapshot of the descriptive human AC16 analysis used to specify the eight expected positive hypoxia directions before inspecting GSE144424 donor effects. It contains eight markers × four recorded oxygen arms (32 rows). The full earlier exploratory analysis is not required to reproduce these reference values.

The source expression matrix and sample metadata are tracked in this repository:

- [Human gene-symbol TPM matrix](../../../data/processed/tsaolab_o2/tsaolab_o2_human_tpm_symbol.tsv.gz)
- [Human sample metadata](../../../data/processed/tsaolab_o2/tsaolab_o2_human_sample_metadata.csv)

[provenance.json](provenance.json) records their repository-relative paths, SHA-256 hashes, source commit, the snapshot hash, sample identities, formula and verification tolerances. No external local paths are needed.

## Profiles and calculation

Select the 12 profiles with `species=human` and `cell_line=AC16`. All have `quality_flag=ok`; all are retained.

| Recorded oxygen arm | Sample IDs | Profiles |
|---|---|---:|
| 20% throughout | 1030390-13, 1030390-14, 1030390-15 | 3 |
| 20%→1% | 1030390-16, 1030390-17, 1030390-18 | 3 |
| 1% throughout | 1030390-19, 1030390-20, 1030390-21 | 3 |
| 1%→20% | 1030390-22, 1030390-23, 1030390-24 | 3 |

The markers, in snapshot order, are CA9, EGLN3, SLC2A1, LDHA, PDK1, VEGFA, BNIP3 and NDRG1. `mean_tpm` is the arithmetic mean over the three profiles for that marker and arm. `log2_ratio_mean_tpm_plus1_vs_20pct` is:

```text
log2((arm mean TPM + 1) / (20% reference mean TPM + 1))
```

The reference arm is zero by construction. Every marker is present in the source matrix. All 32 mean/contrast pairs were recomputed from the tracked inputs and agree with the snapshot at `rtol=1e-12`, `atol=1e-12` (maximum absolute differences: 2.85×10⁻¹⁴ TPM and 1.12×10⁻¹⁶ in the log2 contrast).

## Interpretation limits

These 12 RNA profiles represent one cell-line identity, not 12 independent donors. Independent culture-preparation identities, exposure durations and matched baseline/treated pairing are not established in the inspected metadata. The 20% arm is an assigned incubator condition, not measured intracellular oxygen.

The table gives descriptive expression contrasts, without differential-expression significance tests or claims of cardiac injury, functional recovery or causality. GSE144424 uses donor-derived human iPSC cardiomyocytes, a 10% reference and normalized counts; only the marker directions are compared across the two analyses. Their effect-size scales are different. The snapshot does not pool species, train or score a model, or define or change a frozen benchmark or split.
