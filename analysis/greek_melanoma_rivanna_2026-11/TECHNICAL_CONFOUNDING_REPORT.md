# Technical confounding across expression representations

All Rivanna results below were run and **reported by the user**. The
analysis is exploratory. No batch correction, feature removal or supervised
modelling was performed.

## Design

| Representation | What it is |
|---|---|
| original33 | The committed collaborator matrix |
| mixed42 | 33 collaborator columns plus 9 r3 columns |
| harmonized42 | All 42 from the qualified r2/r3 pipeline |
| harmonized33 | harmonized42 restricted to the same 33 specimens as original33 |

**PCA recipe** (the frozen exploratory recipe):
- genes: protein-coding (the biotype of the gene's first transcript, which is
  the collaborator's convention), with TPM ≥ 1 in at least 50 % of the
  analysed specimens;
- transform: log1p(TPM), per-gene mean centring, no scaling.

**Separating effects by design:**
- **original33 vs harmonized33** (same specimens, different quantification)
  tests the quantification contrast for the same measured specimens under
  this analysis recipe.
- **mixed42 vs harmonized42** (same specimens) tests how the mixed and uniformly
  quantified matrices compare under this recipe.
- **Specimen set** (original 33 vs expansion 9) is a covariate in every
  representation. Uniform quantification removes variation in
  pipeline origin across these specimens. Residual workflow and specimen
  differences can still contribute to associations.

**Associations, stability and clustering:**
- covariates: sequencing run, lane, insert size, duplication, mapping rate,
  adapter content, depth, inferred sex, quality flag and specimen set
  (eta-squared or Spearman; permutation p; BH q);
- stability: leave-one-out and 100 random 80 % subsamples;
- clustering: k-means with k = 2–4 on PC1–5.

## Quantification pipeline effects (user-reported)

| Comparison | PC1–3 score correlation | Clustering |
|---|---|---|
| original33 vs harmonized33 | All above 0.99998 | Identical assignments at every tested k |
| mixed42 vs harmonized42 | All above 0.99999 | Identical assignments at every tested k |

Only the reported correlation bounds are available in this record.
Exact PC1–3 correlations for both comparisons: **pending transcription — not yet verified**.
The source is the historical `confounding_r3/cross_representation.tsv` output.

**Interpretation.** The principal axes and tested clusterings agree closely
for these particular comparisons. PCA stability is conditional on the
protein-coding filter, expression threshold, feature matching, axis
matching and resampling procedures used. Other feature sets or leading
subspaces may behave differently.

## The leading axis (harmonized42, user-reported)

| Measure | Value |
|---|---|
| PC1 share of PCA input variance | 22.74 % |
| PC1 R² on marker scores (melanocytic, keratinocyte, immune, stromal) | 0.9346 |
| PC1 R² on technical covariates | 0.2062 |
| Unique marker R² | 0.7568 |
| Unique technical R² | 0.02838 |
| PC1 subsample stability: median \|r\|, 5th percentile | 0.9994, 0.9978 |
| PC1 after excluding the six suspect specimens | Stable |
| PC2 and PC3 | More sensitive to which specimens are included |

**Interpretation.** Marker scores explain more of the fitted PC1 scores
than the recorded technical covariates in the historical in-sample
regressions. Both use the same 42 specimens on which the axis was defined.
Marker genes also contribute to PCA, and the measured scores can reflect
tissue composition, co-expression and shared technical effects. These
associations do not establish biological causality. The technical design's
effective rank, aliasing and residual degrees of freedom need direct
examination before comparing unique R² values as evidence about separate
sources of variation.

## Clustering (harmonized42, user-reported)

| k | Result |
|---|---|
| 2 | Clusters of 24 and 18. Median subsample ARI 1.0, 5th percentile 0.6672; silhouette 0.3241 |
| 3, 4 | Produce singleton clusters |

The k = 2 partition describes cohort structure for the specified PCA
inputs and clustering settings. Stable k-means partitions can arise from
continuous variation and do not establish discrete molecular subtypes.
The singleton clusters at k = 3 and 4 remain an exploratory observation.

## Flagged specimens and sensitivity

The prespecified quality rule flagged six specimens with mapping rates
below 30 percent. Four were flagged in the collaborator metadata and two
in the expansion. The primary representation includes all six. PC1 remained
stable in the sensitivity analysis without them. PC2 and PC3 changed more.

## Confounding among technical factors

- **In the 33 committed specimens, each lane belongs to a single
  sequencing run** in the committed metadata. Lane is nested in run
  for this subset.
- **For the nine expansion specimens, run and lane were read from the first
  read name** of each R1 file. File-name lanes can differ from read-header
  lanes, so only read-header values are used.
- **How run, lane and specimen set relate** across all 42 specimens is
  tabulated on Rivanna and has not been transcribed here.
- Until that review, associations with run, lane and specimen set are
  reported **jointly** and are not attributed to any one of them.
- Two tables on Rivanna list every result, including the q-values:
  `pc_covariate_associations.tsv` (PC × covariate) and
  `covariate_relations.tsv` (covariate × covariate, including Cramér's V).
  The full-cohort design rank and aliasing have yet to be reviewed.
- Patient-level dependence is unknown. Specimen-level permutation results
  and q-values remain exploratory because exchangeability and the effective
  patient sample size have not been established.

## Independent M1 scientific review

**BLOCK for cohort execution.** The independent review assessed the separate
M1 proposal and its synthetic qualification. GRK-M1 re-review finding
M1-01 identified the historical PC1/marker recovery result as a false pass
because the original implementation passed the recovery gate even though
the marker-excluded PC1 failed to recover the original PC1 in the marker-only
synthetic control. No M1 cohort results or corrective reruns are part of
this report. The historical estimates above
are unchanged. Any future M1 release requires revision and independent
re-review, including marker-disjoint and group-aware evaluation where
feasible, rank-aware technical models and stability analyses, and explicit
failure handling. M2 remains prospective and requires separate review.

## Three kinds of evidence, kept apart

| Kind | Status |
|---|---|
| **Computational reproducibility** | Close agreement in the tested comparisons and resampling analyses |
| **Biological association** | Exploratory: PC1 is associated with marker programmes in-sample. Independent and out-of-sample tests remain unexecuted |
| **Independent biological validation** | None: no histology, orthogonal assay or external cohort has been examined |
| **Clinical validity** | Not addressed: no clinical endpoint was analysed and no classifier was trained |
