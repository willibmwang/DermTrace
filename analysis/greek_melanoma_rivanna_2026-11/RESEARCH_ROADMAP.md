# Research roadmap and limitations

## Completed

| Milestone | Kind of evidence |
|---|---|
| r2 real-tool pilot qualification | Computational reproducibility |
| r3 parser-specific qualification (48 eligibility decisions) | Computational reproducibility |
| Cross-specimen validation (four specimens) | Computational reproducibility |
| Nine-specimen expansion: mixed 42 | Engineering: complete provenance |
| Harmonized 42: 42/42 compatible; 33/33 concordance; 462/462 gates | Computational reproducibility |
| Technical-confounding analysis | Reproducibility, plus **exploratory** association |
| Independent M1 protocol review and local synthetic verification | Review complete. The procedure-20261128 bundle is blocked for cohort execution. See the [independent M1 review summary](TECHNICAL_CONFOUNDING_REPORT.md#independent-m1-scientific-review). |

## Prospective (designed, not run)

| Item | Status |
|---|---|
| **M1:** interpretation and robustness of the existing PCA axes | **Blocked pending revision and independent re-review.** The cohort analysis is unrun. The original thresholds were informed by this cohort. |
| **M2:** NMF programme discovery | The draft design requires separate review and execution authorisation. |

## Limitations

1. **No independent biological validation.** Every biological statement so
   far is an association within these 42 specimens. Patient linkage and
   the effective independent sample size remain unverified. No histology, orthogonal
   assay or external cohort has been examined.
2. **No clinical validity.** No clinical endpoint has been analysed, and no
   classifier has been trained or evaluated. Driver-variant calls are
   annotations and are not used as labels.
3. **Provisional acceptance rule.** Provisional-v1 (transcript Pearson ≥ 0.99;
   mapping-rate difference ≤ 1 percentage point; identity and integrity gates)
   is an engineering threshold, not a validated equivalence margin.
4. **Unknown collaborator settings.** The fastp version and the Salmon k are
   unverified. Our trimming keeps slightly more reads than the collaborator's
   (pilot: +0.37 %).
5. **Small cohort, partly aliased batch structure.** n = 42. Sequencing run,
   lane and specimen set are not independent, so associations involving them
   are reported jointly and treated as hypothesis-generating.
6. **Quality.** FFPE material gives low transcriptome mapping rates. Six
   specimens are flagged suspect (below 30 % mapping). They are kept in the
   primary representation and handled through sensitivity analyses.
7. **Unresolved specimens.** One specimen is excluded pending provenance adjudication.
   Five pairs failed gzip validation.
8. **Annotation convention.** The collaborator's `gene_type` is the biotype
   of the gene's first transcript in the FASTA, not the GENCODE gene biotype.
   No chromosome coordinates are attached, because the GTF is not deployed.
9. **Expression-derived associations.** The existing PC1 marker and
   technical R² values are in-sample. Tissue composition and technical
   variation remain plausible shared explanations. A future M1 proposal
   should evaluate a marker-disjoint axis with patient-grouped splits where
   linkage and design support them.
10. **Full-precision cross-representation PCA correlations: pending transcription — not yet verified.**
    The technical confounding report records the reported bounds.

## Next decisions (for review)

1. **Revise M1 in a new version.** Correct PC1 recovery and failure handling.
   Prespecify a marker-disjoint primary association, group-aware evaluation,
   rank-aware stability and technical models, and descriptive loading
   enrichment. Qualify a pinned numerical environment with synthetic
   negative controls and a production-entry-point rehearsal.
2. **Re-review and freeze M1 before cohort execution.** Bind the revised
   protocol, criteria, input and implementation hashes, and rehearsal
   evidence in the release record. A scientifically inconclusive result is
   an acceptable outcome. Review M2 separately and carry forward the
   applicable provenance and inference decisions.
3. **Choose an independent validation path:** histology or orthogonal data
   for the same specimens, or an external cohort under a protocol fixed in
   advance.
4. **Resolve the outstanding duplicate-object provenance**, or document
   permanent exclusion.
