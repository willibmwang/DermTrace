# Greek melanoma cohort: eligibility and the nine-specimen expansion

All Rivanna results below were run and **reported by the user**.
The committed sample sheet documents the delivered specimen identifiers.
This publication reports cohort-level eligibility and quality results.

## Delivered specimens and eligibility

| Count | Meaning |
|---:|---|
| 48 | Specimens in the delivered sample sheet |
| 43 | Paired FASTQ files that passed full gzip validation (`validated_pairs.tsv`) |
| **42** | **Eligible specimens:** 43 validated pairs minus one pending provenance adjudication |
| 33 | Specimens with committed collaborator quantification |
| 9 | Eligible specimens without collaborator quantification |

### The six excluded specimens

| Count | Reason | Status |
|---:|---|---|
| 1 | Duplicate source objects | Excluded by a built-in rule pending recorded provenance adjudication |
| 5 | Corrupt FASTQ pairs | Excluded after gzip validation failed |

**How eligibility was checked.** On the real `validated_pairs.tsv`, the
parser-qualified r3 implementation made the same eligibility decision as r2,
with byte-identical resolved inputs, for all 48 sample-sheet specimens.
Exactly 42 were eligible, and exactly the six above were blocked. All 42
resolved records match the committed sample sheet (core ID, file names and
both MD5s).

## Nine-specimen production expansion

**How the nine were processed:**
- code: the parser-qualified r3 implementation, with the r2 configuration
  that passed real-tool qualification;
- jobs: nine independent Slurm jobs, each taking a lock for its specimen;
- the specimen list was recomputed on Rivanna from the real files before
  submission.

**Acceptance (user-reported):**
- **All nine specimens** reached `QUANT_OK`, passed structural validation,
  and have complete r3 provenance (status, job ID, tools SHA-256, no lifted
  exclusion, configuration recorded).
- **Concordance** with collaborator output was not required, because none
  exists for these nine. Their evidence is the qualified pipeline plus
  per-specimen validation and provenance.

**Mixed 42-specimen matrices** (33 collaborator columns plus 9 appended r3
columns):

| Matrix | Rows × specimens |
|---|---|
| Ensembl genes (TPM, counts) | 47,671 × 42 |
| Gene symbols (TPM, counts) | 46,896 × 42 |

- **The original 33 columns were preserved exactly** (byte-identical text),
  and all four expanded matrices passed the integrity checks.
- **1,079 Ensembl gene rows were added.** Each is zero in all 33 committed
  specimens and non-zero in at least one new specimen, which is exactly what
  the collaborator's builder would produce. No gene was discarded, and no
  gene set was intersected.
- All nine new TPM columns sum to about 10⁶.
- **The aggregation method was validated before use.** It rebuilt all 33
  committed specimens string for string in all four matrices (about 6.1
  million values) and reproduced the committed gene table. This check ran
  locally on committed data, and on Rivanna for the five qualification specimens.

**Quality flags.** Two expansion specimens are flagged `suspect` under the inherited
collaborator rule (Salmon mapping rate below 30 %). They are flagged, not
excluded. Per-specimen values for the new specimens are not published.

## Exploratory mixed-matrix findings (user-reported)

These results motivated the harmonization in `HARMONIZATION_REPORT.md`:
- **Frozen exploratory PCA:** the nine new specimens sit, on average,
  −0.881 cohort SD along PC1 from the original 33.
- **Refitted PCA:** PC1 is associated with pipeline origin (η² = 0.147,
  exploratory permutation p = 0.0121). PC2 is associated with sequencing run
  (η² = 0.228, p = 0.0010) and with lane (η² = 0.241, p = 0.0013).

In the mixed matrix, pipeline origin is completely confounded with
specimen set (original versus expansion). These associations could not
separate a quantification effect from specimen or batch differences. The
earlier paired-pipeline maximum PCA shift of 0.0049 SD came from 5 specimens
and is empirical, not a guaranteed bound.
