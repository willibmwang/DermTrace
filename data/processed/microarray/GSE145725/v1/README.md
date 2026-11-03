# GSE145725: recoverable processed fibroblast expression

This package makes the **19 deposited profiles** directly available: **9 keloid
and 10 normal**, grouped into **10 recorded source lines**. It restores a study
already listed in DermTrace; it is not a new cohort or a new validation set.
The source labels are not a verified donor crosswalk.

Source: [GEO GSE145725](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145725),
[GPL16043](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL16043), and
[Kang et al., 2020](https://doi.org/10.1016/j.jid.2020.01.036).

## What the profiles measure

These are cultured primary fibroblasts, passage 3. Every sample-specific
treatment field says **20% oxygen and 5% CO2**. The growth medium includes
**0.4 micrograms/ml hydrocortisone**. `treatment=none` means no separate
experimental treatment in these deposited arrays; it does not mean a
hydrocortisone-free medium. Despite the paper's title, this deposited expression
matrix is a baseline normal-versus-keloid comparison, not an oxygen perturbation.

GEO describes 20 intended arrays and one keloid array removed after QC. The
19 deposited GSM records are retained. Keloid source labels are **1, 3, 4, 5, 6**;
source 3 has only `rep1`. Normal source labels are **1–5**, with two cultures
each. Do not invent source 2, renumber the lines, infer a rejected GSM, or pair
keloid and normal sources with the same number. `patient_id` retains the
repository's source-group convention solely to keep repeat cultures together;
`patient_id_basis` explains this limit. Donor overlap across studies is unknown.

## Expression and mapping

The GEO Series Matrix has 49,495 probes and 19 samples. Its values are already
**log2 RMA expression**, with background correction, joint quantile normalization
and median polish performed by the depositor. A deposited processed CHP header
explicitly records that scale and 19 jointly normalized arrays. Independent
parsing of its 49,495 values reproduced the Series Matrix's six-significant-digit
export for GSM4331585. Keep this scale; **do not log-transform it again**.

The pinned platform has:

| Probe disposition | Count |
|---|---:|
| Retained, one usable symbol | 47,239 |
| Multiple gene symbols: excluded | 1,457 |
| Unmapped: excluded | 688 |
| Invalid symbols: excluded | 111 |

Invalid labels include numeric spreadsheet-like strings and symbols with `@`
or `/`. They are not guessed or repaired. Single symbols are uppercased to
match the existing vocabulary convention; original symbols, Ensembl and Entrez
fields remain in `probe_mapping.tsv.gz`. Historical aliases are not updated.
Multiple retained probes for one symbol are averaged on the deposited log2
scale, following the existing microarray aggregation convention. This yields
**19,011 unique genes**. Missing genes are not supplied as zeros.

The older generic mapper can choose the first gene of a multi-symbol probe and
strip some punctuation. Therefore this conservative package differs numerically
and in gene coverage from the historical 19,420-feature representation. It has
the same artifact format; it is **not** an exact replacement of that older
matrix. The archive and frozen results are unchanged.

## Files and compatibility

| File | Meaning |
|---|---|
| `microarray_expression_native.parquet` | Canonical metadata plus sample × gene values on the deposited scale; use for new biological analyses |
| `microarray_expression_wide.parquet` | Same canonical metadata, existing module columns, and sorted gene columns; accession-wise gene z-scores |
| `microarray_X.npy` | Float32 standardized values in the exact sample/vocabulary order |
| `microarray_gene_vocab.txt` | Sorted gene columns, one per line |
| `microarray_sample_metadata.parquet` / `sample_metadata.tsv` | Standard metadata plus original GEO fields, source grouping, culture conditions and replicate IDs |
| `microarray_encoder_labels.csv` | Existing encoder label schema |
| `microarray_encoder_decoder_samples.jsonl` | Existing top-2,048-variable-gene encoder view; not the full matrix |
| `probe_mapping.tsv.gz` | Every original probe and its retained/excluded mapping status |
| `microarray_summary.json`, `validation.json`, `manifest.json` | Scale, counts, checks, code/source hashes and output hashes |
| `source/` | Pinned processed matrix, platform annotation and one processed CHP scale reference; no CEL/raw-array inputs |

The standard exports call the existing `write_expression_artifacts()` unchanged.
Its z-score uses each gene's mean and sample SD (`ddof=1`) across this accession.
That scale is for compatibility; it is not counts, TPM, or a direct fold-change
measurement. Both the depositor's RMA normalization and this compatibility
standardization use all 19 arrays. New prediction benchmarks must handle this
history explicitly and fit any new transformations only on permitted training
data. Use the native table for a new analysis rather than treating the
compatibility export as prospective test normalization.

```python
from pathlib import Path
import pandas as pd

study = Path("data/processed/microarray/GSE145725/v1")
genes = (study / "microarray_gene_vocab.txt").read_text().splitlines()
wide = pd.read_parquet(study / "microarray_expression_native.parquet")
expression = wide.set_index("sample_id")[genes]
```

The package passes a smoke test through `build_training_corpus.ModalityArtifact`
and the actual `load_artifacts()` reader. It keeps the existing
`fibroblast_sensitivity_only` role and exclusion from the primary tissue-skin
endpoint. The active registry is deliberately unchanged: the same GSM records
already exist in legacy microarray artifacts. **Select one representation;
do not concatenate this package with the old GSE145725 rows.** It is not
automatically added to training and makes no new prediction claim.

## Rebuild and verify

From the repository root, install the small pinned environment in
`keloidbench/requirements-gse145725.txt`, then run:

```sh
python keloidbench/scripts/ingest_gse145725.py --out-dir /tmp/gse145725-rebuilt
python -m unittest discover -s keloidbench/tests -p test_gse145725.py -v
```

The rebuild uses committed processed inputs and needs no network. It copies
pinned inputs into a fresh output's `source/`, checks them before and after
processing, validates all representations before publishing, and refuses an
existing release. Output paths in the summary are relative for relocation.
Parquet byte identity can depend on the pinned library versions; numerical
identity is tested independently. `manifest.json` records the actual versions.

Checks include independent aggregation of every native gene/sample value,
metadata and matrix order, labels and replicate groups, vocabulary, standard
writer arithmetic, float32 export, JSONL values, actual corpus-reader behavior,
changed-input rejection, overwrite protection, and a clean relocated rebuild.

Source URLs and scale evidence are recorded in `source/sources.json`. Cite the
original study when reusing these public data.
