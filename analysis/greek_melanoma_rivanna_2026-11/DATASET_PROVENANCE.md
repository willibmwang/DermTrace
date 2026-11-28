# Dataset provenance: Greek melanoma (project 20088-11)

This document records where every input came from, which code produced every
derived dataset, and which identifiers prove it. **No data are published with
it.**
- Raw FASTQ, the expanded, mixed and harmonized matrices, per-specimen
  reports, containers and references stay on the cluster.
- Cluster paths are written as `$MELANOMA_ROOT/...`.
- All Rivanna results were run and **reported by the user**.

## Inputs

| Input | Identity |
|---|---|
| Sample sheet | `data/processed/tsaolab_melanoma/sample_sheets/samples_20088-11.csv` (committed): 48 specimens with R1/R2 file names and MD5s |
| Validated pairs | `$MELANOMA_ROOT/validated_pairs.tsv` on Rivanna. The user checked on Rivanna that it has the header `sample/r1_path/r2_path` and 43 rows (see `scripts/tsaolab_melanoma/rivanna_qc/README.md`). The r3 qualification report confirms core-ID labels and that every one of the 42 resolved file names and MD5s matches the sample sheet |
| Raw FASTQ | On Rivanna only. Each run checks both MD5s before reading |
| Collaborator quantification (33 specimens) | Committed: `data/processed/tsaolab_melanoma/salmon_quant/<GMnn>/` (quant.sf.gz, meta_info.json, lib_format_counts.json, cmd_info.json). All 33 share one set of index fingerprints |
| Collaborator matrices | Committed `tsaolab_melanoma_{tpm,counts}_{ensembl,symbol}.tsv.gz` and `tsaolab_melanoma_genes.tsv.gz`, built by `scripts/tsaolab_melanoma/build_dermtrace.py`. 46,592 genes and 45,873 symbols × 33 specimens. `tpm_ensembl` SHA-256 `e9442fc4…7958` |
| Transcriptome | GENCODE v46 `gencode.v46.transcripts.fa.gz`, MD5 `f1549843717b0fecfa4336fe3d5cc85f` |
| Salmon index | Built with `--gencode -k 31`; transcriptome only, no decoys. All six fingerprints (SeqHash, NameHash, SeqHash512, NameHash512, DecoySeqHash, DecoyNameHash) equal the committed collaborator values. k = 31 here; the collaborator's k is **not recorded** |
| fastp | Container `quay.io/biocontainers/fastp:0.23.4--h125f33a_5`, SIF SHA-256 `ef059319ad789e51ab2b40bae48f110b32f2675deeff67052b0ba2f8a21c3e01`. A **candidate**: the collaborator's fastp version is unknown |
| Salmon | 1.10.3 (`gcc/14.2.0 salmon/1.10.3`), `-l A --validateMappings --gcBias --seqBias`, 6 threads. fastp and Salmon run with separate module sets |

## Code identities

| Code | Identity | Role |
|---|---|---|
| r2 quantification | Bundle `rqc-src-20261124-r2-08dd537f`, manifest SHA-256 `08dd537f…b4b4` | Real-tool pilot and four-specimen validation |
| r3 quantification | Bundle `rqc-src-20261125-r3-4432a801`, manifest SHA-256 `4432a801…8f83`. Committed in `scripts/tsaolab_melanoma/rivanna_qc/` (commit `9e9d814`); `rivanna_qc_tools.py` SHA-256 `d5cb8122…2fefb` | Parser-qualified; expansion (9) and harmonization (28). Its quantification code is byte-identical to r2; only the manifest parser differs |
| Expansion tools | `procedure-20261126` (`scripts/tsaolab_melanoma/rivanna_qc/campaigns/`) | Matrix expansion, QC manifest, mixed analysis |
| Harmonization and confounding | `procedure-20261127` (same folder) | Compatibility check, harmonized matrices, confounding analysis |
| Collection recovery | `procedure-20261127r1` (same folder) | Step 17 replacement |

The [deployment hash ledger](../../scripts/tsaolab_melanoma/rivanna_qc/campaigns/deployed_hashes.tsv)
records the full deployed file hashes. Historical tarball SHA-256 identifiers
were recorded in abbreviated form as `dcda5e7d…960a` (expansion),
`50309921…780a` (harmonization and confounding), and `8a36b965…ae83` (recovery).
Published files marked `identical` match their deployed bytes. The six shell
templates replace only the site project root and absolute Python path with
`@MELANOMA_ROOT@` and `@PYTHON_BIN@` respectively. The campaign renderer can
verify the deployed hashes given the original site values. Historical plan
documents, specimen-level PCA scores and the superseded Step 17 script are
omitted from this publication.

## Rivanna jobs (user-reported)

| Purpose | Outcome |
|---|---|
| Salmon index build | Six fingerprints match |
| Historical pilot candidate run (no provenance) | Audited read-only and concordant. Remained a candidate |
| Production pilot qualification (r2) | 14/14 gates |
| Four-specimen validation (r2) | 56/56 gates |
| Nine expansion jobs (r3) | All `QUANT_OK` with complete provenance |
| Harmonization array for the remaining 28 (r3) | All tasks 0:0; 42/42 compatible; 33/33 concordance |

## Derived datasets (on Rivanna only)

| Dataset | Location | Description |
|---|---|---|
| Mixed 42 | `$MELANOMA_ROOT/expansion_r3/` | 33 collaborator columns kept byte-identical, plus 9 appended r3 columns. 47,671 Ensembl genes and 46,896 symbols. Expansion manifest, QC and provenance manifest, mixed-matrix analysis |
| **Harmonized 42 (primary)** | `$MELANOMA_ROOT/harmonized_r3/` | All 42 from r2/r3. 47,647 Ensembl genes and 46,872 symbols. 33 concordance reports, compatibility report, QC and provenance manifest |
| Confounding analysis | `$MELANOMA_ROOT/confounding_r3/` | Associations, covariate relations, technical-versus-marker R², stability, clustering, cross-representation comparisons |

Each derived folder is read-only and carries a `MANIFEST.sha256`. Each
specimen's quantification folder holds a `provenance.json` that records:
- the code SHA-256s (script, library, tools);
- the Slurm job ID;
- the input MD5s;
- the container SHA-256;
- the module sets of each step;
- the index fingerprints;
- the configuration, with each value labelled `verified-collaborator`,
  `candidate`, `inferred` or `observed-this-run`.

## Not reproduced, or unknown

- **Two collaborator settings are unknown.** The fastp version, and the Salmon
  k-mer length, are not recorded in the collaborator outputs.
- **Trimming differs slightly.** We keep a few more reads after trimming than
  the collaborator did (pilot: +0.37 %). Raw read counts are identical.
- **Salmon is not bit-reproducible across runs** (`--gcBias --seqBias` and
  multithreading). For example, two runs on the pilot specimen differ by 12 mapped fragments.
