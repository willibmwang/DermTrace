# GM01 real-tool qualification of the r2 quantification scripts

> **Supersession note (added 2026-11-25; the original record below is unchanged).**
> This record describes the r2 bundle that ran on Rivanna. The repository now
> holds r3, which differs from r2 in two files. See "Deployment identity
> versus repository identity" at the end of section 1.

This is a record of evidence only. It does not change any script.

The scripts it covers were frozen in bundle `rqc-src-20261124-r2-08dd537f`.
That bundle ran on Rivanna and produced the results below.

All Rivanna results here were run and **reported by the user**. The
developer had no direct access to Rivanna. The tests run locally used stub
tools only.

## 1. What was qualified

| Item | Value |
|---|---|
| Source bundle | `rqc-src-20261124-r2-08dd537f.tar.gz`, SHA-256 `3db3cf2b66fd3e012b3c572eeb6d8b3eb63698050a732c01b8e49d69d2f44ac8` |
| Bundle manifest | `MANIFEST.sha256`, SHA-256 `08dd537fadbb3b70815db42487f843883a12b07d54762bda764688882c92b4b4` |
| Comparison rule | `provisional-v1`. Transcript-level Pearson of log1p TPM ≥ 0.99. Absolute mapping-rate difference ≤ 1 percentage point. Plus the identity and integrity gates. The thresholds were frozen before GM01 was examined and are unchanged |
| Reference | Collaborator GM01 Salmon output committed in this repository (`data/processed/tsaolab_melanoma/salmon_quant/GM01`), compared with `qc/run_qc.tsv` |

The bundle manifest lists these files with their SHA-256:

| File | SHA-256 |
|---|---|
| `build_salmon_index.sbatch` | `c33d269a94c6cd604f021dcd54821ac599920ae2fbbac21ec59aaa4dc28966a7` |
| `config.example.env` | `6abcfa07fa41f0cae224597a7e2afed3b96b2bbf315bc3b3c21e48374c5fd89f` |
| `fetch_gencode_v46.sh` | `9a31c7b2f225d72d69eee6adcc62f8ffccfbd2d1bc8dcbfeff5a8b239236c064` |
| `lib_rivanna_qc.sh` | `0f99c671520492d3113a52c53d6041de98082ebf467c0e94475a5a8b7c4bc65d` |
| `provision_fastp_container.sh` | `1be53de7e4efaafb9ce30a6f5ef450bf23f31c7c05f5255d438c8cc034061c51` |
| `quantify_specimen.sbatch` | `bc119684368b3f0d5d3cc2d27337b7d1e667d0f17a09a8e0bd91040bd8a206e9` |
| `rivanna_qc_tools.py` | `61888cb45011b4e3fd31c7a2d7c09fc654eb5362e2530dc278d64ed5b2c9aaee` |
| `tests/run_local_tests.sh` | `8c1b68ca2c854bfe5d065b13e0f2b1e2f82d5d00ff33c7227a535877da1128f7` |
| `verify_environment.sh` | `df838d6edd423ead9f4cfd7d02ddabce99b38dfbc03906b2c26b58513c390f48` |
| `verify_reference_fingerprint.sh` | `56ff67da9991ea7e31dfc8349c0d7fbff72fde47679b053e9d21d3dd7abafb91` |

These ten files are committed byte-identical to the bundle.

`README.md` is the one exception. The bundle's `README.md` (SHA-256
`c7e6e91c…3195`) was the repository README at commit 7784ce3 plus the
quantification section. The README has since gained the legacy QC sections
in commit 412ad24. So the committed README is the current README plus the
same quantification section, unchanged byte for byte. Its whole-file hash
therefore differs from the bundle's.

### Deployment identity versus repository identity (supersession note, 2026-11-25)

The sentence "These ten files are committed byte-identical to the bundle" was
true when r2 was staged. It is not true of the integrated commit, which holds
r3 (`rqc-src-20261125-r3-4432a801`, manifest SHA-256
`4432a80128f6b090266866e0c472614e0717d68d6db04bbb9fc3eec7c1fb8f83`).

| File | r2 (deployed for job 23402243) | r3 (in the repository) |
|---|---|---|
| `rivanna_qc_tools.py` | `61888cb4…aaee` | `d5cb8122…2fefb`. Only `read_validated_pairs` (with its helper and constants) changed |
| `tests/run_local_tests.sh` | `8c1b68ca…28f7` | `4778ecaa…e7f2`. One fixture adjusted; r3 parser tests added |
| The other eight files listed above | as listed | byte-identical to r2 |

- **Deployment identity.** The r2 hashes in the table above identify exactly
  what produced the GM01 results in this record. They remain the reference
  for that evidence.
- **Repository identity.** The r3 manifest identifies the committed code.
- **How r3 inherits the GM01 result.** r3 passed its own parser-specific
  qualification. It showed that, on the real manifest, every eligibility
  decision and every resolved input is identical to r2, so the quantification
  inputs and arguments are unchanged. See `QUALIFICATION_STATUS.md` and
  `R3_PARSER_CANDIDATE.md`.
- **The parity limitation (limitation 7) is resolved in r3.**

## 2. Software environment (user-reported)

| Component | Value |
|---|---|
| Salmon | 1.10.3 (`gcc/14.2.0 salmon/1.10.3`) |
| Apptainer | 1.4.5 |
| fastp | 0.23.4, as reported by the container. This is a **candidate**: the collaborator's fastp version is unknown |
| fastp container | `quay.io/biocontainers/fastp:0.23.4--h125f33a_5`. SIF SHA-256 `ef059319ad789e51ab2b40bae48f110b32f2675deeff67052b0ba2f8a21c3e01`. This hash identifies the image file, not the software version |
| Python | 3.13.14 |
| Transcriptome | GENCODE v46 transcripts, MD5 `f1549843717b0fecfa4336fe3d5cc85f` |
| Salmon index | Built with `--gencode -k 31`, transcriptome only, no decoys. All six fingerprints (SeqHash, NameHash, SeqHash512, NameHash512, DecoySeqHash, DecoyNameHash) match the committed collaborator values. This index has k=31; the collaborator's k is **not recorded and remains unverified** |
| Module isolation | fastp ran with the Apptainer modules only, and Salmon with the Salmon modules only. Loading Apptainer would otherwise switch bzip2 from 1.0.6 to 1.0.8 under Salmon |

## 3. Results (user-reported)

### 3.1 Production qualification run: job 23402243

The run used `quantify_specimen.sbatch` from r2, wrote into a fresh output
folder, and produced complete provenance.

| Measure | Value |
|---|---|
| Slurm | COMPLETED, exit 0:0, 32 min 20 s, MaxRSS 14,126,644 KB |
| Adoption | `QUANT_OK`. `validate-quant-dir`: VALID |
| provisional-v1 gates | **All 14 PASS** |
| Specimen identity, raw read count | PASS. Raw read count is identical to the reference |
| Transcript IDs, target count | Identical. 253,181 targets |
| Library orientation | ISR on both sides |
| Index fingerprints | All six identical to the reference |
| Transcript-level Pearson, log1p TPM | 0.998858 (threshold ≥ 0.99) |
| Transcript-level Spearman | 0.997946 |
| Gene-level Pearson, log1p TPM | 0.999893. Reported, not gated |
| Gene-level Spearman | 0.999421. Reported, not gated |
| Expressed transcripts with \|log2FC\| > 1 | 129 of 41,896 |
| Expressed genes with \|log2FC\| > 1 | 12 of 18,535 |
| Reads kept after trimming | 75,719,724 (reference 75,441,284) |
| Salmon fragments processed | 37,859,862 (reference 37,720,642) |
| Mapped fragments | 13,232,614 (reference 13,211,561) |
| Mapping rate | 34.95156 % (reference 35.02475 %). Difference −0.0732 percentage points (threshold ±1) |

Definitions:
- **Expressed:** TPM ≥ 1 in the candidate or the reference.
- **log2FC:** log2((candidate + 1) / (reference + 1)).
- **Gene TPM:** the sum of that gene's transcript TPMs.

### 3.2 Historical candidate run: job 23390924

This run predates the provenance schema. It was audited read-only with
`audit-historical`.

| Measure | Value |
|---|---|
| Verdict | `HISTORICAL_CONCORDANT`, exit 0. The run was not adopted, and provenance is recorded as absent, not reconstructed |
| Gates | All 12 PASS. All six fingerprints match |
| Transcript-level Pearson / Spearman | 0.998885 / 0.997995 |
| Gene-level Pearson / Spearman | 0.999895 / 0.999379 |
| Discordant expressed transcripts | 120 of 41,888 |
| Discordant expressed genes | 12 of 18,534 |
| Mapped fragments | 13,232,602 |
| Reads kept after trimming | 75,719,724: identical to production job 23402243 |

The two candidate runs differ by 12 mapped fragments out of about 13.2
million. They have the same inputs and the same configuration.

## 4. Limitations of this qualification

1. **The acceptance rule is provisional.** The provisional-v1 thresholds are
   engineering thresholds, not a validated equivalence margin. Passing them
   shows close agreement for one specimen. It does not show that the
   pipeline as a whole is reproduced.
2. **Only one specimen was tested.** Agreement across specimens, sequencing
   runs and QC profiles has not been tested yet. The cross-specimen plan
   covers GM02, GM08, GM11 and GM31.
3. **fastp version and k are unknown.** We reproduce the collaborator's
   *output* closely using a *candidate* fastp (0.23.4) and an index built
   with k=31. Neither is shown to be what the collaborator actually used.
   - The trimmed read count differs by 278,440 reads, about 0.37 % more
     reads kept than the collaborator's run.
   - So the trimming step is not reproduced exactly.
4. **Runs are not bit-for-bit reproducible.** Salmon's bias correction
   (`--gcBias --seqBias`) and multithreaded inference vary slightly between
   runs. The two candidate runs differ by 12 mapped fragments.
   - Gene-level agreement (Pearson 0.99989) is higher than transcript-level
     agreement (0.99886).
   - Many more transcripts than genes are discordant (129 versus 12).
   - Both point to expression measured per transcript being less stable
     than per gene. Which isoforms and genes are involved is reported by the
     separate discordance analysis, not by this record.
5. **Local tests prove control flow, not numbers.**
   `tests/run_local_tests.sh` uses stub tools. It shows the scripts'
   checks, isolation and provenance work, but not the numerical behaviour
   of fastp or Salmon.
6. **The historical run is not an adopted reproduction.** Job 23390924 is
   evidence of concordance only. It has no provenance record, and none was
   invented for it.
7. **Integration with the legacy QC parser parity test is unresolved.**
   - Commit 412ad24 added `tests/test_legacy_qc_scripts.py` from a separate
     workstream. Its parity cases expect `rivanna_qc_tools.read_validated_pairs`
     to follow the stricter manifest rules of `run_fastqc_all.sbatch`.
   - The qualified r2 parser is more lenient. With r2 committed, 16 of those
     parity sub-tests fail; for example, r2 accepts padded labels.
   - The stricter parser exists only as an unqualified working-tree change.
   - Adopting it changes which specimens are eligible. It should be
     committed separately and re-checked on Rivanna: run the regression suite
     and the eligibility checks of the qualification procedure, step 2.
   - It is not part of this qualification.

## 5. How to reproduce on Rivanna

The procedure is in `qualification-r2`, kept outside this repository:

1. Write the complete configuration and verify the environment.
2. Verify the reference, the container and specimen eligibility.
3. Submit one fresh GM01 job into a separate output folder.
4. Compare the result with the collaborator reference using `compare-quant`.

Each step refuses to overwrite existing configuration, outputs or reports.
Paths are set in the configuration under `$MELANOMA_ROOT`.
