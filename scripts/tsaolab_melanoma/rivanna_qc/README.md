# Greek Melanoma RNA-seq Analysis

Reproducible preprocessing, quality control, and validation
workflows for a Greek melanoma sequencing cohort.

## Project status

Preprocessing validation in progress. No final production
trimming configuration or downstream molecular model has
been approved.

## Dataset

- 48 apparent paired-end sample sets in the original inventory
- 43 complete, gzip-valid R1/R2 pairs
- 8 FASTQ files failed integrity validation
- Illumina paired-end reads, 2 x 151 bp
- FFPE bulk RNA-seq reported by collaborator; source metadata
  and library protocol require independent confirmation

The original sequencing data are not included in this repository.

## Completed workflows

- FASTQ file inventory and gzip validation
- FastQC assessment of 43 valid pairs
- Cutadapt adapter-trimming pilots with minimum overlaps
  of 5 and 10 bases
- Post-trimming FastQC for the O5 pilot
- Investigation of adapter-like and poly-G artifacts

## Current findings

On sample 11-01, the O5 pilot removed approximately 14.96%
of sequenced bases and retained all read-pair records.

A subset of reads showed nearly empty trimmed R1 sequences
and predominantly poly-G R2 mates.

These observations support a technical artifact hypothesis,
but do not independently validate all trimming boundaries.

## Software

The five legacy QC scripts use only:

- FastQC 0.12.1 (`module load fastqc/0.12.1`; the scripts stop
  unless `fastqc --version` prints exactly `FastQC v0.12.1`)
- Cutadapt 5.2 (`module load cutadapt/5.2`; the scripts stop
  unless `cutadapt --version` prints exactly `5.2`)
- GNU coreutils, `gzip`, `awk` and `find`

Samtools, Bowtie2 and STAR were used in separate exploratory
work. No script in this directory calls them.

These are historical research-validation workflows with
configurable data paths. The executed versions had the
Rivanna paths written into them. The committed versions take
`MELANOMA_ROOT` and write only to new output directories (see
below). They have been tested locally with stub tools but not
run on Rivanna.

## Reproducibility and data access

Raw human sequencing data and patient-level metadata are
excluded. Data redistribution permissions have not been
established.

The repository will eventually provide configurable
pipelines, software requirements, validation procedures,
and authorized data-access instructions.

## Next steps

1. Validate trimming boundaries and paired-end filtering.
2. Confirm library protocol and sample metadata.
3. Review existing collaborator RNA-seq quantification.
4. Establish a validated production preprocessing pipeline.

## Running on an HPC system

These scripts require a Slurm-based HPC environment with the
appropriate bioinformatics software installed.

Set MELANOMA_ROOT to the directory containing the sequencing
data and the validated-pairs manifest.

Example on UVA Rivanna:

    export MELANOMA_ROOT=/scratch/USERNAME/greek_melanoma

Submit a pilot job:

    sbatch --export=ALL,MELANOMA_ROOT="$MELANOMA_ROOT" \
      scripts/tsaolab_melanoma/rivanna_qc/run_cutadapt_pilot.sbatch

If your account has no default allocation, add
`--account=<allocation>` to the `sbatch` command. The scripts
do not set an account.

Slurm logs are written relative to the job submission directory.

For the FastQC array, the historical configuration is
0-42 (43 valid sample pairs). Other datasets require adjusting
the array range to match their validated-pairs manifest. A
task whose ID, or an array whose upper bound, falls outside
the manifest's data rows stops before doing any work.

### Output protection (legacy QC scripts)

New runs never write to the historical outputs:

- `gzip_full_validation.tsv`
- `gzip_validation_input.txt`
- `trimmed_pilot_11_01/`
- `trimmed_pilot_11_01_O10/`
- `fastqc_all/`
- `fastqc_trimmed_pilot/`

Each run writes only to a new directory:

| Script | New output directory under `$MELANOMA_ROOT` |
|---|---|
| `validate_gzip.sbatch` | `qc_runs/gzip_validation/job_<SLURM_JOB_ID>/` |
| `run_cutadapt_pilot.sbatch` | `qc_runs/cutadapt_pilot_O5/job_<SLURM_JOB_ID>/` |
| `run_cutadapt_pilot_O10.sbatch` | `qc_runs/cutadapt_pilot_O10/job_<SLURM_JOB_ID>/` |
| `run_fastqc_all.sbatch` | `qc_runs/fastqc_all/job_<SLURM_ARRAY_JOB_ID>/task_<TASK_ID>/` |
| `run_fastqc_trimmed_pilot.sbatch` | `qc_runs/fastqc_trimmed_pilot/job_<SLURM_JOB_ID>/` |

How the protection works:

- **Checks run first.** All input checks finish before any
  output is created.
- **The last directory is created atomically.** It is made with
  one `mkdir`, with no `-p`. If it already exists, for example
  because of a repeated job ID or a concurrent launch, the job
  stops and changes nothing. Shared parent directories can be
  created concurrently without conflict. A parent that is a
  symlink is refused.
- **There is no overwrite option.** Redirections inside the run
  directory use `noclobber`, and jobs are submitted with
  `--no-requeue`.
- **Completion is marked.** Each run directory gets
  `run_info.tsv` (script SHA-256, host, Slurm IDs, tool
  version, input paths and sizes, manifest SHA-256), then
  `COMPLETE` on success. A run that fails after its directory
  exists leaves an `INCOMPLETE` marker. That directory is never
  reused.
- **FASTQs are barely read during validation.** Checks read
  the gzip magic bytes and the first decompressed line of each
  file only. `validate_gzip.sbatch` is the only script that
  reads files in full, because that is its purpose.

**Input checks:**

- `MELANOMA_ROOT` must be an absolute directory.
- The Slurm job ID must be numeric.
- R1/R2 must be in the same directory with mate names
  (`_R1_`/`_R2_` for raw files, `_R1.`/`_R2.` for trimmed
  ones), be non-empty gzip files, and have the same first read
  name.
- `validated_pairs.tsv`: see **validated_pairs.tsv format** below.
- Cutadapt parameters are exactly those of the executed
  pilots.
- `run_fastqc_trimmed_pilot.sbatch` reads the historical O5
  pilot by default. Set `PILOT_DIR` to a completed
  `qc_runs/cutadapt_pilot_O5/job_<id>` to use a new pilot; a
  run without `COMPLETE` is refused.

### validated_pairs.tsv format

**The historical Rivanna file** (`$MELANOMA_ROOT/validated_pairs.tsv`)
has 44 lines: the header `sample<TAB>r1_path<TAB>r2_path`, then
43 specimen rows. Every line has three tab-separated fields,
there are no CRLF line endings, all specimen paths are absolute,
and there are no mate-pairing or duplicate-label errors. This
was checked on Rivanna by the user; these scripts have not yet
parsed the real file. Historical FastQC has 86 ZIP
reports (43 pairs), including both reports for 20088XR-11-01.

`run_fastqc_all.sbatch` applies the rules below. Any defect
rejects the whole manifest, and no row is ever skipped. The
separately maintained quantification helper
(`rivanna_qc_tools.py`, `read_validated_pairs`) applies the same
rules. When that file is present, `ParserParity` in the legacy
tests checks that the two parsers agree. Otherwise those tests
are skipped with a stated reason.

- **Header.** It is optional and may only be line 1. It must
  have three fields, the first must be `sample`, `sample_id`,
  `specimen` or `label` (any case), and fields 2–3 must contain
  no `/` and not end in `.fastq`/`.fq` (with or without `.gz`).
  Any other first line is validated as data, so a malformed first
  row is reported, not taken for a header. Headerless manifests
  are still accepted.
- **Rows.** `SAMPLE<TAB>R1<TAB>R2`: exactly three fields, no
  empty or whitespace-padded fields, no extra columns (a
  trailing tab counts as one), LF line endings (any CR rejects
  the file), and no blank lines.
- **Mates.** R1 is named `*_R1_NNN.fastq.gz`. R2 must be the same
  name with its single `_R1_` replaced by `_R2_`, in the same
  directory. Swapped mates, mates from two specimens and mates
  in two directories are rejected.
- **Duplicates.** Sample labels, file paths and specimens must
  each be unique. The specimen is the `20088XR-NN-NN` prefix of
  the R1 file name, so two different file pairs for the same
  specimen (for example, a second copy of 11-37) reject the
  manifest. A label written as a core ID (`20088XR-NN-NN` or
  `NN-NN`) must name the same specimen as its files. Other
  labels (e.g. `GM01`) are not cross-checked here. The
  quantification helper checks them against the sample sheet.
- **Paths.** R1 and R2 must be absolute. They are used
  verbatim: not normalised, not resolved through symlinks, and
  not checked against `MELANOMA_ROOT`. `run_fastqc_all.sbatch`
  reads the listed files directly. (The quantification helper
  matches rows to the sample sheet by **file name**, then
  requires each listed path to be the same file as
  `RAW_DIR/<file name>`.)
- **Array mapping.** Task *i* reads data row *i* (0-based, header
  excluded), so `--array=0-42` covers the 43 historical rows. A
  task ID or array upper bound beyond the data rows stops the job
  before any output is written.

### Tests (legacy QC scripts)

    python3 -m unittest discover -s scripts/tsaolab_melanoma/rivanna_qc/tests \
      -p 'test_legacy_*.py' -v

The tests use the standard library only, with synthetic
FASTQs and stub `module`, `cutadapt` and `fastqc`. They cover:

- `bash -n` syntax
- identical guard blocks across the scripts
- the exact Cutadapt arguments
- the pinned FastQC version
- every refusal path, with the whole fixture tree byte-identical
  afterwards
- concurrent launches with the same job ID: exactly one
  succeeds
- historical outputs unchanged after successful runs

They do not show how the real Cutadapt, FastQC or Slurm
behave on Rivanna.

## Reference, fastp container and specimen quantification

Scripts that rebuild the collaborator's expression pipeline on
Rivanna, one specimen at a time:

- fastp `--detect_adapter_for_pe -l 25`
- `salmon quant -l A --validateMappings --gcBias --seqBias`
  against GENCODE v46

These are the collaborator's parameters, taken from
`scripts/tsaolab_melanoma/run_quant.sh`. That script is not
modified. The driver-reference minimap2 step is not included.

| Script | Purpose |
|---|---|
| `config.example.env` | Every path, module, version and checksum. Copy it outside the repository and edit it |
| `verify_environment.sh` | Checks module versions, Python ≥ 3.9, curl/gzip, and the fastp SIF if present. Writes a TSV report |
| `fetch_gencode_v46.sh` | Downloads `gencode.v46.transcripts.fa.gz` and verifies it against the release `MD5SUMS` and the pinned MD5 |
| `build_salmon_index.sbatch` | Builds the index with `--gencode -k $SALMON_KMER`, transcriptome only, no decoys. Adopts it only if all six fingerprints match |
| `verify_reference_fingerprint.sh` | Reports MATCH / MISMATCH / UNAVAILABLE for SeqHash, NameHash, SeqHash512, NameHash512, DecoySeqHash and DecoyNameHash of an index or a `meta_info.json` |
| `provision_fastp_container.sh` | Pulls the pinned fastp image and verifies its SIF SHA-256 and `fastp --version`. A pull whose bytes differ is adopted only with `FASTP_ACCEPT_NEW_SIF=1` |
| `quantify_specimen.sbatch` | Runs fastp and Salmon for one eligible specimen, identified by `SPECIMEN` (sample_id, core_id or short id such as `11-01`) |
| `rivanna_qc_tools.py` | Standard-library helpers: `provenance`, `resolve-specimen`, `verify-fingerprint`, `validate-quant-dir`, `compare-quant`, `audit-historical` |
| `tests/run_local_tests.sh` | Local tests with synthetic fixtures, stub tools and an Lmod-like module stub. Re-runs itself under `env -i` so a login shell's Lmod cannot leak in |

**Specimen eligibility.** A specimen is quantified only if all
of these hold:

- It is in the sample sheet.
- It has exactly one pair in `validated_pairs.tsv` with the same
  R1/R2 filenames and a consistent sample label. If that file
  gives absolute paths, they must be the same files as in
  `RAW_DIR`.
- It is not excluded. Built-in exclusions are 20088XR-11-37
  (provenance adjudication pending) and the corrupt pairs
  11-43, 44, 46, 47 and 48. A built-in exclusion can be lifted
  only with `LIFT_EXCLUSIONS="CORE_ID=<adjudication reference>"`,
  and the reference is recorded in provenance.

From r3, `validated_pairs.tsv` is read with the strict rules in
**validated_pairs.tsv format** above. One defective row rejects the
whole manifest, so no specimen is eligible until it is fixed. The
GM01-qualified r2 used a more lenient parser that judged rows one at
a time. For the historical 43-row layout, both give the same 42
eligible specimens and identical resolved inputs. The r3 parser passed
its qualification on Rivanna; see `QUALIFICATION_STATUS.md` and
`R3_PARSER_CANDIDATE.md`.

**Safeguards (all scripts fail closed):**

- Required variables, input files and checksums are checked
  before any work starts.
- Existing outputs are never overwritten. Directories are
  renamed under a lock, and the destination is re-checked
  immediately before the rename. Files are adopted through a
  hard link, which fails if the destination exists.
- Locks (`.<id>.lock` and `<index>.lock`) block concurrent jobs
  on the same specimen or index. A stale lock left by a killed
  job must be removed by hand, after checking it.
- Raw FASTQs are only read. `RAW_DIR` and `QUANT_ROOT` must be
  real directories (not symlinks), and neither may sit inside
  the other. Filenames and sample IDs containing `/` or a
  leading `.` are rejected.
- Work happens in `*.partial*` paths. Before adoption, a
  specimen output must pass structural checks:
  - quant.sf rows = `num_valid_targets` = `EXPECT_NUM_TARGETS`
  - library format = `EXPECT_LIBTYPE`
  - Salmon's processed fragments = fastp's retained reads / 2
  - all six fingerprints match
  - provenance has `status=complete`
- Partial or incomplete directories are rejected by
  `validate-quant-dir` and `compare-quant`.
- Each step purges its modules and loads only its own:
  - fastp runs with `APPTAINER_MODULES`
  - Salmon runs with `SALMON_MODULES`

  On Rivanna, loading Apptainer swaps bzip2/1.0.6 for 1.0.8, so
  co-loading would change Salmon's runtime. Each step's module
  set is recorded, and `verify_environment.sh` reports any swaps
  caused by co-loading.
- Provenance labels every setting in a `configuration` block as
  one of:
  - `verified-collaborator` (from committed outputs or scripts)
  - `candidate`
  - `inferred`
  - `observed-this-run`
  - `user-reported`

  The record also includes script, library and tool hashes, the
  git commit (or the bundle manifest hash), the Slurm job ID,
  input checksums and parameters.
- A SIF SHA-256 identifies one built artifact, not the software.
  Independent pulls of the same tag may differ byte for byte.
  Each run therefore pins the exact SIF it used.

### Typical use

    cp "$RQC_DIR/config.example.env" "$MELANOMA_ROOT/rqc_config.env"   # then edit
    source "$MELANOMA_ROOT/rqc_config.env"
    bash   "$RQC_DIR/verify_environment.sh"
    bash   "$RQC_DIR/fetch_gencode_v46.sh"
    sbatch --export=ALL "$RQC_DIR/build_salmon_index.sbatch"
    bash   "$RQC_DIR/provision_fastp_container.sh"
    sbatch --export=ALL,SPECIMEN=GM01 "$RQC_DIR/quantify_specimen.sbatch"

To compare a reproduced specimen with the committed
collaborator output (files extracted with `git show`):

    python3 "$RQC_DIR/rivanna_qc_tools.py" compare-quant --specimen GM01 \
      --test-dir "$QUANT_ROOT/GM01" \
      --reference-quant GM01_committed_quant.sf.gz \
      --reference-meta GM01_committed_meta_info.json \
      --reference-lib-format GM01_committed_lib_format_counts.json \
      --reference-run-qc run_qc.tsv \
      --gene-map-fasta "$REF_DIR/gencode.v46.transcripts.fa.gz" \
      --json-out GM01_reproduction.json

**Acceptance rule `provisional-v1`.** These are engineering
thresholds, frozen before any GM01 result, and not a validated
equivalence margin. A comparison passes only if all of these
hold:

- **Identity and integrity:**
  - the test output is complete and valid
  - the specimen matches the test provenance and the reference
    `run_qc.tsv`
  - raw input read counts are identical
  - transcript IDs are identical
  - target counts are equal (253,181)
  - the library format is ISR on both sides
  - all six index fingerprints are identical
- **Quantitative:**
  - transcript-level Pearson of log1p TPM is ≥ 0.99
  - the mapping-rate difference is ≤ 1 percentage point

The report also includes, separately and without gating:
gene-level stability (Pearson, Spearman, fraction of genes with
|log2FC| > 1), the same transcript-level difference metrics,
and mapping yield relative to raw input. Any other thresholds
need a new `--rule-version` and a `--rule-justification`, and
the report carries that label.

### Historical runs (before the provenance schema)

A run made before these scripts existed has no `provenance.json`, so
`validate-quant-dir` and `compare-quant` reject it. That is intended
and is not relaxed. `audit-historical` audits such a run instead:

- it reads the directory in place, finds `quant.sf`, `meta_info.json`,
  `lib_format_counts.json`, `cmd_info.json` and the fastp JSON by name,
  and refuses ambiguous layouts unless a `--test-*` option picks a file;
- it never writes into the audited or reference directory;
- it hashes every input it reads and records provenance as absent;
- it applies the same provisional-v1 computation as `compare-quant`;
- every check is PASS, FAIL or UNAVAILABLE. The verdict is
  `HISTORICAL_CONCORDANT` (exit 0), `HISTORICAL_DISCORDANT` (exit 1) or
  `HISTORICAL_INCONCLUSIVE` (exit 2). None of them adopts the run or
  counts as an accepted reproduction.

Specimen identity comes from the input paths recorded in `cmd_info.json`
and the fastp command line, plus the raw read count.

### Test isolation

`tests/run_local_tests.sh` re-executes itself under
`env -i ... bash --noprofile --norc` with a minimal PATH. On an HPC login
node, Lmod exports a `module` function and sets `BASH_ENV`, which every
child `bash` re-sources. That replaced the stub `module` with the real
one and caused 13 failures on Rivanna. Production scripts still use the
site's real Lmod and are unchanged by this.

### Evidence status

**Local tests (developer machine, stub tools):** all pass.
Run `bash tests/run_local_tests.sh`. They cover:

- shell syntax
- every fail-closed path
- protection of raw inputs from modification
- provenance content
- comparison gates

They do **not** show how the real fastp or Salmon binaries
behave.

**Reported by the user from Rivanna, not re-verified by these
scripts:**

- **Software:** Cutadapt 5.2, FastQC 0.12.1, BBMap 39.01,
  Python 3.13.14, gcc/14.2.0 with Salmon 1.10.3, Apptainer 1.4.5.
- **GENCODE v46 transcripts:** MD5 `f1549843717b0fecfa4336fe3d5cc85f`.
- **Index:** built with `--gencode -k 31` (job 23390691). Its
  `info.json` reports k=31. All six fingerprints match the 33
  committed collaborator runs.
- **fastp:** image `quay.io/biocontainers/fastp:0.23.4--h125f33a_5`,
  SIF SHA-256
  `ef059319ad789e51ab2b40bae48f110b32f2675deeff67052b0ba2f8a21c3e01`;
  the executable reports `fastp 0.23.4`.
- **Modules:** loading Apptainer changes bzip2/1.0.6 to 1.0.8.
- **GM01 reproduction:** job 23390924 was still running at the
  time of writing. Reproduction is **not** accepted.

**Unverified:**

- **The collaborator's k-mer length.** k is not part of the
  fingerprints. Our index uses k=31, which is a candidate match.
- **The collaborator's fastp version.** 0.23.4 is a candidate.
- **Image reproducibility.** A fresh `apptainer pull` may not
  produce the same SIF bytes. The provisioning script fails
  closed unless `FASTP_ACCEPT_NEW_SIF=1` is set.
