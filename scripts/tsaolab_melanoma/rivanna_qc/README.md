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
