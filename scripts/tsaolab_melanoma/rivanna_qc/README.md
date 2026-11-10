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

- FastQC 0.12.1
- Cutadapt 5.2
- Samtools 1.21
- Bowtie2 2.5.4
- STAR 2.7.11b

These are historical research-validation workflows with
configurable data paths. They have passed Bash syntax checks
but have not yet been execution-tested after refactoring.

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

Slurm logs are written relative to the job submission directory.

For the FastQC array, the historical configuration is
0-42 (43 valid sample pairs). Other datasets require adjusting
the array range to match their validated-pairs manifest.

These scripts have passed Bash syntax checks but have not
yet been execution-tested following portability refactoring.
