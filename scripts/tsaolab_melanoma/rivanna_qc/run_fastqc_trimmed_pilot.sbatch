#!/bin/bash
#SBATCH --job-name=melanoma_posttrim_qc
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=fastqc_trimmed_%j.out
#SBATCH --error=fastqc_trimmed_%j.err

set -euo pipefail

module load fastqc

ROOT="${MELANOMA_ROOT:?Set MELANOMA_ROOT to the data directory}"

mkdir -p "$ROOT/fastqc_trimmed_pilot"

fastqc \
  --threads 2 \
  --outdir "$ROOT/fastqc_trimmed_pilot" \
  "$ROOT/trimmed_pilot_11_01/11-01_R1.trimmed.fastq.gz" \
  "$ROOT/trimmed_pilot_11_01/11-01_R2.trimmed.fastq.gz"
