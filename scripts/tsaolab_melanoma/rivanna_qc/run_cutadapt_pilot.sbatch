#!/bin/bash
#SBATCH --job-name=melanoma_trim_pilot
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=cutadapt_pilot_%j.out
#SBATCH --error=cutadapt_pilot_%j.err

set -euo pipefail

module load cutadapt/5.2

ROOT="${MELANOMA_ROOT:?Set MELANOMA_ROOT to the data directory}"

mkdir -p "$ROOT/trimmed_pilot_11_01"
IN="$ROOT/raw_fastq"
OUT="$ROOT/trimmed_pilot_11_01"

R1="$IN/20088XR-11-01_S0_L001_R1_001.fastq.gz"
R2="$IN/20088XR-11-01_S0_L001_R2_001.fastq.gz"

test -s "$R1"
test -s "$R2"

cutadapt \
  -j 2 \
  -O 5 \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  --json "$OUT/cutadapt_report.json" \
  -o "$OUT/11-01_R1.trimmed.fastq.gz" \
  -p "$OUT/11-01_R2.trimmed.fastq.gz" \
  "$R1" "$R2"
