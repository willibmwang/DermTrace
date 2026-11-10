#!/bin/bash
#SBATCH --job-name=melanoma_trim_O10
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=03:00:00
#SBATCH --output=cutadapt_O10_%j.out
#SBATCH --error=cutadapt_O10_%j.err

set -euo pipefail

module load cutadapt/5.2

ROOT="${MELANOMA_ROOT:?Set MELANOMA_ROOT to the data directory}"

mkdir -p "$ROOT/trimmed_pilot_11_01_O10"
IN="$ROOT/raw_fastq"
OUT="$ROOT/trimmed_pilot_11_01_O10"

cutadapt \
  -j 2 \
  -O 10 \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  --json "$OUT/cutadapt_report.json" \
  -o "$OUT/11-01_R1.trimmed.fastq.gz" \
  -p "$OUT/11-01_R2.trimmed.fastq.gz" \
  "$IN/20088XR-11-01_S0_L001_R1_001.fastq.gz" \
  "$IN/20088XR-11-01_S0_L001_R2_001.fastq.gz"
