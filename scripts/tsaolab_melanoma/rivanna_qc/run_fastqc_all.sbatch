#!/bin/bash
#SBATCH --job-name=melanoma_fastqc
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --array=0-42%8
#SBATCH --output=fastqc_%A_%a.out
#SBATCH --error=fastqc_%A_%a.err

set -euo pipefail

module load fastqc/0.12.1

ROOT="${MELANOMA_ROOT:?Set MELANOMA_ROOT to the data directory}"

mkdir -p "$ROOT/fastqc_all"
MANIFEST="$ROOT/validated_pairs.tsv"

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 2))p" "$MANIFEST")

IFS=$'\t' read -r SAMPLE R1 R2 <<< "$LINE"

test -n "$SAMPLE"
test -s "$R1"
test -s "$R2"

echo "Sample: $SAMPLE"
echo "R1: $R1"
echo "R2: $R2"

fastqc \
  --threads 2 \
  --outdir "$ROOT/fastqc_all" \
  "$R1" "$R2"

echo "FASTQC_COMPLETE: $SAMPLE"
