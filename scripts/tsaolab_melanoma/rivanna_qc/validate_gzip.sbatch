#!/bin/bash
#SBATCH --job-name=melanoma_gzip_qc
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=gzip_job_%j.out

set -euo pipefail

ROOT="${MELANOMA_ROOT:?Set MELANOMA_ROOT to the data directory}"
DATA="$ROOT/raw_fastq"
RESULTS="$ROOT/gzip_full_validation.tsv"

export DATA

python3 - <<'PY'
from pathlib import Path
root = Path(__import__("os").environ["MELANOMA_ROOT"])
files = sorted((root / "raw_fastq").glob("*.fastq.gz"))
with (root / "gzip_validation_input.txt").open("w") as f:
    for path in files:
        f.write(path.name + "\n")
print("Files scheduled:", len(files))
PY

: > "$RESULTS"

cat "$ROOT/gzip_validation_input.txt" |
  xargs -P 4 -I '{}' bash -c '
    file="$1"
    if gzip -t "$DATA/$file" 2>/dev/null; then
      printf "PASS\t%s\n" "$file"
    else
      printf "FAIL\t%s\n" "$file"
    fi
  ' _ '{}' >> "$RESULTS"

sort -o "$RESULTS" "$RESULTS"

echo "=== FULL GZIP VALIDATION ==="
echo "PASS: $(grep -c $'^PASS\t' "$RESULTS" || true)"
echo "FAIL: $(grep -c $'^FAIL\t' "$RESULTS" || true)"

echo
echo "=== FAILED FILES ==="
grep $'^FAIL\t' "$RESULTS" || true

echo
echo "VALIDATION_COMPLETE"
