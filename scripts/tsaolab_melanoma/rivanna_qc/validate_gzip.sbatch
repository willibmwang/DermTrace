#!/bin/bash
#SBATCH --job-name=melanoma_gzip_qc
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --output=gzip_job_%j.out
#SBATCH --error=gzip_job_%j.err
#
# Full gzip integrity test (gzip -t) of every *.fastq.gz in $MELANOMA_ROOT/raw_fastq. This reads every file
# in full. New runs write gzip_validation_input.txt and gzip_full_validation.tsv (PASS|FAIL<TAB>name, sorted)
# to $MELANOMA_ROOT/qc_runs/gzip_validation/job_<SLURM_JOB_ID>/. The historical
# $MELANOMA_ROOT/gzip_full_validation.tsv and gzip_validation_input.txt are never written.
# Submit: sbatch --export=ALL,MELANOMA_ROOT=/abs/path validate_gzip.sbatch

set -euo pipefail

# >>> legacy-qc-guard v1 >>>
# Identical in all five legacy QC scripts; tests/test_legacy_qc_scripts.py fails if the copies drift.
# Fail closed: every input check runs before any output exists. New outputs go only to a fresh directory,
# $MELANOMA_ROOT/qc_runs/<step>/job_<id>[/task_<n>], whose last component is created with a single atomic
# mkdir; if it already exists the job stops without writing anything. There is no overwrite option.
# Historical outputs (gzip_full_validation.tsv, gzip_validation_input.txt, trimmed_pilot_11_01*,
# fastqc_all, fastqc_trimmed_pilot) are never written. Raw FASTQs are only read; at most one decompressed
# line per file is read during validation.
RUN_DIR=
die() { echo "ERROR: $*" >&2; exit 1; }
sha256_of() { if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | cut -d' ' -f1; else shasum -a 256 "$1" | cut -d' ' -f1; fi; }
size_of() { stat -c %s "$1" 2>/dev/null || stat -f %z "$1"; }
require_root() {
  ROOT="${MELANOMA_ROOT:-}"
  [ -n "$ROOT" ] || die "set MELANOMA_ROOT to the data directory"
  case "$ROOT" in /*) ;; *) die "MELANOMA_ROOT must be an absolute path: $ROOT" ;; esac
  [ -d "$ROOT" ] || die "MELANOMA_ROOT is not a directory: $ROOT"
}
require_slurm_int() {  # require_slurm_int VAR: set, digits only
  local v="${!1:-}"
  [ -n "$v" ] || die "$1 is not set; submit this script with sbatch"
  case "$v" in *[!0-9]*) die "$1 is not a non-negative integer: $v" ;; esac
}
require_gz() {  # regular, readable, non-empty, gzip magic bytes
  [ -f "$1" ] && [ -r "$1" ] && [ -s "$1" ] || die "missing, empty or unreadable input: $1"
  [ "$(head -c 2 "$1" | od -An -tx1 | tr -d ' \n')" = 1f8b ] || die "not gzip-compressed: $1"
}
first_read_name() { { gzip -dc "$1" 2>/dev/null || true; } | head -n 1 | awk '{print $1}'; }
require_mates() {  # require_mates R1 R2 R1_TOKEN R2_TOKEN: names pair up and the first records agree
  local b1 b2 n1 n2
  b1=$(basename "$1"); b2=$(basename "$2")
  [ "$(dirname "$1")" = "$(dirname "$2")" ] || die "R1 and R2 are in different directories: $1 | $2"
  case "$b1" in *"$3"*"$3"*) die "ambiguous R1 name ($3 occurs twice): $b1" ;; *"$3"*) ;; *) die "R1 name lacks $3: $b1" ;; esac
  [ "${b1/$3/$4}" = "$b2" ] || die "R2 is not the mate of R1: $b1 | $b2"
  require_gz "$1"; require_gz "$2"
  n1=$(first_read_name "$1"); n2=$(first_read_name "$2")
  case "$n1" in @?*) ;; *) die "first record of R1 is not a FASTQ header: $1" ;; esac
  [ "${n1%/1}" = "${n2%/2}" ] || die "first read names differ between mates: $n1 | $n2"
}
load_module() {
  type module >/dev/null 2>&1 || die "Lmod 'module' command unavailable"
  module load "$1" || die "module load $1 failed"
}
require_version() {  # require_version EXPECTED CMD...: first output line must equal EXPECTED exactly
  local got
  got=$("${@:2}" 2>&1 | head -n 1) || true
  [ "$got" = "$1" ] || die "$2 version mismatch: expected '$1', got '$got'"
}
make_run_dir() {  # make_run_dir STEP REL: parents created race-safely; last component atomic, refused if present
  local d="$ROOT" comp leaf IFS=/
  for comp in qc_runs $1 $(dirname "$2"); do
    [ "$comp" = . ] && continue
    d="$d/$comp"
    mkdir "$d" 2>/dev/null || true
    [ -d "$d" ] && [ ! -L "$d" ] || die "output parent is missing, not a directory, or a symlink: $d"
  done
  leaf="$d/$(basename "$2")"
  mkdir "$leaf" 2>/dev/null || die "output directory already exists or cannot be created (refusing to overwrite): $leaf"
  RUN_DIR=$leaf
  set -o noclobber
  run_info script "$(basename "$0")" script_sha256 "$(sha256_of "$0")" \
    started_utc "$(date -u +%Y-%m-%dT%H:%M:%SZ)" host "$(hostname)" \
    slurm_job_id "${SLURM_JOB_ID:-}" slurm_array_job_id "${SLURM_ARRAY_JOB_ID:-}" \
    slurm_array_task_id "${SLURM_ARRAY_TASK_ID:-}" melanoma_root "$ROOT"
}
run_info() {  # run_info KEY VALUE ...: appends tab-separated provenance to $RUN_DIR/run_info.tsv
  while [ $# -ge 2 ]; do printf '%s\t%s\n' "$1" "$2" >> "$RUN_DIR/run_info.tsv"; shift 2; done
}
mark_complete() { date -u +%Y-%m-%dT%H:%M:%SZ > "$RUN_DIR/COMPLETE"; echo "Outputs: $RUN_DIR"; }
on_exit() {
  local rc=$?
  if [ -n "$RUN_DIR" ] && [ ! -e "$RUN_DIR/COMPLETE" ]; then printf 'exit_status\t%s\n' "$rc" >> "$RUN_DIR/INCOMPLETE"; fi
}
trap on_exit EXIT
# <<< legacy-qc-guard v1 <<<

require_root
require_slurm_int SLURM_JOB_ID
DATA="$ROOT/raw_fastq"
[ -d "$DATA" ] || die "raw FASTQ directory missing: $DATA"
NFILES=$(find "$DATA" -maxdepth 1 \( -type f -o -type l \) -name '*.fastq.gz' | wc -l | tr -d ' ')
[ "$NFILES" -gt 0 ] || die "no *.fastq.gz files in $DATA"

make_run_dir gzip_validation "job_${SLURM_JOB_ID}"
LIST="$RUN_DIR/gzip_validation_input.txt"
PARTIAL="$RUN_DIR/gzip_full_validation.partial.tsv"
RESULTS="$RUN_DIR/gzip_full_validation.tsv"
find "$DATA" -maxdepth 1 \( -type f -o -type l \) -name '*.fastq.gz' | sed 's#.*/##' | LC_ALL=C sort > "$LIST"
run_info raw_dir "$DATA" files_scheduled "$(wc -l < "$LIST" | tr -d ' ')" gzip_version "$(gzip --version 2>&1 | head -n 1)"
echo "Files scheduled: $(wc -l < "$LIST" | tr -d ' ')"

export DATA
tr '\n' '\0' < "$LIST" |
  xargs -0 -P 4 -I '{}' bash -c '
    file="$1"
    if gzip -t "$DATA/$file" 2>/dev/null; then
      printf "PASS\t%s\n" "$file"
    else
      printf "FAIL\t%s\n" "$file"
    fi
  ' _ '{}' >> "$PARTIAL"

[ "$(wc -l < "$PARTIAL")" -eq "$(wc -l < "$LIST")" ] || die "result count does not match the file list; see $PARTIAL"
LC_ALL=C sort "$PARTIAL" > "$RESULTS"
rm -f "$PARTIAL"

echo "=== FULL GZIP VALIDATION ==="
echo "PASS: $(grep -c $'^PASS\t' "$RESULTS" || true)"
echo "FAIL: $(grep -c $'^FAIL\t' "$RESULTS" || true)"
echo
echo "=== FAILED FILES ==="
grep $'^FAIL\t' "$RESULTS" || true
echo
mark_complete
echo "VALIDATION_COMPLETE"
