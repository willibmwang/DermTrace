#!/bin/bash
#SBATCH --job-name=melanoma_fastqc
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --array=0-42%8
#SBATCH --no-requeue
#SBATCH --output=fastqc_%A_%a.out
#SBATCH --error=fastqc_%A_%a.err
#
# FastQC 0.12.1 on every pair in $MELANOMA_ROOT/validated_pairs.tsv: one array task per data row.
# Manifest format (historical Rivanna file: header sample<TAB>r1_path<TAB>r2_path + 43 rows). Same rules as
# read_validated_pairs in rivanna_qc_tools.py when present (tests check parity); any defect fails every task:
#  - optional header: line 1 only, 3 fields, field 1 one of sample|sample_id|specimen|label (any case), and
#    fields 2-3 contain no '/' and do not end in .fastq/.fq[.gz]. Any other line 1 is validated as data.
#  - rows SAMPLE<TAB>R1<TAB>R2: exactly 3 fields, no surrounding whitespace, LF only, no blank lines
#  - R1/R2 absolute paths in the same directory; R1 *_R1_NNN.fastq.gz, R2 = R1 name with _R1_ -> _R2_
#  - unique labels, files and filename specimens (20088XR-NN-NN); a core-ID label (20088XR-NN-NN or NN-NN)
#    must name the same specimen as its files
# Array task i processes data row i (0-based, header excluded), so --array=0-42 covers the 43 rows.
# The whole manifest is validated in every task; the --array range must lie within its data rows.
# New runs write to $MELANOMA_ROOT/qc_runs/fastqc_all/job_<ARRAY_JOB_ID>/task_<TASK_ID>/;
# the historical $MELANOMA_ROOT/fastqc_all/ is never written.
# Submit: sbatch --export=ALL,MELANOMA_ROOT=/abs/path run_fastqc_all.sbatch   (adjust --array to the manifest)

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
require_slurm_int SLURM_ARRAY_JOB_ID
require_slurm_int SLURM_ARRAY_TASK_ID
MANIFEST="$ROOT/validated_pairs.tsv"
[ -f "$MANIFEST" ] && [ -r "$MANIFEST" ] && [ -s "$MANIFEST" ] || die "missing, empty or unreadable manifest: $MANIFEST"
if grep -q "$(printf '\r')" "$MANIFEST"; then die "manifest has CR (Windows) line endings: $MANIFEST"; fi

HDR=$(awk -F'\t' 'NR == 1 {
    h = (NF == 3 && tolower($1) ~ /^(sample|sample_id|specimen|label)$/ && ($2 $3) !~ /\// \
         && tolower($2) !~ /\.(fastq|fq)(\.gz)?$/ && tolower($3) !~ /\.(fastq|fq)(\.gz)?$/)
    print (h ? 1 : 0); exit }' "$MANIFEST")
ERRS=$(awk -F'\t' -v h="$HDR" '
  NR == 1 && h == 1 { next }   # optional header (rule at top of file)
  $0 == "" { print "line " NR ": blank line"; next }
  NF != 3  { print "line " NR ": expected 3 tab-separated fields, found " NF; next }
  {
    bad = 0
    for (i = 1; i <= 3; i++) if ($i == "" || $i ~ /^[[:space:]]/ || $i ~ /[[:space:]]$/) bad = 1
    if (bad) { print "line " NR ": empty field or surrounding whitespace"; next }
    if ($2 !~ /^\//) print "line " NR ": R1 path is not absolute: " $2
    if ($3 !~ /^\//) print "line " NR ": R2 path is not absolute: " $3
    b1 = $2; sub(/.*\//, "", b1); b2 = $3; sub(/.*\//, "", b2)
    d1 = substr($2, 1, length($2) - length(b1)); d2 = substr($3, 1, length($3) - length(b2))
    if (b1 !~ /_R1_[0-9][0-9][0-9]\.fastq\.gz$/) print "line " NR ": R1 name does not end in _R1_NNN.fastq.gz: " b1
    m = b1; k = gsub(/_R1_/, "_R2_", m)
    if (k != 1 || m != b2 || d1 != d2) print "line " NR ": R2 is not the mate of R1: " $2 " | " $3
    if (seen_s[$1]++) print "line " NR ": duplicate sample label: " $1
    if (match(b1, /^20088XR-[0-9][0-9]-[0-9][0-9]_/)) {
      core = substr(b1, 1, RLENGTH - 1)
      if (seen_c[core]++) print "line " NR ": duplicate specimen " core
      if (($1 ~ /^20088XR-[0-9][0-9]-[0-9][0-9]$/ || $1 ~ /^[0-9][0-9]-[0-9][0-9]$/) && $1 != core && $1 != substr(core, 9))
        print "line " NR ": label " $1 " is inconsistent with file specimen " core
    }
    if (seen_f[$2]++) print "line " NR ": file listed twice: " $2
    if (seen_f[$3]++) print "line " NR ": file listed twice: " $3
    n++
  }
  END { if (n == 0) print "manifest has no data rows" }' "$MANIFEST")
[ -z "$ERRS" ] || die "invalid manifest $MANIFEST:
$(printf '%s\n' "$ERRS" | head -n 20)"

N=$(awk -v h="$HDR" 'END { print NR - h }' "$MANIFEST")
ID=$((10#$SLURM_ARRAY_TASK_ID))
[ "$ID" -lt "$N" ] || die "SLURM_ARRAY_TASK_ID $ID is out of range: manifest has $N data rows (valid 0-$((N - 1)))"
if [ -n "${SLURM_ARRAY_TASK_MAX:-}" ]; then
  require_slurm_int SLURM_ARRAY_TASK_MAX
  [ "$((10#$SLURM_ARRAY_TASK_MAX))" -lt "$N" ] || die "--array range ends at $SLURM_ARRAY_TASK_MAX but the manifest has only $N data rows (max index $((N - 1)))"
fi

LINE=$(awk -v r=$((ID + 1 + HDR)) 'NR == r' "$MANIFEST")
IFS=$'\t' read -r SAMPLE R1 R2 <<< "$LINE"
require_mates "$R1" "$R2" _R1_ _R2_
load_module fastqc/0.12.1
require_version "FastQC v0.12.1" fastqc --version

echo "Sample: $SAMPLE"
echo "R1: $R1"
echo "R2: $R2"

make_run_dir fastqc_all "job_${SLURM_ARRAY_JOB_ID}/task_${ID}"
run_info tool_version "FastQC v0.12.1" manifest "$MANIFEST" manifest_sha256 "$(sha256_of "$MANIFEST")" \
  manifest_data_rows "$N" manifest_header_line "$HDR" sample "$SAMPLE" input_r1 "$R1" input_r1_bytes "$(size_of "$R1")" \
  input_r2 "$R2" input_r2_bytes "$(size_of "$R2")"

fastqc \
  --threads 2 \
  --outdir "$RUN_DIR" \
  "$R1" "$R2"

mark_complete
echo "FASTQC_COMPLETE: $SAMPLE"
