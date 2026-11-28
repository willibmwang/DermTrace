#!/bin/bash
# Shared helpers for the Greek melanoma Rivanna quantification scripts. Source, do not execute.
# All helpers fail closed: any unmet precondition calls die() and the caller exits non-zero.

RQC_LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RQC_TOOLS="$RQC_LIB_DIR/rivanna_qc_tools.py"
PYTHON_BIN="${PYTHON_BIN:-python3}"

die() { echo "ERROR: $*" >&2; exit 1; }
log() { echo "[$(date '+%F %T')] $*" >&2; }

require_var() {  # require_var NAME...: each variable must be set and non-empty
  local v
  for v in "$@"; do [ -n "${!v:-}" ] || die "required variable $v is not set"; done
}

require_file() {  # require_file PATH: exists, regular, non-empty, readable
  [ -f "$1" ] && [ -s "$1" ] && [ -r "$1" ] || die "missing, empty or unreadable input: $1"
}

require_absent() {  # require_absent PATH: refuse to overwrite anything (including dangling symlinks)
  [ ! -e "$1" ] && [ ! -L "$1" ] || die "output already exists (refusing to overwrite): $1"
}

require_real_dir() {  # require_real_dir PATH: an existing directory that is not itself a symlink
  [ -d "$1" ] && [ ! -L "$1" ] || die "not a real directory (missing or symlink): $1"
}

md5_of() {
  if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | cut -d' ' -f1
  else md5 -q "$1"; fi
}

sha256_of() {
  if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | cut -d' ' -f1
  else shasum -a 256 "$1" | cut -d' ' -f1; fi
}

verify_md5() {  # verify_md5 PATH EXPECTED
  local got; got=$(md5_of "$1") || die "cannot hash $1"
  [ "$got" = "$2" ] || die "MD5 mismatch for $1: expected $2, got $got"
  log "MD5 OK $1"
}

verify_sha256() {  # verify_sha256 PATH EXPECTED
  local got; got=$(sha256_of "$1") || die "cannot hash $1"
  [ "$got" = "$2" ] || die "SHA-256 mismatch for $1: expected $2, got $got"
  log "SHA-256 OK $1"
}

load_modules() {  # load_modules "mods": purge (RQC_MODULE_PURGE=1, default) then load; no-op if RQC_SKIP_MODULE_LOAD=1
  [ "${RQC_SKIP_MODULE_LOAD:-0}" = 1 ] && return 0
  type module >/dev/null 2>&1 || die "Lmod 'module' command unavailable (set RQC_SKIP_MODULE_LOAD=1 only for local tests)"
  if [ "${RQC_MODULE_PURGE:-1}" = 1 ]; then module purge >/dev/null 2>&1 || die "module purge failed"; fi
  [ -n "${1// /}" ] || return 0
  # shellcheck disable=SC2086
  module load $1 || die "module load failed: $1"
}

run_with_modules() {  # run_with_modules "mods" CMD...: CMD runs in a subshell holding ONLY these modules
  # The module load is first tried on its own so that a load failure is reported with Lmod's own message
  # (callers often capture 2>&1 of CMD; the message is then part of the captured text they report).
  local mods=$1 load_out; shift
  load_out=$( (load_modules "$mods") 2>&1 ) \
    || die "module load failed for step modules '$mods': $(printf '%s' "$load_out" | grep -v '^[[:space:]]*$' | tail -3 | tr '\n' ' ')"
  ( load_modules "$mods" >/dev/null 2>&1; "$@" )
}

run_capture() {  # run_capture VAR "what" "mods" CMD...: VAR=first line of CMD output (stdout+stderr) under mods;
  # on a non-zero exit dies with the exit status and the last lines of output instead of exiting silently.
  local __var=$1 what=$2 mods=$3 out rc resolved; shift 3
  out=$(run_with_modules "$mods" "$@" 2>&1) && rc=0 || rc=$?
  out=$(printf '%s\n' "$out" | tr -d '\r')
  if [ "$rc" -ne 0 ]; then
    resolved=$( (load_modules "$mods" >/dev/null 2>&1; command -v "$1") 2>/dev/null || echo "not found")
    die "$what failed (exit $rc) under step modules '${mods:-none}' running '$*' [resolved=$resolved; module=$(type -t module 2>/dev/null || echo none)]: $(printf '%s' "$out" | grep -v '^[[:space:]]*$' | tail -3 | tr '\n' ' ')"
  fi
  printf -v "$__var" '%s' "$(printf '%s\n' "$out" | grep -v '^[[:space:]]*$' | head -1)"
}

module_list() {  # module_list "mods": the module set a step actually runs with (one line, ';'-separated)
  # Never fails silently under set -e/pipefail: an empty set is reported as "none"; a load failure is fatal.
  if [ "${RQC_SKIP_MODULE_LOAD:-0}" = 1 ]; then echo "module-loading-skipped"; return 0; fi
  local raw load_out
  load_out=$( (load_modules "$1") 2>&1 ) \
    || die "module load failed for step modules '$1': $(printf '%s' "$load_out" | grep -v '^[[:space:]]*$' | tail -3 | tr '\n' ' ')"
  raw=$( (load_modules "$1" >/dev/null 2>&1; module -t list 2>&1) ) || die "module_list: 'module -t list' failed for '$1'"
  raw=$(printf '%s\n' "$raw" | grep -v -e '^Currently' -e ':$' -e '^[[:space:]]*$' -e '^No modules' | sort | paste -sd';' - || true)
  echo "${raw:-none}"
}

adopt_file() {  # adopt_file SRC DEST: atomic no-clobber adoption (hard link fails if DEST exists), then drop SRC
  ln "$1" "$2" 2>/dev/null || die "cannot adopt $1 -> $2 (destination exists or hard link failed)"
  rm -f "$1"
}

acquire_lock() {  # acquire_lock DIR: atomic mkdir; a stale lock must be removed by hand after inspection
  mkdir "$1" 2>/dev/null || die "lock held or stale: $1 (another job may be running; remove only after checking)"
  printf 'job=%s host=%s pid=%s time=%s\n' "${SLURM_JOB_ID:-none}" "$(hostname)" "$$" "$(date '+%F %T')" > "$1/owner"
  RQC_LOCK="$1"
}

release_lock() {
  if [ -n "${RQC_LOCK:-}" ] && [ -d "$RQC_LOCK" ]; then rm -f "$RQC_LOCK/owner"; rmdir "$RQC_LOCK" 2>/dev/null || true; fi
}

path_is_under() {  # path_is_under CHILD_DIR PARENT_DIR: true if CHILD resolves inside (or equal to) PARENT
  local c p
  c=$(cd "$1" 2>/dev/null && pwd -P) || return 1
  p=$(cd "$2" 2>/dev/null && pwd -P) || return 1
  case "$c/" in "$p"/*) return 0 ;; *) return 1 ;; esac
}

provenance() {  # provenance OUT_JSON args...: script/git/host/slurm context is added automatically
  local out=$1; shift
  "$PYTHON_BIN" "$RQC_TOOLS" provenance --out "$out" \
    --set "script=${RQC_SCRIPT_NAME:-$(basename "$0")}" \
    --set "script_sha256=$(sha256_of "$0")" \
    --set "lib_sha256=$(sha256_of "$RQC_LIB_DIR/lib_rivanna_qc.sh")" \
    --set "tools_sha256=$(sha256_of "$RQC_TOOLS")" \
    --set "git_commit=$(git -C "$RQC_LIB_DIR" rev-parse HEAD 2>/dev/null || echo unavailable)" \
    --set "git_dirty=$(git -C "$RQC_LIB_DIR" status --porcelain -- . 2>/dev/null | grep -q . && echo yes || echo no)" \
    --set "bundle_manifest_sha256=$( [ -f "$RQC_LIB_DIR/MANIFEST.sha256" ] && sha256_of "$RQC_LIB_DIR/MANIFEST.sha256" || echo none)" \
    --set "slurm_job_id=${SLURM_JOB_ID:-none}" \
    "$@" || die "failed to write provenance $out"
  log "provenance -> $out"
}

verify_index_fingerprints() {  # verify_index_fingerprints OBSERVED JSON_OUT [extra verify-fingerprint args]
  local obs=$1 out=$2; shift 2
  require_var EXPECT_SEQ_HASH EXPECT_NAME_HASH EXPECT_SEQ_HASH512 EXPECT_NAME_HASH512 EXPECT_DECOY_HASH
  "$PYTHON_BIN" "$RQC_TOOLS" verify-fingerprint --observed "$obs" \
    --expect-seq-hash "$EXPECT_SEQ_HASH" --expect-name-hash "$EXPECT_NAME_HASH" \
    --expect-seq-hash512 "$EXPECT_SEQ_HASH512" --expect-name-hash512 "$EXPECT_NAME_HASH512" \
    --expect-decoy-seq-hash "$EXPECT_DECOY_HASH" --expect-decoy-name-hash "$EXPECT_DECOY_HASH" \
    --json-out "$out" "$@"
}
