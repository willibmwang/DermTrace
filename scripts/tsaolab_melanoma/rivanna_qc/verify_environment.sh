#!/bin/bash
# Verify the Rivanna software environment used by the Greek melanoma QC and quantification scripts.
# Each tool is checked in its own subshell after loading its modules, so module conflicts cannot mask failures.
# Writes a TSV report (and provenance JSON when PROVENANCE_DIR is set); exits non-zero if any check fails.
#
# Usage:  source config.env && bash verify_environment.sh [report.tsv]
# Override expected versions with EXPECT_CUTADAPT, EXPECT_FASTQC, EXPECT_BBMAP, EXPECT_SALMON, EXPECT_APPTAINER,
# EXPECT_PYTHON_MIN (minimum major.minor), and module names with CUTADAPT_MODULES, FASTQC_MODULES, BBMAP_MODULES,
# SALMON_MODULES, APPTAINER_MODULES. BBMAP_VERSION_CMD selects the BBMap version command (default bbversion.sh).
set -euo pipefail
. "$(dirname "$0")/lib_rivanna_qc.sh"

EXPECT_CUTADAPT=${EXPECT_CUTADAPT:-5.2}
EXPECT_FASTQC=${EXPECT_FASTQC:-0.12.1}
EXPECT_BBMAP=${EXPECT_BBMAP:-39.01}
EXPECT_SALMON=${EXPECT_SALMON:-1.10.3}
EXPECT_APPTAINER=${EXPECT_APPTAINER:-1.4.5}
EXPECT_PYTHON_MIN=${EXPECT_PYTHON_MIN:-3.9}
CUTADAPT_MODULES=${CUTADAPT_MODULES:-cutadapt/5.2}
FASTQC_MODULES=${FASTQC_MODULES:-fastqc/0.12.1}
BBMAP_MODULES=${BBMAP_MODULES:-bbmap/39.01}
SALMON_MODULES=${SALMON_MODULES:-gcc/14.2.0 salmon/1.10.3}
APPTAINER_MODULES=${APPTAINER_MODULES:-apptainer}
BBMAP_VERSION_CMD=${BBMAP_VERSION_CMD:-bbversion.sh}

if [ -n "${1:-}" ]; then REPORT=$1
else require_var PROVENANCE_DIR; mkdir -p "$PROVENANCE_DIR"; REPORT="$PROVENANCE_DIR/environment_$(date +%Y%m%dT%H%M%S)_$$.tsv"; fi
require_absent "$REPORT"
printf "check\tmodules\texpected\tobserved\tstatus\n" > "$REPORT"
nfail=0

check() {  # check LABEL MODULES EXPECTED_SUBSTRING COMMAND...
  local label=$1 mods=$2 want=$3; shift 3
  local out status
  if ! ( load_modules "$mods" ) >/dev/null 2>&1; then
    out="module load failed: $mods"
  else  # module chatter is discarded; only the tool's own first non-empty output line is compared
    out=$( (load_modules "$mods" >/dev/null 2>&1; "$@") 2>&1 | tr -d '\r' | grep -v '^[[:space:]]*$' | head -1 ) || true
  fi
  if [ -n "$out" ] && [[ "$out" == *"$want"* ]]; then status=PASS; else status=FAIL; nfail=$((nfail + 1)); fi
  printf "%s\t%s\t%s\t%s\t%s\n" "$label" "${mods:-none}" "$want" "${out:-<no output>}" "$status" >> "$REPORT"
  echo "$status  $label: ${out:-<no output>} (expected $want)" >&2
}

check cutadapt  "$CUTADAPT_MODULES"  "$EXPECT_CUTADAPT"  cutadapt --version
check fastqc    "$FASTQC_MODULES"    "v$EXPECT_FASTQC"   fastqc --version
check bbmap     "$BBMAP_MODULES"     "$EXPECT_BBMAP"     "$BBMAP_VERSION_CMD"
check salmon    "$SALMON_MODULES"    "salmon $EXPECT_SALMON" "${SALMON_BIN:-salmon}" --version
check apptainer "$APPTAINER_MODULES" "$EXPECT_APPTAINER" "${APPTAINER_BIN:-apptainer}" --version

pyv=$("$PYTHON_BIN" -c 'import sys; print("%d.%d.%d" % sys.version_info[:3])' 2>/dev/null || echo none)
if [ "$pyv" != none ] && "$PYTHON_BIN" -c "import sys; sys.exit(0 if sys.version_info[:2] >= tuple(map(int, '$EXPECT_PYTHON_MIN'.split('.'))) else 1)"; then
  pst=PASS; else pst=FAIL; nfail=$((nfail + 1)); fi
printf "python\tnone\t>=%s\t%s\t%s\n" "$EXPECT_PYTHON_MIN" "$pyv" "$pst" >> "$REPORT"
echo "$pst  python: $pyv (expected >= $EXPECT_PYTHON_MIN)" >&2

for tool in curl gzip; do
  if command -v "$tool" >/dev/null; then st=PASS; else st=FAIL; nfail=$((nfail + 1)); fi
  printf "%s\tnone\tpresent\t%s\t%s\n" "$tool" "$(command -v "$tool" || echo missing)" "$st" >> "$REPORT"
done

if [ -n "${FASTP_SIF:-}" ] && [ -e "$FASTP_SIF" ]; then
  require_var FASTP_SIF_SHA256 FASTP_EXPECTED_VERSION
  got=$(sha256_of "$FASTP_SIF")
  if [ "$got" = "$FASTP_SIF_SHA256" ]; then st=PASS; else st=FAIL; nfail=$((nfail + 1)); fi
  printf "fastp_sif_sha256\tnone\t%s\t%s\t%s\n" "$FASTP_SIF_SHA256" "$got" "$st" >> "$REPORT"
  check fastp_container "$APPTAINER_MODULES" "fastp $FASTP_EXPECTED_VERSION" "${APPTAINER_BIN:-apptainer}" exec "$FASTP_SIF" fastp --version
else
  printf "fastp_container\tnone\t%s\tnot provisioned (FASTP_SIF unset or absent)\tSKIP\n" "${FASTP_EXPECTED_VERSION:-unset}" >> "$REPORT"
  echo "SKIP  fastp container not provisioned" >&2
fi

# Module-set isolation: record what each step runs with, and flag any version swap caused by co-loading
# (observed on Rivanna: loading Apptainer changes bzip2/1.0.6 to bzip2/1.0.8). The quantification scripts never
# co-load; this is reported as INFO so the swap is documented, and as FAIL only if a step's own set is unloadable.
if [ "${RQC_SKIP_MODULE_LOAD:-0}" = 1 ]; then
  printf "module_isolation\tnone\tn/a\tmodule loading skipped\tSKIP\n" >> "$REPORT"
else
  s_set=$(module_list "$SALMON_MODULES"); a_set=$(module_list "$APPTAINER_MODULES"); both=$(module_list "$SALMON_MODULES $APPTAINER_MODULES")
  printf "salmon_step_modules\t%s\tloadable\t%s\t%s\n" "$SALMON_MODULES" "${s_set:-<none>}" "$([ -n "$s_set" ] && echo PASS || echo FAIL)" >> "$REPORT"
  printf "fastp_step_modules\t%s\tloadable\t%s\t%s\n" "$APPTAINER_MODULES" "${a_set:-<none>}" "$([ -n "$a_set" ] && echo PASS || echo FAIL)" >> "$REPORT"
  [ -n "$s_set" ] && [ -n "$a_set" ] || nfail=$((nfail + 1))
  swaps=$(comm -13 <(tr ';' '\n' <<< "$s_set" | sort) <(tr ';' '\n' <<< "$both" | sort) | while read -r m; do
            base=${m%%/*}; old=$(tr ';' '\n' <<< "$s_set" | grep "^$base/" || true)
            [ -n "$old" ] && printf '%s->%s ' "$old" "$m"; done)
  printf "coload_swaps\t%s\tnone expected under isolation\t%s\tINFO\n" "$SALMON_MODULES + $APPTAINER_MODULES" "${swaps:-none}" >> "$REPORT"
  echo "INFO  co-loading swaps (steps are isolated, so these do not affect Salmon): ${swaps:-none}" >&2
fi

if [ -n "${PROVENANCE_DIR:-}" ]; then
  provenance "${REPORT%.tsv}.provenance.json" --set "report=$REPORT" --set "n_fail=$nfail"
fi
echo "environment report: $REPORT ($nfail failed)"
[ "$nfail" -eq 0 ]
