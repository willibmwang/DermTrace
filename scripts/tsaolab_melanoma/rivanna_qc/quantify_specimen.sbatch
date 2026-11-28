#!/bin/bash
#SBATCH --job-name=melanoma_quant
#SBATCH --partition=standard
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=melanoma_quant_%j.out
#SBATCH --error=melanoma_quant_%j.err
#
# One specimen: fastp (CANDIDATE version, pinned container) -> Salmon 1.10.3, with the collaborator's parameters
# from scripts/tsaolab_melanoma/run_quant.sh (fastp --detect_adapter_for_pe -l 25; salmon quant -l A
# --validateMappings --gcBias --seqBias). The driver-reference minimap2 step is not included.
#
# Eligibility: the specimen must be in the sample sheet AND have exactly one gzip-validated pair in
# validated_pairs.tsv with the same R1/R2 files; built-in exclusions (11-37 pending provenance adjudication;
# corrupt 11-43/44/46/47/48) can be lifted only with a recorded adjudication reference (LIFT_EXCLUSIONS).
#
# Isolation: fastp runs with APPTAINER_MODULES only and Salmon with SALMON_MODULES only (separate purged
# subshells), so loading Apptainer cannot swap libraries (e.g. bzip2) under Salmon. Module sets are recorded.
#
# Safety: a per-specimen lock blocks concurrent jobs; work happens in $QUANT_ROOT/.<id>.partial.<job>; outputs
# are structurally validated, provenance is written with status=complete, and only then is the directory
# renamed to $QUANT_ROOT/<id> (re-checked absent under the lock). Raw inputs are only read.
#
# Submit:  source config.env && sbatch --export=ALL,SPECIMEN=GM01 "$RQC_DIR/quantify_specimen.sbatch"
set -euo pipefail
RQC_SCRIPT_NAME=quantify_specimen.sbatch
RQC_DIR=${RQC_DIR:-$(cd "$(dirname "$0")" && pwd)}
[ -f "$RQC_DIR/lib_rivanna_qc.sh" ] || { echo "ERROR: set RQC_DIR to the rivanna_qc script directory" >&2; exit 1; }
. "$RQC_DIR/lib_rivanna_qc.sh"

require_var SPECIMEN SAMPLE_SHEET VALIDATED_PAIRS RAW_DIR QUANT_ROOT INDEX_DIR FASTP_SIF FASTP_SIF_SHA256 \
  FASTP_EXPECTED_VERSION EXPECT_SEQ_HASH EXPECT_NAME_HASH EXPECT_SEQ_HASH512 EXPECT_NAME_HASH512 EXPECT_DECOY_HASH \
  EXPECT_NUM_TARGETS EXPECT_LIBTYPE PROVENANCE_DIR
require_file "$SAMPLE_SHEET"
require_file "$VALIDATED_PAIRS"
require_real_dir "$RAW_DIR"
require_real_dir "$INDEX_DIR"
case "$(cd "$INDEX_DIR" && pwd -P)" in *.partial*) die "INDEX_DIR is a partial (unadopted) index: $INDEX_DIR" ;; esac
require_file "$FASTP_SIF"

LIFT_ARGS=()
for x in ${LIFT_EXCLUSIONS:-}; do LIFT_ARGS+=(--lift-exclusion "$x"); done
RES=$("$PYTHON_BIN" "$RQC_TOOLS" resolve-specimen --sheet "$SAMPLE_SHEET" --validated-pairs "$VALIDATED_PAIRS" \
      --specimen "$SPECIMEN" --raw-dir "$RAW_DIR" --exclude "${EXCLUDED_SPECIMENS:-}" ${LIFT_ARGS[@]+"${LIFT_ARGS[@]}"}) \
  || die "specimen $SPECIMEN is not eligible (see message above)"
IFS=$'\t' read -r SID CORE R1N R2N MD5_R1 MD5_R2 VP_LABEL LIFTED <<< "$RES"
[ -n "${LIFTED:-}" ] || die "incomplete specimen record for $SPECIMEN"
R1="$RAW_DIR/$R1N"; R2="$RAW_DIR/$R2N"
require_file "$R1"; require_file "$R2"

mkdir -p "$QUANT_ROOT" "$PROVENANCE_DIR"
require_real_dir "$QUANT_ROOT"
path_is_under "$QUANT_ROOT" "$RAW_DIR" && die "QUANT_ROOT must not be inside RAW_DIR"
path_is_under "$RAW_DIR" "$QUANT_ROOT" && die "RAW_DIR must not be inside QUANT_ROOT"
FINAL="$QUANT_ROOT/$SID"
WORK="$QUANT_ROOT/.$SID.partial.${SLURM_JOB_ID:-$$}"
acquire_lock "$QUANT_ROOT/.$SID.lock"
trap 'rm -f "$WORK/tmp/"*.fq.gz 2>/dev/null || true; release_lock' EXIT
require_absent "$FINAL"
require_absent "$WORK"

SALMON_BIN=${SALMON_BIN:-salmon}; APPTAINER_BIN=${APPTAINER_BIN:-apptainer}
SALMON_MODS="${SALMON_MODULES:-}"; APPTAINER_MODS="${APPTAINER_MODULES:-}"
run_capture sv "salmon version check" "$SALMON_MODS" "$SALMON_BIN" --version
[ "$sv" = "salmon ${EXPECT_SALMON:-1.10.3}" ] || die "unexpected salmon version '$sv'"
salmon_mods_seen=$(module_list "$SALMON_MODS") || die "cannot record salmon step modules"
fastp_mods_seen=$(module_list "$APPTAINER_MODS") || die "cannot record fastp step modules"

log "specimen $SID ($CORE); validated_pairs label '$VP_LABEL'; lifted exclusions: $LIFTED"
verify_md5 "$R1" "$MD5_R1"
verify_md5 "$R2" "$MD5_R2"
verify_sha256 "$FASTP_SIF" "$FASTP_SIF_SHA256"
run_capture fv "fastp version check" "$APPTAINER_MODS" "$APPTAINER_BIN" exec "$FASTP_SIF" fastp --version
[ "$fv" = "fastp $FASTP_EXPECTED_VERSION" ] || die "unexpected fastp version '$fv'"

mkdir "$WORK"
mkdir "$WORK/fastp" "$WORK/tmp"
verify_index_fingerprints "$INDEX_DIR" "$WORK/index_dir_fingerprint.json" > "$WORK/index_dir_fingerprint.txt" \
  || die "index directory fingerprints did not all MATCH; nothing quantified (see $WORK/index_dir_fingerprint.txt)"
k_obs=$("$PYTHON_BIN" -c 'import json,sys; print(json.load(open(sys.argv[1]))["kmer"]["observed_k"])' "$WORK/index_dir_fingerprint.json")
T1="$WORK/tmp/${SID}_R1.fq.gz"; T2="$WORK/tmp/${SID}_R2.fq.gz"

log "fastp ($fv, CANDIDATE) with modules: $fastp_mods_seen"
run_with_modules "$APPTAINER_MODS" "$APPTAINER_BIN" exec "$FASTP_SIF" fastp -i "$R1" -I "$R2" -o "$T1" -O "$T2" \
  --detect_adapter_for_pe -l "${FASTP_MIN_LENGTH:-25}" -w "${FASTP_THREADS:-4}" \
  -j "$WORK/fastp/$SID.json" -h "$WORK/fastp/$SID.html" 2> "$WORK/fastp/$SID.log" < /dev/null \
  || die "fastp failed (see $WORK/fastp/$SID.log)"

log "salmon quant ($sv) with modules: $salmon_mods_seen"
run_with_modules "$SALMON_MODS" "$SALMON_BIN" quant -i "$INDEX_DIR" -l A -1 "$T1" -2 "$T2" -p "${SALMON_THREADS:-6}" \
  --validateMappings --gcBias --seqBias -o "$WORK/quant" > "$WORK/salmon.log" 2>&1 < /dev/null \
  || die "salmon failed (see $WORK/salmon.log)"
rm -f "$T1" "$T2"; rmdir "$WORK/tmp"

"$PYTHON_BIN" "$RQC_TOOLS" validate-quant-dir --dir "$WORK" --stage --expect-targets "$EXPECT_NUM_TARGETS" \
  --expect-libtype "$EXPECT_LIBTYPE" > "$WORK/structure_check.txt" \
  || die "staged output failed structural validation (see $WORK/structure_check.txt)"
verify_index_fingerprints "$WORK/quant/aux_info/meta_info.json" "$WORK/index_fingerprint.json" > /dev/null \
  || die "quantification index fingerprints did not all MATCH; output left at $WORK"

C=verified-collaborator
provenance "$WORK/provenance.json" --set "status=complete" \
  --set "specimen=$SID" --set "core_id=$CORE" --set "validated_pairs_label=$VP_LABEL" --set "lifted_exclusions=$LIFTED" \
  --set "r1=$R1" --set "r2=$R2" --set "md5_r1=$MD5_R1" --set "md5_r2=$MD5_R2" \
  --set "fastp_sif=$FASTP_SIF" --set "fastp_sif_sha256=$FASTP_SIF_SHA256" \
  --set "salmon_step_modules=$salmon_mods_seen" --set "fastp_step_modules=$fastp_mods_seen" \
  --config "salmon_version=$sv::$C::committed cmd_info.json / meta_info.json" \
  --config "salmon_quant_args=-l A --validateMappings --gcBias --seqBias::$C::committed cmd_info.json" \
  --config "salmon_threads=${SALMON_THREADS:-6}::$C::committed cmd_info.json threads=6" \
  --config "index_fingerprints=six SHA-256/SHA-512 hashes::observed-this-run::matched EXPECT_* (committed meta_info.json values)" \
  --config "index_kmer=$k_obs::observed-this-run::this index info.json; collaborator k NOT recorded, remains unverified" \
  --config "transcriptome=GENCODE v46 transcripts::$C::committed index name/sequence hashes" \
  --config "fastp_args=--detect_adapter_for_pe -l ${FASTP_MIN_LENGTH:-25}::$C::scripts/tsaolab_melanoma/run_quant.sh" \
  --config "fastp_version=$fv::candidate::collaborator fastp version unknown" \
  --config "fastp_container=$FASTP_SIF_SHA256::candidate::pinned SIF used for this run" \
  --config "fastp_threads=${FASTP_THREADS:-4}::inferred::collaborator THREADS inferred as salmon threads minus 2" \
  --json "index_fingerprint=$WORK/index_fingerprint.json"
"$PYTHON_BIN" "$RQC_TOOLS" validate-quant-dir --dir "$WORK" --stage --specimen "$SID" \
  --expect-targets "$EXPECT_NUM_TARGETS" --expect-libtype "$EXPECT_LIBTYPE" > /dev/null || die "post-provenance validation failed"

require_absent "$FINAL"   # re-check under the lock immediately before adoption
mv "$WORK" "$FINAL"
[ -f "$FINAL/provenance.json" ] && [ ! -e "$WORK" ] || die "adoption of $WORK as $FINAL did not complete as expected"
"$PYTHON_BIN" "$RQC_TOOLS" validate-quant-dir --dir "$FINAL" --specimen "$SID" \
  --expect-targets "$EXPECT_NUM_TARGETS" --expect-libtype "$EXPECT_LIBTYPE" > /dev/null || die "adopted output failed validation: $FINAL"
echo "QUANT_OK $FINAL"
