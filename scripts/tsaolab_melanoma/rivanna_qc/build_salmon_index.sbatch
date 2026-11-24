#!/bin/bash
#SBATCH --job-name=salmon_index_v46
#SBATCH --partition=standard
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=salmon_index_v46_%j.out
#SBATCH --error=salmon_index_v46_%j.err
#
# Build the GENCODE v46 transcriptome-only Salmon index (no decoys) and verify its fingerprints against the
# collaborator's committed runs before adopting it. Candidate configuration: --gencode -k $SALMON_KMER.
# The collaborator's original k is unverified (k is not part of the sequence/name fingerprints).
#
# A lock ($INDEX_DIR.lock) blocks concurrent builds. Builds into $INDEX_DIR.partial.<job>; renamed to $INDEX_DIR
# only if all six fingerprints MATCH and the observed k equals SALMON_KMER. Refuses if $INDEX_DIR exists.
# Salmon runs in a purged subshell holding only SALMON_MODULES.
#
# Submit:  source config.env && sbatch --export=ALL "$RQC_DIR/build_salmon_index.sbatch"
set -euo pipefail
RQC_SCRIPT_NAME=build_salmon_index.sbatch
RQC_DIR=${RQC_DIR:-$(cd "$(dirname "$0")" && pwd)}
[ -f "$RQC_DIR/lib_rivanna_qc.sh" ] || { echo "ERROR: set RQC_DIR to the rivanna_qc script directory" >&2; exit 1; }
. "$RQC_DIR/lib_rivanna_qc.sh"

require_var REF_DIR GENCODE_TRANSCRIPTS GENCODE_TRANSCRIPTS_MD5 INDEX_DIR SALMON_KMER \
  EXPECT_SEQ_HASH EXPECT_NAME_HASH EXPECT_SEQ_HASH512 EXPECT_NAME_HASH512 EXPECT_DECOY_HASH PROVENANCE_DIR
FA="$REF_DIR/$GENCODE_TRANSCRIPTS"
require_file "$FA"
acquire_lock "$INDEX_DIR.lock"
trap release_lock EXIT
require_absent "$INDEX_DIR"
PART="$INDEX_DIR.partial.${SLURM_JOB_ID:-$$}"
require_absent "$PART"
SALMON_BIN=${SALMON_BIN:-salmon}
SALMON_MODS="${SALMON_MODULES:-}"
run_capture sv "salmon version check" "$SALMON_MODS" "$SALMON_BIN" --version
[ "$sv" = "salmon ${EXPECT_SALMON:-1.10.3}" ] || die "unexpected salmon version '$sv'"
mods_seen=$(module_list "$SALMON_MODS") || die "cannot record salmon step modules"

verify_md5 "$FA" "$GENCODE_TRANSCRIPTS_MD5"
log "$sv: indexing $FA -> $PART (k=$SALMON_KMER, --gencode, no decoys); modules: $mods_seen"
run_with_modules "$SALMON_MODS" "$SALMON_BIN" index -t "$FA" -i "$PART" --gencode -k "$SALMON_KMER" -p "${SALMON_INDEX_THREADS:-8}" \
  || die "salmon index failed; partial index left at $PART"

mkdir -p "$PROVENANCE_DIR"
STAMP=$(date +%Y%m%dT%H%M%S)
FP_JSON="$PROVENANCE_DIR/index_fingerprint_${STAMP}_${SLURM_JOB_ID:-$$}.json"
set +e
verify_index_fingerprints "$PART" "$FP_JSON" --expect-k "$SALMON_KMER"
fp=$?
set -e
[ $fp -eq 0 ] || die "index fingerprints/k did not all MATCH (exit $fp); index left unadopted at $PART"

require_absent "$INDEX_DIR"   # re-check under the lock
mv "$PART" "$INDEX_DIR"
chmod -R a-w "$INDEX_DIR"
C=verified-collaborator
provenance "$PROVENANCE_DIR/build_salmon_index_${STAMP}_${SLURM_JOB_ID:-$$}.json" --set "status=complete" \
  --set "transcriptome=$FA" --set "index_dir=$INDEX_DIR" --set "salmon_step_modules=$mods_seen" \
  --config "salmon_version=$sv::$C::committed meta_info.json" \
  --config "transcriptome_md5=$GENCODE_TRANSCRIPTS_MD5::observed-this-run::verified against pinned MD5" \
  --config "index_fingerprints=six SHA-256/SHA-512 hashes::observed-this-run::matched EXPECT_* (committed meta_info.json values)" \
  --config "decoys=none::$C::committed decoy hashes are SHA-256 of empty input" \
  --config "gencode_flag=--gencode::inferred::committed quant.sf transcript names lack '|' fields; name hash matches" \
  --config "kmer=$SALMON_KMER::candidate::observed in this index info.json; collaborator k NOT recorded, remains unverified" \
  --json "fingerprint_report=$FP_JSON"
echo "INDEX_OK $INDEX_DIR"
