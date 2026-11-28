#!/bin/bash
# Verify Salmon reference fingerprints (SeqHash, NameHash, SeqHash512, NameHash512, DecoySeqHash, DecoyNameHash)
# of an index directory or of a quantification's aux_info/meta_info.json against the collaborator values (EXPECT_*).
# Prints MATCH / MISMATCH / UNAVAILABLE per fingerprint; the k-mer length is reported as observed only, because the
# collaborator's k is not recorded and remains UNVERIFIED.
#
# Usage:  source config.env && bash verify_reference_fingerprint.sh <index_dir | meta_info.json> [report.json]
# Exit:   0 all MATCH, 1 any MISMATCH, 2 a fingerprint UNAVAILABLE
set -euo pipefail
. "$(dirname "$0")/lib_rivanna_qc.sh"

[ $# -ge 1 ] || die "usage: $0 <index_dir | meta_info.json> [report.json]"
[ -e "$1" ] || die "not found: $1"
OUT=${2:-}
if [ -z "$OUT" ]; then
  require_var PROVENANCE_DIR
  mkdir -p "$PROVENANCE_DIR"
  OUT="$PROVENANCE_DIR/fingerprint_$(date +%Y%m%dT%H%M%S)_$$.json"
fi
require_absent "$OUT"
set +e
verify_index_fingerprints "$1" "$OUT" ${SALMON_KMER:+--expect-k "$SALMON_KMER"}
rc=$?
set -e
log "report -> $OUT"
exit $rc
