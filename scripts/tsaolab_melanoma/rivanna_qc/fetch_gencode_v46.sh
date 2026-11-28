#!/bin/bash
# Download the GENCODE v46 transcriptome FASTA and verify it against BOTH the release MD5SUMS file and a pinned MD5.
# Idempotent and fail-closed:
#   - final file absent  -> download to *.partial, verify, then rename
#   - final file present -> verify only; never re-download or overwrite
#   - any mismatch       -> exit non-zero; the partial download is removed and no final file is created
#
# Usage:  source config.env && bash fetch_gencode_v46.sh
# Env:    REF_DIR GENCODE_BASE_URL GENCODE_TRANSCRIPTS GENCODE_TRANSCRIPTS_MD5 [PROVENANCE_DIR]
set -euo pipefail
. "$(dirname "$0")/lib_rivanna_qc.sh"

require_var REF_DIR GENCODE_BASE_URL GENCODE_TRANSCRIPTS GENCODE_TRANSCRIPTS_MD5
[[ "$GENCODE_TRANSCRIPTS_MD5" =~ ^[0-9a-f]{32}$ ]] || die "GENCODE_TRANSCRIPTS_MD5 is not a 32-hex MD5"
command -v curl >/dev/null || die "curl not found"

mkdir -p "$REF_DIR"
FA="$REF_DIR/$GENCODE_TRANSCRIPTS"
SUMS="$REF_DIR/MD5SUMS.release_46"
mode=verify_existing

if [ ! -e "$SUMS" ]; then
  require_absent "$SUMS.partial"
  curl -fsSL --retry 3 -o "$SUMS.partial" "$GENCODE_BASE_URL/MD5SUMS" || { rm -f "$SUMS.partial"; die "MD5SUMS download failed"; }
  adopt_file "$SUMS.partial" "$SUMS"
fi
listed=$(awk -v f="$GENCODE_TRANSCRIPTS" '$2==f || $2=="*"f {print $1}' "$SUMS")
[ -n "$listed" ] || die "$GENCODE_TRANSCRIPTS not listed in $SUMS"
[ "$listed" = "$GENCODE_TRANSCRIPTS_MD5" ] || die "release MD5SUMS ($listed) disagrees with pinned GENCODE_TRANSCRIPTS_MD5 ($GENCODE_TRANSCRIPTS_MD5)"

if [ ! -e "$FA" ]; then
  mode=downloaded
  require_absent "$FA.partial"
  log "downloading $GENCODE_BASE_URL/$GENCODE_TRANSCRIPTS"
  curl -fsSL --retry 3 -o "$FA.partial" "$GENCODE_BASE_URL/$GENCODE_TRANSCRIPTS" || { rm -f "$FA.partial"; die "download failed"; }
  got=$(md5_of "$FA.partial")
  if [ "$got" != "$GENCODE_TRANSCRIPTS_MD5" ]; then rm -f "$FA.partial"; die "MD5 mismatch on download: got $got"; fi
  gzip -t "$FA.partial" || { rm -f "$FA.partial"; die "gzip integrity test failed"; }
  chmod a-w "$FA.partial"
  adopt_file "$FA.partial" "$FA"
fi

verify_md5 "$FA" "$GENCODE_TRANSCRIPTS_MD5"
n_tx=$(gzip -dc "$FA" | grep -c '^>')
n_par=$(gzip -dc "$FA" | grep '^>' | grep -c '_PAR_Y' || true)
log "transcripts=$n_tx PAR_Y=$n_par mode=$mode"

if [ -n "${PROVENANCE_DIR:-}" ]; then
  provenance "$PROVENANCE_DIR/fetch_gencode_v46_$(date +%Y%m%dT%H%M%S)_$$.json" --set "status=complete" \
    --set "mode=$mode" --set "url=$GENCODE_BASE_URL/$GENCODE_TRANSCRIPTS" --set "path=$FA" \
    --set "md5=$GENCODE_TRANSCRIPTS_MD5" --set "sha256=$(sha256_of "$FA")" \
    --set "n_transcripts=$n_tx" --set "n_par_y=$n_par"
fi
echo "GENCODE_OK $FA"
