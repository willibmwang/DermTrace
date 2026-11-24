#!/bin/bash
# Provision the version-pinned fastp Apptainer image and verify its SHA-256 and reported version.
# fastp 0.23.4 is a CANDIDATE: the collaborator's original fastp version is unknown.
#
# A SIF's SHA-256 identifies one built ARTIFACT, not the software: independently pulled SIFs of the same image tag
# (or even the same registry digest) are not guaranteed to be byte-identical, because SIF assembly can embed
# build metadata. Therefore:
#   - SIF present                 -> verify FASTP_SIF_SHA256 + version only; never re-pull or overwrite
#   - SIF absent, hash matches    -> adopt (no-clobber)
#   - SIF absent, hash differs    -> fail closed, keep the pull at *.partial for inspection, unless
#                                    FASTP_ACCEPT_NEW_SIF=1, which adopts it only if `fastp --version` matches,
#                                    and records the new hash and the reason. Update FASTP_SIF_SHA256 afterwards:
#                                    every run pins the exact artifact it used.
# Optional FASTP_IMAGE_DIGEST (sha256:...) pins the registry manifest; the URI must then use @sha256.
#
# Usage:  source config.env && bash provision_fastp_container.sh
set -euo pipefail
. "$(dirname "$0")/lib_rivanna_qc.sh"

require_var FASTP_IMAGE_URI FASTP_SIF FASTP_SIF_SHA256 FASTP_EXPECTED_VERSION
[[ "$FASTP_SIF_SHA256" =~ ^[0-9a-f]{64}$ ]] || die "FASTP_SIF_SHA256 is not a 64-hex SHA-256"
if [ -n "${FASTP_IMAGE_DIGEST:-}" ]; then
  [[ "$FASTP_IMAGE_DIGEST" =~ ^sha256:[0-9a-f]{64}$ ]] || die "FASTP_IMAGE_DIGEST must be sha256:<64 hex>"
  [[ "$FASTP_IMAGE_URI" == *"@$FASTP_IMAGE_DIGEST" ]] || die "FASTP_IMAGE_URI must pin @$FASTP_IMAGE_DIGEST"
else
  [[ "$FASTP_IMAGE_URI" == *":$FASTP_EXPECTED_VERSION"--* ]] || die "image tag does not pin fastp $FASTP_EXPECTED_VERSION: $FASTP_IMAGE_URI"
fi
APPTAINER_BIN=${APPTAINER_BIN:-apptainer}
APPTAINER_MODS="${APPTAINER_MODULES:-}"
app_mods_seen=$(module_list "$APPTAINER_MODS") || die "cannot record apptainer step modules"
run_capture app_ver "apptainer version check" "$APPTAINER_MODS" "$APPTAINER_BIN" --version
fastp_version_of() { local _v; run_capture _v "fastp version check" "$APPTAINER_MODS" "$APPTAINER_BIN" exec "$1" fastp --version; printf '%s\n' "$_v"; }

mode=verify_existing; sif_hash=$FASTP_SIF_SHA256; note="matches pinned artifact"
if [ ! -e "$FASTP_SIF" ] && [ ! -L "$FASTP_SIF" ]; then
  mkdir -p "$(dirname "$FASTP_SIF")"
  require_absent "$FASTP_SIF.partial"
  run_with_modules "$APPTAINER_MODS" "$APPTAINER_BIN" pull "$FASTP_SIF.partial" "$FASTP_IMAGE_URI" \
    || { rm -f "$FASTP_SIF.partial"; die "apptainer pull failed"; }
  got=$(sha256_of "$FASTP_SIF.partial")
  if [ "$got" = "$FASTP_SIF_SHA256" ]; then
    mode=pulled_matching
  elif [ "${FASTP_ACCEPT_NEW_SIF:-0}" = 1 ]; then
    v=$(fastp_version_of "$FASTP_SIF.partial") || die "fastp version check failed for $FASTP_SIF.partial (kept; see message above)"
    [ "$v" = "fastp $FASTP_EXPECTED_VERSION" ] || die "new SIF reports '$v'; not adopted ($FASTP_SIF.partial kept)"
    mode=pulled_new_artifact_accepted; sif_hash=$got
    note="new artifact; differs from previously pinned $FASTP_SIF_SHA256; adopted by FASTP_ACCEPT_NEW_SIF=1 after version check"
  else
    echo "Pulled image kept for inspection at $FASTP_SIF.partial (not adopted)." >&2
    echo "Independently pulled SIFs can differ byte-wise. To adopt it deliberately: FASTP_ACCEPT_NEW_SIF=1, then set FASTP_SIF_SHA256=$got" >&2
    die "SIF SHA-256 mismatch: expected $FASTP_SIF_SHA256, got $got"
  fi
  chmod a-w "$FASTP_SIF.partial"
  adopt_file "$FASTP_SIF.partial" "$FASTP_SIF"
fi

verify_sha256 "$FASTP_SIF" "$sif_hash"
ver=$(fastp_version_of "$FASTP_SIF") || die "fastp version check failed for $FASTP_SIF (see message above)"
[ "$ver" = "fastp $FASTP_EXPECTED_VERSION" ] || die "unexpected fastp version string: '$ver' (want 'fastp $FASTP_EXPECTED_VERSION')"
log "$ver in $FASTP_SIF ($app_ver); $note"

if [ -n "${PROVENANCE_DIR:-}" ]; then
  provenance "$PROVENANCE_DIR/provision_fastp_$(date +%Y%m%dT%H%M%S)_$$.json" --set "status=complete" \
    --set "mode=$mode" --set "image_uri=$FASTP_IMAGE_URI" --set "image_digest=${FASTP_IMAGE_DIGEST:-not pinned}" \
    --set "sif=$FASTP_SIF" --set "sif_sha256=$sif_hash" --set "sif_note=$note" --set "apptainer_version=$app_ver" \
    --set "apptainer_step_modules=$app_mods_seen" \
    --config "fastp_version=$ver::candidate::collaborator fastp version unknown" \
    --config "fastp_sif_sha256=$sif_hash::observed-this-run::artifact identity only; not a cross-pull software identity"
fi
[ "$sif_hash" = "$FASTP_SIF_SHA256" ] || echo "ACTION: set FASTP_SIF_SHA256=$sif_hash in your config before quantifying" >&2
echo "FASTP_OK $FASTP_SIF ($ver, candidate; sha256 $sif_hash)"
