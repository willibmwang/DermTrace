#!/bin/bash
# Step 15 (login node; writes only a new config and new provenance reports): harmonization preflight.
#  - manifests (procedures 26/27, r3, committed matrices, 33-specimen reference bundle)
#  - config = r2-qualified configuration + r3 code + new output roots; no exclusion lifted
#  - inventory: eligible (r3, REAL validated_pairs.tsv) = 42; the 14 existing outputs are eligible; the remainder is
#    exactly the 28 expected specimens
#  - the 14 existing outputs are compatible (structure, provenance, identical configuration, fingerprints): reuse
#    only if this PASSES
set -euo pipefail
. "$(dirname "$0")/common.sh"
verify_manifests "$PD" "$PD26" "$R3" "$CM" "$REF33"
echo "OK    manifests"
if [ ! -e "$HC" ]; then
  cat > "$HC" <<CFG
# Harmonized quantification of the 28 remaining eligible specimens: r3 code, r2-qualified configuration unchanged.
source $QC_CONF
export RQC_DIR=$R3
export QUANT_ROOT=$M/quant_harmonize_r3
export PROVENANCE_DIR=$M/provenance_harmonize_r3
export EXCLUDED_SPECIMENS=""
export LIFT_EXCLUSIONS=""
CFG
  chmod a-w "$HC"; echo "OK    wrote $HC"
else echo "OK    reusing $HC (read-only)"; fi
# shellcheck disable=SC1090
source "$HC"
[ "$RQC_DIR" = "$R3" ] && [ -z "$LIFT_EXCLUSIONS" ] && [ -z "$EXCLUDED_SPECIMENS" ] || { echo "STOP config is not r3 / exclusions altered"; exit 1; }
if [ -e "$QUANT_ROOT" ]; then echo "STOP $QUANT_ROOT exists"; exit 1; fi
mkdir -p "$PROVENANCE_DIR"; RUN=$PROVENANCE_DIR/preflight_$(date +%Y%m%dT%H%M%S); mkdir "$RUN"
eligible=""
for s in $(awk -F, 'NR>1 {print $1}' "$SAMPLE_SHEET"); do
  if "$PYTHON_BIN" -I "$RQC_DIR/rivanna_qc_tools.py" resolve-specimen --sheet "$SAMPLE_SHEET" --validated-pairs "$VALIDATED_PAIRS" \
       --specimen "$s" --raw-dir "$RAW_DIR" >> "$RUN/resolve_all.tsv" 2>> "$RUN/resolve_blocked.txt"; then eligible="$eligible $s"; fi
done
existing_ids=$(for e in $EXISTING; do echo "${e%%=*}"; done | sort)
n_el=$(printf '%s\n' $eligible | grep -c .)
remaining=$(comm -23 <(printf '%s\n' $eligible | sort) <(printf '%s\n' "$existing_ids") | tr '\n' ' ' | sed 's/ $//')
echo "eligible: $n_el; existing: $(printf '%s\n' "$existing_ids" | grep -c .); remaining: $(printf '%s\n' $remaining | grep -c .)"
[ "$n_el" -eq 42 ] || { echo "STOP expected 42 eligible specimens, found $n_el"; exit 1; }
[ -z "$(comm -13 <(printf '%s\n' $eligible | sort) <(printf '%s\n' "$existing_ids"))" ] || { echo "STOP an existing output is not eligible"; exit 1; }
[ "$remaining" = "$(printf '%s\n' $REMAINING | sort | tr '\n' ' ' | sed 's/ $//')" ] || { echo "STOP remaining set differs from the expected 28: $remaining"; exit 1; }
for b in GM37 GM43 GM44 GM46 GM47 GM48; do case " $eligible " in *" $b "*) echo "STOP $b eligible"; exit 1 ;; esac; done
echo "OK    inventory: 42 eligible = 14 existing + exactly the 28 expected; GM37 and corrupt specimens blocked"
args=(); for e in $EXISTING; do args+=(--specimen "$e"); done
"$PYTHON_BIN" -I "$PD/campaign_tools.py" verify-compat --rqc-tools "$RQC_DIR/rivanna_qc_tools.py" "${args[@]}" --reference-specimen GM01 \
  --expect-targets "$EXPECT_NUM_TARGETS" --expect-libtype "$EXPECT_LIBTYPE" --expect-seq-hash "$EXPECT_SEQ_HASH" --expect-name-hash "$EXPECT_NAME_HASH" \
  --expect-seq-hash512 "$EXPECT_SEQ_HASH512" --expect-name-hash512 "$EXPECT_NAME_HASH512" --expect-decoy-hash "$EXPECT_DECOY_HASH" \
  --expect-sif-sha256 "$FASTP_SIF_SHA256" --json-out "$RUN/compat_existing14.json"
printf '%s\n' $REMAINING > "$RUN/specimens_28.txt"
sha256sum "$RUN/specimens_28.txt" | cut -d' ' -f1 > "$RUN/PASS"
chmod -R a-w "$RUN"
echo "STEP 15 PASS ($RUN)"
