#!/usr/bin/env bash
# Unattended run for ONE Quintara order: finish Drive download -> verify MD5 -> salmon quant -> per-order tables.
# When every order is quantified, the last one to finish builds the combined DermTrace dataset.
# Safe to re-run: rclone skips finished files, run_quant.sh skips quantified samples.
#
# Usage: DATA_DIR=... REF_DIR=... [OUT_DIR=...] [RNASEQ_ENV=...] [PYTHON=...] run_all.sh <order_id>
#        <order_id> is 0013831748 or 1030390
# Progress: datasets/processed/STATUS_<order_id> (one line) and run_all_<order_id>.log
set -uo pipefail

ID=$1
HERE="$(cd "$(dirname "$0")" && pwd)"
DS="${DATA_DIR:?set DATA_DIR to the folder holding the <order> Raw_Data folders}"
OUT="${OUT_DIR:-$DS/processed}"
export REF_DIR="${REF_DIR:?set REF_DIR to the folder holding gencode_human_v46/ and gencode_mouse_vM35/}"
ENV="${RNASEQ_ENV:-/opt/anaconda3/envs/rnaseq}/bin"  # conda env with salmon 1.10.3, fastp, rclone
PY="${PYTHON:-python3}"  # needs pandas
RC="$ENV/rclone"
export PATH="$ENV:$PATH"
SHEETS="$(cd "$HERE/../../data/processed/tsaolab_o2/sample_sheets" && pwd)"
mkdir -p "$OUT"
status() { echo "$(date '+%F %T') $*" | tee "$OUT/STATUS_$ID"; }
fail() { status "FAILED: $*"; exit 1; }

# order: raw dir | Drive folder id | lanes | how to verify
case $ID in
  0013831748) RAW="0013831748-2024 Raw_Data"; FID=1Wf0_DocioVgJOXxGWpdgESjGNxJm8mUC; LANES="L001 L002" ;;
  1030390)    RAW="1030390-2025 Raw_Data";    FID=1sjsuRjayueMCZmLZFhiYsbZUeqttWQ42; LANES="L001" ;;
  *) echo "unknown order $ID"; exit 1 ;;
esac
D="$DS/$RAW"
ALL_ORDERS=(0013831748 1030390)

# 1. let a download of this order started earlier finish, then (re)run the copy until complete
while pgrep -f "rclone copy yale: $D " >/dev/null; do status "waiting for running download"; sleep 120; done
for attempt in 1 2 3; do
  status "download (attempt $attempt)"
  "$RC" copy yale: "$D" --drive-root-folder-id "$FID" --exclude md5sum.txt \
    --transfers 4 --checkers 8 --retries 5 --low-level-retries 20 --stats 10m --stats-one-line -v && break
done

# 2. verify: Quintara's md5sum.txt when the order has one, otherwise Drive's own MD5s
status "verifying checksums"
if "$RC" lsf yale: --drive-root-folder-id "$FID" --files-only --max-depth 1 | grep -qx md5sum.txt; then
  "$RC" copyto yale:md5sum.txt "$D/md5sum.txt" --drive-root-folder-id "$FID"
  ( cd "$D" && "$RC" md5sum . --exclude md5sum.txt | awk '{print $1, $2}' | sort -k2 ) > "$OUT/md5_local_$ID.txt"
  awk '{print $1, $2}' "$D/md5sum.txt" | sort -k2 > "$OUT/md5_expected_$ID.txt"
  diff "$OUT/md5_expected_$ID.txt" "$OUT/md5_local_$ID.txt" || fail "MD5 mismatch vs md5sum.txt (see $OUT/md5_*_$ID.txt)"
else
  "$RC" check yale: "$D" --drive-root-folder-id "$FID" --one-way || fail "MD5 check against Drive failed"
fi
status "checksums OK ($(find "$D" -name '*.fastq.gz' | wc -l | tr -d ' ') FASTQs)"

# 3. quantify + per-order tables
status "salmon"
LANES="$LANES" bash "$HERE/run_quant.sh" "$SHEETS/samples_$ID.csv" "$D" "$OUT/quant_$ID" || fail "run_quant"
grep -q incomplete "$OUT/quant_$ID/sample_status.tsv" && fail "incomplete samples (see quant_$ID/sample_status.tsv)"
"$PY" -I "$HERE/aggregate.py" "$SHEETS/samples_$ID.csv" "$OUT/quant_$ID" "$OUT/per_order" || fail "aggregate"
touch "$OUT/quant_$ID/.complete"
status "quantified"

# 4. combined DermTrace dataset, built once by whichever order finishes last
for o in "${ALL_ORDERS[@]}"; do [ -e "$OUT/quant_$o/.complete" ] || { status "DONE (waiting on order $o for combined build)"; exit 0; }; done
mkdir "$OUT/.build_lock" 2>/dev/null || { status "DONE (other order is building combined dataset)"; exit 0; }
status "building DermTrace dataset"
args=(); for o in "${ALL_ORDERS[@]}"; do args+=("$SHEETS/samples_$o.csv=$OUT/quant_$o"); done
"$PY" -I "$HERE/build_dermtrace.py" "$OUT/tsaolab_o2_dermtrace" "${args[@]}" || { rmdir "$OUT/.build_lock"; fail "build_dermtrace"; }
rmdir "$OUT/.build_lock"
status "DONE -> $OUT/tsaolab_o2_dermtrace"
