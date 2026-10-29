#!/usr/bin/env bash
# Trim (fastp) + quantify (salmon, GENCODE v46) + align to the melanoma driver reference (minimap2) for every
# sample of the 20088-11 sheet whose two FASTQs are downloaded and match md5.txt. One lane per sample, flat folder.
# Re-runnable: finished samples are skipped. While an rclone download is still running, samples that are not
# ready yet are retried every 5 minutes; otherwise they are reported as incomplete.
#
# Usage: run_quant.sh <sample_sheet.csv> <raw_data_dir> <out_dir>
#   env: THREADS=8  RNASEQ_ENV=/opt/anaconda3/envs/rnaseq  VARCALL_ENV=/opt/anaconda3/envs/varcall
#        DRIVER_REF=<references>/melanoma_driver_ref/driver_ref.fa
set -uo pipefail
SHEET=$1; RAW=$2; OUT=$3
HERE="$(cd "$(dirname "$0")" && pwd)"
REFS="$(cd "$HERE/../../../references" && pwd)"
THREADS=${THREADS:-8}
RNASEQ="${RNASEQ_ENV:-/opt/anaconda3/envs/rnaseq}/bin"
VARCALL="${VARCALL_ENV:-/opt/anaconda3/envs/varcall}/bin"
DRIVER_REF=${DRIVER_REF:-$REFS/melanoma_driver_ref/driver_ref.fa}
IDX="$REFS/gencode_human_v46/salmon_index"
mkdir -p "$OUT/quant" "$OUT/fastp" "$OUT/driver_bam" "$OUT/tmp"
STATUS="$OUT/sample_status.tsv"; MD5OK="$OUT/md5_verified.tsv"
touch "$MD5OK"
[ -s "$OUT/driver_ref.mmi" ] || "$VARCALL/minimap2" -d "$OUT/driver_ref.mmi" "$DRIVER_REF" 2>/dev/null
log() { echo "$(date '+%F %T') $*"; }

verified() {  # <file> <md5>: true once the file's md5 matched (cached)
  grep -q "^$1	ok" "$MD5OK" && return 0
  [ -s "$RAW/$1" ] || return 1
  local h; h=$(md5 -q "$RAW/$1" 2>/dev/null || md5sum "$RAW/$1" | cut -d' ' -f1)
  if [ "$h" = "$2" ]; then printf "%s\tok\n" "$1" >> "$MD5OK"; return 0; fi
  printf "%s\tMISMATCH %s\n" "$1" "$h" >> "$MD5OK"; log "MD5 MISMATCH $1"; return 1
}

while :; do
  pending=0
  printf "sample_id\tcore_id\tstatus\tdetail\n" > "$STATUS"
  tail -n +2 "$SHEET" | while IFS=, read -r sid core _project _species _tissue _spec _dis _lane _s r1 r2 m1 m2; do
    if [ -s "$OUT/quant/$sid/quant.sf" ] && [ -s "$OUT/driver_bam/$sid.bam.bai" ]; then
      printf "%s\t%s\tdone\t\n" "$sid" "$core" >> "$STATUS"; continue
    fi
    if ! verified "$r1" "$m1" || ! verified "$r2" "$m2"; then
      printf "%s\t%s\tincomplete\tFASTQ not downloaded or md5 mismatch\n" "$sid" "$core" >> "$STATUS"; continue
    fi
    log "[$sid] fastp"
    t1="$OUT/tmp/${sid}_R1.fq.gz"; t2="$OUT/tmp/${sid}_R2.fq.gz"
    "$RNASEQ/fastp" -i "$RAW/$r1" -I "$RAW/$r2" -o "$t1" -O "$t2" --detect_adapter_for_pe -l 25 -w "$THREADS" \
      -j "$OUT/fastp/$sid.json" -h "$OUT/fastp/$sid.html" 2> "$OUT/fastp/$sid.log" < /dev/null \
      || { log "[$sid] fastp FAILED"; printf "%s\t%s\tfailed\tfastp\n" "$sid" "$core" >> "$STATUS"; continue; }
    log "[$sid] salmon + minimap2"
    "$RNASEQ/salmon" quant -i "$IDX" -l A -1 "$t1" -2 "$t2" -p $((THREADS - 2)) --validateMappings --gcBias --seqBias \
      -o "$OUT/quant/$sid" > "$OUT/quant/$sid.log" 2>&1 < /dev/null &
    sp=$!
    "$VARCALL/minimap2" -t 3 -ax sr --secondary=no -R "@RG\tID:$sid\tSM:$sid" "$OUT/driver_ref.mmi" "$t1" "$t2" 2> "$OUT/driver_bam/$sid.minimap2.log" \
      | "$VARCALL/samtools" view -u -q 20 -F 0x904 - | "$VARCALL/samtools" sort -@ 2 -o "$OUT/driver_bam/$sid.bam" - 2>/dev/null \
      && "$VARCALL/samtools" index "$OUT/driver_bam/$sid.bam"
    mm=$?
    wait $sp; sm=$?
    rm -f "$t1" "$t2"
    if [ $sm -ne 0 ] || [ $mm -ne 0 ] || [ ! -s "$OUT/quant/$sid/quant.sf" ]; then
      log "[$sid] FAILED (salmon=$sm minimap2=$mm)"; printf "%s\t%s\tfailed\tsalmon=%s minimap2=%s\n" "$sid" "$core" "$sm" "$mm" >> "$STATUS"
      rm -rf "$OUT/quant/$sid" "$OUT/driver_bam/$sid.bam" "$OUT/driver_bam/$sid.bam.bai"; continue
    fi
    log "[$sid] done ($(grep -o '"percent_mapped": [0-9.]*' "$OUT/quant/$sid/aux_info/meta_info.json"))"
    printf "%s\t%s\tdone\t\n" "$sid" "$core" >> "$STATUS"
  done
  pending=$(grep -c incomplete "$STATUS")
  if [ "$pending" -gt 0 ] && pgrep -f "rclone copy" >/dev/null; then
    log "$pending samples not downloaded yet; waiting 5 min"; sleep 300; continue
  fi
  break
done
rmdir "$OUT/tmp" 2>/dev/null || true
log "finished: $(grep -c $'\tdone\t' "$STATUS") done, $(grep -c $'\tfailed\t' "$STATUS") failed, $(grep -c incomplete "$STATUS") incomplete"
