#!/usr/bin/env bash
# Trim (fastp) + quantify (salmon) every sample in a sample sheet whose 4 FASTQs
# (each lane in $LANES x R1/R2) are all present. Re-runnable: finished samples are skipped,
# incomplete samples are reported and skipped.
#
# Usage: [LANES=L001] run_quant.sh <sample_sheet.csv> <raw_data_dir> <out_dir>
# Requires the conda env "rnaseq" (salmon 1.10.3, fastp).
set -euo pipefail

SHEET=$1; RAW=$2; OUT=$3
REF="${REF_DIR:?set REF_DIR to the folder holding gencode_human_v46/ and gencode_mouse_vM35/}"
THREADS=${THREADS:-8}
LANES=${LANES:-"L001 L002"}  # lanes every sample must have; 1030390 was sequenced on L001 only
mkdir -p "$OUT/quant" "$OUT/fastp" "$OUT/tmp"
STATUS="$OUT/sample_status.tsv"
printf "sample_id\tqb_id\tstatus\tdetail\n" > "$STATUS"

tail -n +2 "$SHEET" | while IFS=, read -r sid qb name project species _rest; do
  dir="$RAW/$qb"
  missing=()
  for lane in $LANES; do for rd in R1 R2; do
    ls "$dir"/"${qb}"_S*_"${lane}"_"${rd}"_001.fastq.gz >/dev/null 2>&1 || missing+=("${lane}_${rd}")
  done; done
  if [ ${#missing[@]} -gt 0 ]; then
    printf "%s\t%s\tincomplete\tmissing %s\n" "$sid" "$qb" "${missing[*]}" >> "$STATUS"
    echo "[$sid] SKIP incomplete: missing ${missing[*]}"; continue
  fi
  if [ -s "$OUT/quant/$sid/quant.sf" ]; then
    printf "%s\t%s\tdone\tprevious run\n" "$sid" "$qb" >> "$STATUS"
    echo "[$sid] already quantified"; continue
  fi
  case $species in
    human) idx="$REF/gencode_human_v46/salmon_index" ;;
    mouse) idx="$REF/gencode_mouse_vM35/salmon_index" ;;
    *) echo "unknown species $species"; exit 1 ;;
  esac

  echo "[$sid] fastp + salmon ($species)"
  t1=(); t2=()
  for lane in $LANES; do
    r1=$(ls "$dir"/"${qb}"_S*_"${lane}"_R1_001.fastq.gz); r2=${r1/_R1_001/_R2_001}
    o1="$OUT/tmp/${sid}_${lane}_R1.fq.gz"; o2="$OUT/tmp/${sid}_${lane}_R2.fq.gz"
    fastp -i "$r1" -I "$r2" -o "$o1" -O "$o2" --detect_adapter_for_pe -w "$THREADS" \
      -j "$OUT/fastp/${sid}_${lane}.json" -h "$OUT/fastp/${sid}_${lane}.html" 2> "$OUT/fastp/${sid}_${lane}.log" < /dev/null
    t1+=("$o1"); t2+=("$o2")
  done
  salmon quant -i "$idx" -l A -1 "${t1[@]}" -2 "${t2[@]}" -p "$THREADS" \
    --validateMappings --gcBias --seqBias -o "$OUT/quant/$sid" > "$OUT/quant/$sid.log" 2>&1 < /dev/null
  rm -f "${t1[@]}" "${t2[@]}"
  printf "%s\t%s\tdone\t\n" "$sid" "$qb" >> "$STATUS"
done
rmdir "$OUT/tmp" 2>/dev/null || true
