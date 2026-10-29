#!/usr/bin/env bash
# From the driver-reference BAMs of run_quant.sh: per sample, mark and drop PCR duplicates (fixmate + markdup; FFPE
# libraries here are 60-95% duplicate, and one over-amplified molecule otherwise shows up as a 100% "variant"), and keep
# only reads that align over >= 60 nt with <= 4% mismatches (the 3' UTRs of several driver transcripts hold repeat
# elements that attract reads from all over the genome). Then one joint bcftools mpileup over the deduplicated BAMs and
# one over the duplicate-containing BAMs, and variants.py call turns both into variant tables.
#
# Usage: call_variants.sh <quant_out_dir> <sample_sheet.csv> <out_dir>
#   env: VARCALL_ENV=/opt/anaconda3/envs/varcall  DRIVER_REF=<references>/melanoma_driver_ref/driver_ref.fa  JOBS=4
set -euo pipefail
QDIR=$1; SHEET=$2; OUT=$3
HERE="$(cd "$(dirname "$0")" && pwd)"
REFS="$(cd "$HERE/../../../references" && pwd)"
BIN="${VARCALL_ENV:-/opt/anaconda3/envs/varcall}/bin"
DRIVER_REF=${DRIVER_REF:-$REFS/melanoma_driver_ref/driver_ref.fa}
PY=${PY:-/opt/anaconda3/bin/python3}
FILTER='(qlen-sclen)>=60 && [NM]<=0.04*(qlen-sclen)'
mkdir -p "$OUT/bam"
"$BIN/samtools" faidx "$DRIVER_REF"

one() {  # <bam>: writes bam/<sid>.dedup.bam and bam/<sid>.raw.bam (both identity-filtered)
  local bam=$1 sid; sid=$(basename "$bam" .bam)
  [ -s "$OUT/bam/$sid.dedup.bam.bai" ] && [ "$OUT/bam/$sid.dedup.bam.bai" -nt "$bam" ] && return 0
  "$BIN/samtools" sort -n -@2 "$bam" 2>/dev/null | "$BIN/samtools" fixmate -m - - | "$BIN/samtools" sort -@2 - 2>/dev/null \
    | "$BIN/samtools" markdup -r -s - - 2> "$OUT/bam/$sid.markdup.log" \
    | "$BIN/samtools" view -b -e "$FILTER" -o "$OUT/bam/$sid.dedup.bam" -
  "$BIN/samtools" index "$OUT/bam/$sid.dedup.bam"
  "$BIN/samtools" view -b -e "$FILTER" -o "$OUT/bam/$sid.raw.bam" "$bam"
  "$BIN/samtools" index "$OUT/bam/$sid.raw.bam"
}
export -f one; export BIN OUT FILTER
ls "$QDIR"/driver_bam/*.bam | sort -V | tr '\n' '\0' | xargs -0 -P "${JOBS:-4}" -n 1 bash -c 'one "$0"'

for kind in dedup raw; do
  ls "$OUT"/bam/*.$kind.bam | sort -V > "$OUT/bams_$kind.txt"
  # -d: per-file depth cap; -B: no BAQ (it would down-weight mismatches near the ends of these short inserts)
  "$BIN/bcftools" mpileup -f "$DRIVER_REF" -b "$OUT/bams_$kind.txt" -a AD,DP -d 100000 -q 1 -Q 20 -B -Ou 2>/dev/null \
    | "$BIN/bcftools" query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%INDEL[\t%AD]\n' > "$OUT/sites_$kind.tsv"
  echo "$kind: $(wc -l < "$OUT/sites_$kind.tsv") pileup lines over $(wc -l < "$OUT/bams_$kind.txt" | tr -d ' ') BAMs"
done
"$PY" -I "$HERE/variants.py" call "$(dirname "$DRIVER_REF")" "$OUT/sites_dedup.tsv" "$OUT/sites_raw.tsv" "$SHEET" "$OUT"
