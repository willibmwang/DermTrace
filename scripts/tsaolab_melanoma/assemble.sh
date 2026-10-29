#!/usr/bin/env bash
# Assemble the DermTrace dataset from whatever run_quant.sh has finished so far: read-level QC, driver variants,
# expression matrices + metadata, and (once every sample is done) the donor/genotype check. Idempotent.
#
# Usage: assemble.sh <quant_out_dir> <sample_sheet.csv> <raw_data_dir> <dest_dir (data/processed/tsaolab_melanoma)>
set -euo pipefail
QDIR=$1; SHEET=$2; RAW=$3; DEST=$4
HERE="$(cd "$(dirname "$0")" && pwd)"
PY=${PY:-/opt/anaconda3/bin/python3}
mkdir -p "$DEST/qc" "$DEST/variants" "$DEST/salmon_quant"
"$PY" -I "$HERE/run_qc_summary.py" "$DEST/qc/run_qc.tsv" "$SHEET" "$QDIR" "$RAW" > /dev/null
JOBS=${JOBS:-3} bash "$HERE/call_variants.sh" "$QDIR" "$SHEET" "$QDIR/variants" > "$QDIR/variants.log" 2>&1
cp "$QDIR"/variants/{driver_variants_all,driver_variants_by_sample,driver_variants_somatic_candidates,hotspot_summary}.tsv "$DEST/variants/"
"$PY" -I "$HERE/build_dermtrace.py" "$DEST" "$SHEET" "$QDIR" --runqc "$DEST/qc/run_qc.tsv" --hotspots "$DEST/variants/hotspot_summary.tsv" > "$QDIR/build.log"
for d in "$QDIR"/quant/GM*/; do  # gzipped transcript-level salmon output, local paths scrubbed
  sid=$(basename "$d"); [ -s "$d/quant.sf" ] || continue
  mkdir -p "$DEST/salmon_quant/$sid/aux_info"
  [ -s "$DEST/salmon_quant/$sid/quant.sf.gz" ] || gzip -c "$d/quant.sf" > "$DEST/salmon_quant/$sid/quant.sf.gz"
  sed -E 's#/Users/[^"]*/HST Tsao Lab#<LAB>#g' "$d/cmd_info.json" > "$DEST/salmon_quant/$sid/cmd_info.json"
  cp "$d/lib_format_counts.json" "$DEST/salmon_quant/$sid/"; cp "$d/aux_info/meta_info.json" "$DEST/salmon_quant/$sid/aux_info/"
done
echo "assembled $(grep -c . "$DEST/qc/run_qc.tsv" | awk '{print $1-1}') samples -> $DEST"
