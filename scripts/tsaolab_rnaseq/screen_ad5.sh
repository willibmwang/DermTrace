#!/usr/bin/env bash
# Cell-identity screen: count reads from human adenovirus 5 (AC_000008.1) in the first N read pairs
# of each sample's L001 FASTQs. HEK293 carries an integrated Ad5 E1A/E1B fragment, so E1 reads mark
# HEK293(-derived) cells; other lines should have ~0.
#
# Usage: screen_ad5.sh <sample_sheet.csv> <raw_data_dir> <out.tsv>   (env: REF_DIR, N_PAIRS=1000000, JOBS=4)
set -euo pipefail
SHEET=$1; RAW=$2; OUT=$3
IDX="${REF_DIR:?}/adenovirus5_AC_000008.1/salmon_index"
N=${N_PAIRS:-1000000}
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT

one() {  # <sid> <qb>
  local sid=$1 qb=$2 r1 r2 t
  r1=$(ls "$RAW/$qb"/"${qb}"_S*_L001_R1_001.fastq.gz); r2=${r1/_R1_001/_R2_001}
  t="$TMP/$sid"; mkdir -p "$t"
  # gzcat exits on SIGPIPE once head has enough reads
  { gzcat "$r1" || true; } | head -n $((N * 4)) > "$t/1.fq"
  { gzcat "$r2" || true; } | head -n $((N * 4)) > "$t/2.fq"
  salmon quant -i "$IDX" -l A -1 "$t/1.fq" -2 "$t/2.fq" -p 1 --writeMappings="$t/m.sam" -o "$t/q" >/dev/null 2>&1
  # E1A/E1B = Ad5 nt ~500-4100
  awk -v s="$sid" -v n="$(( $(wc -l < "$t/1.fq") / 4 ))" \
    '!/^@/ {tot++; if ($4 >= 500 && $4 <= 4100) e1++} END {printf "%s\t%d\t%d\t%d\n", s, n, tot/2, e1/2}' "$t/m.sam"
  rm -rf "$t"
}
export -f one; export RAW TMP IDX N

printf "sample_id\tpairs_screened\tad5_pairs\tad5_e1_pairs\n" > "$OUT"
tail -n +2 "$SHEET" | cut -d, -f1,2 | tr , ' ' | xargs -P "${JOBS:-4}" -n 2 bash -c 'one "$0" "$1"' | sort -V >> "$OUT"
