#!/usr/bin/env bash
# Cell-line genotype check: align the first N read pairs of each human sample's L001 FASTQs to a small
# transcript reference (genotype.py ref), call SNPs jointly with bcftools, and write <out_dir>/snps.tsv
# (per-sample genotype and allele depths) for genotype.py compare.
# Extra samples, e.g. public RNA-seq of the expected cell line, are given as <name>=<R1 FASTQ path or URL>;
# R2 is the same path with _1.fastq / _R1_ replaced by _2.fastq / _R2_.
#
# Usage: genotype_check.sh <mini_ref.fa> <out_dir> <sample_sheet.csv>=<raw_data_dir> [...] [<name>=<R1>] [...]
#        env: GENO_ENV (conda env with minimap2, samtools, bcftools), N_PAIRS=2000000, JOBS=3
set -euo pipefail
REF=$1; OUT=$2; shift 2
BIN="${GENO_ENV:?set GENO_ENV to a conda env with minimap2, samtools and bcftools}/bin"
N=${N_PAIRS:-2000000}
mkdir -p "$OUT/bam"
"$BIN/samtools" faidx "$REF"
"$BIN/minimap2" -d "$OUT/ref.mmi" "$REF" 2>/dev/null

one() {  # <name> <R1>
  local name=$1 r1=$2 r2
  r2=$(sed -e 's/_R1_/_R2_/' -e 's/_1\.fastq/_2.fastq/' <<<"$r1")
  # streams stop once head has enough reads
  read_n() { case $1 in http*) curl -s "$1" ;; *) cat "$1" ;; esac | { gzip -dc || true; } | head -n $((N * 4)); }
  "$BIN/minimap2" -t 4 -ax sr -R "@RG\tID:$name\tSM:$name" "$OUT/ref.mmi" <(read_n "$r1") <(read_n "$r2") 2>/dev/null \
    | "$BIN/samtools" view -u -q 20 -F 0x904 - | "$BIN/samtools" sort -o "$OUT/bam/$name.bam" - 2>/dev/null
  "$BIN/samtools" index "$OUT/bam/$name.bam"
}
export -f one; export BIN OUT N

jobs=()
for a in "$@"; do
  if [[ $a == *.csv=* ]]; then  # human samples of a sample sheet
    sheet=${a%%=*}; raw=${a#*=}
    while IFS=, read -r sid qb _name _project species _rest; do
      [ "$species" = human ] && jobs+=("$sid" "$(ls "$raw/$qb"/"${qb}"_S*_L001_R1_001.fastq.gz)")
    done < <(tail -n +2 "$sheet")
  else
    jobs+=("${a%%=*}" "${a#*=}")
  fi
done
printf '%s\0' "${jobs[@]}" | xargs -0 -P "${JOBS:-3}" -n 2 bash -c 'one "$0" "$1"'

ls "$OUT"/bam/*.bam | sort -V > "$OUT/bams.txt"
"$BIN/bcftools" mpileup -f "$REF" -b "$OUT/bams.txt" -a AD,DP -d 2000 -q 20 -Q 20 2>/dev/null \
  | "$BIN/bcftools" call -mv 2>/dev/null | "$BIN/bcftools" view -v snps -i 'QUAL>50' \
  | "$BIN/bcftools" query -H -f '%CHROM\t%POS[\t%GT:%AD]\n' > "$OUT/snps.tsv"
echo "$(($(wc -l < "$OUT/snps.tsv") - 1)) SNPs in $(wc -l < "$OUT/bams.txt" | tr -d ' ') samples -> $OUT/snps.tsv"
