#!/bin/bash
# Render one campaign procedure into a NEW deployment directory.
#   - files committed byte-identical are copied unchanged
#   - *.template files are rendered by substituting exactly two site tokens:
#       @MELANOMA_ROOT@  (project root on the cluster, e.g. /scratch/<user>/greek_melanoma)
#       @PYTHON_BIN@     (absolute Python used by the deployed jobs)
#   - with --verify, every rendered/copied file is compared with the SHA-256 recorded for the deployed, qualified
#     bundle in deployed_hashes.tsv; any difference fails (exit 1). Excluded deployed files are listed, not created.
# Usage: MELANOMA_ROOT=... PYTHON_BIN=... bash render_templates.sh <procedure-dir-name> <new-output-dir> [--verify]
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
P=${1:?procedure directory name, e.g. procedure-20261127}; OUT=${2:?new output directory}; VERIFY=${3:-}
: "${MELANOMA_ROOT:?set MELANOMA_ROOT}" "${PYTHON_BIN:?set PYTHON_BIN}"
case "$MELANOMA_ROOT$PYTHON_BIN" in *"#"*|*"@"*) echo "ERROR: values must not contain '#' or '@'" >&2; exit 2 ;; esac
[ -d "$HERE/$P" ] || { echo "ERROR: no procedure $P under $HERE" >&2; exit 2; }
[ -z "$VERIFY" ] || [ "$VERIFY" = --verify ] || { echo "ERROR: third argument must be --verify" >&2; exit 2; }
test ! -e "$OUT" || { echo "ERROR: $OUT exists (refusing to overwrite)" >&2; exit 1; }
sha() { if command -v sha256sum >/dev/null; then sha256sum "$1" | cut -d' ' -f1; else shasum -a 256 "$1" | cut -d' ' -f1; fi; }
mkdir -p "$OUT/$P"
(cd "$HERE/$P" && find . -type f ! -name '*.pyc' | sed 's#^\./##' | sort) > "$OUT/.files"
while read -r f; do
  case "$f" in
    *.template)
      dst=${f%.template}; mkdir -p "$OUT/$P/$(dirname "$dst")"
      sed -e "s#@MELANOMA_ROOT@#$MELANOMA_ROOT#g" -e "s#@PYTHON_BIN@#$PYTHON_BIN#g" "$HERE/$P/$f" > "$OUT/$P/$dst"
      if grep -q '@[A-Z_]*@' "$OUT/$P/$dst"; then echo "ERROR: unrendered token in $dst" >&2; exit 1; fi ;;
    *) mkdir -p "$OUT/$P/$(dirname "$f")"; cp "$HERE/$P/$f" "$OUT/$P/$f" ;;
  esac
done < "$OUT/.files"
rm -f "$OUT/.files"
find "$OUT/$P" -type f \( -name '*.sh' -o -name '*.sbatch' \) -exec chmod 755 {} +
bad=0; n=0
if [ "$VERIFY" = --verify ]; then
  while IFS=$'\t' read -r proc path want form; do
    [ "$proc" = "$P" ] || continue
    if [ "$form" = excluded ]; then echo "EXCLUDED  $path (omitted in this checkout; see deployed_hashes.tsv)"; continue; fi
    n=$((n + 1))
    if [ ! -f "$OUT/$P/$path" ]; then echo "MISSING   $path"; bad=1
    elif [ "$(sha "$OUT/$P/$path")" = "$want" ]; then echo "MATCH     $path ($form)"
    else echo "DIFFERS   $path ($form)"; bad=1; fi
  done < <(tail -n +2 "$HERE/deployed_hashes.tsv")
  [ "$bad" -eq 0 ] && echo "VERIFY PASS: $n included files match their deployed SHA-256" || { echo "VERIFY FAIL"; exit 1; }
fi
(cd "$OUT/$P" && find . -type f | sed 's#^\./##' | sort | while read -r f; do printf '%s  %s\n' "$(sha "$f")" "$f"; done) > "$OUT/$P/RENDERED.sha256"
echo "rendered $P -> $OUT/$P"
