#!/bin/bash
# Run every campaign regression suite (synthetic fixtures and stub tools only; no cluster, no real data).
# Usage: bash run_campaign_tests.sh [python3]
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; PY=${1:-python3}
export RQC_TOOLS="$HERE/../rivanna_qc_tools.py"   # the committed, qualified r3 tools (required by test_campaign)
fail=0
for d in procedure-20261126 procedure-20261127 procedure-20261127r1 .; do
  echo "== $d"
  (cd "$HERE/$d" && "$PY" -I -m unittest discover -s tests -p 'test_*.py' 2>&1 | grep -v -e ResourceWarning -e 'csv.DictReader(open' | tail -4) || true
  (cd "$HERE/$d" && "$PY" -I -m unittest discover -s tests -p 'test_*.py' >/dev/null 2>&1) || fail=1
done
[ "$fail" -eq 0 ] && echo "ALL CAMPAIGN SUITES PASS" || { echo "SOME CAMPAIGN SUITE FAILED"; exit 1; }
