#!/bin/bash
# Step 16: one bounded Slurm array for the 28 remaining specimens (MAX_PARALLEL concurrent tasks, default 8).
set -euo pipefail
. "$(dirname "$0")/common.sh"
[ -s "$HC" ] || { echo "STOP $HC missing: run step 15 first"; exit 1; }
# shellcheck disable=SC1090
source "$HC"
RUN=$(ls -d "$PROVENANCE_DIR"/preflight_* 2>/dev/null | sort | tail -1 || true)
LIST=$RUN/specimens_28.txt
[ -n "$RUN" ] && [ -s "$RUN/PASS" ] && [ "$(sha256sum "$LIST" | cut -d' ' -f1)" = "$(cat "$RUN/PASS")" ] || { echo "STOP latest step 15 preflight did not PASS (${RUN:-none})"; exit 1; }
[ "$(grep -c . "$LIST")" -eq 28 ] || { echo "STOP specimen list is not 28 lines"; exit 1; }
test ! -e "$QUANT_ROOT" || { echo "STOP $QUANT_ROOT exists"; exit 1; }
test ! -e "$LOG/array_job.txt" || { echo "STOP already submitted: $(cat "$LOG/array_job.txt")"; exit 1; }
mkdir -p "$LOG"; cd "$LOG"
ACCT=$(sacct -j 23402243 -X -n -o Account%40 2>/dev/null | awk 'NF {print $1; exit}')
J=$(sbatch --parsable ${ACCT:+--account="$ACCT"} --array=0-27%"${MAX_PARALLEL:-8}" \
      --export=ALL,SPECIMEN_LIST="$LIST",SPECIMEN_LIST_SHA256="$(cat "$RUN/PASS")" "$PD/harmonize_array.sbatch")
( set -C; printf '%s\t%s\n' "$J" "$LIST" > "$LOG/array_job.txt" )
echo "submitted array job $J (28 tasks, at most ${MAX_PARALLEL:-8} concurrent)"
echo "monitor: sacct -j $J -o JobID,JobName,State,Elapsed,MaxRSS,ExitCode"
