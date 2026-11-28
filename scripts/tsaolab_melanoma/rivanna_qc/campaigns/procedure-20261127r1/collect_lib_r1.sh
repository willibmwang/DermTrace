# Step 17 recovery library (r1). Sourced; defines functions only.

# show_sacct JOBID [MAX_LINES]: print at most MAX_LINES (default 40) of Slurm accounting for JOBID.
# The complete sacct output is captured first, so sacct can never receive SIGPIPE, and it is displayed through a
# here-string (no pipe, so nothing downstream can raise SIGPIPE either). A genuine sacct failure (non-zero exit)
# is reported with its first lines and returns 1; callers must treat that as fatal.
show_sacct() {
  local job=$1 max=${2:-40} out rc n
  out=$(sacct -j "$job" -o JobID,JobName,State,Elapsed,MaxRSS,ExitCode 2>&1) && rc=0 || rc=$?
  if [ "$rc" -ne 0 ]; then
    echo "STOP sacct failed for job $job (exit $rc): $(sed -n '1,5p' <<< "$out" | tr '\n' ' ')" >&2
    return 1
  fi
  n=$(awk 'END {print NR}' <<< "$out")
  sed -n "1,${max}p" <<< "$out"
  if [ "$n" -gt "$max" ]; then echo "... ($((n - max)) further sacct lines not shown; all $n lines were read)"; fi
  return 0
}

# install_partial_trap: on exit, report an incomplete build left in $PARTIAL_DIR and exit with the ORIGINAL status.
# Set PARTIAL_DIR="" once the build has been adopted.
on_exit_preserve_status() {
  local rc=$?
  trap - EXIT
  if [ -n "${PARTIAL_DIR:-}" ] && [ -d "${PARTIAL_DIR:-}" ]; then
    echo "NOTE: incomplete build left at $PARTIAL_DIR (inspect, then remove before rerunning); original exit status $rc" >&2
  fi
  exit "$rc"
}
install_partial_trap() { PARTIAL_DIR=${1:-}; trap on_exit_preserve_status EXIT; }
