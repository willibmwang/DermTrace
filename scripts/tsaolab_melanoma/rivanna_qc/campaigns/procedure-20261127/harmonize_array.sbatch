#!/bin/bash
#SBATCH --job-name=melanoma_harmonize
#SBATCH --partition=standard
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=harmonize_%A_%a.out
#SBATCH --error=harmonize_%A_%a.err
# Array wrapper: task i quantifies line i+1 of the checksummed specimen list with the UNCHANGED r3
# quantify_specimen.sbatch (resources above are copied from its qualified #SBATCH header). Per-specimen locks,
# fail-closed checks and provenance are those of quantify_specimen.sbatch.
set -euo pipefail
: "${SPECIMEN_LIST:?}" "${SPECIMEN_LIST_SHA256:?}" "${RQC_DIR:?}" "${SLURM_ARRAY_TASK_ID:?}"
[ "$(sha256sum "$SPECIMEN_LIST" | cut -d' ' -f1)" = "$SPECIMEN_LIST_SHA256" ] || { echo "ERROR: specimen list checksum mismatch" >&2; exit 1; }
n=$(grep -c . "$SPECIMEN_LIST"); [ "$SLURM_ARRAY_TASK_ID" -lt "$n" ] || { echo "ERROR: task $SLURM_ARRAY_TASK_ID >= $n specimens" >&2; exit 1; }
SPECIMEN=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SPECIMEN_LIST"); export SPECIMEN
echo "ARRAY ${SLURM_ARRAY_JOB_ID:-?}_${SLURM_ARRAY_TASK_ID} SPECIMEN $SPECIMEN SLURM_JOB_ID ${SLURM_JOB_ID:-?}"
exec bash "$RQC_DIR/quantify_specimen.sbatch"
