#!/usr/bin/env bash
#SBATCH --job-name=keloid-narrative-audit
#SBATCH --partition=gpu
#SBATCH --account=dijk
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=results/slurm/%x-%j.out
#SBATCH --error=results/slurm/%x-%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export TOKENIZERS_PARALLELISM=false
export HF_HUB_OFFLINE=1
KB_PYTHON="${KB_PYTHON:-/gpfs/radev/home/wbw7/.conda/envs/encdec_llm/bin/python}"
KB_MODEL="${KB_MODEL:-/gpfs/radev/home/wbw7/.cache/huggingface/hub/models--Qwen--Qwen3-1.7B/snapshots/70d244cc86ccca08cf5af4e1e306ecf908b1ad5e}"
"${KB_PYTHON}" scripts/audit_copilot_faithfulness.py --arm llm --repeats 3 \
  --local-model "${KB_MODEL}" --batch-size 16 \
  --output results/reviewer_revision/narrative_audit/qwen3_1_7b_v1
