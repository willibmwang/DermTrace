#!/usr/bin/env bash
#SBATCH --job-name=keloid-narrative-v2
#SBATCH --partition=gpu
#SBATCH --account=dijk
#SBATCH --gres=gpu:1
#SBATCH --constraint="a40|a100|l40s|h100|h200"
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=results/slurm/%x-%j.out
#SBATCH --error=results/slurm/%x-%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export TOKENIZERS_PARALLELISM=false
export HF_HUB_OFFLINE=1
KB_PYTHON="${KB_PYTHON:-/gpfs/radev/home/wbw7/.conda/envs/encdec_llm/bin/python}"
KB_MODEL="${KB_MODEL:-/gpfs/radev/project/dijk/wbw7/hf_models/models--Qwen--Qwen3-8B/snapshots/b968826d9c46dd6066d109eabc6255188de91218}"
KB_OUTPUT="${KB_OUTPUT:-results/reviewer_revision/narrative_audit/qwen3_8b_v2}"
"${KB_PYTHON}" scripts/audit_copilot_faithfulness.py --arm llm --repeats "${KB_REPEATS:-3}" ${KB_LIMIT:+--limit "${KB_LIMIT}"} \
  --local-model "${KB_MODEL}" --model-name Qwen/Qwen3-8B --narrative-format prose \
  --max-new-tokens 160 --batch-size 16 --output "${KB_OUTPUT}"
