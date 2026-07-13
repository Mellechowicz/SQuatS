#!/bin/bash
#SBATCH --job-name=exsqs-ladder
#SBATCH --account=project_465002828
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2G
#SBATCH --time=12:00:00
#SBATCH --array=0-9

# Flow 3 -- campaign scan: one array task per config, alphabetical order, each
# with its own in-job resume chain. The default glob covers the 10-config
# WCrTaFePb set; match --array to the glob count when overriding, e.g. the
# 16-cell K=5 ladder:
#     sbatch --array=0-15 --export=ALL,CONFIG_GLOB='configs/hea*_g0.yaml' \
#         scripts/hpc/lumi/run_ladder_array.sbatch
# Task outputs land in $RUNDIR_BASE/<config-stem>/ and resubmission resumes.
# Submit from the repository root.
set -euo pipefail
source scripts/hpc/lumi/env.sh

CONFIG_GLOB=${CONFIG_GLOB:-configs/wcrtafepb_*.yaml}
mapfile -t CFGS < <(compgen -G "$CONFIG_GLOB" | sort)
CFG=${CFGS[$SLURM_ARRAY_TASK_ID]:?array index ${SLURM_ARRAY_TASK_ID} exceeds the ${#CFGS[@]}-config glob}
RUNDIR=${RUNDIR_BASE:-runs}/$(basename "$CFG" .yaml)
echo "array task ${SLURM_ARRAY_TASK_ID}: $CFG -> $RUNDIR"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
rc=0
CFG="$CFG" RUNDIR="$RUNDIR" BIN=./build/exsqs MAX_SEGMENTS=${MAX_SEGMENTS:-10} \
  bash scripts/chain_resume.sh \
    --set parallel.omp_threads="${SLURM_CPUS_PER_TASK}" \
    --set evolution.max_wall_s="${SEGMENT_WALL_S:-3600}" "$@" || rc=$?
[[ $rc -eq 3 ]] && echo "budget exhausted -- resubmit this array index to continue"
exit "$rc"
