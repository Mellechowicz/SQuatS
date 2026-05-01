#!/bin/bash
#SBATCH --job-name=exsqs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
##SBATCH --account=<grant>          # e.g. PLGrid allocation
##SBATCH --partition=<partition>

# Single-node OpenMP run (SPEC 8.1). Thread counts never change results; resume
# is bit-exact (T-K1), so requeue/chain freely. Exit 0 = e_tol reached,
# 3 = budget exhausted (state.ckpt kept for the next segment).
set -euo pipefail
CFG=${CFG:-configs/w70cr30_4x4x4.yaml}
RUNDIR=${RUNDIR:-runs/slurm_${SLURM_JOB_ID}}
BIN=${BIN:-./build/exsqs}

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
RESUME=()
[[ -f "$RUNDIR/state.ckpt" ]] && RESUME=(--resume "$RUNDIR")

"$BIN" "$CFG" \
  --set parallel.omp_threads="${SLURM_CPUS_PER_TASK}" \
  --out "$RUNDIR" "${RESUME[@]}"
