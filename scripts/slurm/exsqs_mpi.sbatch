#!/bin/bash
#SBATCH --job-name=exsqs-mpi
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
##SBATCH --account=<grant>          # e.g. PLGrid allocation
##SBATCH --partition=<partition>

# Rank-distributed islands (SPEC 8.2): one rank per node, OpenMP inside each.
# Results are bit-identical for ANY rank count (T-MPI1); rank 0 writes
# state.ckpt, all ranks read it on resume (shared filesystem). Serial and MPI
# state files are interchangeable, so segments may move between 1 and N nodes.
set -euo pipefail
CFG=${CFG:-configs/w70cr30_4x4x4.yaml}
RUNDIR=${RUNDIR:-runs/slurm_${SLURM_JOB_ID}}
BIN=${BIN:-./build/exsqs_mpi}

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
RESUME=()
[[ -f "$RUNDIR/state.ckpt" ]] && RESUME=(--resume "$RUNDIR")

srun "$BIN" "$CFG" \
  --set parallel.islands="${SLURM_NTASKS}" \
  --set parallel.omp_threads="${SLURM_CPUS_PER_TASK}" \
  --out "$RUNDIR" "${RESUME[@]}"

# Self-resubmitting chain (uncomment; pairs with the exit-code contract):
# rc=$?
# if [[ $rc -eq 3 && ${SEGMENT:-1} -lt ${MAX_SEGMENTS:-10} ]]; then
#   sbatch --export=ALL,RUNDIR="$RUNDIR",SEGMENT=$(( ${SEGMENT:-1} + 1 )) "$0"
# fi
