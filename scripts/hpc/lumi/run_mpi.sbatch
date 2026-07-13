#!/bin/bash
#SBATCH --job-name=exsqs-mpi
#SBATCH --account=project_465002828
#SBATCH --partition=standard
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00

# Flow 2 -- multi-node rank-distributed islands (SPEC 8.2) on whole LUMI-C
# nodes: 8 ranks x 16 threads fills the 128 cores of an EPYC 7763 node along
# the measured scaling law (islands scale, inner threads saturate at ~16).
# Results are bit-identical for ANY rank count (T-MPI1), and serial/MPI
# checkpoints are interchangeable. On budget exhaustion (rc 3) the job
# resubmits itself up to MAX_SEGMENTS times.
#
# Submit from the repository root:
#     sbatch --export=ALL,CFG=configs/w70cr30_4x4x4.yaml,RUNDIR=runs/mpi \
#         scripts/hpc/lumi/run_mpi.sbatch
set -euo pipefail
source scripts/hpc/lumi/env.sh

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
CFG=${CFG:-configs/w70cr30_4x4x4.yaml}
RUNDIR=${RUNDIR:-runs/mpi_${SLURM_JOB_ID}}

RESUME=()
[[ -f "$RUNDIR/state.ckpt" ]] && RESUME=(--resume "$RUNDIR")

rc=0
srun ./build/exsqs_mpi "$CFG" \
    --set parallel.islands="${SLURM_NTASKS}" \
    --set parallel.omp_threads="${SLURM_CPUS_PER_TASK}" \
    --out "$RUNDIR" "${RESUME[@]}" "$@" || rc=$?

if [[ $rc -eq 3 && ${SEGMENT:-1} -lt ${MAX_SEGMENTS:-10} ]]; then
  echo "segment ${SEGMENT:-1}: budget exhausted -- chaining next segment"
  sbatch --export=ALL,CFG="$CFG",RUNDIR="$RUNDIR",SEGMENT=$(( ${SEGMENT:-1} + 1 )) \
      scripts/hpc/lumi/run_mpi.sbatch "$@"
  exit 0
fi
exit "$rc"
