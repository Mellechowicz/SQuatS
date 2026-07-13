#!/bin/bash
#SBATCH --job-name=exsqs-omp
#SBATCH --account=plgtopologyybbi2-cpu
#SBATCH --partition=plgrid
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2G
#SBATCH --time=24:00:00

# Flow 1 -- single-node islands+OpenMP production run with an in-job resume
# chain: run -> exit 3 (budget exhausted) -> --resume, until convergence
# (exit 0) or MAX_SEGMENTS. Each segment gets SEGMENT_WALL_S seconds of fresh
# wall budget; the chained trajectory is bit-identical to an uninterrupted
# run (T-K1). 16 cpus on purpose: the inner OpenMP layer saturates there --
# for more compute use the MPI flow, not more threads.
#
# Submit from the repository root:
#     sbatch --export=ALL,CFG=configs/w70cr30_4x4x4.yaml,RUNDIR=runs/myrun \
#         scripts/hpc/helios/run_omp_chain.sbatch
# Resubmitting with the same RUNDIR resumes where the job ended.
set -euo pipefail
source scripts/hpc/helios/env.sh

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
CFG=${CFG:-configs/w70cr30_4x4x4.yaml}
RUNDIR=${RUNDIR:-runs/omp_${SLURM_JOB_ID}}
SEGMENT_WALL_S=${SEGMENT_WALL_S:-14400}

rc=0
CFG="$CFG" RUNDIR="$RUNDIR" BIN=./build/exsqs MAX_SEGMENTS=${MAX_SEGMENTS:-10} \
  bash scripts/chain_resume.sh \
    --set parallel.omp_threads="${SLURM_CPUS_PER_TASK}" \
    --set evolution.max_wall_s="$SEGMENT_WALL_S" "$@" || rc=$?
[[ $rc -eq 3 ]] && echo "budget exhausted -- resubmit with RUNDIR=$RUNDIR to continue"
exit "$rc"
