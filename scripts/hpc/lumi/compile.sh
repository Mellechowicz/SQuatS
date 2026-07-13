#!/bin/bash
#SBATCH --job-name=exsqs-build
#SBATCH --account=project_465002828
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:30:00

# Configure + build + smoke-run on LUMI-C.
#
# FetchContent downloads the pinned third-party sources AT CONFIGURE TIME, so
# the FIRST run must happen where internet exists -- the login node:
#     bash scripts/hpc/lumi/compile.sh
# After that, build/_deps/ is populated and clean rebuilds may go through the
# queue instead:  sbatch scripts/hpc/lumi/compile.sh
# Always submit/run from the repository root.
set -euo pipefail
source scripts/hpc/lumi/env.sh

case "$(cmake --version | head -1)" in
  *" 4."*) echo "cmake 4.x cannot configure the bundled yaml-cpp 0.8.0 -- use a 3.x cmake" >&2; exit 1 ;;
esac

cmake -B build -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=CC -DMPI_CXX_COMPILER=CC
cmake --build build -j "${SLURM_CPUS_PER_TASK:-16}"

rc=0
./build/exsqs configs/smoke_sc27.yaml --set evolution.max_generations=4 \
    --out "${TMPDIR:-/tmp}/exsqs_smoke_$$" || rc=$?
if [[ $rc -eq 0 || $rc -eq 3 ]]; then
  echo "LUMI build OK (smoke rc=$rc; 3 = budget exhausted, normal)"
else
  echo "smoke run FAILED rc=$rc" >&2
  exit 1
fi
