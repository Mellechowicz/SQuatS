# Helios (Cyfronet) build/run environment. Source from the repository root:
#     source scripts/hpc/helios/env.sh
# lmod is hierarchical: compiler first, then MPI, then the rest.
# Check exact versions after system upgrades with `module spider CMake` etc.

module load GCC/13.2.0
module load OpenMPI/5.0.3
module load CMake/3.27.6 2>/dev/null || module load CMake   # must stay >= 3.20 and < 4.0

export CC=gcc CXX=g++
export OMP_PLACES=cores OMP_PROC_BIND=close
