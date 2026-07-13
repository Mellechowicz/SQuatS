# LUMI-C build/run environment. Source from the repository root:
#     source scripts/hpc/lumi/env.sh
# Stack versions move with system upgrades -- check `module avail LUMI`.

module load LUMI/24.03 partition/C
module load cpeGNU/24.03        # GNU compilers; cray-mpich linked by the cc/CC wrappers
module load buildtools/24.03    # cmake 3.x + ninja (cmake must stay < 4.0, see docs/INSTALL.md)

export CC=cc CXX=CC             # Cray wrappers -- MPI comes for free
export OMP_PLACES=cores OMP_PROC_BIND=close
