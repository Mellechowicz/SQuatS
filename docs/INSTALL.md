# Installing and running SQuatS / EXSQS

Everything the engine needs at build time is either a compiler feature or a
pinned third-party source fetched by CMake itself. There is nothing to install
system-wide beyond the toolchain.

## Requirements

Build (mandatory):

| what              | constraint                                   | why |
|-------------------|----------------------------------------------|-----|
| CMake             | >= 3.20 **and < 4.0**                        | cmake 4.x rejects the `cmake_minimum_required` of the bundled yaml-cpp 0.8.0 |
| C++17 compiler    | GCC >= 9 (tested 12.2, 14.2); OpenMP support | engine + inner parallel layer |
| git + network     | **at configure time only**                   | `FetchContent` pins and downloads spglib v2.7.0, yaml-cpp 0.8.0, Catch2 v3.5.4 |

Optional:

| what              | constraint                                   | enables |
|-------------------|----------------------------------------------|---------|
| MPI (C++)         | any; tested OpenMPI 4.1/5.0, cray-mpich      | the `exsqs_mpi` rank-distributed driver (built automatically when FindMPI succeeds) |
| Python >= 3.9 with `pymatgen`, `numpy` | venv recommended        | validation gates T-V1/T-X2 and the `tools/py/` interop scripts — never needed to build or run the engine |

The configure step must run on a machine with internet access (login node,
workstation). After the first configure the sources live in `build/_deps/`
and rebuilds are fully offline.

## Build and verify

    cmake -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build -j
    bash tools/run_all_tests.sh          # full 24-gate matrix, ~2 min
    ./build/exsqs configs/smoke_sc27.yaml --set evolution.max_generations=8 --out /tmp/s

`SKIP_MPI=1` / `SKIP_V1=1` / `SKIP_E2E=1` skip the gates that need MPI /
pymatgen / time. Exit codes everywhere: 0 = converged, 3 = budget exhausted
(resumable, **not** an error), 1 = error.

## Running efficiently

- **Determinism is free** — one trajectory per config+seed, bit-identical
  across OpenMP thread counts, MPI rank counts and checkpoint/resume chains
  (T-R1, T-K1, T-MPI1). Chunk long runs into as many budget segments as the
  queue requires; the chained result equals the uninterrupted one.
- **Scale via islands (MPI ranks), not inner threads.** The inner OpenMP
  layer saturates near 16 threads (measured Amdahl fit on N=250: p = 0.787,
  S_inf = 4.7); beyond that, extra cores are better spent on more islands.
  On a whole node prefer e.g. 8 ranks x 16 threads over 1 x 128.
- **Pin threads**: `OMP_PLACES=cores OMP_PROC_BIND=close` (the SLURM scripts
  below set this).
- **Budget by cell size** — per-evaluation cost grows like `|Pi| * N`
  (recorded T-B1 law), so large cells need proportionally larger
  `evolution.max_generations` or wall budgets.
- Resume chains change only budget caps between segments:
  `scripts/chain_resume.sh` loops run -> exit 3 -> `--resume` until exit 0.

## HPC sites

Ready-made scripts live under `scripts/hpc/<site>/`; each site has an
environment file, a compile driver and three run flows. All of them are
submitted **from the repository root**, and all honour the exit-0/3 contract:

| script                    | flow |
|---------------------------|------|
| `env.sh`                  | module loads + pinning (source it, don't run it) |
| `compile.sh`              | configure + build + smoke; first run on a **login node** (FetchContent needs internet), afterwards `sbatch`-able |
| `run_omp_chain.sbatch`    | flow 1: single-node islands+OpenMP with in-job resume chain |
| `run_mpi.sbatch`          | flow 2: multi-node MPI islands, self-resubmitting on budget exhaustion |
| `run_ladder_array.sbatch` | flow 3: campaign scan, one array task per config with per-task resume chains |

Common knobs (via `--export=ALL,...`): `CFG`, `RUNDIR` (same value resumes),
`MAX_SEGMENTS`, `SEGMENT_WALL_S`, `CONFIG_GLOB`/`RUNDIR_BASE` (array flow).

### LUMI-C (`scripts/hpc/lumi/`, account `project_465002828`)

- HPE Cray EX, 2x AMD EPYC 7763 = 128 cores/node, 256 GiB (standard nodes).
- Modules: `LUMI/24.03 partition/C` + `cpeGNU/24.03` (GNU compilers,
  cray-mpich linked automatically by the `cc`/`CC` wrappers) +
  `buildtools/24.03` (cmake 3.x). Versions move with system upgrades —
  `module avail LUMI`.
- Partitions used: `debug` (compile, 30 min), `small` (by-core: OMP chain and
  array flows, 16 cores each, max 3 days), `standard` (whole-node MPI flow,
  8 ranks x 16 threads per node, max 2 days).
- MPI launches with `srun` (no mpirun on LUMI).

### Helios, Cyfronet (`scripts/hpc/helios/`, account `plgtopologyybbi2-cpu`)

- HPE Cray EX, 2x AMD EPYC 9654 = 192 cores/node, 384 GiB (`plgrid`
  partition); PLGrid accounts need the `-cpu` suffix — the bare grant name
  fails.
- Modules (hierarchical lmod, compiler -> MPI -> rest):
  `GCC/13.2.0 OpenMPI/5.0.3 CMake/3.27.6` — check drifted versions with
  `module spider CMake`; cmake must stay >= 3.20 and < 4.0.
- Partitions used: `plgrid-now` (compile, high-priority, 12 h cap),
  `plgrid` (all three flows, 72 h cap); `plgrid-long` exists for >72 h runs
  with special grant permission.
- MPI flow fills nodes with 12 ranks x 16 threads; if `srun` lacks PMIx
  support for OpenMPI, swap in `mpiexec -n $SLURM_NTASKS`.
