# SQuatS — Symmetric QUasirandom ATomic Structures

**SQuatS** (**S**ymmetric **QU**asirandom **AT**omic **S**tructures) are
special quasirandom structures that keep a deliberate residue of crystal
symmetry: instead of scrambling the cell into P1, the generator trades
correlation error `E_pure` against the number `D` of symmetry-inequivalent
displacements a downstream phonon/DFT workflow must compute, minimizing
`E_obj = E_pure * D^gamma`. Classical SQS is the `gamma = 0` member of the
family.

This repository is the SQuatS project. Its engine, **EXSQS**, is a C++17
implementation of the extinction evolutionary algorithm from *"On generating
Special Quasirandom Structures: Optimization for the DFT computational
efficiency"* (arXiv:2602.10872) — at `gamma = 1` on the reference system it
finds cells needing **six-fold fewer displacements** than a random (P1)
decoration at controlled correlation cost.

Design guarantees, all machine-verified (see Testing): one deterministic
trajectory per config+seed — **bitwise identical** across OpenMP thread
counts, MPI rank counts, and checkpoint/resume chains; a provable
quantization floor `E_floor` making convergence targets floor-relative
(SPEC section 4.1); loud-fail configuration; and a self-auditing coherence
gate (spec <-> code <-> tests <-> docs).

Version 1.7.0. `docs/SPEC.md` (v1.7) is the normative specification;
`CHANGELOG.md` summarizes releases.

## Build

    cmake -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build -j

Requires CMake >= 3.20 (< 4.0), a C++17 compiler, OpenMP. Bundled/fetched:
yaml-cpp, spglib, Catch2. The MPI driver builds automatically when FindMPI
succeeds (if several MPI stacks are installed, point CMake at one, e.g.
`-DMPI_CXX_COMPILER=/usr/bin/mpicxx.openmpi`). Full dependency list, tuning
guidance and site instructions: `docs/INSTALL.md`.

## Quickstart

    ./build/exsqs configs/w70cr30_4x4x4.yaml --out runs/demo

Key outputs in `--out`: `best_XX.vasp` (POSCAR structures), `summary.json`
(config echo, per-structure records, generation log, migration ledger),
`checkpoint.json` + `state.ckpt` (resumable state). Exit code 0 = converged
to `e_tol`, 3 = budget exhausted (resumable). Continue a run, optionally
raising budgets:

    ./build/exsqs configs/w70cr30_4x4x4.yaml --resume runs/demo \
        --set evolution.max_generations=200

Reference configs: `configs/w70cr30_4x4x4.yaml` (binary bcc, 128 sites),
`configs/w50mo25cr25_4x4x4.yaml` (ternary, exact-zero floor),
`configs/smoke_sc27.yaml` (seconds-fast smoke).

## Scoring external structures

Compare *raw structures* (paper supplements, ATAT output), never reported
scalars — absolute error scales are not comparable across codes (SPEC 4.1):

    ./build/exsqs score configs/w70cr30_4x4x4.yaml some.vasp [--json out.json]

The scorer is orientation-strict by design. For rotated / reduced / permuted
files, align them onto the config frame first:

    ./build/exsqs geom configs/w70cr30_4x4x4.yaml -o geom.vasp
    python3 tools/py/align_to_config.py geom.vasp their_file.cif aligned.vasp
    ./build/exsqs score configs/w70cr30_4x4x4.yaml aligned.vasp

`tools/py/to_poscar.py` converts anything pymatgen reads to POSCAR.

## HPC

MPI driver (islands sharded over ranks, rank-count invariant):

    mpirun -n 8 ./build/exsqs_mpi configs/w70cr30_4x4x4.yaml --out runs/mpi

SLURM templates in `__scripts/slurm/` (`exsqs_omp.sbatch`, `exsqs_mpi.sbatch`);
`__scripts/chain_resume.sh` chains jobs on the exit-0/3 contract until
convergence. Site-ready script sets for LUMI-C and Cyfronet's Helios —
compile driver plus three run flows each — live in `__scripts/hpc/lumi/` and
`__scripts/hpc/helios/` (see `docs/INSTALL.md`).

## Testing

    bash tools/run_all_tests.sh

runs the full verification matrix (24 gates, ~2 min on a multicore box): unit
suites (27k+ assertions), the phonopy displacement gate, end-to-end
integration, MPI rank invariance, python cross-validation of every correlation
function (`tools/py/validate.py`, |diff| = 0 gates), the T-X2 align/score
round trip, recorded benchmarks — headed by the coherence audit
(`python3 tools/check_coherence.py`), which cross-checks versions, test-id
citations, tag coverage, the trajectory-signature field ledger, executable
spec samples, and README references.

## Repository layout

`src/`, `include/exsqs/` — engine; `tests/` — Catch2 suites; `tools/` —
runner, coherence audit, benchmarks, MPI/alignment gates; `tools/py/` —
validation and interop scripts; `configs/` — reference configs (including the
16-cell K=5 supercell ladder); `_scripts/` — SLURM + chaining; `docs/` — SPEC.md,
TEST_MATRIX.txt, per-step reports, DEV_NOTES.md, SUPERCELL_STUDY.md.

## Citing

See `CITATION.cff` (cite the arXiv paper and this software).

## License

MIT — see `LICENSE`.
