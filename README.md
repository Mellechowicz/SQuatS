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

Multiplet sectors (v1.9, SPEC 4.2): `error.multiplets: {lambda3: 1.0}`
adds cancellation-free triplet (and with `lambda4` quadruplet) L1 sectors
over symmetry-invariant cluster classes, `E = E_2 + l3*E_3 + l4*E_4`;
`configs/mul_smoke_sc27.yaml` is the smoke,
`configs/hea5_bcc_5x5x5_mul.yaml` the production example. Defaults (both
lambdas 0) keep the pair-only engine bit-identical.

## Tags

Every key the parser accepts, with the value used when the key is absent.
Any of them is CLI-overridable: `--set evolution.population=64`. Keys the
parser does not know (typos included) are ignored without a message.

### System — required

| tag | example | default | notes |
|---|---|---|---|
| `lattice.type` | `bcc` | — | `sc`, `bcc`, `fcc`, `hcp` |
| `lattice.a` | `3.165` | — | lattice constant (A), must be > 0 |
| `lattice.c` | `5.21` | ideal `a*sqrt(8/3)` | `hcp` only |
| `lattice.file` | `proto.vasp` | — | VASP5 POSCAR prototype; its geometry is used, species are re-decorated |
| `supercell.diag` | `[4, 4, 4]` | — | diagonal H, three integers |
| `supercell.matrix` | `[[2,1,0],[0,2,0],[0,0,2]]` | — | full 3x3 integer H, `det H > 0` |
| `composition.<El>` | `W: 0.70` | — | >= 2 species, fractions summing to 1 (1e-6); counts by largest remainder [A5] |

### Zones and error

| tag | example | default | notes |
|---|---|---|---|
| `zones.n_shells` | `5` | `7` | coordination shells entering `E_pure` |
| `zones.shell_tol` | `1.0e-3` | `1.0e-3` | relative distance tolerance merging pairs into a shell |
| `error.weights.form` | `inv_n` | `inv_n` | `inv_n`, `inv_n_pow`, `inv_R`, `custom` |
| `error.weights.p` | `2.0` | — | required by `inv_n_pow` |
| `error.weights.values` | `[1, 0.5, 0.33, 0.25, 0.2]` | — | required by `custom`; length must equal `zones.n_shells` |
| `error.mode` | `auto` | `auto` | `auto` = `diagonal` for K = 2, `full_pairs` for K >= 3 [A16] |
| `error.gamma` | `1.0` | `1.0` | `E_obj = E_pure * D^gamma`; `0` = classical SQS |
| `error.multiplets.lambda3` | `1.0` | `0.0` | triplet sector weight (>= 0); `0` keeps the pair-only engine bit-identical |
| `error.multiplets.lambda4` | `0.5` | `0.0` | quadruplet sector weight (>= 0) |
| `error.multiplets.shell3` | `2` | `2` | triplet cutoff in shells, `1..zones.n_shells` |
| `error.multiplets.shell4` | `1` | `1` | quadruplet cutoff in shells, `1..zones.n_shells` |

### Evolution

| tag | example | default | notes |
|---|---|---|---|
| `evolution.population` | `100` | `200` | >= 2 |
| `evolution.outputs` | `6` | `10` | M best inequivalent structures written |
| `evolution.e_tol` | `auto` | `auto` | `auto` = `3.0 * E_floor` (SPEC 4.1), or a positive number |
| `evolution.max_generations` | `300` | `5000` | budget cap; exit code 3 when hit |
| `evolution.max_wall_s` | `3600` | unlimited | wall cap; stops are clock-dependent, so not bit-comparable |
| `evolution.stagnation_stop` | `25` | `3` | consecutive all-fallback generations before stopping; `0` disables [A13] |
| `evolution.survival.mode` | `ratio` | `ratio` | `ratio` or `metropolis` |
| `evolution.survival.beta` | `3000` | `auto` | `metropolis` only; `auto` = `ln2 / median(E_i - E_min)` per generation |
| `evolution.survival.schedule` | `geometric` | `const` | `const` or `geometric` (`beta_g = beta * growth^g`) [A11] |
| `evolution.survival.beta_growth` | `1.02` | `1.0` | `geometric` only, > 0 |
| `evolution.elitism_best` | `1` | `1` | protected best members; `1 <= x < population` [D5] |
| `evolution.p1_elite_quota` | `0` | `0` | P1 structures tolerated in the population [A8] |
| `evolution.mutation.swaps` | `1` | `1` | integer >= 1, or `poisson` |
| `evolution.mutation.lambda` | `1.0` | `1.0` | `swaps: poisson` only: `k ~ 1 + Poisson(lambda)` [A12] |
| `evolution.mutation.symmetry_preserving` | `true` | `true` | cyclic-subgroup orbit swaps [D6] |
| `evolution.mutation.retry_budget` | `100` | `100` | attempts per slot before constructive fallback |
| `evolution.seeding.mode` | `mixed` | `mixed` | `rejection`, `constructive`, `mixed` [D4] |
| `evolution.seeding.sg_pool` | `auto` | `auto` | only `auto` is defined |

### Symmetry, displacements, parallel, rng, output

| tag | example | default | notes |
|---|---|---|---|
| `symmetry.symprec` | `1.0e-5` | `1.0e-5` | spglib tolerance [A6] |
| `symmetry.filter.policy` | `reject_p1` | `reject_p1` | only `reject_p1` is defined [A7] |
| `displacements.convention` | `phonopy_default` | `phonopy_default` | only value defined [A15] |
| `parallel.islands` | `4` | `1` | >= 1; trajectories are island-keyed |
| `parallel.omp_threads` | `auto` | `auto` | `auto` or >= 1; performance only, never changes results [A14] |
| `parallel.migration_every` | `8` | `50` | rounds between ring migrations; `0` disables |
| `parallel.migrants` | `2` | `2` | pool-best copies sent per island; `0` disables |
| `rng.seed` | `20260712` | `42` | one seed + config = one trajectory |
| `output.dir` | `./runs/demo` | `./run1` | `--out DIR` is shorthand for this |
| `output.formats` | `[poscar]` | `[poscar]` | other formats warn and are skipped |
| `output.checkpoint_every` | `20` | `100` | rounds between `state.ckpt` and mid-run `checkpoint.json` writes |
| `output.log_level` | `info` | `info` | `info`/`debug` print per-generation rows; anything else is quiet |

### Conflicts and precedence

| combination | outcome |
|---|---|
| `lattice.file` + `lattice.type`/`a` | file wins, the others are unused |
| `supercell.diag` + `supercell.matrix` | `diag` wins, `matrix` is unused |
| `survival.schedule: geometric` + `mode: ratio` | error: requires mode metropolis |
| `survival.schedule: geometric` + `beta: auto` | error: requires a numeric beta |
| `survival.schedule: geometric` + `beta_growth <= 0` | error: beta_growth must be > 0 |
| `survival.beta` + `mode: ratio` | ignored silently (beta is read only under metropolis) |
| `weights.form: inv_n_pow` without `p` / `custom` without `values` | error at load |
| `weights.values` length != `zones.n_shells` | error at context build: custom weight size mismatch |
| `multiplets.shell3`/`shell4` > `zones.n_shells` | error: must be in [1, zones.n_shells] |
| `multiplets.lambda3`/`lambda4` > 0 | `E_floor` gains the multiplet sectors, so `e_tol: auto` rises with them |
| `elitism_best` >= `population` | error: elitism_best out of range |
| composition deviating more than 2% from the target (fixed `comp_tol`) | error; pin representable fractions on small cells |
| composition where a species rounds to 0 atoms | error: species rounds to zero atoms |
| `e_tol: auto` on a commensurate composition (`E_floor = 0`) | target degenerates to 0 and is unreachable — set a numeric `e_tol` |
| numeric `e_tol` below `E_floor` | warning; the run can only end on caps |
| `islands: 1`, or `migrants: 0`, or `migration_every: 0` | migration inactive (ledger stays zero) |
| `mutation.lambda` without `swaps: poisson` | ignored silently (`lambda` is read only in the poisson branch) |
| `lattice.c` on non-`hcp` types | ignored |
| `zones.n_shells` reaching past the cell on small supercells | error: `build_zones: sites do not form a single orbit`; lower `n_shells` (8-site cells take 3) |

### Changing tags on `--resume`

Only `evolution.max_generations`, `evolution.max_wall_s`, `parallel.omp_threads`
and the `output.*` keys may differ from the checkpointed run. Every other tag
enters the trajectory signature — including `parallel.islands`, `rng.seed`,
`error.*` and the whole system block — and a change is refused with
`state: trajectory signature mismatch`. The `multiplets.shell3`/`shell4` keys
are signed only while a lambda is non-zero, so they are free to change on a
pair-only run (where they have no effect anyway).

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

SLURM templates in `scripts/slurm/` (`exsqs_omp.sbatch`, `exsqs_mpi.sbatch`);
`scripts/chain_resume.sh` chains jobs on the exit-0/3 contract until
convergence. Site-ready script sets for LUMI-C and Cyfronet's Helios —
compile driver plus three run flows each — live in `scripts/hpc/lumi/` and
`scripts/hpc/helios/` (see `docs/INSTALL.md`).

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
16-cell K=5 supercell ladder); `scripts/` — SLURM + chaining; `docs/` — SPEC.md,
TEST_MATRIX.txt, per-step reports, DEV_NOTES.md, SUPERCELL_STUDY.md.

## Citing

Please cite the Acta Physica Polonica B paper (DOI
[10.5506/APhysPolB.57.5-A15](https://doi.org/10.5506/APhysPolB.57.5-A15))
and this software (github.com/Mellechowicz/SQuatS). Ready-made BibTeX
entries live in `CITATION.bib`; citation metadata in `CITATION.cff`
(GitHub's "Cite this repository" button uses it).

## License

MIT — see `LICENSE`.

## Disclaimer
For the coding and version control, we used assist by [Claude Code](https://claude.com/claude-code).
