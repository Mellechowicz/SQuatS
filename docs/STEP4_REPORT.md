# EXSQS — Step 4 Report: HPC Layer — Checkpoint/Restart + MPI + SLURM (v1.3.0, SPEC v1.3)

Date 2026-04-30. Scope per plan: bit-exact checkpoint/restart with a `--resume` CLI and
resubmission-chain pattern, the rank-distributed `exsqs_mpi` driver reusing the §8.1 migration
protocol, SLURM workflow templates, and validation demos on the reference system. Everything from
the previous milestones is unchanged and green.

Suite status: fast suite **29 test cases, 27,043 assertions, all passing** (all earlier gates
included); e2e **T-E1 and T-E2 passing** — the reference trajectories survive the v1.3 driver
unchanged; new **T-K1** (resume bit-exactness), **T-K2** (signature guard), **T-MPI1**
(rank invariance, ctest-registered) all passing.

## 1. The headline: one trajectory, any execution topology, any interruption

v1.3 closes the equivalence class the project has been building since [A14]: for a given
(seed, config), the following all produce **bit-identical** outputs, logs, and migration
ledgers — the 1.1.0 sequential engine, the lockstep OpenMP engine at any thread count, the MPI
driver at any rank count, and any chain of budget-limited segments stitched with `--resume`.
Concretely verified: T-K1 splits a 12-generation two-island sc27 run at generation 5 (a
migration round falls at 6, *after* the split) and matches the uninterrupted run field-for-field;
the W₉₀Cr₃₈ CLI demo does the same at 30 generations split 15+15 with migration every 10 rounds
(ledger `[2, 2]` in both, best 1.2701e-2 = 3.91× floor); T-MPI1 matches serial vs `mpirun -n 1`
vs `-n 3` including full per-generation logs; and the W-system MPI demo reproduces the serial
best (1.820570e-2) exactly from 3 ranks. Wall-capped stops remain the one deliberate exclusion —
they depend on the clock, as documented since v1.2.

## 2. Checkpoint/restart design

A new little-endian binary encoding (`serialize.hpp`, format v1; doubles as IEEE-754 bit
patterns, so round trips are exact) serves both the state file and the MPI wire. The engine was
already resumable in memory because the counter-keyed RNG is stateless; v1.3 serializes the full
per-island state — population, dedup archive (sorted for deterministic files), pool, log,
counters, stop status — behind the new public `EngineHandle` facade. The driver writes
`state.ckpt` (atomic tmp + rename) into the output directory at every `checkpoint_every` round
and at run end, so a budget-exhausted run (exit code 3) always leaves a resumable state.

`--resume DIR` restores everything and re-arms islands stopped on *raisable* caps
(`max_generations`, `wall_time`); `e_tol` success and `stagnation` stops are terminal. A
*trajectory signature* embedded in the file — every config field that influences the trajectory:
geometry, composition, zones, error model, evolution dynamics, seed, islands, migration —
guards the resume: T-K2 shows a changed seed is refused with a clear message while raising
`max_generations` continues the run. Budget caps, output settings and thread counts are
deliberately outside the signature, which is exactly the resubmission-chain contract.

## 3. MPI driver

`exsqs_mpi` distributes islands round-robin (island *i* → rank *i* mod *R*) and runs the same
lockstep rounds as the serial driver: per-round activity flags travel by `MPI_Allreduce`, §8.1
migration payloads as `(island id, blob)` records over `MPI_Allgatherv` — every rank sees every
active island's emigrants and applies only its owned receivers, so the ring is identical to the
single-process one by construction. Rank 0 gathers per-island blobs to write the *same*
`state.ckpt` format, meaning serial and MPI segments are interchangeable: start on a laptop,
resume on 8 nodes, finish serially — one trajectory throughout. OpenMP nests inside each rank
(outer over the rank's islands, remainder inner), and a homogeneous little-endian cluster is
assumed and guarded at every file/wire boundary. Under MPI, per-generation stdout rows come from
rank 0's islands and the mid-run `checkpoint.json` snapshots are superseded by `state.ckpt`.

Environment note: local verification uses a user-local OpenMPI 4.1.6; the
`--allow-run-as-root --oversubscribe` flags in `tools/test_mpi_invariance.sh` are
container/CI conveniences and harmless elsewhere. On clusters, `srun` with the site MPI is the
intended launcher.

## 4. SLURM workflow

`scripts/slurm/exsqs_omp.sbatch` (single node, threads from `SLURM_CPUS_PER_TASK`) and
`scripts/slurm/exsqs_mpi.sbatch` (one rank per node, OpenMP inside, `srun`) both auto-resume when
`$RUNDIR/state.ckpt` exists, so `sbatch` requeues and dependency chains work out of the box; the
MPI template includes a commented self-resubmission block keyed to the exit-code contract
(0 = `e_tol` reached, 3 = budget exhausted with state kept). `scripts/chain_resume.sh` is the
login-node/local variant: it loops run → resume until success or `MAX_SEGMENTS`, and because of
T-K1 the chain is provably the same trajectory as one long run. The floor-relative metric
(E/E_floor, logged at startup and in `summary.json`) is the scale-free convergence measure for
those long campaigns.

## 5. Evidence summary

| check | system | result |
|---|---|---|
| T-K1 (test) | sc27, 2 islands, migration every 3, split 5/12 | bit-exact vs straight run |
| T-K2 (test) | sc27 | seed change refused ("signature"); a cap raise resumes the run |
| T-MPI1 (ctest) | sc27, 3 islands | serial ≡ 1 rank ≡ 3 ranks (outputs, logs, ledger) |
| chain demo | W₉₀Cr₃₈, 2 islands, β=3000, 30 gens split 15+15 | bit-exact incl. ledger `[2,2]`; best 3.91× floor |
| MPI demo | W₉₀Cr₃₈, 3 islands, 6 gens | serial ≡ 3 ranks (best 1.820570e-2) |
| e2e regression | T-E1/T-E2 | pass — reference trajectories unchanged through the v1.3 driver |

## 6. Deferrals and Step-5 handoff

Deliberately out of v1.3 scope: heterogeneous or big-endian clusters (guarded, not supported),
mid-run JSON checkpoints under MPI, and topology-comparability of wall-capped stops. What remains
of the overall plan: the long-budget validation campaign — multi-node, multi-hour chained runs on
W₇₀Cr₃₀ and larger/K≥3 systems, tracking E/E_floor and D against the paper's claims, with
`validate.py` as the scale-proof cross-check — plus release polish (CI, packaging, docs) to make
the repository public-ready. The machinery for all of it now exists; step 5 will validate the
multicomponent path and consolidate the verification matrix.
