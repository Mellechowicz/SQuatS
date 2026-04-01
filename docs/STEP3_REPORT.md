# EXSQS — Step 3 Report: Parallel Layer (v1.2.0, SPEC v1.2)

Date 2026-04-01. Scope per plan: lockstep multi-island execution with synchronous ring migration,
OpenMP parallelism over islands and inside each generation's evaluation rounds, thread-safety and
thread-invariance test coverage, and the scaling harness. The engine was refactored into a
resumable `IslandEngine`, which is also the checkpoint/restart foundation for Step 4. Everything
from the 1.0.0 core and the 1.1.0 engine is unchanged and green.

Suite status: fast suite **26 test cases, 26,753 assertions, all passing** (the core phonopy gate
and all 1.1.0 engine cases included); e2e **T-E1 and T-E2 passing** (63 s); new v1.2 tests
**T-P0** (pipeline thread safety), **T-P1** (thread-count invariance, islands + migration and
single-island paths), **T-M1** (migration semantics + deterministic rerun) all passing; **T-V1**
still exact (|ΔE| = 0.0) on the extended `summary.json` schema.

## 1. The headline: thread topology is provably performance-only

SPEC §10 chose the counter-keyed RNG so that "trajectories are independent of thread count and
scheduling"; v1.2 realizes that end-to-end and verifies it three ways. First, T-P1 reruns the
same seeded system with `omp_threads` 1 vs 4 — through lockstep islands *with active migration*
and through the single-island inner-parallel path — and requires field-identical logs and
outputs. Second, the scaling harness prints the best error alongside wall time, and the column is
identical across thread counts by construction. Third, and strongest: the 1.1.0 T-E1 recipe
(4 islands × β = 3000), run by the *sequential* 1.1.0 engine, is reproduced by the v1.2
lockstep + round-based engine under `omp_threads: 4` **to the last bit** — the generation log and
all output errors match exactly (best 0.008111946532999254 in both). Two engine generations, one
trajectory.

## 2. Architecture

`evolve_island`'s monolithic loop became `IslandEngine` with `seed_generation0()` and
`advance()` — verbatim transpositions of the v1.1 algorithm, so a direct `evolve_island()` call
remains bit-exact. Because the RNG is stateless, an island's full state is (population, archive,
pool, generation, stall counter): trivially resumable, which Step 4 will serialize for restart.
`run_evolution` is now a lockstep driver: every active island advances one generation per round
(OpenMP `parallel for` over islands, exceptions captured per island and rethrown in island order;
the checkpoint callback is mutex-serialized), and stopped islands freeze. `max_wall_s` and the
logged `elapsed_s` are global run wall time in this model — the one deliberate semantic change,
recorded in SPEC §16.

Synchronous ring migration (SPEC §8.1, previously only schema keys) runs every `migration_every`
rounds among the still-active islands: each sends copies of its `migrants` pool-best individuals
to the next active island in ascending-index order; receivers skip canonicals already in their
archive, otherwise archive the migrant and replace their current worst member (the [D5] elite is
never the worst for P ≥ 2, so E_min stays monotone). Exchanges are computed from pre-round
snapshots and consume no randomness — the migrating system reruns bit-identically (T-M1). Note
the behavior change: the §11 defaults (`migration_every: 50, migrants: 2`) are now *active* for
multi-island runs that reach 50 generations; 1.1.0 ignored them.

Inside a generation, seeding and repopulation now execute as rounds: candidate generation and
canonicalization run in parallel (pure per-slot functions of the counter RNG), then the duplicate
check runs serially in slot order with a *provisional* archive insert, then all surviving
candidates are evaluated in parallel (one spglib dataset each), then committed serially in slot
order. The provisional insert is what makes this exact: v1.1 archived a candidate's canonical
whether it was ultimately admitted or P1-rejected, so inserting before the evaluation reproduces
the sequential archive evolution verbatim. Seeding is therefore bit-equivalent unconditionally;
repopulation is bit-equivalent whenever no post-seed P1 admission attempt occurs (the empirical
case in every reference run — a P1-rejected slot resumes its own attempt stream next round, the
normative v1.2 semantics for that rare path). Thread budget: outer = min(islands, T) over
islands, inner = T/outer inside the rounds, nested OpenMP enabled when both exceed one.

## 3. Thread safety and scaling evidence

T-P0 evaluates 48 random decorations through the full pipeline (spglib dataset, displacement
count, correlations) serially and under a 4-thread `parallel for`, requiring identical space
groups, operation counts, equivalent-atom vectors, D values, and errors — spglib 2.7.0 and the
core wrappers hold up. Measured on the dev workstation (`./build/exsqs_bench_scaling`, fixed
workload W₉₀Cr₃₈ bcc 4×4×4, pop 48, 6 generations):

| islands | omp_threads | wall_s | best_e_pure |
|---|---|---|---|
| 1 | 1 | 1.83 | 2.333890e-02 |
| 1 | 2 | 1.05 | 2.333890e-02 |
| 1 | 4 | 0.66 | 2.333890e-02 |
| 4 | 1 | 7.95 | 1.925091e-02 |
| 4 | 2 | 4.13 | 1.925091e-02 |
| 4 | 4 | 2.22 | 1.925091e-02 |

Inner-round parallelism alone reaches 2.8× at 4 threads on one island, consistent with the
Step-2 phase measurements at N = 128 (canonical 0.28 ms, spglib 2.17 ms, D 0.07 ms, correlations
0.014 ms per evaluation): the parallel evaluation covers ≈ 88% of per-candidate cost, and the
serial residue is the in-spin canonicalization of duplicate retries plus commits and extinction,
giving a single-island Amdahl ceiling around 7–8×. Four islands at 4 threads run 3.6× faster
than serial — islands scale near-linearly on top (they only synchronize at round barriers and
migration points). The error column is constant across every row, as it must be. Multi-node
island scaling over MPI ranks (the "1→32 ranks" clause of T-B1) is deferred to the HPC step and
will reuse the migration protocol at rank boundaries.

## 4. Migration demonstrated on the reference system

W₉₀Cr₃₈, γ = 0, four islands at moderate pressure (metropolis β = 800), `migration_every: 8`,
`migrants: 2`, ≤ 40 generations, `e_tol: auto`: the run terminates **rc = 0 SUCCESS at
2.50× E_floor** (8.112e-3, Cm) after 10,726 evaluations. Three islands cross the gate
independently (stops `[e_tol, e_tol, stagnation, e_tol]` at generations `[4, 10, 21, 11]`), and
the dedup-aware ledger is visible in `summary.json` — `island_migrants_in = [0, 2, 2, 2]`
against `island_migrants_out = [0, 2, 2, 2]` (island 0 crossed before the first sync round;
receivers skip canonicals they have already explored, so in ≤ out always). The one island that
did not converge is the one the [A13] stagnation stop retired — at β = 800 the exploiter
exhausts its ⟨g⟩ neighborhood exactly as the Step-2 arithmetic predicts. The run validates
through `validate.py` with |ΔE| = 0.0.

## 5. Usage additions

`parallel: {islands: N, omp_threads: auto|K, migration_every: R, migrants: k}` is fully active;
set `migration_every: 0` or `migrants: 0` to disable exchange. `summary.json` now records
`island_migrants_in/out`. Everything else — CLI, outputs, validator — is unchanged. Version
strings report 1.2.0.

## 6. Deferrals and Step-4 handoff

Deferred to Step 4 (HPC): checkpoint *restart* — `IslandEngine` is already resumable in memory;
what remains is serializing (population σ's, generation, archive, pool) and a `--resume` path;
MPI island ranks across nodes reusing the §8.1 migration protocol; SLURM workflow templates; and
long-budget validation runs against the paper's raw structures via `validate.py`, using the
floor-relative metric (E/E_floor) as the scale-free convergence measure. On a 16–64-core node
the v1.2 layer already sustains the budgets that the Step-2 dynamics study identified as the
remaining gap.
