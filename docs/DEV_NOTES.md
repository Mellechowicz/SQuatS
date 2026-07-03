# EXSQS development notes (accreted per-release README sections)

Moved here verbatim at the 1.7.0 release. The authoritative history is
`docs/SPEC.md` section 16 plus the dated `docs/STEP*_REPORT.md` snapshots.

## v1.1 notes (Step 2)

`e_tol: auto` (default) resolves to 3.0 x E_floor, the provable L1 quantization lower bound
(SPEC 4.1) that the engine logs at startup; the paper's absolute error scale lies below this
bound for its own reference system, so compare on raw structures (tools/py/validate.py), never
on reported scalars. Serial exploitation recipe (seed-reproducible):
`--set parallel.islands=4 --set evolution.survival.mode=metropolis --set evolution.survival.beta=3000`.
New in v1.1: constructive space-group seeding [D4] behind the P1 rejection filter,
cyclic-subgroup [D6] moves, mutation `swaps: poisson` (+ `lambda`), `stagnation_stop`.

## v1.2 notes (Step 3)

Islands now advance in lockstep with synchronous ring migration (`migration_every`, `migrants`);
generations execute as parallel rounds under OpenMP (`parallel.omp_threads`). Thread counts never
change results: runs are bit-identical for 1 vs N threads and match 1.1.0 sequential artifacts
exactly. Record single-node scaling with `./build/exsqs_bench_scaling`.

## v1.3 notes (Step 4)

Checkpoint/restart: every run writes `state.ckpt` (at `checkpoint_every` rounds and at exit);
`--resume DIR` continues bit-exactly — only budget caps may change (signature-guarded), so
wall-limited segments chain into one deterministic trajectory (`scripts/chain_resume.sh`,
SLURM templates in `scripts/slurm/`). Multi-node: `mpirun -n R ./build/exsqs_mpi ...`
distributes islands over ranks with results bit-identical to the serial driver for any R;
serial and MPI segments are interchangeable. Wire/state format assumes little-endian hosts.

## v1.4 notes (Step 5) and testing

K >= 3 is validated end-to-end (`configs/w50mo25cr25_4x4x4.yaml`; `error.mode: auto` resolves to
full_pairs [A16]); note that commensurate compositions can have E_floor = 0 exactly -- use a
numeric `e_tol` there. Run the complete verification matrix with:

    bash tools/run_all_tests.sh          # everything
    SKIP_E2E=1 bash tools/run_all_tests.sh   # fast gates only

Individual suites: `./build/tests/exsqs_tests "~[e2e]"` (fast), `"[e2e]"` (integration),
`tools/test_mpi_invariance.sh` (T-MPI1), `tools/py/validate.py <rundir>` (T-V1),
`./build/exsqs_bench` and `./build/exsqs_bench_scaling` (T-B1, recorded). CI stub in
`.github/workflows/ci.yml`.

## v1.5 notes (Step 6)

Score any external structure (paper supplementary files, ATAT output, hand-built cells) under a
config -- the raw-structure comparison the SPEC 4.1 scale caveat demands:

    ./build/exsqs score configs/w70cr30_4x4x4.yaml runs/*/best_00.vasp [--json scores.json]

Input site order is free; lattice, species and composition must match the config (convert CIF
etc. with `tools/py/to_poscar.py`). Also new: the [A11] geometric beta schedule
(`survival: {mode: metropolis, beta: B, schedule: geometric, beta_growth: g}`); enabling it
invalidates earlier state files under `--resume` (signature-guarded, by design).

## v1.6 notes (Step 7)

`python3 tools/check_coherence.py` audits spec<->code<->tests<->configs<->docs coherence
(versions, section-12 test-id citations, tag coverage, the trajectory-signature field ledger,
executable section-11 sample and reference configs, README references). It runs as the S0 row of
`tools/run_all_tests.sh`.
