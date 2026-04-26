# exsqs — extinction-based evolutionary SQS generator (minimal)

Minimal C++17 implementation of the extinction evolutionary algorithm for
Special Quasirandom Structures (arXiv:2602.10872). The engine minimizes
E_obj = E_pure * D^gamma over decorations of a supercell: correlation error
against the ideal random alloy, traded against the number D of
symmetry-inequivalent displacements a phonon workflow must compute.

Build:  cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
Run:    ./build/exsqs configs/w70cr30_4x4x4.yaml --out runs/demo
Exit codes: 0 converged to e_tol, 3 budget exhausted, 1 error.

Outputs: best_XX.vasp (POSCAR) and summary.json in --out.

Testing: `./build/tests/exsqs_tests` runs the Catch2 suites (geometry, correlations,
symmetry/dedup, the T-D1 phonopy displacement gate, RNG, config); cross-validate any
run directory independently with `python3 tools/py/validate.py <dir>` (T-V1).

## v1.1 notes

`e_tol: auto` (default) resolves to 3.0 x E_floor, the provable L1 quantization lower bound
(SPEC 4.1) that the engine logs at startup; the paper's absolute error scale lies below this
bound for its own reference system, so compare on raw structures (tools/py/validate.py), never
on reported scalars. Serial exploitation recipe (seed-reproducible):
`--set parallel.islands=4 --set evolution.survival.mode=metropolis --set evolution.survival.beta=3000`.
New in v1.1: constructive space-group seeding [D4] behind the P1 rejection filter,
cyclic-subgroup [D6] moves, mutation `swaps: poisson` (+ `lambda`), `stagnation_stop`.

## v1.2 notes

Islands now advance in lockstep with synchronous ring migration (`migration_every`, `migrants`);
generations execute as parallel rounds under OpenMP (`parallel.omp_threads`). Thread counts never
change results: runs are bit-identical for 1 vs N threads and match 1.1.0 sequential artifacts
exactly. Record single-node scaling with `./build/exsqs_bench_scaling`.

## v1.3 notes

Checkpoint/restart: every run writes `state.ckpt` (at `checkpoint_every` rounds and at exit);
`--resume DIR` continues bit-exactly - only budget caps may change (signature-guarded), so
wall-limited segments chain into one deterministic trajectory (`scripts/chain_resume.sh`).
Multi-node: `mpirun -n R ./build/exsqs_mpi ...` distributes islands over ranks with results
bit-identical to the serial driver for any R; serial and MPI segments are interchangeable.
Wire/state format assumes little-endian hosts.

License: MIT (see LICENSE).
