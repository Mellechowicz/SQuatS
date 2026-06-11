# EXSQS — Step 6 Report: Interoperability + Spec Completion (v1.5.0, SPEC v1.5)

Date 2026-06-11. The original plan ended at Step 5, so Step 6 was defined from the project's own
deferral ledger and one long-standing tooling gap: SPEC §4.1 established that absolute error
scalars are not comparable across codes and that comparisons must happen on raw structures —
but until now only this engine's own runs could be evaluated. Step 6 delivers that tool,
completes the one genuinely deferred spec feature ([A11]'s geometric β schedule), closes two
coverage gaps, and re-runs the full verification matrix.

Verification: **all matrix gates green**; fast suite **37 test cases / 27,111 assertions**; e2e
T-E1/T-E2/T-E3 pass; T-MPI1 and all T-V1 instances hold at |diff| = 0.

## 1. `exsqs score` — raw-structure comparison, operational

The new subcommand evaluates any POSCAR on a config's geometry with the config's zones, weights,
error mode and γ: `exsqs score <config> [--set ...] <files...> [--json PATH]` prints E_pure,
E/E_floor, D, D(P1)/D, SG and E_obj per structure. Input site *order* is free (sites are matched
by wrapped fractional coordinate); lattice orientation, species names and the [A5] composition
must match the config and fail loudly otherwise. The API (`score_structure` in `score.hpp`) is
gated by **T-X1**: scoring the engine's own outputs reproduces every recorded value bit-exactly,
including through a full POSCAR write/read round trip, and including under site permutation.

The demonstration is the comparison table this project was built to enable — one command, three
structures on the reference system:

| structure | E_pure | E/E_floor | D | D(P1)/D | SG |
|---|---|---|---|---|---|
| random decoration (baseline) | 5.105e-2 | 15.73× | 768 | 1.00 | P1 |
| campaign champion | 8.112e-3 | **2.50×** | 416 | 1.85 | Cm |
| γ = 1 knee | 2.788e-2 | 8.59× | 128 | **6.00** | P4mm |

To score the paper's supplementary structures or ATAT output, convert anything pymatgen reads
with the new `tools/py/to_poscar.py` (which also closes the CIF-prototype deferral
pragmatically) and point `exsqs score` at the result — the same cell and orientation are
required, by design.

## 2. [A11] geometric β schedule — the last deferred spec feature

`survival: {mode: metropolis, beta: β₀, schedule: geometric, beta_growth: g}` now applies
β_g = β₀·gᵍᵉⁿ; `beta: auto` recomputes per generation and *rejects* the schedule, per the
project's loud-fail config policy (a first draft parsed the keys only inside the metropolis
branch, which would have silently ignored `ratio` + `geometric` — the test caught it and the
validation was hoisted). The schedule tests pin two properties: growth = 1.0 is **bit-identical**
to the const schedule, and growth > 1 diverges deterministically. The two fields have been part
of the trajectory signature since v1.3, so enabling the schedule invalidates earlier state files
under `--resume` — the dynamics vocabulary changed, and the signature exists precisely to notice
that.

## 3. Coverage closure

Non-diagonal `supercell.matrix` turned out to be implemented since the 1.0.0 schema but never
exercised; it now has a test (sheared 8-site cell: |det H| counts, zone/context build, sane
empty-cell symmetry). The runner gained the three S6 rows and reports v1.5.

## 4. Ledger status after Step 6

On the software side the deferral ledger is now effectively empty: `seeding.sg_pool` variants
and non-`reject_p1` filter policies were never *specified* beyond their names, so the
implementation is spec-complete as written (both remain loud-fail placeholders); the remaining
open items — heterogeneous clusters, wall-cap topology comparability, mid-run JSON checkpoints
under MPI — are documented architectural choices rather than gaps. What remains of the overall
program is unchanged from Step 5: the multi-node, long-budget campaigns, for which every needed
piece (bit-exact chains, rank-invariant MPI, SLURM templates, floor-relative convergence metric,
and now external-structure scoring for the final paper comparison) is in place.
