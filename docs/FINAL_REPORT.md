# EXSQS 1.7.0 — Final Report (Releases 1.1.0 – 1.7.0)

2026-02-16 → 2026-07-03. This is the consolidated record of the full-program development on top
of the minimal 1.0.0 release: a C++17 implementation of the extinction-based evolutionary SQS
algorithm from *"On generating Special Quasirandom Structures: Optimization for the DFT
computational efficiency"* (arXiv:2602.10872), built spec-first, verified to the bit, released
as 1.7.0 with a 24-gate machine-checked verification matrix. Every claim below is reproducible
from the repository with the commands in §5.

The one-table summary of what the algorithm buys, produced by one `exsqs score` call on the
reference system (bcc W₉₀Cr₃₈, 4×4×4, 128 sites, E_floor = 3.246078e-3):

| structure | E_pure | E/E_floor | D | D(P1)/D | SG |
|---|---|---|---|---|---|
| random decoration | 5.105e-2 | 15.73× | 768 | 1.00 | P1 |
| γ = 0 champion | 8.112e-3 | **2.50×** | 416 | 1.85 | Cm |
| γ = 1 knee | 2.788e-2 | 8.59× | 128 | **6.00** | P4mm |

Read: against a random special-quasirandom start, the engine cuts correlation error six-fold at
γ = 0; at γ = 1 it instead buys a **six-fold reduction in the symmetry-inequivalent
displacements** a downstream phonon workflow must compute, at a controlled error cost — the
paper's central trade, demonstrated end to end.

## 1. What was built, release by release

The 1.0.0 core (released separately, February) delivered geometry, coordination zones, pair
correlations, the spglib symmetry wrapper, canonical deduplication, the counter-based RNG [A14]
and a minimal sequential engine; `feat/tests-core` then gated it bit-exactly against phonopy's
displacement convention (T-D1). **1.1.0** completed the engine per SPEC §8 — constructive
seeding behind the P1 filter, [D6] mutation, stagnation stops — and produced the project's
first scientific result, the quantization floor doctrine (§2). **1.2.0** made islands lockstep
and OpenMP-parallel with **bitwise** thread-count invariance. **1.3.0** added the HPC plumbing:
signature-guarded checkpoint/resume with raisable budgets, an MPI driver with rank-count
invariance (T-MPI1: serial ≡ n1 ≡ n3, including logs and migration ledgers), SLURM templates
and a chaining script on the exit-0/3 contract. **1.4.0** extended to K ≥ 3 (ternary
W₆₄Mo₃₂Cr₃₂, `full_pairs` [A16]), ran the validation campaigns, and unified everything under
`tools/run_all_tests.sh`. **1.5.0** delivered interoperability: `exsqs score` for external
structures (T-X1 bit-exact against engine records) and the [A11] geometric β schedule.
**1.6.0** mechanized project coherence (T-CO1, 17 spec↔code↔tests↔docs cross-checks as the
matrix's "step-0 test"). **1.7.0** is this release: `exsqs geom` + `tools/py/align_to_config.py`
close the scorer's orientation-strictness (T-X2: rotate + shuffle + align + score reproduces
engine records bitwise), plus CITATION.cff, CHANGELOG.md, a restructured README, and the
`v1.7.0` tag. Dated per-step reports live in `docs/STEP*_REPORT.md`.

## 2. Scientific findings

**The quantization floor.** On a finite supercell the achievable correlation error is
quantized: pair counts are integers, so `E_pure` has a provable lower bound `E_floor`
computable in closed form from the composition and zone structure (SPEC §4.1). On the paper's
own reference system the floor is 3.246078e-3 — yet the paper's Figure 1 reports converged
errors near 2.92e-4, *below* any value the defined objective can attain on that cell. The
inescapable conclusion is a hidden normalization in the plotted quantity (a per-term mean,
≈ E/14, would land exactly in its range), which is why this project treats absolute error
scales as non-comparable across codes and compares **raw structures only** (`exsqs score`,
T-V1/T-X1/T-X2); the discrepancy is flagged for an erratum note in any revision of the paper.
The e_tol defaults are therefore floor-relative (auto = 3×E_floor), which also handles the
K = 3 discovery that commensurate (dyadic) compositions make the floor **exactly zero**, where
only numeric tolerances are meaningful.

**Dynamics.** The γ-knee above is robust across seeds; γ-extreme runs reach R-3m cells with
D = 40 (19.2× fewer displacements). Exploitation regimes (fixed high β) descend fast and then
exhaust their ⟨g⟩-move neighborhoods against the uniqueness archive — the stagnation stop
[A13] retires them, and resume correctly refuses to re-arm the terminal state — an honest
negative result that motivated the migration defaults and the [A11] geometric schedule now
available for cluster-scale sweeps.

## 3. Engineering guarantees

One config + seed defines **one trajectory**, and the release proves it survives every
execution topology: same bits sequential vs. lockstep-OpenMP (any thread count), vs. MPI (any
rank count), vs. checkpoint/resume chains crossing migration rounds (T-K1) — outputs,
generation logs, and migration ledgers all compared. Resume is signature-guarded: every
trajectory-relevant `RunConfig` field is either serialized into the signature or on a
documented exclusion list, and the coherence audit fails if a future field goes unclassified —
so "it resumed" can never silently mean "it resumed a different run". Configuration follows a
loud-fail policy (unknown or inapplicable keys reject, with version-neutral messages). And the
project audits itself: T-CO1's seventeen checks run live (the §11 spec sample and every
reference config are *executed*, versions are checked against a real `summary.json`), and its
own history is instructive — its first runs caught drift (an uncited test id) *and a bug in
its own extractor*.

## 4. Verification status

`docs/TEST_MATRIX.txt` (committed, refreshed at release) records **24 gates, all green**: the
coherence audit; four core unit groups (including T-D1 and 20,826 RNG assertions); config,
engine, floor, and bitwise-rerun groups; parallel invariance; serializer + T-K1/T-K2; ternary;
scoring T-X1; schedule T-A11; non-diagonal H; the T-X2 round trip; T-E1/T-E2/T-E3 integration
(~1–3 min, load-dependent); T-MPI1; three T-V1 python cross-validations at |diff| = 0; and
two recorded benchmarks. Fast suite: 37 Catch2 cases, 27k+ assertions, ~10 s.

## 5. Reproducing everything

    cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
    bash tools/run_all_tests.sh          # the 24-gate matrix
    python3 tools/check_coherence.py     # the audit alone, ~2 s

The campaign evidence referenced throughout lives under `runs/` (untracked outputs; regenerate
with the commands in the step reports), so the T-V1-on-campaign rows and T-X2 reproduce as-is;
fresh campaigns are one `./build/exsqs configs/w70cr30_4x4x4.yaml --out runs/mine` away.

## 6. Comparing against the paper (or ATAT)

For each external structure: convert if needed, align, score —

    python3 tools/py/to_poscar.py their_file.cif their.vasp
    ./build/exsqs geom configs/w70cr30_4x4x4.yaml -o geom.vasp
    python3 tools/py/align_to_config.py geom.vasp their.vasp aligned.vasp
    ./build/exsqs score configs/w70cr30_4x4x4.yaml aligned.vasp --json theirs.json

then compare E_pure/E_floor, D, and SG directly against the engine's `summary.json` — same
geometry, same zones, same objective, no scale ambiguity.

## 7. What requires cluster hardware (the remaining work)

Multi-node scaling curves (`scripts/slurm/exsqs_mpi.sbatch`; islands shard rank-invariantly, so
results are comparable across node counts by construction); the long-budget frontier — γ sweeps
at β-schedule settings without the exhaustion limit seen at small budgets
(`scripts/chain_resume.sh` chains to convergence on the exit-0/3 contract); and K ≥ 3 physics
campaigns on the ternary config. Every such run remains bit-reproducible and scoreable with the
tools above.

## 8. Honest limitations

Single sublattice and pair correlations only (the paper's scope); `diagonal` and `full_pairs`
error modes; POSCAR as the only writer; the scorer is orientation-strict by design (the
aligner handles frames; sites relaxed beyond 1e-4 fractional are rejected because exsqs scores
ideal decorations); wall-clock-capped runs are not comparable across machine topologies
(generation-capped runs are); absolute error values are not comparable across codes — that one
is not a limitation of this implementation but of the field's reporting conventions, and it is
exactly what `exsqs score` exists to bypass.

## 9. Release facts

Version 1.7.0; SPEC v1.7; tag `v1.7.0`. Before publishing externally: set a personal copyright
holder in LICENSE if desired (it currently names the project).
