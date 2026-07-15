# EXSQS — Specification v1.8 (2026-07-15) — v1.0 frozen 2026-02-16; changes tracked in §16

Supersedes the v0.1 draft (`SPEC_step0.md`). All Step-0 blocking decisions are resolved (§14);
Step 0 is complete. Changes vs draft: Q1→[A15], Q2→[A16], Q3→§1 scope; new test T-C4.

Reference algorithm: A. P. Kądzielawa, *On generating Special Quasirandom Structures: Optimization
for the DFT computational efficiency*, arXiv:2602.10872 (2026). Alg. 1 (extinction workflow),
Alg. 2 (error function).

Tags used below:
- **[A#]** — adopted convention (binding for v1)
- **[D#]** — documented deviation from / disambiguation of the paper's pseudocode

---

## 1. Scope of v1

In scope:
- Single active sublattice: prototypes whose decorated sites form one Wyckoff orbit
  (sc, bcc, fcc, hcp built-in; arbitrary POSCAR/CIF prototype accepted if its sites are one orbit).
- **K-component data model from day one; v1 *validated* scope is binary (K = 2).** K-ary code
  paths are compiled and unit-tested but flagged `experimental` at runtime until the K-ary
  validation pass (v1.1). *(Q3 resolution.)*
- Fixed supercell per run: diagonal `n1×n2×n3` or full 3×3 integer matrix H (shape *scan* = outer
  loop over configs, embarrassingly parallel; shape *optimization* is out of scope).
- Pair correlations by coordination zone (paper-faithful). Triplets/CE clusters: out of scope v1.
- Objective = correlation error × displacement-count penalty; pure-SQS mode via `gamma: 0`.

Out of scope v1: multi-sublattice decoration, GPU, on-the-fly cell-shape search, cluster orders > 2,
relaxed/strained geometries (ideal lattice positions only).

## 2. Notation and data model

- Prototype lattice cell `A0` (3×3, columns = vectors), basis sites; supercell `A = A0·H`,
  `N = m·n_basis` sites at ideal fractional coordinates `f_i` (exact rationals by construction).
- Species set `T = {t_1..t_K}`; decoration `σ: {1..N} → T`; counts `N_t`, achieved concentrations
  `x̃_t = N_t/N` (§9).
- Background `B`: the 27 periodic copies of the supercell, offsets `(i,j,k) ∈ {−1,0,1}³` (Alg. 2:2–7).
  Ordered pairs `(i ∈ S, j ∈ B)` with `d_ij > 0` — i.e., the zero-distance self pair is excluded,
  but nonzero self-images **are** legitimate PBC neighbors and are included. **[A1]**
- `D(σ)`: number of irreducible finite-displacement calculations for the phonon run (§6).

## 3. Geometry: coordination zones

- Zones are a property of the **undecorated** geometry: sorted pair distances of `(S,B)` merged into
  shells with tolerance `ε_shell` (relative, default `1e-3`); shell n has radius `R_n` and per-site
  count `Z_n` (constant across sites because active sites are one orbit — v1 scope guard). **[A2]**
- Pair→zone lookup precomputed **once per (lattice, H, n_shells)** and shared by all evaluations;
  stored as per-site neighbor lists up to shell `n_max` (not the full 27N² table).
- `n_max` default 7. If `R_{n_max}` exceeds half the minimal supercell width, emit a warning
  (correlations there are constrained by self-images — allowed, but user should know). **[A3]**
- Sanity identities: `Σ_i (neighbors of i in shell n) = N·Z_n`; bcc coordination sequence
  8, 6, 12, 24, 8; fcc 12, 6, 24, 12, 24 (used in tests, §12).

## 4. Correlation and error function

Counters over ordered pairs within shells `n = 1..n_max`, for center type `t`, neighbor type `t'`:

```
C_{t,t'}(n) = #{ (i,j) : σ_i = t, σ_j = t', zone(i,j) = n }
P_t(n)      = Σ_{t'} C_{t,t'}(n)                     # pairs whose central atom is type t
Π_{t,t'}(n) = C_{t,t'}(n) / P_t(n)                   # P(neighbor is t' | center is t), Σ_{t'} Π = 1
```

**[D1]** Alg. 2:19 `population` is read as `P_t(n)` above (pairs whose *central* atom is type t).
Justification: only then does the random limit give `Π_{t,t}(n) → x_t`, which is what line 20
compares against. Normalizing by all pairs would give `x_t²` and never converge to `x_t`.

**[D2]** Alg. 2:20 `correlation = x[atomType]` is read as `-=` (subtraction). The error term is the
absolute deviation from the random limit.

**[A16] Error modes** *(Q2 resolution)* — config `error.mode: auto | diagonal | full_pairs`,
default `auto` = `diagonal` for K = 2, `full_pairs` for K ≥ 3:

```
diagonal   (paper-faithful):  E_pure = Σ_n A_n Σ_t     | Π_{t,t}(n)  − x̃_t  |
full_pairs (K ≥ 3 default):   E_pure = Σ_n A_n Σ_{t,t'} | Π_{t,t'}(n) − x̃_t' |
E_obj = E_pure · D^γ                                  # γ = 1 default; γ = 0 → pure-SQS mode
```

- Rationale: for K = 2 the off-diagonal deviations are fully determined by the diagonal ones; for
  K ≥ 3 diagonal-only under-constrains the SRO matrix.
- Exact identities (used as tests): `C_{t,t'}(n) = C_{t',t}(n)` (ordered-pair bijection), and for
  **K = 2: `E_full ≡ 2·E_diag` for every decoration** — mode choice never changes rankings for
  binaries, only the scale. Consequence: `e_tol` is interpreted in the *active* mode's scale;
  logs always report both `E_diag` and `E_full` for K = 2.
- Weights `A_n`: `1/n` default; options `n^{−p}`, `R_n^{−1}`, explicit list (paper tested these;
  `1/n` converged best for cheap supercells). **[A4]**
- **[D3]** Exponent γ generalizes Alg. 2:22 (`γ=1` reproduces the paper; `γ=0` gives an
  ATAT/sqsgenerator-comparable pure error for benchmarking).
- Interpretation (documentation + validator): `|Π_{t,t}(n) − x̃_t| = x̃_t·|α_tt(n)|` with the diagonal
  Warren–Cowley parameter `α_tt(n) = 1 − Π_{t,t}(n)/x̃_t`. So E_pure is a weighted sum of
  concentration-scaled |WC SRO| values.
- Targets are the **achieved** concentrations `x̃_t`, not the requested ones (§9). **[A5]**
- Note for tests: a random *canonical* decoration has `E[Π_{t,t}(n)] = (N_t−1)/(N−1)`, not `x̃_t`
  (hypergeometric, offset `(1−x̃_t)/(N−1)`); the statistical test (§12, T-C3) must expect this.

### 4.1 Quantization floor E_floor (v1.1)

`C` counters are integers while `P_t(n) = N_t·Z_n` is fixed by composition and geometry, so every
error term is bounded below: `|Π_{t,t'}(n) − x̃_{t'}| ≥ dist(x̃_{t'}·P_t(n), ℤ)/P_t(n)`. Summing the
independent per-term minima gives a computable lower bound on E_pure:

```
E_floor = Σ_n A_n Σ_t dist(x̃_t·P_t(n), ℤ) / P_t(n)     # diagonal; full_pairs sums t' too
                                                        # (K = 2: exactly 2×E_floor_diag)
```

Reference system W₉₀Cr₃₈ / bcc 4×4×4 / 7 shells / A_n = 1/n: `E_floor = 3.246e-3`. The engine
computes and logs E_floor at startup and warns when the requested `e_tol` lies below it
(infeasible; the run then terminates on caps only). For *commensurate* compositions every
`x̃_t·P_t(n)` can be integral and `E_floor = 0` exactly (e.g. W₆₄Mo₃₂Cr₃₂ on bcc 4×4×4, v1.4):
the bound is then uninformative and `e_tol: auto` degenerates to 0 — set a numeric tolerance.

**Scale caveat (v1.1):** the paper's reported absolute errors (Fig. 1: 𝔈 = 2.01–2.92·10⁻⁴ for this
very system) lie *below* this provable bound, so the paper's reported scalar carries an unstated
normalization. Absolute 𝔈 values are not comparable across codes; comparisons must recompute the
error on raw structures (the T-V1 route). v1.0's absolute defaults (`e_tol: 3e-4`, T-E1's `1e-3`),
inherited from the paper's scale, were infeasible and are replaced by floor-relative values
(§11, §12).

## 5. Symmetry and filtration

- Backend: spglib on ideal decorated cells (positions exact ⇒ tolerance-robust);
  `symprec` default `1e-5` Å, logged. **[A6]**
- Filter policy (Alg. 1:4, 1:21), configurable: default `reject_p1` = reject iff space-group
  number 1 (P1). Options: `min_sg_index N`, explicit whitelist, `min_pointgroup_order N`. **[A7]**
- P1-elite quota `q ≥ 0` (paper's "small ≪ P best P1 structures" variant); default `q = 0`. **[A8]**
- Every emitted structure records: SG number + international symbol, symprec used, point-group order,
  number of inequivalent sites.

## 6. Displacement count D

- **[A15] Pinned convention** *(Q1 resolution)*: `D(σ)` := the number of displaced supercells
  phonopy generates for σ used *directly as the phonon supercell* (`supercell_matrix = 1`) with
  **default settings** (`is_diagonal = True`, `is_plusminus = 'auto'`). Displacement distance is
  irrelevant to the count. Any alternative convention in the future is a *new* enum value of
  `displacements.convention`, never a redefinition of `phonopy_default`.
- Implementation: native C++ from the spglib symmetry dataset — inequivalent atoms
  (`equivalent_atoms`) → site-symmetry group per representative → minimal direction set per
  phonopy's reduction, ± doubling when site symmetry cannot reverse the displacement.
- **Hard validation gate:** exact match against `phonopy` on the full test set (§12, T-D1) before
  the optimizer is trusted.
- Reference magnitudes: P1, 128 atoms ⇒ 6·128 = 768 displacements (identity site symmetry, ±, 3 dirs);
  the paper's C2/Amm2-class cells reach ~5× fewer than P1 ATAT cells.

## 7. Uniqueness and canonical form

- Two decorations are duplicates iff related by a symmetry operation of the **undecorated**
  supercell (lattice translations of H-cell × point ops of the empty supercell), acting on site
  indices. **[A9]**
- Canonical form: precompute the permutation group Π of the undecorated supercell once (from
  spglib on the empty cell); canonical label = lexicographic minimum of `π(σ)` over `π ∈ Π`;
  dedup via 128-bit hash of the canonical label, hash set per island + merge at output. Cost is
  `O(|Π|·N)` per new structure — acceptable at these sizes; a cheap invariant pre-filter
  (composition + per-shell `C_{t,t'}(n)` vector) short-circuits most comparisons. **[A10]**
- Rare hash collisions: accepted risk (verify by full compare on collision).

## 8. Evolutionary operators and termination

- **Seeding** (`rejection | constructive | mixed`, default `mixed`): **[D4]**
  - *rejection* (paper-faithful): uniform random canonical decorations, keep non-P1 & unique.
  - *constructive*: pick a target space group G from a pool of subgroups of the empty supercell's
    group; partition sites into G-orbits; assign species to whole orbits so that composition is met
    (subset-sum over orbit sizes, DP solver); guarantees symmetry ⇒ immune to the
    rejection-collapse failure mode at large N.
- **Evaluation**: `E_obj` for every population member (embarrassingly parallel).
- **Extinction** (Alg. 1:14–17): survival probability `E_min/E_i` (`ratio`, default) or
  `min(1, exp(−β(E_i − E_min)))` (`metropolis`). β: fixed value or `auto` =
  `ln 2 / median(E_i − E_min)` recomputed each generation (median structure survives 50%);
  optional geometric schedule. **[A11]** v1.5 makes the schedule normative: with a numeric β,
  `β_g = β·growth^g` (`survival.schedule: geometric`, `beta_growth > 0`); `beta: auto`
  recomputes per generation and rejects the schedule (loud-fail config policy).
- **[D5]** The current best-by-`E_obj` always survives (`elitism_best = 1`) ⇒ `E_min` monotone.
  (Paper does not state this; without it the stochastic extinction can lose the optimum.)
- **Repopulation**: parent = uniform random survivor; mutation = `k` composition-preserving swaps
  of unlike-species site pairs (`k = 1` default; option `k ~ 1 + Poisson(λ)`); accept offspring iff
  passes filter + unique; retry budget per slot (default 100), then fall back to constructive
  seeding for that slot. **[A12]**
  - Optional `symmetry_preserving: true`: swap species assignments between equal-size orbits of a
    subgroup H of the parent's group (keeps H by construction; avoids the
    mutate→P1→reject loop). Recommended default `true`, with plain swaps as fallback. **[D6]**
- **Termination [D7]:** disambiguation of Alg. 1:9–12 (E "includes displacements", success tests
  `E_min/D_min`): selection pressure uses `E_obj`; success is declared when
  `min_i E_pure(σ_i) ≤ E_tol` (the tolerance certifies "is a proper SQS"; D only steers selection).
  Additional caps: `max_generations`, wall-time budget → checkpoint + best-effort output. **[A13]**
  v1.1 adds a *stagnation stop*: if repopulation degenerates to all-fallback seeding
  (`fb ≥ vacancies`) for `stagnation_stop` consecutive generations (default 3), the island
  stops with reason `stagnation`.
- **Output set [D8]:** Alg. 1:2's outer loop is realized as: each island evolves a population of
  size P; the M outputs are the best M *pairwise-inequivalent* structures pooled across
  islands/restarts, ranked by `E_obj` (report `E_pure`, `D`, SG for each).

### 8.1 Parallel execution model (v1.2)

- **Lockstep islands:** all active islands advance one generation per round; a stopped island
  (success/caps/stagnation) freezes and no longer participates. `max_wall_s` and the logged
  `elapsed_s` measure *global run* wall time in this model.
- **Synchronous ring migration** among the still-active islands every `migration_every` rounds:
  island *i* sends copies of its `migrants` pool-best individuals to the next active island
  (ascending-index ring). Receivers skip canonicals already in their archive; otherwise the
  migrant is archived and replaces the receiver's current worst member (the [D5] elite is never
  the worst for P ≥ 2, so `E_min` stays monotone). Exchanges are computed from pre-round
  snapshots and use no RNG — fully deterministic.
- **Round-based generations:** within seeding and repopulation, candidate generation +
  canonicalization and per-candidate evaluation (spglib, D, E) run in parallel; the duplicate
  check with *provisional* archive insert and the quota commit run serially in slot order.
  Provisional insertion is outcome-independent (v1.1 archived admitted and P1-rejected candidates
  alike), so with the default `p1_elite_quota: 0` and no post-seed P1 admissions the rounds are
  bit-equivalent to the v1.1 sequential loop. A P1-rejected slot resumes its own attempt stream
  next round (normative semantics for the rare quota/P1 cases).
- **Thread budget:** `omp_threads` (auto = OpenMP default) split as outer = min(islands, T) over
  islands, inner = T/outer inside each island's evaluation rounds. **Thread counts never change
  results** ([A14]) — verified bitwise for 1 vs 4 threads (T-P1) and against v1.1 sequential run
  artifacts.

### 8.2 Checkpoint/restart and distributed execution (v1.3)

**State files.** The driver writes `state.ckpt` (atomic tmp+rename) into the output directory at
every `checkpoint_every` round and at run end. The file carries a format version and a
*trajectory signature* — every config field that influences the trajectory (system, composition,
zones, error model, evolution dynamics, seed, islands, migration). `--resume DIR` restores the
complete engine state (population, archive, pool, log, counters) bit-exactly; a signature
mismatch is refused, so only budget caps (`max_generations`, `max_wall_s`) and output/threading
fields may change between segments. Islands stopped on a raisable cap re-arm on resume; `e_tol`
success and `stagnation` stops are terminal. Chained segments are bit-identical to an
uninterrupted run (T-K1), including across migration rounds.

**MPI ranks (`exsqs_mpi`).** Islands are distributed round-robin (island *i* → rank *i* mod *R*)
and advance in the same lockstep rounds; activity flags and §8.1 migration payloads travel as
little-endian byte records over collectives, and rank 0 writes the identical `state.ckpt` and
pooled outputs. Because trajectories are island-keyed [A14] and the exchange protocol is
deterministic, **results are bit-identical to the single-process driver for any rank count**
(T-MPI1); serial and MPI state files are interchangeable, so a run may move between 1 and N
nodes across segments. Homogeneous little-endian clusters are assumed (guarded at every
file/wire boundary). Wall-capped stops remain wall-clock-dependent and are excluded from
bit-exactness claims, as in v1.2.

## 9. Composition handling

- `N_t = round(x_t·N)` by largest-remainder so `Σ N_t = N`; achieved `x̃_t` logged and used as the
  correlation target **[A5]**. Hard error if `|x̃_t − x_t| > comp_tol` (default 0.02).
- Consistency check vs paper: 70:30 on 128 sites → W₉₀Cr₃₈ (x̃ = 0.7031/0.2969) — matches the
  paper's stated supercell.
- Constructive seeding feasibility: composition must be a subset-sum of orbit sizes for the chosen
  G; infeasible G are skipped (logged).

## 10. RNG and reproducibility

- Counter-based generator (Philox/PCG) keyed by `(master_seed, island, generation, slot, purpose)`
  ⇒ trajectories are independent of thread count and scheduling; bit-identical reruns from the same
  config+seed. **[A14]**
- v1.2 realizes the thread-invariance guarantee end-to-end (verified bitwise for 1 vs 4 threads
  and against v1.1 sequential artifacts); wall-time semantics are global in the lockstep model.
- `summary.json` records: config echo, achieved composition, seeds, git commit, spglib/phonopy
  versions, timings.

## 11. Config schema (v1.0; YAML, CLI-overridable)

```yaml
lattice:    {type: bcc, a: 3.165, file: null}          # sc|bcc|fcc|hcp|file (POSCAR/CIF prototype)
supercell:  {diag: [4, 4, 4]}                          # or matrix: [[...],[...],[...]]
composition: {W: 0.70, Cr: 0.30}                       # largest-remainder rounding; achieved x logged
zones:      {n_shells: 7, shell_tol: 1.0e-3}
error:
  weights:  {form: inv_n}                              # inv_n | inv_n_pow: p | inv_R | custom: [...]
  mode:     auto                                       # auto = diagonal (K=2) | full_pairs (K>=3) [A16]
  gamma:    1.0                                        # E_obj = E_pure * D^gamma ; 0 = pure-SQS mode
evolution:
  population: 200
  outputs: 10                                          # M
  e_tol: auto                                          # v1.1: auto = 3.0×E_floor (§4.1); numeric override; warn if < E_floor                                        # interpreted in the ACTIVE error mode's scale
  max_generations: 5000
  stagnation_stop: 3                                   # v1.1 (§8 [A13]); 0 disables
  survival: {mode: ratio}                              # ratio | metropolis: {beta: auto, schedule: const}
  elitism_best: 1
  p1_elite_quota: 0
  mutation: {swaps: 1, symmetry_preserving: true, retry_budget: 100}
  # v1.1: swaps also accepts 'poisson' (+ lambda): k ~ 1+Poisson(λ) [A12]
  seeding:  {mode: mixed, sg_pool: auto}
symmetry:   {symprec: 1.0e-5, filter: {policy: reject_p1}}
displacements: {convention: phonopy_default}           # pinned [A15]
parallel:   {islands: 1, omp_threads: auto, migration_every: 50, migrants: 2}
rng:        {seed: 42}
output:     {dir: ./run1, formats: [poscar], checkpoint_every: 100, log_level: info}
```

## 12. Test matrix

| ID | Level | Content | Acceptance |
|----|-------|---------|-----------|
| T-G1 | unit | Shell radii + coordination sequences, sc/bcc/fcc (bcc: 8,6,12,24,8; fcc: 12,6,24,12,24) | exact |
| T-G2 | unit | Pair table: `zone(i,j)=zone(j,i)`, `C_{t,t'}=C_{t',t}`, self excluded, `Σ_i n_i(n)=N·Z_n` | exact |
| T-C1 | unit | B2 (CsCl) decoration: `Π_{t,t}(1)=0`, `Π_{t,t}(2)=1`; hand-computed E_pure | exact |
| T-C2 | unit | L1₀ and L1₂ on fcc: analytic Π per shell | exact |
| T-C3 | statistical | 10⁴ random canonical decorations: mean `Π_{t,t}(n) = (N_t−1)/(N−1)` | within 3σ |
| T-C4 | unit | K=2 mode identity: `E_full = 2·E_diag` on random decorations | ≤ 1e-15 rel |
| T-S1 | unit | SG detection: B2→221, L1₀→123, L1₂→221, pure bcc supercell→229; wrapper vs pymatgen | exact |
| T-S2 | unit | Dedup: translated/rotated copies collide; 10⁴ random distinct pairs don't (fuzz) | exact / no false-pos |
| T-D1 | unit (gate) | D-count vs `phonopy` on ≥12 cells: pure, B2, L1₀, L1₂, and low-symmetry SQS incl. C2/Amm2-class | exact match |
| T-E1 | integration | W₇₀Cr₃₀ bcc 4×4×4 (binary), `gamma=0`: reaches `E_pure ≤ e_tol(auto) = 3.0×E_floor` (§4.1) and ≥ M unique non-P1 outputs within budget | pass |
| T-E2 | integration | same, `gamma=1`: best output has `D ≤ D(P1 reference)/3` at `E_pure ≤ 3×e_tol(auto)` (v1.1 reading of “within ~3× of T-E1 best”) | pass |
| T-R1 | reproducibility | same seed/config, 1 vs 16 threads, 1 vs 4 islands re-run: identical trajectories/outputs | bit-exact |
| T-V1 | cross-validation | `validate.py` (pymatgen/numpy) recomputes E_pure of every emitted structure | ≤ 1e-10 abs diff |
| T-B1 | benchmark | evaluations/s vs N (54…1024 sites); OMP scaling on one node; island scaling 1→32 ranks | recorded, no gate |
| T-P0 | unit (v1.2) | symmetry/displacement/correlation pipeline under 4 concurrent threads vs serial | identical |
| T-P1 | reproducibility (v1.2) | `omp_threads` 1 vs 4, islands+migration and single-island paths | bit-exact |
| T-M1 | integration (v1.2) | ring migration exchanges pool-best; dedup drops only; deterministic rerun | pass |
| T-K1 | reproducibility (v1.3) | budget-stop → `--resume` chain vs uninterrupted run (migration boundary crossed) | bit-exact |
| T-K2 | unit (v1.3) | resume refuses trajectory-signature mismatch; raising budget caps allowed | pass |
| T-MPI1 | reproducibility (v1.3) | serial vs `mpirun -n 1` vs `-n 3`: outputs, logs, migration ledger | bit-exact |
| T-E3 | integration (v1.4) | ternary bcc W₆₄Mo₃₂Cr₃₂ (`full_pairs` [A16]): `E_pure ≤ 1.2e-1`, non-P1, composition preserved | pass |
| T-X1 | unit (v1.5) | `score_structure` reproduces engine records bit-exactly; input site order free; loud mismatch errors | pass |
| T-A11 | unit (v1.5) | geometric β schedule: growth 1.0 ≡ const bitwise; growth > 1 diverges deterministically | pass |
| T-CO1 | audit (v1.6) | spec↔code↔tests↔configs↔docs coherence: versions, §12 id citations, tag coverage, signature field ledger, §11 sample + configs execute, README refs | pass |
| T-X2 | tool (v1.7) | geom → rotate+permute → align_to_config → score reproduces engine records bitwise (cross-tool, file-level) | pass |

K-ary validation (full_pairs end-to-end, ternary integration case) is deferred to v1.1 per Q3;
until then K ≥ 3 runs print an `experimental` banner.

## 13. Deviations from the paper (summary)

D1 normalizer = pairs with central atom of type t; D2 line 20 is `-=`; D3 γ-exponent
generalization of the ×D penalty; D4 constructive (subgroup-orbit) seeding added alongside
rejection; D5 best-structure elitism added; D6 symmetry-preserving mutation added; D7 termination
= `min E_pure ≤ E_tol` with selection on `E_obj`, plus generation/wall-time caps; D8 outer loop
realized as island pooling of M unique outputs.

Each deviation is either paper-faithful under a config switch or strictly additive; defaults are
chosen to reproduce the paper's setup where the pseudocode is unambiguous.

## 14. Resolved decisions (Step-0 sign-off, 2026-02-16)

| # | Question | Resolution |
|---|----------|-----------|
| Q1 | D-count convention | **phonopy defaults** (`is_plusminus='auto'`, `is_diagonal=True`, supercell used as-is) → **[A15]** |
| Q2 | K-ary error mode default | **auto: diagonal for K=2, full_pairs for K≥3** → **[A16]** |
| Q3 | v1 validated scope | **Binary first, K-ary later** → §1, §12 |

Non-blocking assumptions accepted (no veto received): symprec 1e-5 [A6]; filter = reject SG #1
only [A7]; targets = achieved concentrations [A5]; mixed seeding default [D4]; license
Apache-2.0; toolchain GCC ≥ 12 + OpenMPI; spglib vendored via CMake FetchContent (BSD).

## 15. Status and handoff to Step 1

Step-0 exit criteria **met** (Q1–Q3 answered; A/D items accepted). This document is the binding
reference; changes from here require a version bump and a changelog entry.

Steps 0–1 are realized by the released 1.0.0 core library and its suites (`v1.0.0` tag,
`feat/tests-core`): the CMake ≥ 3.20 + FetchContent skeleton, T-G1/T-G2, T-C1–T-C4, T-S1/T-S2
and the T-D1 phonopy gate, over `lattice` → `structure` → `zones` → `correlation` → `symmetry`
→ `displacements` → `dedup`. Step 2 (the full §8 engine) proceeds tests-first from there.

## 16. Changelog

- **v1.8 (2026-07-15)** — housekeeping release: repository-wide `*.tex` ignore (manuscript and
  presentation sources live outside the code history); first public mirror at
  github.com/Mellechowicz/SQuatS; LUMI-C validation campaigns T1/T2 of the two-level parallel
  design executed (single-node probe Amdahl p = 0.479 dominated by the fixed serial floor;
  island throughput saturating at ranks = islands, mean lockstep imbalance 133 ms per
  generation). No engine behaviour changes; trajectories bit-identical to v1.7.1.

- **v1.7 (2026-07-12, post-release note 3)** — 16-cell supercell scan (2×2×2 → 7×7×7, K = 5
  equiatomic bcc, γ ∈ {0, 1}; docs/SUPERCELL_STUDY.md): the floor law **E_floor = 0 ⟺ 25 | N**
  holds on all sixteen cells; equiatomic x = 1/5 is unrepresentable within `comp_tol` below 54
  sites; the γ = 1 ordered-compound collapse (Cmcm, D = 10, E_pure = 7.938889 bit-identical at
  200 and 300 sites) reproduces exactly on the two 5,5-axis cells and nowhere else; at ≥ 504
  sites γ = 1 reaches 4× displacement reductions at correlation error equal to or better than
  γ = 0.
- **v1.7 (2026-07-10, post-release note 2)** — supercell study fix: `--resume` without
  `--out` silently retargeted outputs and subsequent state saves to the config's `output.dir`
  (observed overwriting a sibling run at 6×6×6). Resume now defaults `output.dir` to the state
  directory and says so on stderr; an explicit `--out` overrides as before.
- **v1.7 (2026-07-08, post-release note)** — K = 5 exercise (HEA session): equiatomic bcc
  MoNbTaVW (Mo₂₆Nb₂₆Ta₂₆V₂₅W₂₅, 128 sites) run end-to-end; **T-V1 at K = 5 passes with
  |diff| = 0** (`tools/py/validate.py` independently reproduces every `full_pairs` energy).
  Finding: E_floor = 8.0379e-2 — ~25× the binary floor; the K = 3 dyadic zero-floor phenomenon
  does not recur at K = 5.
- **v1.7 (2026-07-03)** — step-8 release 1.7.0:
  - `exsqs geom` writes the undecorated supercell frame; `tools/py/align_to_config.py`
    (pymatgen StructureMatcher with a species-blind FrameworkComparator) maps rotated /
    permuted / relabelled external files onto it, closing the scorer's deliberate
    orientation-strictness for real paper/ATAT inputs. **T-X2** gates the whole pipeline:
    rotate the recorded champion 90°, shuffle its sites, align, score — bitwise equal to the
    engine's records.
  - Release packaging: LICENSE (MIT; holder to be set by the maintainer), CITATION.cff (paper
    + software), CHANGELOG.md, README restructured with accreted
    step notes moved to docs/DEV_NOTES.md, tagged v1.7.0.
  - Version scheme: plain triples continuing the 1.0.0 line; T-CO1's C1 requires the CMake
    VERSION to equal the driver literal exactly.
- **v1.6 (2026-06-24)** — step-7 coherence audit:
  - **T-CO1 / `tools/check_coherence.py`**: a mechanical audit joins the matrix as its "step-0
    test" — one version literal everywhere (checked live against `summary.json`), every §12 test
    id cited by an actual test/tool, every Catch2 tag exercised by the runner, no version claims
    in user-facing strings, a formal *signature field ledger* (every `RunConfig` field is either
    in the trajectory signature or on the documented exclusion list: budgets, output/logging,
    `omp_threads` [A14], `x_target` [A5]), the §11 sample config and every `configs/*.yaml`
    executed end-to-end, SPEC header == latest changelog version, sections 1–16 present, and all
    README-referenced paths/binaries resolving.
  - Drift found and fixed by the audit's first run: T-A11 was defined in §12 but cited nowhere
    (schedule tests referenced decision [A11] only — now carry the id). The first run also caught
    a bug in the auditor itself: the struct-field extractor stopped at the first `};`, which any
    brace-initialized member produces — fixed with brace counting.
  - Policy, now enforced: step reports are immutable dated snapshots; SPEC, README, configs and
    code must stay mutually current.
- **v1.5 (2026-06-09)** — step-6 interoperability + spec completion:
  - **`exsqs score`** subcommand (score.hpp): evaluates external POSCARs on the config geometry —
    E_pure, E/E_floor, D, D(P1)/D, SG, E_obj — with free input site order (coordinate matching),
    strict lattice/species/composition checks, and optional `--json`. This operationalizes the
    §4.1 doctrine: cross-code comparisons happen on raw structures (paper/ATAT outputs), never on
    reported scalars. T-X1 pins bit-exact agreement with the engine's own records.
  - [A11] geometric β schedule implemented (`survival.schedule: geometric`, `beta_growth`);
    growth = 1.0 is bit-identical to const. The two fields enter the trajectory signature, and
    enabling the schedule invalidates earlier state files (correct: the dynamics vocabulary changed).
  - Non-diagonal `supercell.matrix` (parsed since the 1.0.0 schema) gains first test coverage;
    `tools/py/to_poscar.py` converts CIF/anything-pymatgen-reads to POSCAR, closing the CIF
    prototype deferral pragmatically.
- **v1.4 (2026-05-19)** — step-5 validation + release:
  - **K ≥ 3 validated end-to-end**: [A16] `mode: auto → full_pairs` exercised; the §4.1 floor is
    K-general (K=2's exact `full = 2×diag` identity does not extend); constructive seeding, orbit
    mutation and the canonical dedup are species-count-agnostic; T-V1 reproduces the K=3
    full-pairs error at |diff| = 0. New reference config `w50mo25cr25_4x4x4.yaml` and tests
    (`[ternary]`, incl. T-E3).
  - §4.1 note: commensurate compositions give `E_floor = 0` exactly (dyadic x on bcc 4×4×4);
    use a numeric `e_tol` there.
  - Campaigns recorded (single-core container scale): binary flagship (6 islands + migration,
    β = 3000) reconfirms 2.50×E_floor with rc = 0 in 102 s / 14.8k evals; the β = 800 + migration
    resume chain is stagnation-terminated near generation 50 at 3.61× — exploitation regimes are
    exhaustion-limited by design, so the genuine long-budget frontier belongs to the
    non-exhausting survival modes on multi-core hardware. Ternary campaign best:
    9.55e-2 (P2) at 25 generations.
  - Release: unified verification runner `tools/run_all_tests.sh` (steps 1–5 matrix), CI stub
    (`.github/workflows/ci.yml`), version 1.4.0.
- **v1.3 (2026-04-28)** — step-4 HPC layer:
  - §8.2: checkpoint/restart (signature-guarded `state.ckpt`, atomic writes, raisable budget
    caps, bit-exact chains) and the rank-distributed `exsqs_mpi` driver (round-robin islands,
    collective-based §8.1 migration, rank-invariant to the bit; T-B1's "island scaling over
    ranks" is now runnable). Serial and MPI state files are interchangeable.
  - New binary encoding (`serialize.hpp`): little-endian, IEEE-754 bit patterns; homogeneous
    little-endian clusters assumed and guarded. Engine exposed via `EngineHandle`
    (island_engine.hpp) for external drivers and tests.
  - CLI `--resume DIR`; SLURM templates (`scripts/slurm/`) and a resubmission-chain pattern
    (`scripts/chain_resume.sh`) keyed to the 0/3 exit-code contract.
  - New tests T-K1/T-K2/T-MPI1. Mid-run `checkpoint.json` snapshots are not written under MPI
    (state.ckpt supersedes them there); per-generation stdout rows come from rank 0's islands.
- **v1.2 (2026-03-31)** — step-3 parallel layer:
  - §8.1: lockstep islands, synchronous ring migration (defined; previously only schema keys), and
    round-based generations with provisional archive inserts; bit-equivalent to v1.1 for the
    default quota and empirically preserved on all reference runs (verified against 1.1.0
    artifacts to the last bit, including under `omp_threads: 4`).
  - `parallel.omp_threads` active; nested outer(islands)×inner(evaluation) OpenMP; exceptions
    captured per island and rethrown deterministically; checkpoint callback serialized.
  - `summary.json` gains `island_migrants_in/out`. `max_wall_s`/`elapsed_s` are global run wall
    time. New tests T-P0/T-P1/T-M1; engine refactored into a resumable IslandEngine (the
    checkpoint/restart foundation for step 4). MPI multi-node island scaling ("ranks" in T-B1)
    is deferred to the HPC step; `exsqs_bench_scaling` records single-node numbers.
- **v1.1 (2026-03-10)** — step-2 integration findings:
  - §4.1: quantization floor `E_floor` defined, computed and logged; **the paper's absolute error
    scale is not reproducible from Alg. 2** (reported 2.92e-4 < floor 3.246e-3 for its own
    reference system) — cross-code comparisons go via raw structures only (T-V1).
  - §11: `e_tol: auto` (= 3.0×E_floor) is the new default; numeric values below E_floor trigger an
    infeasibility warning. §12 T-E1/T-E2 acceptance re-pinned floor-relative.
  - [A12]: Poisson swap-count option implemented (`swaps: poisson`, `lambda`; k ~ 1+Poisson(λ)).
  - [A13]: stagnation early-stop (all-fallback repopulation × `stagnation_stop` generations).
  - [D6] refined: H = cyclic ⟨g⟩ of one randomly chosen nontrivial parent op per move, replacing
    the full-stabilizer variant, which ratchets symmetry monotonically upward and traps the search
    in high-symmetry strata (empirical stall at ~5× floor; see STEP2_REPORT).
  - Empirical serial reference (seed 42): 4 islands × metropolis β=3000 reach 2.50×E_floor
    (E_pure = 8.11e-3, SG 8/Cm) on W₉₀Cr₃₈; the productive island crossed the gate at gen 10.
- **v1.0 (2026-02-16)** — frozen after step-0 sign-off (§14).
