# EXSQS — Supercell-Scaling Study: K = 5 equiatomic bcc, 2×2×2 → 7×7×7

Date 2026-07-12; engine 1.7.1 (seed 20260704 throughout). This study runs the full diagonal
supercell ladder 2×2×2, 2×2×3, 2×3×3, 3×3×3, 3×3×4, 3×4×4, 4×4×4, 4×4×5, 4×5×5, 5×5×5, 5×5×6,
5×6×6, 6×6×6, 6×6×7, 6×7×7, 7×7×7 (16 bcc cells, 16 → 686 sites) for the 5-component equiatomic
MoNbTaVW-class alloy, at γ = 0 (pure SQS) and γ = 1 (displacement-penalized), under identical
evolution settings (population 100, 4 islands, ring migration every 8, ratio survival,
60-generation budget, `full_pairs` error [A16], 5 shells). It extends the earlier seven-cell
scan (4×4×4 → 6×6×6) down into the under-resolved small-cell regime and up to 686 sites,
re-tests the two empirical laws found there, and records the failure modes that matter for
production use. Every run is bit-reproducible from its config; `tools/py/validate.py`
re-derives the reported errors independently (spot checks below at |diff| = 0).

## The ladder at γ = 0 (pure SQS)

| cell | N | \|PG\| | \|Π\| | 25\|N | E_floor | best E_pure | E/floor | g95/gmax | evals | seed waste | stop |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2×2×2 | 16 | 48 | 768 | no | 7.287e-1 | 3.3060 | 4.54× | 0/60 | 4,489 | 69% | cap |
| 2×2×3 | 24 | 16 | 384 | no | 3.736e-1 | 1.9811 | 5.30× | 16/60 | 7,844 | 69% | cap |
| 2×3×3 | 36 | 16 | 576 | no | 2.230e-1 | 1.4604 | 6.55× | 45/60 | 8,248 | 53% | cap |
| 3×3×3 | 54 | 48 | 2,592 | no | 1.411e-1 | 1.0808 | 7.66× | 42/60 | 9,268 | 78% | cap |
| 3×3×4 | 72 | 16 | 1,152 | no | 1.327e-1 | 0.9183 | 6.92× | 47/60 | 7,825 | 59% | cap |
| 3×4×4 | 96 | 16 | 1,536 | no | 8.363e-2 | 0.8885 | 10.63× | 37/60 | 7,752 | 59% | cap |
| 4×4×4 | 128 | 48 | 6,144 | no | 8.038e-2 | 0.7667 | 9.54× | 33/60 | 8,034 | 67% | cap |
| 4×4×5 | 160 | 16 | 2,560 | no | 5.708e-2 | 0.7520 | 13.18× | 46/60 | 7,073 | 50% | cap |
| 4×5×5 | 200 | 16 | 3,200 | **yes** | **0** | 0.5364 | — | 37/60 | 8,530 | 50% | cap |
| 5×5×5 | 250 | 48 | 12,000 | **yes** | **0** | 0.7786 | — | 55/60 | 6,286 | 50% | cap |
| 5×5×6 | 300 | 16 | 4,800 | **yes** | **0** | 0.6794 | — | 60/60 | 6,216 | 50% | cap |
| 5×6×6 | 360 | 16 | 5,760 | no | 2.537e-2 | 0.5826 | 22.97× | 51/60 | 6,491 | 50% | cap |
| 6×6×6 | 432 | 48 | 20,736 | no | 1.794e-2 | 0.5506 | 30.69× | 58/60 | 6,406 | 66% | cap |
| 6×6×7 | 504 | 16 | 8,064 | no | 1.557e-2 | 0.8410 | 54.01× | 56/60 | 6,100 | 70% | cap |
| 6×7×7 | 588 | 16 | 9,408 | no | 1.354e-2 | 0.7088 | 52.34× | 57/60 | 5,368 | 50% | cap |
| 7×7×7 | 686 | 48 | 32,928 | no | 1.217e-2 | 1.0006 | 82.24× | 54/60 | 5,501 | 73% | cap |

g95 = first generation whose pooled best is within 5% of the run's final best; "seed waste" =
P1-rejected fraction of generation-0 evaluations; every stop is the 60-generation cap.

**The floor law holds on all sixteen cells: E_floor = 0 ⟺ 25 | N.** Only 200, 250 and 300
sites divide the 25 sites-per-composition-period of the equiatomic quinary; perfect
stoichiometry at any other size leaves a positive, computable floor that decays roughly as 1/N
(0.729 at 16 sites → 0.0122 at 686). Below 54 sites the composition itself is not
representable within the 2% guard (§Vulnerabilities, item 1) — the equiatomic quinary
effectively does not exist at 2×2×2.

**Convergence at fixed budget degrades predictably with N.** The 60-generation budget that
brings 16–128 sites to their plateau (g95 well inside the cap, best-of-run reached by
generation ≈ 33–47) is nowhere near sufficient from 300 sites up: g95 ≈ gmax and the pooled
best is still descending when the cap fires, with the floor-relative error rising monotonically
to 82× at 686 sites. In absolute terms the picture inverts — best E_pure *falls* from 3.31 to
≈ 0.55–1.0 across the ladder — so budget comparisons must be floor-relative (or
random-normalized, see below), never absolute. At the small end the opposite regime appears:
2×2×2 finds its final best already in generation 0 and spends 60 generations re-discovering
archived duplicates — the unique non-P1 space at N = 16 is nearly exhausted by seeding alone.

**Seeding waste tracks the point group, extending the seven-cell law.** The mixed [D4] seeder
loses its rejection slots to the P1 filter at 50% on every two-equal-axes cell (|PG| = 16) and
at 66–78% on the cubic cells (|PG| = 48) — the constructive slots carry the population
everywhere, and the rejection halves are pure waste from 96 sites up.

**Cost.** Aggregate throughput across the 4-island runs spans 6,800 → 125 evaluations/s as
|Π|·N grows by 1,800× — damped relative to the strict 1/(|Π|·N) law because duplicate
rejections (which dominate small cells) skip the spglib call. The clean per-evaluation law
lives in the recorded T-B1 bench: evals/s × |Π| × N is constant within ×1.5 from 128 to 1,024
sites on one worker. 7×7×7 remains entirely practical: a full 60-generation, 4-island run is
under a minute of wall on a 32-core node.

## The ladder at γ = 1 (displacement-penalized)

| cell | D(P1)=6N | best D | D-reduction | E_pure(γ1) | SG | E(γ1)/E(γ0) |
|---|---|---|---|---|---|---|
| 2×2×2 | 96 | 30 | 3.2× | 4.587 | Cmm2 | 1.39 |
| 2×2×3 | 144 | 44 | 3.3× | 3.616 | Pmm2 | 1.83 |
| 2×3×3 | 216 | 108 | 2.0× | 1.764 | C2 | 1.21 |
| 3×3×3 | 324 | 162 | 2.0× | 1.129 | C2 | 1.04 |
| 3×3×4 | 432 | 216 | 2.0× | 1.037 | C2 | 1.13 |
| 3×4×4 | 576 | 147 | 3.9× | 1.175 | P4 | 1.32 |
| 4×4×4 | 768 | 128 | 6.0× | 1.383 | R-3 | 1.80 |
| 4×4×5 | 960 | 80 | 12.0× | 1.816 | P2₁/m | 2.41 |
| 4×5×5 | 1,200 | **10** | **120×** | **7.939** | **Cmcm** | 14.80 |
| 5×5×5 | 1,500 | 150 | 10.0× | 1.046 | Pc | 1.34 |
| 5×5×6 | 1,800 | **10** | **180×** | **7.939** | **Cmcm** | 11.68 |
| 5×6×6 | 2,160 | 120 | 18.0× | 1.601 | Pmc2₁ | 2.75 |
| 6×6×6 | 2,592 | 432 | 6.0× | 0.822 | R-3 | 1.49 |
| 6×6×7 | 3,024 | 763 | 4.0× | 0.810 | P4 | 0.96 |
| 6×7×7 | 3,528 | 883 | 4.0× | 0.743 | P-4 | 1.05 |
| 7×7×7 | 4,116 | 1,036 | 4.0× | 0.887 | P4 | 0.89 |

**The ordered-compound collapse reproduces exactly where the seven-cell study found it — and
nowhere else.** The two cells with a 5,5-axis pair (4×5×5 and 5×5×6) both fall onto the *same*
Cmcm-layered compound at D = 10 and E_pure = 7.938889 — the identical value at 200 and 300
sites, an order of magnitude above their γ = 0 bests (14.8×/11.7×). This is an anti-SQS: a
crystalline decoration winning an SQS objective purely through its D collapse, the clearest
argument for normalizing the objective (§Vulnerabilities, item 3). The other equal-axis cells
(2,2 / 3,3 / 4,4 / 6,6 / 7,7 pairs) do *not* collapse at this budget — the compatible layering
period appears specific to the 5-axis pair at x = 1/5, i.e. it is the commensuration of the
composition with the axis, not the axis equality itself, that opens the trapdoor.

**Away from the trapdoor, γ = 1 is productive at every size.** Typical reductions are 2–6×
(3.2× at 16 sites; 6× at 128 and 432 — both landing on R-3 cells, the K = 5 analogue of the
binary knee), 10–18× on the mid-size non-collapsing cells, and a striking regime appears at
504+ sites: the γ = 1 runs reach D-reductions of 4× at *equal or better* correlation error
than their γ = 0 twins (ratios 0.96, 1.05, 0.89) — at large N the D-penalty acts as a search
regularizer rather than a trade-off, steering the walk through symmetric strata where the
correlation optimum apparently also lives.

## Comparison with the literature

**Objective.** Classic SQS construction, from the original 8–16-atom fcc structures of Zunger,
Wei, Ferreira and Bernard [PRL **65**, 353 (1990)] to the stochastic `mcsqs` of van de Walle et
al. [Calphad **42**, 13 (2013)], minimizes the mismatch of cluster correlation functions against
the ideal random alloy — the analogue of our `E_pure` (γ = 0), usually extended beyond pairs to
triplets. `sqsgenerator` [Gehringer, Friák and Holec, Comput. Phys. Commun. **286**, 108664
(2023)] casts the same target in Warren–Cowley SRO parameters with highly parallel exhaustive or
stochastic search, and `icet` [Ångqvist et al., Adv. Theory Simul. **2**, 1900015 (2019)]
provides enumeration- and MC-based SQS in a cluster-expansion framework. None of these treats
the *downstream phonon cost* as part of the objective: a converged mcsqs cell is generically P1,
so a phonopy workflow on N sites pays the full D = 6N displaced supercells. The E·D^γ objective
of the reference paper (arXiv:2602.10872), implemented here, is — to our knowledge — unique in
buying that cost down explicitly; the γ = 1 column above shows 6N cut to D values one to two
orders of magnitude smaller at quantified correlation cost.

**Deduplication.** Our canonical-label archive over the empty-cell permutation group [A9][A10]
plays the same role as the symmetry reduction in the `supercell` program [Okhotnikov,
Charpentier and Cadars, J. Cheminform. **8**, 17 (2016)] — there used to enumerate distinct
substitutional configurations, here used online to keep the evolutionary search from re-paying
spglib evaluations on symmetry-equivalent decorations.

**Scale conventions.** Absolute correlation-error scalars are not comparable across codes (each
normalizes differently — see the SPEC §4.1 scale caveat and the floor analysis in
STEP2_REPORT); comparisons belong on raw structures, which is what `exsqs score` +
`tools/py/align_to_config.py` exist for. Within this constraint, the qualitative picture
matches the field: at 16–54 sites a 5-component random alloy is barely representable (the
literature's small SQS are 2–4-component for this reason), and from ~128 sites upward the
equiatomic quinary becomes well-conditioned — consistent with the 100–250-site cells typically
used for refractory MoNbTaVW-class HEAs with mcsqs.

**Determinism.** We are not aware of another SQS code offering bit-identical trajectories
across thread counts, MPI rank counts and checkpoint/resume chains; for scaling studies like
this one that property converts every anomaly into a reproducible test case (the resume
foot-gun below was found, reproduced and regression-tested in a day).

## Vulnerabilities and problematic parts, with proposed solutions

1. **Representability guard vs small cells.** Equiatomic x = 0.2 violates the fixed 2%
   `comp_tol` below 54 sites (N = 16/24/36 round to 4/3/3/3/3, 5/5/5/5/4, 8/7/7/7/7), and the
   loud-fail guard rejects the config outright. We pinned the nearest representable fractions
   explicitly in those configs. *Proposal:* make `comp_tol` a config key (loud default 0.02),
   and report the rounded counts in the error message so the fix is copy-pasteable.
2. **`e_tol: auto` degenerates on both ends.** On zero-floor cells (25 | N) auto = 0 is
   unreachable; on small cells the floor is so large (0.73 at N = 16) that auto = 3×floor is
   meaningless. *Proposal:* normalize the target by the analytic random-alloy expectation
   E[random] (hypergeometric, SPEC §4 note) — a scale-free `e_tol: fraction_of_random` that is
   well-defined on every cell, zero-floor or not.
3. **γ = 1 rewards anti-SQS collapse.** On cells with a compatible axis pair the run walks into
   an *ordered* low-D compound whose correlation error is an order of magnitude above the
   random baseline — the objective E·D^γ is unnormalized, so a D collapse pays for any E.
   *Proposal (top-ranked):* divide E_pure by the analytic E[random] in the objective, or keep a
   Pareto archive of (E_pure, D) and select post hoc; either prevents the degenerate corner
   from dominating.
4. **Fixed generation budgets do not scale with N.** Ratio survival never stagnated here (all
   stops are `max_generations`), and on the large cells the best error is still descending at
   the cap (g95 ≈ gmax from 5×5×6 upward): at a fixed 60-generation budget the relative error
   degrades monotonically with N, from 4.5× floor at 16 sites to 82× at 686. The opposite
   failure belongs to exploitation regimes: at fixed high β the [A13] stagnation stop retires
   every island by generation ~35 (1.4.0 chain demo). *Proposal:* budget generations ∝ N (or
   chain wall-capped segments with `scripts/chain_resume.sh`), and report
   generations-to-plateau (g95) rather than the raw cap as the convergence metric.
5. **Canonicalization cost at large N.** The per-evaluation cost is dominated by spglib and the
   |Π|·N canonical-label scan; at N = 686 (|Π| = 32,928) throughput drops to a few evaluations
   per second per worker. *Proposal:* the [A10] invariant pre-filter (composition + per-shell
   pair vector) before the full lex-min scan, and a 128-bit label hash in the archive instead
   of full labels (the spec explicitly allows it) — memory at N = 686 is ~0.7 kB per admission
   otherwise.
6. **Seeding waste tracks the point group.** Mixed seeding [D4] discards its even (rejection)
   slots at large N with probability ≈ 1 (all P1) — measured above at ~50% of generation-0
   evaluations across the ladder. *Proposal:* switch `seeding.mode` default to `constructive`
   above a size threshold, or demote rejection slots to a fixed small fraction.
7. **The `--resume` output foot-gun (fixed here as 1.7.1).** `--resume` without `--out`
   retargeted outputs and subsequent state saves to the config's `output.dir`, overwriting a
   sibling run (observed at 6×6×6, reproduced deterministically, regression noted in SPEC §16).
   Resume now defaults `output.dir` to the state directory.

## Reproducing

    cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
    for c in configs/hea5_bcc_*.yaml; do ./build/exsqs "$c"; done              # gamma = 0
    for c in configs/hea5_bcc_*.yaml; do ./build/exsqs "$c" \
        --set error.gamma=1 --out runs/$(basename "$c" .yaml)_g1; done         # gamma = 1
    python3 tools/py/validate.py runs/hea777_g0                               # T-V1 spot check

Every table value re-derives from `runs/*/summary.json` and `runs/*/generations.csv`.
