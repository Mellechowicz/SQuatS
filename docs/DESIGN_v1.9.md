# Design v1.9 — multiplet correlations and restricted cell-shape search

Status: PLANNING (feat/v1.9-design, 2026-07-17). Two missing utilities
identified by the CPC-manuscript survey: (i) triplet and quadruplet
correlation terms in the objective, (ii) a supercell-shape search under
user restrictions. Structure of this document: theory -> algorithms ->
development & deployment.

## Part I — theoretical background

### I.1 Multiplet correlations (triplets, quadruplets)

The v1.8 objective compares the conditional zone occupancies
Pi(t->t'; n) = C_tt'^(n) / (c_t z_n) with the random-alloy target x_t'
(SPEC section 4). Pair statistics do not pin down three- and four-site
short-range order; ATAT and icet match multiplet correlations for this
reason, and reviewers will ask.

**Cluster classes.** Enumerate site triples (and quadruples) of the
undecorated lattice with all pair separations inside a cutoff, and
classify them into orbit classes c under the empty-lattice group Pi (the
same site_permutations the dedup archive uses). A class is labelled by
its sorted distance multiset; the orbit classification makes the count
tables symmetry-exact rather than distance-binned. Practical default
(ATAT-like): triplets from the first two zones, quadruplets from the
first zone only.

**Objective extension.** For class c with per-site multiplicity m_c,
the joint occupancy Pi_c(t1..tk) = C_c(t1..tk) / (N m_c) is compared
with the ideal x_t1 ... x_tk times the tuple-degeneracy factor; the L1
form is kept (no cancellation), giving

    E = E_pair + lambda_3 sum_c A_c sum_tuple |Pi_c - ideal|  (+ lambda_4 ...)

with per-class weights A_c decreasing with the class diameter, and
lambda_3, lambda_4 user knobs (0 by default: v1.8 behaviour and
bit-compatibility preserved).

**Quantization floor.** Counts are integers, so every term is bounded
below by dist(ideal * N m_c, Z) / (N m_c). Summing gives an analytic
lower bound on the multiplet sector. Unlike the pair sector this bound
need not be attained (pair/triplet/quad counts are coupled through the
decoration); it is reported as a bound, and the exact pair floor keeps
its own line in the summary.

**Symmetry interaction.** Multiplet counts are invariant under
Pi-relabelling, so the canonical-label archive and the
symmetry-preserving mutation [D6] carry over untouched. The irreducible-
problem measure D is unaffected by definition.

**Cost model.** Pair evaluation is O(N zbar); triplets are O(N zbar^2)
and quadruplets O(N zbar^3) with zbar the within-cutoff coordination
(8-12) — roughly 10x and 100x the pair cost if evaluated naively. See
II.1 for the mitigation ladder. The v1.8 throughput law
(evals/s x |Pi| x N ~ const) will not survive as-is; a new law is to be
measured (T-B1 extension).

### I.2 Restricted cell-shape search

The engine decorates a FIXED user cell; mcsqs's central advantage is
searching cell shapes. The search space is the superlattice set
{H in Z^(3x3), |det H| = M} enumerated by Hermite Normal Forms and
reduced by the parent point group (Hart & Forcade, PRB 77, 224115) —
the same enumeration the dedup theory already leans on.

**Restrictions** (user-facing; the requested "fit better within some
restrictions"):
- size window: N in [N_min, N_max];
- commensuration: exact integer counts, and optionally the pair-floor-0
  test (x_t c_t z_n integer for all t, n) evaluated with the analytic
  floor formula — SQuatS can RANK CELLS ANALYTICALLY before any search,
  a lever no competitor has;
- geometry: Minkowski-reduced basis, minimal slab width
  w_min >= 2 R_shell(n_s) (eliminates the [A3] half-width warning),
  aspect-ratio cap;
- symmetry headroom: larger |Pi(H)| means richer retained-symmetry
  strata and cheaper D plateaus — used as a ranking term.

**Scoring.** Primary ranking is analytic and instantaneous:
(E_floor, -|Pi|, skewness) lexicographically or weighted. An optional
budgeted tournament then runs short seeded searches on the top-k cells
and picks by E_obj; each mini-run is bit-reproducible, so the whole
selection is.

## Part II — workflow (algorithms)

### II.1 Multiplets

- A-M1 cluster tables: new clusters.{hpp,cpp}; RunContext gains a
  ClusterTable built once (pair zones reused; triples/quads enumerated
  per site, orbit-classified under Pi; per-site tuple lists stored).
- A-M2 counting: CorrData extension holding C_c(t1..tk); read-only
  tables, OpenMP-safe per-candidate loops.
- A-M3 objective: e_pure gains the lambda_3/lambda_4 sectors; config
  block error.clusters {triplet_shells, quadruplet_shells, lambda3,
  lambda4, weights}; all defaults preserve v1.8 trajectories bitwise.
- A-M4 floor: per-class bound added to e_floor reporting (bound, not
  exact — labelled as such in summary.json and the banner).
- A-M5 parity: score subcommand, summary.json fields, and — critical —
  trajectory_signature extended with every new knob (resume guard).
- Optimization milestone (M2, after the MVP): delta-evaluation on swap
  mutations — a swap touches two sites, so only their tuple lists need
  recounting: O(zbar^2) per swap instead of O(N zbar^2) per candidate.

### II.2 Cell search

- A-C1 `exsqs cells` subcommand: parent lattice + composition +
  restrictions in, ranked cell table out (csv/json + ready config
  snippets). HNF enumeration, point-group reduction, Minkowski
  reduction, filters, analytic floors.
- A-C2 tournament driver (`--tournament <budget>`): deterministic
  seeded mini-runs over the top-k cells, winner by E_obj.
- A-C3 config integration: supercell: {search: {sites: [min,max],
  commensurate: true, min_width_shells: 5, max_aspect: 2.0, top: 10}}
  resolves to a concrete H BEFORE RunContext::build — the engine core
  and the reproducibility contract stay untouched.

### II.3 Test workflow

- T-M1 golden counts: brute-force multiplet counting vs ClusterTable on
  tiny cells (<= 16 sites), all lattices.
- T-M2 floor bound: exhaustive decoration enumeration on <= 16 sites
  confirms bound <= true minimum, and equality where attainable.
- T-M3 interop: triplet correlations cross-checked against icet on
  shared structures (tools/py, T-V1 pattern).
- T-C1' HNF census: enumeration counts vs the published Hart-Forcade
  sequences; T-C2' floor-0 detection reproduces the 25|N law of the
  supercell study.
- Determinism: T-K1/T-MPI1 re-run with multiplets enabled; signature
  guard tested against every new knob.

## Part III — development and deployment

Implementation follows the repository workflow (feature branches,
>= 13 commits for significant branches, every commit configures/builds/
runs, per-file bullets for 3+-file commits):

1. **feat/multiplet-correlations** — commit sketch: cluster enumeration;
   orbit classification; tuple storage; counting kernel; objective
   sectors; config keys + validation; floor bound; score parity;
   summary fields; signature extension; T-M1; T-M2; T-M3; SPEC section
   4/16 updates; T-B1 cost-law re-measurement.
2. **feat/cell-search** — commit sketch: HNF enumerator; point-group
   reduction; Minkowski + geometry filters; analytic-floor ranking;
   `cells` subcommand; tournament driver; config integration; T-C1';
   T-C2'; docs; ladder demo (does the search rediscover 5x5x5 on the
   W-Cr family, or beat it?).

Release ritual: SPEC v1.9 header + changelog entry, the four version
literals to 1.9.0, coherence C1..C9, full 24+ gate suite, merge --no-ff,
tag v1.9.0, GitHub release; LUMI/Helios rebuild note (see the T3
post-mortem: stale binaries abort under upgraded stacks). Packaging
groundwork for conda-forge follows the survey's criticism list.

Risks: (i) evaluation-cost regression — mitigate with shell-restricted
defaults and the M2 delta-evaluation; (ii) multiplet floor is a bound,
not a target — communicate clearly; (iii) HNF growth with N — cap the
window and rely on point-group reduction; (iv) config-surface bloat —
defaults must reproduce v1.8 bitwise (guarded by T-K1).

Milestones: M1 multiplet MVP -> M2 delta-eval -> C1 cells MVP ->
C2 tournament -> R: release 1.9.0 with re-run benchmark subset and a
manuscript addendum if the results warrant one.
