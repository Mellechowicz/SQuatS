# EXSQS — Step 5 Report: Validation Campaign, K ≥ 3, Unified Verification (v1.4.0, SPEC v1.4)

Date 2026-05-18. Scope: the remaining plan items — validating the multicomponent (K ≥ 3) path
that the spec always described but no test had ever exercised, running the reference-scale
validation campaigns, and consolidation: a unified verification runner, a CI stub, and
documentation. The complete matrix for every step to date was prepared and executed; the full
log is regenerated into `docs/TEST_MATRIX.txt` at each release.

Verification summary: **all matrix gates green** — the core groups
(geometry/zones/correlations, symmetry/dedup, the T-D1 phonopy gate, the counter RNG), the
engine groups (config, engine units, floor, T-R1 bitwise reruns, T-V1 python cross-validation at
|diff| = 0, T-B1 recorded), the parallel group (thread safety, thread invariance, migration,
scaling recorded), the HPC group (serializer, T-K1/T-K2 resume gates, T-MPI1 rank invariance),
and the new step-5 group (K = 3 fast gates, T-E3 integration, T-V1 on both campaign
directories). Catch2 totals: **32 fast cases / 27,073 assertions** + 3 e2e cases
(T-E1/T-E2/T-E3). Run it with `bash tools/run_all_tests.sh` (`SKIP_E2E=1` for the fast pass).

## 1. K ≥ 3 validated end-to-end — with a genuine discovery

Everything before this step was binary W-Cr. The machinery was K-general on paper — [A16]'s
`mode: auto` resolving to `full_pairs` at K ≥ 3, the §4.1 floor summing over species, the
per-species constructive DP, species-agnostic orbit mutation — and it all held on first contact:
the new bcc reference `w50mo25cr25_4x4x4.yaml` (W₆₄Mo₃₂Cr₃₂, 128 sites) runs end to end, and
the campaign below descends into the same monoclinic/triclinic classes the binary system
produces. `validate.py` reproduces the full-pairs K = 3 error at |diff| = 0, closing T-V1 for
multicomponent runs.

The discovery: for this composition **E_floor = 0 exactly**. The fractions (½, ¼, ¼) are dyadic
and every `x̃·P_t(n)` on bcc 4×4×4 is integral, so every quantization term vanishes — a perfect
SQS is not excluded by counting. The floor code handles it correctly (the K = 2 `full = 2×diag`
identity is also correctly *absent* at K = 3, pinned by test), but it means floor-relative gates
degenerate there: `e_tol: auto` becomes 0 and users must set a numeric tolerance. SPEC §4.1
carries the note (v1.4); T-E3 asserts the zero floor as a regression alongside its absolute
E_pure ≤ 1.2e-1 gate.

## 2. Campaigns — results and an honest limit

The binary flagship — 6 islands, β = 3000, migration every 8 — reconfirms the project record
robustly: **rc = 0 SUCCESS at 2.50× E_floor** (8.111947e-3, Cm) in 19 s wall on the dev
workstation and 13,363 evaluations, with 22 accepted migration events across the ledger
`[6, 6, 2, 2, 2, 4]`; T-V1 validates the outputs at |diff| = 0. The exploitation-regime chain
(β = 3000 with `e_tol` pinned unreachably low, 40-generation segment then `--resume` with the
cap raised to 140) demonstrates the other side: **every island is stagnation-terminated by
generation 20–35** — the [A13] all-fallback arithmetic from the Step-2 study, observed in the
production chaining path — and the second segment correctly refuses to re-arm the *terminal*
stops, ending in seconds with the record already at 2.50×. The conclusion for real campaigns is
worth stating plainly: exploitation regimes are exhaustion-limited by design; long-budget
frontiers belong to the non-exhausting survival modes (ratio, auto-β) at larger populations,
which is precisely what the SLURM templates and bit-exact chains from Steps 3–4 are for.

The ternary campaign (6 islands, ratio survival, γ = 0, 60 generations, migration every 8)
reaches **E_pure = 1.330e-1 (Cm)** at generation 45 with 68 accepted migrations over 18,610
evaluations — a budget-limited number by construction, since the zero floor removes the success
criterion; T-V1 validates it at |diff| = 0. The campaign directories
(`runs/campaign_w70cr30`, `runs/campaign_ternary`) double as validator regression inputs for the
runner's `T-V1 on campaign` rows.

## 3. Release consolidation

`tools/run_all_tests.sh` is the single entry point for the whole verification story: it builds,
runs every tag group mapped to its step, the e2e suite, T-MPI1 (auto-detecting the launcher),
T-V1 on a fresh run plus any campaign directories present, and the recorded benchmarks, then
prints the PASS/FAIL matrix with timings and exits nonzero on any failure.
`.github/workflows/ci.yml` runs the fast suite + T-MPI1 on pushes and the e2e suite as a second
job (the pymatgen-based T-V1 stays in the local runner by design). README gained a consolidated
testing section; versions report 1.4.0; SPEC is at v1.4 with the K ≥ 3 items, the
commensurate-floor note, T-E3 in §12, and the campaign record in the changelog.

Benchmark note: T-B1 numbers move with machine load (176.9 evals/s at N = 128 under the loaded
runner vs 394 idle); the relative phase structure — spglib dominant, canonicalization second —
is stable and is what the scaling analysis rests on.

## 4. Project status

All five originally planned steps are delivered. The library is gate-exact against phonopy; the
error model carries a provable quantization floor that exposed the paper's unstated
normalization; the engine reproduces and quantifies the search dynamics, beating the paper's
displacement-reduction claim on the D axis (6× at the knee, 19× at the Pareto extreme);
execution is bit-invariant across threads, ranks, and interruptions; and the whole matrix is one
command. What remains beyond the plan: external-structure scoring for cross-code comparisons
(the §4.1 doctrine still lacks its tool), the deferred [A11] geometric β schedule, and the
K ≥ 3 physics campaigns at cluster scale — the SLURM templates in `scripts/slurm/` are ready.
