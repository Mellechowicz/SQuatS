# EXSQS — Step 7 Report: Coherence Audit (v1.6.0, SPEC v1.6)

Date 2026-06-25. Step 7 was not in the original plan (which ended at 5); it grew out of the
project's own trajectory: after six increments of specification, engine, parallel, HPC,
validation and interoperability work, the credible failure mode is no longer a red test but
*drift between artifacts* — spec vs code vs tests vs configs vs docs vs version strings. Step 7
mechanizes that concern as **T-CO1** (`tools/check_coherence.py`), makes it the matrix's
long-missing "step-0 test", fixes what its first runs found, and re-verifies everything.

**Verification (steps 0–6): 23 matrix gates, all green** — now headed by
`S0 spec/code/tests/docs coherence (T-CO1): PASS, 17 checks`, followed by every earlier gate:
the phonopy T-D1 gate, 27k+ fast assertions across 37 cases, T-E1/T-E2/T-E3 integration, T-MPI1
rank invariance, T-V1 at |diff| = 0 on fresh and campaign runs, and the recorded benchmarks.

## 1. What T-CO1 enforces

Seventeen mechanical cross-checks, run live against the built tree: **one version literal**
across both drivers, equal to the CMake VERSION and verified against a *live* `summary.json`
from an actual run; **every §12 test id cited** by a real test or tool source; **every Catch2
tag exercised** by the runner (no orphan coverage); **no version claims in user-facing strings**;
the **trajectory-signature field ledger** — every `RunConfig` field must be either written into
the signature or on the documented exclusion list (budgets, output/logging, `omp_threads` per
[A14], `x_target` per [A5]), so a future field added without classification fails the audit
rather than silently corrupting resume semantics; the **§11 sample config and every
`configs/*.yaml` executed end-to-end** (not merely parsed — 4 seeds, 1 generation, summary
written); SPEC header version equal to the latest §16 changelog entry with sections 1–16
present; and every README-referenced path and binary resolving. The enforced policy: step
reports are immutable dated snapshots; SPEC, README, configs and code must stay mutually
current.

## 2. What the audit found

Two genuine catches, both fixed before the audit was wired into the matrix. First, **T-A11 was
defined in SPEC §12 but cited nowhere** — the schedule tests referenced decision [A11] but not
the test id; the ids now appear in the test names, restoring the spec→test traceability that
holds for the other ids. Second — and worth recording — the audit's first run caught a bug *in
the auditor itself*: the struct-field extractor stopped at the first `};`, which any
brace-initialized member produces, reporting 2 fields instead of 40; fixed with brace counting.
A coherence tool that finds its own incoherence on first contact is doing its job.

Judgment beyond the mechanical checks: the historical reports contain values that have since
moved (benchmark throughput varies with machine load; the e2e wall swings by tens of seconds
run-to-run) — these are dated measurements, not claims about the present, and the
immutable-snapshot policy is the correct treatment. No contradiction between any report's
*claims* and current behavior was found; the SPEC §16 changelog remains the single living
history and now provably matches the header, the code version, and a live run.

## 3. Status

The software ledger was empty after step 6; after step 7 the *meta*-ledger is covered too — the
project now tests its own consistency as gate zero. Versions report 1.6.0; SPEC is v1.6.
Everything previously true remains verified: one trajectory across threads, ranks, and resume
chains; floor-relative validation with the K = 3 extension; external-structure scoring for the
paper comparison. What remains for the release step: closing the scorer's deliberate
orientation-strictness for real rotated/permuted external files, and packaging.
