# Changelog

Software versions below; the authoritative, fully detailed history (including
every decision id [A*] and test id T-*) is `docs/SPEC.md` section 16, whose
spec versions are noted in parentheses.

## 1.7.1 (2026-07-10) — supercell-study fix
- `--resume` without `--out` silently retargeted outputs and subsequent state
  saves to the config's `output.dir` (observed overwriting a sibling run at
  6x6x6 in the K=5 study). Resume now defaults `output.dir` to the state
  directory and says so on stderr; an explicit `--out` overrides as before.

## 1.7.0 (2026-07-03) — release (SPEC v1.7)
- `exsqs geom` subcommand: writes the undecorated supercell frame of a config.
- `tools/py/align_to_config.py`: maps rotated/permuted/relabelled external
  structures (paper supplements, ATAT output) onto the config frame so that
  `exsqs score` accepts them; gated by T-X2, a cross-tool file-level round trip
  that reproduces engine records bitwise.
- Release packaging: CITATION.cff, CHANGELOG.md, restructured README (accreted
  step notes moved to `docs/DEV_NOTES.md`), the recorded verification matrix in
  `docs/TEST_MATRIX.txt`, and the `v1.7.0` tag.

## 1.6.0 (2026-06-26) — coherence audit (SPEC v1.6)
- T-CO1 `tools/check_coherence.py`: 17 mechanical spec/code/tests/configs/docs
  cross-checks, wired into the matrix as its "step-0 test". Its first runs found
  and fixed real drift (uncited T-A11) and a bug in the auditor's own field
  extractor (fixed with brace counting).

## 1.5.0 (2026-06-12) — interoperability + spec completion (SPEC v1.5)
- `exsqs score`: evaluate any POSCAR on a config's geometry (E_pure, E/E_floor,
  D, D(P1)/D, SG, E_obj), site-order-free, loud mismatch errors; T-X1 pins
  bit-exact agreement with engine records.
- [A11] geometric beta schedule (`schedule: geometric`, `beta_growth`); the two
  fields have been in the trajectory signature since 1.3.0, so enabling the
  schedule invalidates earlier state files under `--resume` (by design).
- Non-diagonal `supercell.matrix` test coverage; `tools/py/to_poscar.py`.

## 1.4.0 (2026-05-20) — K >= 3 + campaigns + unified runner (SPEC v1.4)
- Ternary support validated end to end (W64Mo32Cr32, `full_pairs` [A16]);
  exact-zero floor for commensurate compositions handled and documented.
- `tools/run_all_tests.sh` verification matrix; CI workflow.

## 1.3.0 (2026-05-01) — HPC: checkpoint/resume + MPI + SLURM (SPEC v1.3)
- Bit-exact checkpoint/resume (`--resume`, signature-guarded, raisable budget
  caps); `exsqs_mpi` with rank-count invariance (T-MPI1: serial == n1 == n3);
  SLURM templates and `scripts/chain_resume.sh` (exit 0/3 contract).

## 1.2.0 (2026-04-03) — lockstep islands + OpenMP (SPEC v1.2)
- Round-based generations, ring migration; thread-count invariance proven
  bitwise against the 1.1.0 sequential artifacts.

## 1.1.0 (2026-03-12) — engine completion (SPEC v1.1)
- Constructive space-group seeding [D4] behind the P1 rejection filter [A7][A8],
  symmetry-preserving [D6] mutation with Poisson swap counts [A12], stagnation
  stop [A13]; floor-relative acceptance (`e_tol: auto` = 3.0 x E_floor,
  SPEC 4.1) promoted to the acceptance criterion.

## 1.0.0 (2026-02-13) — the minimal release
- Sequential extinction engine over the core library: coordination zones, pair
  correlations, spglib symmetry + phonopy-convention displacement counts,
  canonical-label deduplication, counter-keyed RNG [A14], quantization floor
  `E_floor`, YAML config with CLI overrides, POSCAR + summary outputs.
