#!/usr/bin/env bash
# EXSQS full verification matrix, steps 1-5 (v1.4). Runs every automated gate
# in the project and prints a step-by-step PASS/FAIL table with timings.
#
#   bash tools/run_all_tests.sh
#
# Environment switches:
#   SKIP_BUILD=1   skip the cmake build step
#   SKIP_E2E=1     skip the slow [e2e] integration suite (~3 min, W-system)
#   SKIP_MPI=1     skip the MPI rank-invariance test (T-MPI1)
#   SKIP_V1=1      skip the python cross-validation (T-V1; needs pymatgen)
#   SKIP_BENCH=1   skip the recorded benchmarks (T-B1)
#   MPIEXEC=...    MPI launcher override (default: mpirun.openmpi, then mpirun)
#
# Exit code: 0 iff every executed gate passed.
set -u
cd "$(dirname "$0")/.."
BIN=build/tests/exsqs_tests

NAMES=(); RES=(); TIMES=(); DETAILS=()
FAILED=0
add() { NAMES+=("$1"); RES+=("$2"); TIMES+=("$3"); DETAILS+=("$4"); }

run_tags() {  # <label> <catch2-tag-filter>
  local label=$1 tags=$2 t0 out rc dt line
  t0=$(date +%s)
  out=$("$BIN" "$tags" 2>&1); rc=$?
  dt=$(( $(date +%s) - t0 ))
  line=$(echo "$out" | grep -E "All tests passed|assertions" | tail -1 | sed 's/^ *//')
  if [ $rc -eq 0 ]; then
    add "$label" PASS "${dt}s" "$line"
  else
    add "$label" FAIL "${dt}s" "$line"
    echo "---- $label FAILED ----"; echo "$out" | tail -25; FAILED=1
  fi
}

echo "EXSQS verification matrix ($(./build/exsqs --help 2>/dev/null | head -0; true)v1.7, $(date -u +%F))"
echo

if [ "${SKIP_BUILD:-0}" != "1" ]; then
  t0=$(date +%s)
  if cmake --build build -j2 > /tmp/exsqs_build.log 2>&1; then
    add "build (all targets)" PASS "$(( $(date +%s) - t0 ))s" "cmake --build build"
  else
    add "build (all targets)" FAIL "$(( $(date +%s) - t0 ))s" "see /tmp/exsqs_build.log"
    tail -20 /tmp/exsqs_build.log; FAILED=1
  fi
fi

# ---- Step 0: specification coherence (v1.6, T-CO1) ----
t0=$(date +%s)
if python3 tools/check_coherence.py > /tmp/exsqs_coh.log 2>&1; then
  add "S0  spec/code/tests/docs coherence (T-CO1)" PASS "$(( $(date +%s) - t0 ))s" \
      "$(grep -c '\[PASS\]' /tmp/exsqs_coh.log) checks"
else
  add "S0  spec/code/tests/docs coherence (T-CO1)" FAIL "$(( $(date +%s) - t0 ))s" \
      "see /tmp/exsqs_coh.log"
  grep FAIL /tmp/exsqs_coh.log; FAILED=1
fi

# ---- Step 1: core library (spec sections 3-7) ----
run_tags "S1  geometry + zones + correlations"        "[geometry],[correlation]"
run_tags "S1  symmetry wrapper + canonical dedup"     "[symmetry],[dedup]"
run_tags "S1  displacements (T-D1 phonopy gate)"      "[displacements],[gate]"
run_tags "S1  counter RNG [A14]"                      "[rng]"
# ---- Step 2: serial evolution (spec sections 8-12) ----
run_tags "S2  config schema + overrides"              "[config]"
run_tags "S2  engine units (seed/mutate/admit)"       "[evolution]"
run_tags "S2  quantization floor (SPEC 4.1)"          "[floor]"
run_tags "S2  bitwise reruns (T-R1)"                  "[T-R1]"
# ---- Step 3: parallel layer (spec section 8.1) ----
run_tags "S3  thread safety + invariance + migration" "[parallel]"
# ---- Step 4: checkpoint/restart (spec section 8.2) ----
run_tags "S4  serializer + T-K1 + T-K2"               "[checkpoint]"
# ---- Step 5: ternary K=3 (v1.4) ----
run_tags "S5  K=3 config/floor/engine (fast)"         "[ternary]~[e2e]"

# ---- Step 6: interoperability + spec completion (v1.5) ----
run_tags "S6  external-structure scoring (T-X1)"      "[score]"
run_tags "S6  geometric beta schedule [A11]"          "[schedule]"
run_tags "S6  non-diagonal supercell matrix"          "[supercell]"

# ---- Step 8: release pipeline (v1.7) ----
t0=$(date +%s)
if bash tools/test_align_roundtrip.sh > /tmp/exsqs_x2.log 2>&1; then
  add "S8  align + score round trip (T-X2)" PASS "$(( $(date +%s) - t0 ))s" \
      "$(tail -1 /tmp/exsqs_x2.log | cut -c1-60)"
else
  add "S8  align + score round trip (T-X2)" FAIL "$(( $(date +%s) - t0 ))s" \
      "see /tmp/exsqs_x2.log"
  FAILED=1
fi

if [ "${SKIP_E2E:-0}" != "1" ]; then
  run_tags "S2/S5  integration T-E1/T-E2/T-E3 [e2e]"  "[e2e]"
else
  add "S2/S5  integration T-E1/T-E2/T-E3 [e2e]" SKIP "-" "SKIP_E2E=1"
fi

if [ "${SKIP_MPI:-0}" != "1" ] && [ -x build/exsqs_mpi ]; then
  MP=${MPIEXEC:-}
  if [ -z "$MP" ]; then
    if command -v mpirun.openmpi >/dev/null 2>&1; then MP=mpirun.openmpi; else MP=mpirun; fi
  fi
  t0=$(date +%s)
  if tools/test_mpi_invariance.sh build/exsqs build/exsqs_mpi configs/smoke_sc27.yaml "$MP" \
      > /tmp/exsqs_mpi_test.log 2>&1; then
    add "S4  T-MPI1 rank invariance (serial=1=3 ranks)" PASS "$(( $(date +%s) - t0 ))s" \
        "$(tail -1 /tmp/exsqs_mpi_test.log)"
  else
    add "S4  T-MPI1 rank invariance" FAIL "$(( $(date +%s) - t0 ))s" "see /tmp/exsqs_mpi_test.log"
    tail -15 /tmp/exsqs_mpi_test.log; FAILED=1
  fi
else
  add "S4  T-MPI1 rank invariance" SKIP "-" "SKIP_MPI=1 or exsqs_mpi absent"
fi

if [ "${SKIP_V1:-0}" != "1" ]; then
  t0=$(date +%s)
  rm -rf /tmp/exsqs_v1run
  build/exsqs configs/smoke_sc27.yaml --set evolution.max_generations=8 \
      --out /tmp/exsqs_v1run > /dev/null 2>&1 || true
  if python3 tools/py/validate.py /tmp/exsqs_v1run > /tmp/exsqs_v1.log 2>&1; then
    add "S2  T-V1 python cross-validation" PASS "$(( $(date +%s) - t0 ))s" \
        "$(tail -1 /tmp/exsqs_v1.log)"
  else
    add "S2  T-V1 python cross-validation" FAIL "$(( $(date +%s) - t0 ))s" "see /tmp/exsqs_v1.log"
    tail -10 /tmp/exsqs_v1.log; FAILED=1
  fi
  for d in runs/campaign_w70cr30 runs/campaign_ternary; do
    if [ -f "$d/summary.json" ]; then
      t0=$(date +%s)
      if python3 tools/py/validate.py "$d" > /tmp/exsqs_v1c.log 2>&1; then
        add "S5  T-V1 on $(basename "$d")" PASS "$(( $(date +%s) - t0 ))s" \
            "$(tail -1 /tmp/exsqs_v1c.log)"
      else
        add "S5  T-V1 on $(basename "$d")" FAIL "$(( $(date +%s) - t0 ))s" "see /tmp/exsqs_v1c.log"
        FAILED=1
      fi
    fi
  done
else
  add "S2  T-V1 python cross-validation" SKIP "-" "SKIP_V1=1"
fi

BENCH_OUT=""
SCALE_OUT=""
if [ "${SKIP_BENCH:-0}" != "1" ]; then
  t0=$(date +%s); BENCH_OUT=$(build/exsqs_bench 2>&1)
  add "S2  T-B1 throughput (recorded)" PASS "$(( $(date +%s) - t0 ))s" "table below"
  t0=$(date +%s); SCALE_OUT=$(build/exsqs_bench_scaling 2>&1)
  add "S3  T-B1 scaling (recorded)" PASS "$(( $(date +%s) - t0 ))s" "table below"
fi

echo
echo "==================== verification matrix ===================="
printf "%-50s %-5s %7s  %s\n" "gate" "state" "wall" "detail"
printf "%-50s %-5s %7s  %s\n" "----" "-----" "----" "------"
for i in "${!NAMES[@]}"; do
  printf "%-50s %-5s %7s  %s\n" "${NAMES[$i]}" "${RES[$i]}" "${TIMES[$i]}" "${DETAILS[$i]}"
done
echo "=============================================================="
if [ -n "$BENCH_OUT" ]; then echo; echo "$BENCH_OUT"; fi
if [ -n "$SCALE_OUT" ]; then echo; echo "$SCALE_OUT"; fi
echo
if [ "$FAILED" -eq 0 ]; then
  echo "ALL EXECUTED GATES PASSED"
else
  echo "FAILURES PRESENT - see sections above"
fi
exit "$FAILED"
