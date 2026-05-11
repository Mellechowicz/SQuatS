#!/usr/bin/env bash
# EXSQS full verification matrix, steps 1-5 (v1.4). Runs every automated gate
# in the project and prints a step-by-step PASS/FAIL table with timings.
#
#   bash tools/run_all_tests.sh
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

echo "EXSQS verification matrix ($(./build/exsqs --help 2>/dev/null | head -0; true)v1.4, $(date -u +%F))"
echo

t0=$(date +%s)
if cmake --build build -j2 > /tmp/exsqs_build.log 2>&1; then
  add "build (all targets)" PASS "$(( $(date +%s) - t0 ))s" "cmake --build build"
else
  add "build (all targets)" FAIL "$(( $(date +%s) - t0 ))s" "see /tmp/exsqs_build.log"
  tail -20 /tmp/exsqs_build.log; FAILED=1
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

echo
echo "==================== verification matrix ===================="
printf "%-50s %-5s %7s  %s\n" "gate" "state" "wall" "detail"
printf "%-50s %-5s %7s  %s\n" "----" "-----" "----" "------"
for i in "${!NAMES[@]}"; do
  printf "%-50s %-5s %7s  %s\n" "${NAMES[$i]}" "${RES[$i]}" "${TIMES[$i]}" "${DETAILS[$i]}"
done
echo "=============================================================="
echo
if [ "$FAILED" -eq 0 ]; then
  echo "ALL EXECUTED GATES PASSED"
else
  echo "FAILURES PRESENT - see sections above"
fi
exit "$FAILED"
