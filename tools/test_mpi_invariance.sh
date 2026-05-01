#!/usr/bin/env bash
# T-MPI1 (v1.3, SPEC 8.2): results are independent of MPI rank count and
# identical to the single-process exsqs. Compares pooled outputs, per-island
# generation counts, the migration ledger and the full per-generation logs
# (elapsed_s, the wall-clock column, excluded).
set -euo pipefail
EXSQS=$1; EXSQS_MPI=$2; CFG=$3; MPIEXEC=${4:-mpirun}

FLAGS=""
if "$MPIEXEC" --version 2>/dev/null | head -1 | grep -qi "open"; then
  FLAGS="--allow-run-as-root --oversubscribe"   # container/CI artifacts, harmless elsewhere
fi

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
COMMON=(--set parallel.islands=3 --set parallel.migration_every=4 --set parallel.migrants=1
        --set evolution.population=12 --set evolution.max_generations=12
        --set evolution.e_tol=1e-12 --set output.log_level=warn)

"$EXSQS" "$CFG" "${COMMON[@]}" --out "$TMP/serial" > /dev/null 2>&1 || true
"$MPIEXEC" $FLAGS -n 1 "$EXSQS_MPI" "$CFG" "${COMMON[@]}" --out "$TMP/r1" > /dev/null 2>&1 || true
"$MPIEXEC" $FLAGS -n 3 "$EXSQS_MPI" "$CFG" "${COMMON[@]}" --out "$TMP/r3" > /dev/null 2>&1 || true

python3 - "$TMP" << 'PY'
import csv, json, sys
tmp = sys.argv[1]

def load(d):
    s = json.load(open(f"{tmp}/{d}/summary.json"))
    outs = [(o["e_pure"], o["D"], o["sg"]) for o in s["outputs"]]
    rows = []
    with open(f"{tmp}/{d}/generations.csv") as f:
        for r in csv.reader(f):
            rows.append(tuple(r[:-1]))  # drop elapsed_s (wall clock)
    return outs, rows, s["island_generations"], s["island_migrants_in"]

a, b, c = load("serial"), load("r1"), load("r3")
assert a == b, "serial vs 1-rank MISMATCH"
assert a == c, "serial vs 3-rank MISMATCH"
print("T-MPI1 PASS: serial == mpirun -n 1 == mpirun -n 3 "
      "(outputs, logs, generations, migration ledger)")
PY
