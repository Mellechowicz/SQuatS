#!/usr/bin/env bash
# Local/login-node resubmission chain: run -> exit 3 (budget exhausted) ->
# --resume until success (exit 0) or MAX_SEGMENTS. The chain is bit-identical
# to an uninterrupted run (T-K1) because only budget caps change between
# segments. Give each segment more budget via extra args, e.g.:
#   MAX_SEGMENTS=3 scripts/chain_resume.sh --set evolution.max_wall_s=3600
# (each segment then gets a fresh hour), or raise max_generations per segment.
set -euo pipefail
CFG=${CFG:-configs/w70cr30_4x4x4.yaml}
RUNDIR=${RUNDIR:-runs/chain}
BIN=${BIN:-./build/exsqs}
MAX_SEGMENTS=${MAX_SEGMENTS:-10}

for seg in $(seq 1 "$MAX_SEGMENTS"); do
  RESUME=()
  [[ -f "$RUNDIR/state.ckpt" ]] && RESUME=(--resume "$RUNDIR")
  echo "== chain segment $seg =="
  rc=0
  "$BIN" "$CFG" --out "$RUNDIR" "$@" "${RESUME[@]}" || rc=$?
  if [[ $rc -eq 0 ]]; then
    echo "chain: SUCCESS after $seg segment(s)"
    exit 0
  fi
  if [[ $rc -ne 3 ]]; then
    echo "chain: error rc=$rc"
    exit "$rc"
  fi
done
echo "chain: budget exhausted after $MAX_SEGMENTS segments (state kept in $RUNDIR)"
exit 3
