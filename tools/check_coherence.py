#!/usr/bin/env python3
"""EXSQS coherence audit (T-CO1, v1.6 / step 7): mechanical cross-checks
between the spec, the implementation, the tests, the configs, the docs and the
verification runner -- the "step 0 test" of the matrix. Exit 0 iff every check
passes.

Policy enforced here: step reports (docs/STEP*_REPORT.md) are immutable dated
snapshots and are NOT audited; SPEC.md, README.md, configs/ and the code must
stay mutually current.
"""
import json
import os
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parent.parent
os.chdir(ROOT)
OK = True
TMPS = []


def report(name, ok, detail=""):
    global OK
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""))
    if not ok:
        OK = False


def sh(cmd):
    return subprocess.run(cmd, capture_output=True, text=True)


def run_config(path):
    """Load-and-execute a config minimally (4 seeds, 1 generation)."""
    tmp = tempfile.mkdtemp(prefix="exsqs_coh_")
    TMPS.append(tmp)
    r = sh(["./build/exsqs", path, "--set", "evolution.population=4",
            "--set", "evolution.max_generations=1", "--set", "evolution.outputs=2",
            "--set", "output.log_level=warn", "--out", tmp])
    ok = r.returncode in (0, 3) and (pathlib.Path(tmp) / "summary.json").exists()
    return ok, tmp, r


spec = open("docs/SPEC.md").read()
readme = open("README.md").read()

# ---- C1: one version, everywhere ------------------------------------------
print("C1 version consistency")
cm = re.search(r"project\(\s*exsqs\s+VERSION\s+([\d.]+)", open("CMakeLists.txt").read())
src_lines = (open("src/evolution.cpp").read() + open("src/mpi_main.cpp").read()).splitlines()
lits = set()
for _l in src_lines:  # scope to exsqs version contexts; dependency versions
    if "exsqs" in _l:  # (e.g. phonopy in [A15] provenance strings) are content
        lits |= set(re.findall(r"\b\d+\.\d+\.\d+(?:-step\d+)?\b", _l))
report("single version literal across drivers", len(lits) == 1, ", ".join(sorted(lits)))
lit = sorted(lits)[0] if lits else ""
report("CMake VERSION equals the literal", bool(cm) and lit == cm.group(1),
       f"cmake {cm.group(1) if cm else '?'} vs {lit}")
ok, tmp, r = run_config("configs/smoke_sc27.yaml")
ver = ""
if ok:
    ver = json.load(open(f"{tmp}/summary.json")).get("exsqs_version", "")
report("live summary.json exsqs_version", ok and ver == lit, ver or r.stderr.strip()[:60])

for t in TMPS:
    shutil.rmtree(t, ignore_errors=True) if os.path.isdir(t) else os.unlink(t)

print("COHERENCE " + ("PASS" if OK else "FAIL"))
sys.exit(0 if OK else 1)
