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

# ---- C2: every SPEC section-12 test id is cited by a test/tool -------------
print("C2 SPEC 12 test-id coverage")
spec_ids = set(re.findall(r"^\| (T-[A-Z0-9]+) ", spec, re.M))
cited = set()
for d in ("tests", "tools"):
    for p in pathlib.Path(d).rglob("*"):
        if p.suffix in (".cpp", ".sh", ".py"):
            cited |= set(re.findall(r"T-[A-Z0-9]+", p.read_text(errors="ignore")))
cited |= set(re.findall(r"T-[A-Z0-9]+", open("tools/run_all_tests.sh").read()))
missing = sorted(spec_ids - cited)
report("all spec ids cited in tests/tools", not missing, ",".join(missing) or f"{len(spec_ids)} ids")

# ---- C3: every Catch2 tag is exercised by the runner -----------------------
print("C3 tag coverage in run_all_tests.sh")
r = sh(["./build/tests/exsqs_tests", "--list-tags"])
tags = set(re.findall(r"\[(\w[\w-]*)\]", r.stdout))
runner = open("tools/run_all_tests.sh").read()
orphans = sorted(t for t in tags if f"[{t}]" not in runner)
report("no orphan tags", not orphans, ",".join(orphans) or f"{len(tags)} tags")

# ---- C4: no stale version claims in user-facing strings --------------------
print("C4 user-facing strings")
bad = []
for p in list(pathlib.Path("src").glob("*.cpp")):
    for i, line in enumerate(p.read_text().splitlines(), 1):
        if ("fail(" in line or "fprintf" in line or "printf" in line) and re.search(r"v0\.\d", line):
            bad.append(f"{p.name}:{i}")
report("no vX.Y claims in errors/prints", not bad, ",".join(bad))

# ---- C5: trajectory signature covers every RunConfig field -----------------
print("C5 signature field ledger")
# Documented exclusions (rationale): budgets are raisable on resume
# (max_generations, max_wall_s); output/logging never affect the trajectory
# (outdir, checkpoint_every, log_info, config_echo); thread counts are
# performance-only [A14] (omp_threads); x_target is redundant once counts and
# x_achieved are signed [A5].
EXCLUDED = {"max_generations", "max_wall_s", "outdir", "checkpoint_every",
            "log_info", "config_echo", "omp_threads", "x_target"}
hpp = open("include/exsqs/config.hpp").read()
i = hpp.index("struct RunConfig")
j = hpp.index("{", i)
depth = 0
k = j
for k in range(j, len(hpp)):  # brace-count: member initializers contain "};"
    if hpp[k] == "{":
        depth += 1
    elif hpp[k] == "}":
        depth -= 1
        if depth == 0:
            break
block = hpp[j:k]
fields = set(re.findall(r"^\s+[A-Za-z_][\w:<>, ]*?\s+(\w+)\s*[=;{]", block, re.M))
fields -= {"RunConfig"}
signed = set(re.findall(r"c\.(\w+)", open("src/serialize.cpp").read()))
unclassified = sorted(fields - signed - EXCLUDED)
phantom = sorted(signed - fields)
report("every field signed or excluded", not unclassified, ",".join(unclassified) or f"{len(fields)} fields")
report("signature references real fields", not phantom, ",".join(phantom))
report("no excluded field is also signed", not (EXCLUDED & signed), ",".join(sorted(EXCLUDED & signed)))

# ---- C6: the SPEC section-11 sample config actually runs -------------------
print("C6 SPEC 11 sample config")
m = re.search(r"^## 11\..*?```(?:yaml)?\n(.*?)```", spec, re.S | re.M)
if not m:
    report("section-11 yaml block found", False)
else:
    tmpf = tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False)
    tmpf.write(m.group(1))
    tmpf.close()
    TMPS.append(tmpf.name)
    ok, _, r = run_config(tmpf.name)
    report("section-11 sample loads and runs", ok, (r.stderr.strip().splitlines() or [""])[-1][:70])

# ---- C7: reference configs run --------------------------------------------
print("C7 reference configs")
for cfg in sorted(pathlib.Path("configs").glob("*.yaml")):
    ok, _, r = run_config(str(cfg))
    report(f"{cfg.name} loads and runs", ok, (r.stderr.strip().splitlines() or [""])[-1][:70] if not ok else "")

# ---- C8: SPEC structure -----------------------------------------------------
print("C8 SPEC structure")
head = re.search(r"Specification v([\d.]+)", spec)
chlog = re.search(r"^- \*\*v([\d.]+)", spec[spec.index("## 16."):], re.M)
report("header version == latest changelog entry",
       bool(head and chlog) and head.group(1) == chlog.group(1),
       f"{head.group(1) if head else '?'} vs {chlog.group(1) if chlog else '?'}")
missing_sec = [str(i) for i in range(1, 17) if not re.search(rf"^## {i}\. ", spec, re.M)]
report("sections 1..16 present", not missing_sec, ",".join(missing_sec))

# ---- C9: README references resolve -----------------------------------------
print("C9 README references")
refs = set(re.findall(r"(?<!\w)(?:tools|configs)/[\w./-]+", readme))  # _scripts/ is untracked by design
dead = sorted(x for x in refs if not pathlib.Path(x.rstrip(".,)")).exists())
report("referenced paths exist", not dead, ",".join(dead))
bins = set(re.findall(r"\./build/(\w+)", readme))
deadb = sorted(b for b in bins if not pathlib.Path(f"build/{b}").exists())
report("referenced binaries built", not deadb, ",".join(deadb))

for t in TMPS:
    shutil.rmtree(t, ignore_errors=True) if os.path.isdir(t) else os.unlink(t)

print("COHERENCE " + ("PASS" if OK else "FAIL"))
sys.exit(0 if OK else 1)
