#!/bin/bash
# T-X2 (v1.7): cross-tool, file-level closure of the scoring pipeline.
# Take the recorded campaign champion, rotate its cell 90 degrees and shuffle
# its site order (pymatgen), align it back onto the config frame with
# tools/py/align_to_config.py, score it with `exsqs score`, and require the
# result to reproduce the engine's recorded values BITWISE.
set -e
cd "$(dirname "$0")/.."
TMP=$(mktemp -d /tmp/exsqs_x2.XXXX)
trap 'rm -rf "$TMP"' EXIT

./build/exsqs geom configs/w70cr30_4x4x4.yaml -o "$TMP/geom.vasp" > /dev/null

python3 - "$TMP" << 'PY'
import random
import sys

from pymatgen.core import Structure
from pymatgen.core.operations import SymmOp

tmp = sys.argv[1]
s = Structure.from_file('runs/campaign_w70cr30/best_00.vasp')
rot = SymmOp.from_rotation_and_translation(((0, -1, 0), (1, 0, 0), (0, 0, 1)), (0, 0, 0))
s.apply_operation(rot)  # 90-degree rotation about z: new cell setting, same crystal
idx = list(range(len(s)))
random.seed(3)
random.shuffle(idx)
s2 = Structure(s.lattice, [s[i].specie for i in idx], [s[i].frac_coords for i in idx])
s2.to(fmt='poscar', filename=f'{tmp}/rotated.vasp')
PY

python3 tools/py/align_to_config.py "$TMP/geom.vasp" "$TMP/rotated.vasp" "$TMP/aligned.vasp" > /dev/null
./build/exsqs score configs/w70cr30_4x4x4.yaml "$TMP/aligned.vasp" --json "$TMP/sc.json" > /dev/null

python3 - "$TMP" << 'PY'
import json
import sys

tmp = sys.argv[1]
r = json.load(open(f'{tmp}/sc.json'))['scores'][0]
e = json.load(open('runs/campaign_w70cr30/summary.json'))['outputs'][0]
ok = r['e_pure'] == e['e_pure'] and r['D'] == e['D'] and r['sg'] == e['sg']
if ok:
    print('T-X2 PASS: rotated+shuffled champion scores bitwise-equal after alignment '
          f"(E_pure={r['e_pure']:.6e}, D={r['D']}, SG={r['sg_symbol']})")
else:
    print('T-X2 FAIL:', r, 'vs', e)
sys.exit(0 if ok else 1)
PY
