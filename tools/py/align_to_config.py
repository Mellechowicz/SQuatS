#!/usr/bin/env python3
"""Align an arbitrary structure file onto an exsqs config geometry (v1.7, T-X2).

Rotated, translated, permuted or relabelled cells (paper supplements, ATAT
output, Niggli-reduced files) cannot be scored directly, because `exsqs score`
is deliberately orientation-strict. This tool maps the input onto the config's
supercell frame with pymatgen's StructureMatcher using a species-blind
FrameworkComparator -- the *frame* is matched, the input's species decoration
is transferred -- and writes a POSCAR that `exsqs score` accepts.

Workflow:
    ./build/exsqs geom config.yaml -o geom.vasp
    tools/py/align_to_config.py geom.vasp input.(vasp|cif|...) aligned.vasp
    ./build/exsqs score config.yaml aligned.vasp

Note: the aligner fixes the frame, not the coordinates. If the input's sites
are relaxed away from the ideal lattice by more than the scorer's matching
tolerance (1e-4 fractional), `exsqs score` will still reject it -- by design,
since exsqs scores ideal decorations.
"""
import sys

from pymatgen.analysis.structure_matcher import FrameworkComparator, StructureMatcher
from pymatgen.core import Structure

if len(sys.argv) != 4:
    print(__doc__)
    sys.exit(2)

geom = Structure.from_file(sys.argv[1])
inp = Structure.from_file(sys.argv[2])
m = StructureMatcher(ltol=0.2, stol=0.3, angle_tol=5.0, primitive_cell=False,
                     scale=False, attempt_supercell=False,
                     comparator=FrameworkComparator())
aligned = m.get_s2_like_s1(geom, inp)
if aligned is None:
    sys.exit("align: input does not match the config geometry as a framework "
             "(different supercell, composition count, or heavily relaxed sites)")
aligned.to(fmt="poscar", filename=sys.argv[3])
print(f"aligned structure written to {sys.argv[3]} ({len(aligned)} sites)")
