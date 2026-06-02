#!/usr/bin/env python3
"""Convert anything pymatgen reads (CIF, POSCAR, xyz, ...) to a VASP5 POSCAR
for exsqs (`lattice.file` prototypes or `exsqs score` inputs). v1.5 utility.

Usage: to_poscar.py <input> <output.vasp>
"""
import sys

from pymatgen.core import Structure

if len(sys.argv) != 3:
    print(__doc__)
    sys.exit(2)
Structure.from_file(sys.argv[1]).to(fmt="poscar", filename=sys.argv[2])
print(f"wrote {sys.argv[2]}")
