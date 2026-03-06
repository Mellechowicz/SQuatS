#!/usr/bin/env python3
"""Generate reference data for the exsqs step-1 gate tests.

Outputs (into <data>/):
  * dcounts_ref.txt : displaced-supercell counts from phonopy (T-D1 oracle) [A15]
  * sg_ref.txt      : space-group numbers from pymatgen (T-S1 cross-check),
                      falling back to python-spglib if pymatgen is unavailable.

Usage: python3 tools/py/gen_reference.py [tests/data]
"""
import glob
import os
import sys

import numpy as np
from phonopy import Phonopy
from phonopy.interface.vasp import read_vasp


def dcount(path: str) -> int:
    import warnings

    cell = read_vasp(path)
    with warnings.catch_warnings():
        # primitive_matrix only affects band-path bookkeeping, not the
        # displacement dataset, which is generated from the supercell symmetry.
        warnings.simplefilter("ignore")
        ph = Phonopy(cell, supercell_matrix=np.eye(3, dtype=int), symprec=1e-5, log_level=0)
    # pinned convention [A15]: defaults is_plusminus='auto', is_diagonal=True
    ph.generate_displacements(distance=0.01)
    return len(ph.supercells_with_displacements)


def sg_pymatgen(path: str) -> int:
    from pymatgen.core import Structure
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    return SpacegroupAnalyzer(Structure.from_file(path), symprec=1e-5).get_space_group_number()


def sg_spglib(path: str) -> int:
    import spglib

    cell = read_vasp(path)
    ds = spglib.get_symmetry_dataset(
        (cell.cell, cell.scaled_positions, cell.numbers), symprec=1e-5
    )
    return int(ds.number)


def main() -> None:
    data = sys.argv[1] if len(sys.argv) > 1 else "tests/data"
    cells = sorted(glob.glob(os.path.join(data, "cells", "*.vasp")))
    if not cells:
        sys.exit(f"no POSCARs under {data}/cells - run exsqs_gentests first")

    with open(os.path.join(data, "dcounts_ref.txt"), "w") as f:
        for p in cells:
            name = os.path.splitext(os.path.basename(p))[0]
            f.write(f"{name} {dcount(p)}\n")
    print(f"wrote {data}/dcounts_ref.txt ({len(cells)} cells)")

    try:
        sg_fn, src = sg_pymatgen, "pymatgen"
        sg_fn(cells[0])
    except Exception:
        sg_fn, src = sg_spglib, "python-spglib (pymatgen unavailable)"
    with open(os.path.join(data, "sg_ref.txt"), "w") as f:
        for p in cells:
            name = os.path.splitext(os.path.basename(p))[0]
            f.write(f"{name} {sg_fn(p)}\n")
    print(f"wrote {data}/sg_ref.txt via {src}")


if __name__ == "__main__":
    main()
