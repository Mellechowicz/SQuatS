#!/usr/bin/env python3
"""T-V1 [SPEC section 12]: independent recomputation of E_pure for every emitted
structure. Reads <rundir>/summary.json + best_*.vasp, rebuilds the 27-image pair
table in numpy using the recorded shell radii / tolerance / weights, and compares
E_pure against the engine's values. Acceptance: max |diff| <= 1e-10.

Usage: python3 validate.py <rundir>
"""
import json
import sys

import numpy as np
from pymatgen.core import Structure


def e_pure_of(path, species, x, radii, tol, w, mode):
    st = Structure.from_file(path)
    L = np.array(st.lattice.matrix)
    f = np.array([s.frac_coords for s in st.sites])
    t = np.array([species.index(str(s.specie)) for s in st.sites])
    N, K, S = len(f), len(species), len(radii)
    imgs = np.array(
        [[a, b, c] for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)], dtype=float
    )
    ub = np.array(radii) * (1.0 + tol)  # shell upper bounds (first-match rule)
    C = np.zeros((S, K, K), dtype=np.int64)
    P = np.zeros((S, K), dtype=np.int64)
    for i in range(N):
        d = f[None, :, :] + imgs[:, None, :] - f[i][None, None, :]  # (27, N, 3)
        r = np.linalg.norm(d @ L, axis=2).ravel()  # (27*N,)
        jj = np.tile(np.arange(N), 27)
        keep = (r > 1e-9) & (r <= ub[-1])  # exclude only the exact self pair [A1]
        rr, jk = r[keep], jj[keep]
        sh = np.searchsorted(ub, rr, side="left")  # first n with r <= radii[n]*(1+tol)
        for n, jx in zip(sh, jk):
            C[n, t[i], t[jx]] += 1
            P[n, t[i]] += 1
    E = 0.0
    for n in range(S):
        for a in range(K):
            if P[n, a] == 0:
                continue
            if mode == "diagonal":
                E += w[n] * abs(C[n, a, a] / P[n, a] - x[a])
            else:
                for b in range(K):
                    E += w[n] * abs(C[n, a, b] / P[n, a] - x[b])
    return E


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    rundir = sys.argv[1]
    s = json.load(open(f"{rundir}/summary.json"))
    species, x = s["species"], s["x_achieved"]
    radii, tol = s["shell_radii"], s["shell_tol"]
    w, mode = s["weights"], s["error_mode"]
    worst = 0.0
    for o in s["outputs"]:
        E = e_pure_of(f"{rundir}/{o['file']}", species, x, radii, tol, w, mode)
        diff = abs(E - o["e_pure"])
        worst = max(worst, diff)
        print(f"{o['file']}: E_ref={E:.12e}  E_run={o['e_pure']:.12e}  |diff|={diff:.3e}")
    ok = worst <= 1e-10
    print(("T-V1 PASS" if ok else "T-V1 FAIL") + f"  (worst |diff| = {worst:.3e}, gate 1e-10)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
