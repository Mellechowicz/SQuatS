# exsqs — extinction-based evolutionary SQS generator (minimal)

Minimal C++17 implementation of the extinction evolutionary algorithm for
Special Quasirandom Structures (arXiv:2602.10872). The engine minimizes
E_obj = E_pure * D^gamma over decorations of a supercell: correlation error
against the ideal random alloy, traded against the number D of
symmetry-inequivalent displacements a phonon workflow must compute.

Build:  cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
Run:    ./build/exsqs configs/w70cr30_4x4x4.yaml --out runs/demo
Exit codes: 0 converged to e_tol, 3 budget exhausted, 1 error.

Outputs: best_XX.vasp (POSCAR) and summary.json in --out.
License: MIT (see LICENSE).
