# exsqs

Extinction-based evolutionary generator of Special Quasirandom
Structures (SQS): minimizes E_obj = E_pure * D^gamma over supercell
decorations of a single lattice.

Build:  cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
Run:    ./build/exsqs configs/w70cr30_4x4x4.yaml --out runs/demo
Exit codes: 0 converged to e_tol, 3 budget exhausted, 1 error.
