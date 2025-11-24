#pragma once
#include <array>
#include <vector>

#include "exsqs/structure.hpp"
#include "exsqs/symmetry.hpp"

namespace exsqs {

// Exact port of phonopy 4.3.1 harmonic/displacement.py get_least_displacements
// under the pinned convention [A15]: is_plusminus='auto', is_diagonal=True,
// is_trigonal=False, supercell used as-is (supercell_matrix = identity).
struct DispEntry {
  int atom;
  std::array<int, 3> dir;  // lattice-vector units, as in phonopy
};

std::vector<DispEntry> least_displacements(const Structure& s, double symprec = 1e-5);
// Overloads reusing a precomputed dataset (one spglib call per candidate).
std::vector<DispEntry> least_displacements(const Structure& s, const SymmetryInfo& info,
                                           double symprec = 1e-5);

// D(sigma): number of displaced supercells phonopy would generate [A15].
int displacement_count(const Structure& s, double symprec = 1e-5);
int displacement_count(const Structure& s, const SymmetryInfo& info, double symprec = 1e-5);

}  // namespace exsqs
