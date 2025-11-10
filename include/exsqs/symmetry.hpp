#pragma once
#include <string>
#include <vector>

#include "exsqs/structure.hpp"

namespace exsqs {

struct SymOp {
  Mat3i R;  // fractional rotation, acts as f' = R f + t
  Vec3 t;
};

struct SymmetryInfo {
  int sg_number = 0;
  std::string sg_symbol;
  std::string pointgroup;
  std::vector<SymOp> ops;
  std::vector<int> equivalent_atoms;
  std::vector<int> independent_atoms() const;  // phonopy: i with eq[i] == i, ascending
};

// spglib wrapper; ideal positions => symprec-robust. [A6] default 1e-5.
// Species are mapped to atomic numbers where known (phonopy parity); unknown
// names fall back to 200 + species index.
SymmetryInfo get_symmetry(const Structure& s, double symprec = 1e-5);

// Exact port of phonopy Symmetry._get_site_symmetry (4.3.1): rotations R of
// ops (R,t) with || (f_a - (R f_a + t) - rint(...)) * L ||_2 < symprec.
std::vector<Mat3i> site_symmetry(const Structure& s, const SymmetryInfo& info, int atom,
                                 double symprec);

inline bool is_p1(const SymmetryInfo& i) { return i.sg_number == 1; }

// Runtime spglib version string, e.g. "2.7.0".
std::string spglib_version();

// Point-group order of a dataset: n_operations / n_pure_translations (R == I).
int pointgroup_order(const SymmetryInfo& info);


}  // namespace exsqs
