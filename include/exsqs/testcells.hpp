#pragma once
#include <string>
#include <vector>

#include "exsqs/structure.hpp"

namespace exsqs {

// Gate-cell set for T-D1 / T-S1 [SPEC section 12]. Every decorated cell is
// species-grouped so the in-memory atom order equals the POSCAR order (the
// T-D1 oracle compares against phonopy reading those POSCARs, and phonopy's
// choice of orbit representatives depends on atom enumeration).
struct NamedCell {
  std::string name;
  Structure s;
  int expected_sg = -1;  // asserted only for ordered/pure prototypes
};

std::vector<NamedCell> gate_cells();

}  // namespace exsqs
