#pragma once
#include <vector>

#include "exsqs/structure.hpp"

namespace exsqs {

// Coordination-zone table of an (undecorated) geometry, built from the 27-image
// background [A1]: ordered pairs (i in S, j in B), the zero-distance self pair
// excluded, nonzero self-images included as legitimate PBC neighbours.
struct ZoneTable {
  int n_shells = 0;
  std::vector<double> radii;    // representative (smallest) distance per shell
  std::vector<int> coord_num;   // Z_n, identical across sites (single-orbit guard)
  struct Nbr {
    int j;
    int shell;
  };
  std::vector<std::vector<Nbr>> nbrs;  // per site i: neighbours up to n_shells
  double half_min_width = 0.0;
  bool shells_exceed_half_width = false;  // [A3] warning flag
};

ZoneTable build_zones(const Structure& geom, int n_shells, double shell_rel_tol = 1e-3);

}  // namespace exsqs
