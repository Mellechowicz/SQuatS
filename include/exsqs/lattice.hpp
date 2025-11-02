#pragma once
#include <string>
#include <vector>

#include "exsqs/linalg.hpp"

namespace exsqs {

// Undecorated prototype: cell (rows = basis vectors) + fractional basis sites.
// v1 scope: all sites belong to a single Wyckoff orbit (guard lives in zones).
struct Prototype {
  std::string name;
  Mat3 cell;
  std::vector<Vec3> sites;
};

Prototype make_sc(double a);
Prototype make_bcc(double a);                    // conventional cubic cell, 2 sites
Prototype make_fcc(double a);                    // conventional cubic cell, 4 sites
Prototype make_hcp(double a, double c = -1.0);   // c < 0 -> ideal c/a = sqrt(8/3)

}  // namespace exsqs
