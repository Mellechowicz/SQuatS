#include "exsqs/displacements.hpp"

#include <vector>

namespace exsqs {
namespace {

using Dir = std::array<int, 3>;

// phonopy directions_diag; the exact order matters for the greedy search.
const std::vector<Dir> kDiag = {{1, 0, 0},  {0, 1, 0},  {0, 0, 1},  {1, 1, 0},  {1, 0, 1},
                                {0, 1, 1},  {1, -1, 0}, {1, 0, -1}, {0, 1, -1}, {1, 1, 1},
                                {1, 1, -1}, {1, -1, 1}, {-1, 1, 1}};

// phonopy _get_displacement_one: first direction whose site-symmetry orbit
// contains two images that, together with the direction, span 3D.
bool disp_one(const std::vector<Mat3i>& ss, const std::vector<Dir>& dirs, Dir& out) {
  for (const Dir& d : dirs) {
    std::vector<Dir> rd;
    rd.reserve(ss.size());
    for (const auto& r : ss) rd.push_back(apply_rot(r, d));
    const int n = static_cast<int>(ss.size());
    for (int i = 0; i < n; ++i)
      for (int j = i + 1; j < n; ++j)
        if (det_rows(d, rd[i], rd[j]) != 0) {
          out = d;
          return true;
        }
  }
  return false;
}

// phonopy _get_displacement_two
bool disp_two(const std::vector<Mat3i>& ss, const std::vector<Dir>& dirs,
              std::array<Dir, 2>& out) {
  for (const Dir& d1 : dirs) {
    std::vector<Dir> rd;
    rd.reserve(ss.size());
    for (const auto& r : ss) rd.push_back(apply_rot(r, d1));
    const int n = static_cast<int>(ss.size());
    for (int i = 0; i < n; ++i)
      for (const Dir& d2 : dirs)
        if (det_rows(d1, rd[i], d2) != 0) {
          out = {d1, d2};
          return true;
        }
  }
  return false;
}

// phonopy is_minus_displacement: the minus displacement is NOT needed iff some
// site-symmetry rotation maps d exactly to -d.
bool minus_needed(const Dir& d, const std::vector<Mat3i>& ss) {
  for (const auto& r : ss) {
    const Dir rd = apply_rot(r, d);
    if (rd[0] + d[0] == 0 && rd[1] + d[1] == 0 && rd[2] + d[2] == 0) return false;
  }
  return true;
}

// phonopy get_displacement with is_trigonal = False
std::vector<Dir> get_displacement(const std::vector<Mat3i>& ss) {
  Dir one;
  if (disp_one(ss, kDiag, one)) return {one};
  std::array<Dir, 2> two;
  if (disp_two(ss, kDiag, two)) return {two[0], two[1]};
  return {kDiag[0], kDiag[1], kDiag[2]};
}

}  // namespace

std::vector<DispEntry> least_displacements(const Structure& s, double symprec) {
  return least_displacements(s, get_symmetry(s, symprec), symprec);
}

std::vector<DispEntry> least_displacements(const Structure& s, const SymmetryInfo& info,
                                           double symprec) {
  std::vector<DispEntry> out;
  for (int a : info.independent_atoms()) {
    const auto ss = site_symmetry(s, info, a, symprec);
    for (const Dir& d : get_displacement(ss)) {
      out.push_back({a, d});
      if (minus_needed(d, ss)) out.push_back({a, {-d[0], -d[1], -d[2]}});
    }
  }
  return out;
}

int displacement_count(const Structure& s, double symprec) {
  return static_cast<int>(least_displacements(s, symprec).size());
}

int displacement_count(const Structure& s, const SymmetryInfo& info, double symprec) {
  return static_cast<int>(least_displacements(s, info, symprec).size());
}

}  // namespace exsqs
