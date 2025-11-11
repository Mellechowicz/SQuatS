#include "exsqs/symmetry.hpp"

#include <cmath>
#include <map>
#include <stdexcept>

extern "C" {
#include <spglib.h>
}

namespace exsqs {
namespace {

// Minimal symbol -> atomic number table (phonopy passes atomic numbers to
// spglib; matching the labels removes one axis of potential divergence).
int atomic_number(const std::string& sym, int fallback_index) {
  static const std::map<std::string, int> z = {
      {"H", 1},   {"He", 2},  {"Li", 3},  {"Be", 4},  {"B", 5},   {"C", 6},   {"N", 7},
      {"O", 8},   {"Mg", 12}, {"Al", 13}, {"Si", 14}, {"Ti", 22}, {"V", 23},  {"Cr", 24},
      {"Mn", 25}, {"Fe", 26}, {"Co", 27}, {"Ni", 28}, {"Cu", 29}, {"Zn", 30}, {"Nb", 41},
      {"Mo", 42}, {"Ta", 73}, {"W", 74},  {"Re", 75}, {"Pt", 78}, {"Au", 79}};
  const auto it = z.find(sym);
  if (it != z.end()) return it->second;
  return 200 + fallback_index;
}

}  // namespace

std::vector<int> SymmetryInfo::independent_atoms() const {
  std::vector<int> out;
  for (size_t i = 0; i < equivalent_atoms.size(); ++i)
    if (static_cast<int>(i) == equivalent_atoms[i]) out.push_back(static_cast<int>(i));
  return out;
}

SymmetryInfo get_symmetry(const Structure& s, double symprec) {
  const int N = s.natoms();
  if (N <= 0) throw std::runtime_error("get_symmetry: empty structure");

  // spglib C API stores basis vectors as COLUMNS: lat[i][j] = i-th cartesian
  // component of the j-th basis vector. Our Mat3 rows are basis vectors.
  double lat[3][3];
  for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 3; ++c) lat[r][c] = s.cell[c][r];

  std::vector<double> pos(static_cast<size_t>(N) * 3);
  for (int i = 0; i < N; ++i)
    for (int k = 0; k < 3; ++k) pos[static_cast<size_t>(i) * 3 + k] = s.frac[i][k];
  std::vector<int> types(N);
  for (int i = 0; i < N; ++i) types[i] = atomic_number(s.names[s.species[i]], s.species[i]);

  SpglibDataset* d = spg_get_dataset(lat, reinterpret_cast<double(*)[3]>(pos.data()),
                                     types.data(), N, symprec);
  if (d == nullptr) throw std::runtime_error("get_symmetry: spg_get_dataset failed");

  SymmetryInfo out;
  out.sg_number = d->spacegroup_number;
  out.sg_symbol = d->international_symbol;
  out.pointgroup = d->pointgroup_symbol;
  auto trim = [](std::string& x) {
    while (!x.empty() && (x.back() == ' ' || x.back() == '\0')) x.pop_back();
    while (!x.empty() && x.front() == ' ') x.erase(x.begin());
  };
  trim(out.sg_symbol);
  trim(out.pointgroup);

  // operation extraction lands with the permutation work
  out.equivalent_atoms.assign(d->equivalent_atoms, d->equivalent_atoms + N);
  spg_free_dataset(d);
  return out;
}

std::string spglib_version() {
  return std::to_string(spg_get_major_version()) + "." + std::to_string(spg_get_minor_version()) +
         "." + std::to_string(spg_get_micro_version());
}

}  // namespace exsqs

