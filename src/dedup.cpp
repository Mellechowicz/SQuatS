#include "exsqs/dedup.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>

#include "exsqs/symmetry.hpp"

namespace exsqs {
namespace {

// Quantize a fractional coordinate wrapped into [0,1) to a 1e-6 grid. Ideal
// lattice positions are exact rationals, so 1e-6 quantization with the 1e-8
// wrap guard is collision-free and rounding-stable.
std::array<long long, 3> qkey(const Vec3& f) {
  std::array<long long, 3> k;
  for (int c = 0; c < 3; ++c) {
    double y = f[c] - std::floor(f[c]);
    if (y > 1.0 - 1e-8) y = 0.0;
    k[c] = std::llround(y * 1e6);
  }
  return k;
}

std::map<std::array<long long, 3>, int> site_lut(const Structure& geom) {
  std::map<std::array<long long, 3>, int> lut;
  for (int i = 0; i < geom.natoms(); ++i) {
    if (!lut.emplace(qkey(geom.frac[static_cast<size_t>(i)]), i).second)
      throw std::runtime_error("dedup: duplicate quantized site position");
  }
  return lut;
}

std::vector<int> perm_from_op(const std::map<std::array<long long, 3>, int>& lut,
                              const Structure& geom, const Mat3i& R, const Vec3& t) {
  std::vector<int> perm(static_cast<size_t>(geom.natoms()));
  for (int i = 0; i < geom.natoms(); ++i) {
    const Vec3 fp = apply_rot(R, geom.frac[static_cast<size_t>(i)]) + t;
    const auto it = lut.find(qkey(fp));
    if (it == lut.end())
      throw std::runtime_error("dedup: op does not map site set onto itself");
    perm[static_cast<size_t>(i)] = it->second;
  }
  return perm;
}

}  // namespace

std::vector<int> permutation_of_op(const Structure& geom, const Mat3i& R, const Vec3& t) {
  return perm_from_op(site_lut(geom), geom, R, t);
}

std::vector<std::vector<int>> site_permutations(const Structure& geom, double symprec) {
  Structure g = geom;
  std::fill(g.species.begin(), g.species.end(), 0);
  g.names = {"X"};
  const SymmetryInfo info = get_symmetry(g, symprec);
  const auto lut = site_lut(g);

  std::set<std::vector<int>> uniq;
  for (const auto& op : info.ops) uniq.insert(perm_from_op(lut, g, op.R, op.t));
  return std::vector<std::vector<int>>(uniq.begin(), uniq.end());
}

std::vector<uint8_t> canonical_labels(const std::vector<int>& sigma,
                                      const std::vector<std::vector<int>>& perms) {
  const size_t N = sigma.size();
  std::vector<uint8_t> best, cur(N);
  for (const auto& p : perms) {
    for (size_t i = 0; i < N; ++i) cur[static_cast<size_t>(p[i])] = static_cast<uint8_t>(sigma[i]);
    if (best.empty() || cur < best) best = cur;
  }
  if (best.empty()) {  // no perms supplied: identity only
    best.resize(N);
    for (size_t i = 0; i < N; ++i) best[i] = static_cast<uint8_t>(sigma[i]);
  }
  return best;
}

uint64_t hash_labels(const std::vector<uint8_t>& labels) {
  uint64_t h = 1469598103934665603ULL;
  for (uint8_t b : labels) {
    h ^= b;
    h *= 1099511628211ULL;
  }
  return h;
}

}  // namespace exsqs
