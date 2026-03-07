#include "exsqs/testcells.hpp"

#include <algorithm>
#include <random>
#include <stdexcept>

#include "exsqs/dedup.hpp"
#include "exsqs/lattice.hpp"
#include "exsqs/symmetry.hpp"

namespace exsqs {
namespace {

const std::vector<std::string> kWCr = {"W", "Cr"};

std::vector<int> random_sigma(int n, int n_b, uint64_t seed) {
  std::vector<int> s(static_cast<size_t>(n), 0);
  for (int i = 0; i < n_b; ++i) s[static_cast<size_t>(i)] = 1;
  std::mt19937_64 rng(seed);
  std::shuffle(s.begin(), s.end(), rng);
  return s;
}

// Site permutations of the empty geometry induced by involutive ops.
// kind 0: inversion-like (R == -I, t ~ integer); kind 1: proper two-fold
// rotation (det R = +1, R^2 = I, R != I).
std::vector<std::vector<int>> involution_perms(const Structure& geom, int kind) {
  Structure g = geom;
  std::fill(g.species.begin(), g.species.end(), 0);
  g.names = {"X"};
  const SymmetryInfo info = get_symmetry(g, 1e-5);

  std::vector<std::vector<int>> out;
  for (const auto& op : info.ops) {
    bool is_minus_i = true, is_i = true;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        const int diag = (i == j) ? 1 : 0;
        if (op.R[i][j] != -diag) is_minus_i = false;
        if (op.R[i][j] != diag) is_i = false;
      }
    if (kind == 0) {
      if (!is_minus_i) continue;
      bool t_int = true;
      for (int c = 0; c < 3; ++c)
        if (std::abs(op.t[c] - std::nearbyint(op.t[c])) > 1e-8) t_int = false;
      if (!t_int) continue;
    } else {
      if (is_i || is_minus_i) continue;
      if (det(op.R) != 1) continue;
      Mat3i R2{};
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
          int v = 0;
          for (int k = 0; k < 3; ++k) v += op.R[i][k] * op.R[k][j];
          R2[i][j] = v;
        }
      bool ident = true;
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          if (R2[i][j] != ((i == j) ? 1 : 0)) ident = false;
      if (!ident) continue;
    }
    auto perm = permutation_of_op(g, op.R, op.t);
    bool is_id = true, invol = true;
    for (int i = 0; i < geom.natoms(); ++i) {
      if (perm[static_cast<size_t>(i)] != i) is_id = false;
      if (perm[static_cast<size_t>(perm[static_cast<size_t>(i)])] != i) invol = false;
    }
    if (is_id || !invol) continue;
    out.push_back(std::move(perm));
  }
  return out;
}

// Decorate whole {id, g}-orbits so exactly n_b sites carry species 1.
bool try_symmetrized(const std::vector<int>& perm, int n_b, uint64_t seed,
                     std::vector<int>& out) {
  const int n = static_cast<int>(perm.size());
  std::vector<int> fixed, pair_rep;
  std::vector<char> seen(static_cast<size_t>(n), 0);
  for (int i = 0; i < n; ++i) {
    if (seen[static_cast<size_t>(i)]) continue;
    seen[static_cast<size_t>(i)] = 1;
    if (perm[static_cast<size_t>(i)] == i) {
      fixed.push_back(i);
    } else {
      seen[static_cast<size_t>(perm[static_cast<size_t>(i)])] = 1;
      pair_rep.push_back(i);
    }
  }
  int p = std::min<int>(static_cast<int>(pair_rep.size()), n_b / 2);
  const int f = n_b - 2 * p;
  if (f > static_cast<int>(fixed.size())) return false;

  std::mt19937_64 rng(seed);
  std::shuffle(fixed.begin(), fixed.end(), rng);
  std::shuffle(pair_rep.begin(), pair_rep.end(), rng);

  out.assign(static_cast<size_t>(n), 0);
  for (int k = 0; k < p; ++k) {
    out[static_cast<size_t>(pair_rep[static_cast<size_t>(k)])] = 1;
    out[static_cast<size_t>(perm[static_cast<size_t>(pair_rep[static_cast<size_t>(k)])])] = 1;
  }
  for (int k = 0; k < f; ++k) out[static_cast<size_t>(fixed[static_cast<size_t>(k)])] = 1;
  return true;
}

std::vector<int> symmetric_sigma(const Structure& geom, int kind, int n_b, uint64_t seed) {
  for (const auto& perm : involution_perms(geom, kind)) {
    std::vector<int> s;
    if (try_symmetrized(perm, n_b, seed, s)) return s;
  }
  throw std::runtime_error("symmetric_sigma: no feasible involutive op found");
}

NamedCell cell(const std::string& name, const Structure& s, int sg = -1) {
  return {name, grouped_by_species(s), sg};
}

}  // namespace

std::vector<NamedCell> gate_cells() {
  std::vector<NamedCell> out;

  // --- pure prototypes ---
  out.push_back(cell("sc_pure_222", make_supercell_diag(make_sc(3.0), 2, 2, 2, {"W"}), 221));
  out.push_back(cell("bcc_pure_222", make_supercell_diag(make_bcc(3.165), 2, 2, 2, {"W"}), 229));
  out.push_back(cell("fcc_pure_222", make_supercell_diag(make_fcc(4.05), 2, 2, 2, {"W"}), 225));
  out.push_back(cell("hcp_pure_221", make_supercell_diag(make_hcp(3.2), 2, 2, 1, {"Mg"}), 194));

  // --- ordered prototypes (analytic space groups) ---
  {
    auto g = make_supercell_diag(make_bcc(3.165), 1, 1, 1, kWCr);
    out.push_back(cell("b2_111", decorate(g, {0, 1}, kWCr), 221));
  }
  {
    auto g = make_supercell_diag(make_bcc(3.165), 2, 2, 2, kWCr);
    std::vector<int> sig(static_cast<size_t>(g.natoms()));
    for (int i = 0; i < g.natoms(); ++i) sig[static_cast<size_t>(i)] = i % 2;  // basis parity
    out.push_back(cell("b2_222", decorate(g, sig, kWCr), 221));
  }
  {
    auto g = make_supercell_diag(make_fcc(4.05), 1, 1, 1, kWCr);
    out.push_back(cell("l10_111", decorate(g, {0, 0, 1, 1}, kWCr), 123));
  }
  {
    auto g = make_supercell_diag(make_fcc(4.05), 1, 1, 1, kWCr);
    out.push_back(cell("l12_111", decorate(g, {0, 1, 1, 1}, kWCr), 221));
  }
  {
    auto g = make_supercell_diag(make_fcc(4.05), 2, 2, 2, kWCr);
    std::vector<int> sig(static_cast<size_t>(g.natoms()));
    for (int i = 0; i < g.natoms(); ++i) sig[static_cast<size_t>(i)] = (i % 4 == 0) ? 0 : 1;
    out.push_back(cell("l12_222", decorate(g, sig, kWCr), 221));
  }

  // --- symmetrized random decorations (residual symmetry; SG recorded, not asserted) ---
  {
    auto g = make_supercell_diag(make_bcc(3.165), 2, 2, 2, kWCr);
    out.push_back(cell("bcc16_inv_a", decorate(g, symmetric_sigma(g, 0, 5, 101), kWCr)));
    out.push_back(cell("bcc16_inv_b", decorate(g, symmetric_sigma(g, 0, 8, 202), kWCr)));
    out.push_back(cell("bcc16_c2_a", decorate(g, symmetric_sigma(g, 1, 5, 303), kWCr)));
    out.push_back(cell("bcc16_rand_a", decorate(g, random_sigma(16, 5, 404), kWCr)));
    out.push_back(cell("bcc16_rand_b", decorate(g, random_sigma(16, 8, 505), kWCr)));
  }
  {
    auto g = make_supercell_diag(make_fcc(4.05), 2, 2, 2, kWCr);
    out.push_back(cell("fcc32_rand_a", decorate(g, random_sigma(32, 10, 606), kWCr)));
  }
  {
    auto g = make_supercell_diag(make_sc(3.0), 3, 3, 3, kWCr);
    out.push_back(cell("sc27_rand_a", decorate(g, random_sigma(27, 13, 707), kWCr)));
    out.push_back(cell("sc27_inv_a", decorate(g, symmetric_sigma(g, 0, 13, 808), kWCr)));
  }

  return out;
}

}  // namespace exsqs
