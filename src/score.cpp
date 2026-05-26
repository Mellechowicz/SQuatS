#include "exsqs/score.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "exsqs/correlation.hpp"
#include "exsqs/displacements.hpp"
#include "exsqs/symmetry.hpp"

namespace exsqs {

namespace {

double wrap_half(double d) {  // to (-0.5, 0.5]
  d -= std::floor(d);
  if (d > 0.5) d -= 1.0;
  return d;
}

std::string vec_str(const std::vector<int>& v) {
  std::ostringstream o;
  o << "[";
  for (size_t i = 0; i < v.size(); ++i) o << (i ? "," : "") << v[i];
  o << "]";
  return o.str();
}

}  // namespace

std::vector<int> sigma_on_geometry(const RunConfig& cfg, const RunContext& ctx,
                                   const Structure& s) {
  const Structure& g = ctx.geom;
  // lattice: elementwise, orientation-strict (pre-transform externally if needed)
  double mx = 0.0;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) mx = std::max(mx, std::abs(g.cell[i][j]));
  const double ltol = 1e-6 * mx + 1e-9;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      if (std::abs(s.cell[i][j] - g.cell[i][j]) > ltol)
        throw std::runtime_error(
            "score: lattice does not match the config geometry (same cell and orientation "
            "required; convert/standardize the input first)");
  if (s.natoms() != g.natoms())
    throw std::runtime_error("score: site count differs from the config geometry (" +
                             std::to_string(s.natoms()) + " vs " + std::to_string(g.natoms()) +
                             ")");
  // species name -> cfg index
  std::vector<int> name_map(static_cast<size_t>(s.nspecies()), -1);
  for (int t = 0; t < s.nspecies(); ++t) {
    for (size_t k = 0; k < cfg.species.size(); ++k)
      if (s.names[static_cast<size_t>(t)] == cfg.species[k]) name_map[static_cast<size_t>(t)] = static_cast<int>(k);
    if (name_map[static_cast<size_t>(t)] < 0)
      throw std::runtime_error("score: species '" + s.names[static_cast<size_t>(t)] +
                               "' is not in the config composition");
  }
  // geometry site i -> unique input site j by wrapped fractional distance
  const double stol = 1e-4;
  std::vector<char> used(static_cast<size_t>(s.natoms()), 0);
  std::vector<int> sigma(static_cast<size_t>(g.natoms()), -1);
  for (int i = 0; i < g.natoms(); ++i) {
    int hit = -1;
    for (int j = 0; j < s.natoms(); ++j) {
      if (used[static_cast<size_t>(j)]) continue;
      bool ok = true;
      for (int c = 0; c < 3 && ok; ++c)
        ok = std::abs(wrap_half(s.frac[static_cast<size_t>(j)][c] -
                                g.frac[static_cast<size_t>(i)][c])) <= stol;
      if (!ok) continue;
      if (hit >= 0)
        throw std::runtime_error("score: ambiguous site match (two input sites within tolerance "
                                 "of one geometry site)");
      hit = j;
    }
    if (hit < 0)
      throw std::runtime_error("score: geometry site " + std::to_string(i) +
                               " has no matching input site (same supercell required)");
    used[static_cast<size_t>(hit)] = 1;
    sigma[static_cast<size_t>(i)] = name_map[static_cast<size_t>(s.species[static_cast<size_t>(hit)])];
  }
  // composition must equal the config's [A5] counts
  std::vector<int> have(cfg.species.size(), 0);
  for (int t : sigma) have[static_cast<size_t>(t)]++;
  if (have != cfg.counts)
    throw std::runtime_error("score: composition " + vec_str(have) +
                             " differs from the config counts " + vec_str(cfg.counts));
  return sigma;
}

ScoreResult score_structure(const RunConfig& cfg, const RunContext& ctx, const Structure& s) {
  const std::vector<int> sigma = sigma_on_geometry(cfg, ctx, s);
  const Structure dec = decorate(ctx.geom, sigma, cfg.species);
  ScoreResult r;
  const SymmetryInfo info = get_symmetry(dec, cfg.symprec);
  r.sg = info.sg_number;
  r.sg_symbol = info.sg_symbol;
  r.pg_order = pointgroup_order(info);
  r.D = displacement_count(dec, info, cfg.symprec);
  const CorrData cd = count_pairs(dec, ctx.zones);
  r.e_pure = cfg.full_pairs ? e_pure_full(cd, cfg.x_achieved, ctx.weights)
                            : e_pure_diagonal(cd, cfg.x_achieved, ctx.weights);
  r.e_obj = r.e_pure * std::pow(static_cast<double>(r.D), cfg.gamma);
  return r;
}

}  // namespace exsqs
