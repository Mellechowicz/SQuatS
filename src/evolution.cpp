#include "exsqs/evolution.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <set>
#include <stdexcept>

#include "exsqs/correlation.hpp"
#include "exsqs/dedup.hpp"
#include "exsqs/displacements.hpp"
#include "exsqs/lattice.hpp"
#include "exsqs/rng.hpp"

namespace exsqs {

RunContext RunContext::build(const RunConfig& cfg) {
  RunContext ctx;
  ctx.geom = make_supercell(cfg.proto, cfg.H, cfg.species, 0);
  const int N = ctx.geom.natoms();
  int csum = 0;
  for (int c : cfg.counts) csum += c;
  if (csum != N) throw std::runtime_error("evolution: composition counts do not sum to N");
  ctx.zones = build_zones(ctx.geom, cfg.n_shells, cfg.shell_tol);
  ctx.weights = make_weights(cfg.wform, ctx.zones, cfg.wpow, cfg.wcustom);
  ctx.e_floor = cfg.full_pairs
                    ? e_floor_full(ctx.zones, cfg.counts, cfg.x_achieved, ctx.weights)
                    : e_floor_diagonal(ctx.zones, cfg.counts, cfg.x_achieved, ctx.weights);
  ctx.perms = site_permutations(ctx.geom, cfg.symprec);
  return ctx;
}

double effective_e_tol(const RunConfig& cfg, const RunContext& ctx) {
  return cfg.e_tol < 0 ? 3.0 * ctx.e_floor : cfg.e_tol;
}

std::vector<int> composition_sigma(const RunConfig& cfg) {
  std::vector<int> s;
  for (size_t t = 0; t < cfg.counts.size(); ++t)
    for (int k = 0; k < cfg.counts[t]; ++k) s.push_back(static_cast<int>(t));
  return s;
}

Individual evaluate(const RunConfig& cfg, const RunContext& ctx, std::vector<int> sigma) {
  Individual I;
  const Structure dec = decorate(ctx.geom, sigma, cfg.species);
  const CorrData cd = count_pairs(dec, ctx.zones);
  I.e_pure = cfg.full_pairs ? e_pure_full(cd, cfg.x_achieved, ctx.weights)
                            : e_pure_diagonal(cd, cfg.x_achieved, ctx.weights);
  const SymmetryInfo info = get_symmetry(dec, cfg.symprec);
  I.sg = info.sg_number;
  I.sg_symbol = info.sg_symbol;
  I.D = displacement_count(dec, info, cfg.symprec);
  // scalarized objective: correlation error times D^gamma
  I.e_obj = I.e_pure * std::pow(static_cast<double>(I.D), cfg.gamma);
  I.sigma = std::move(sigma);
  return I;
}

namespace {

void consider_output(std::vector<Individual>& top, const Individual& I, int M) {
  top.push_back(I);
  std::sort(top.begin(), top.end(),
            [](const Individual& a, const Individual& b) { return a.e_obj < b.e_obj; });
  if (static_cast<int>(top.size()) > M) top.resize(static_cast<size_t>(M));
}

}  // namespace

}  // namespace exsqs
