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

}  // namespace exsqs
