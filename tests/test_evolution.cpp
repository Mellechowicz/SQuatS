#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <set>
#include <vector>

#include "exsqs/dedup.hpp"
#include <cmath>

#include "exsqs/correlation.hpp"
#include "exsqs/evolution.hpp"
#include "exsqs/lattice.hpp"
#include "exsqs/structure.hpp"
#include "exsqs/symmetry.hpp"

using namespace exsqs;
using Catch::Approx;

namespace {
RunConfig sc27_cfg() {
  RunConfig c;
  c.proto = make_sc(3.0);
  c.H = {{{3, 0, 0}, {0, 3, 0}, {0, 0, 3}}};
  c.species = {"W", "Cr"};
  c.counts = {14, 13};
  c.x_target = {14.0 / 27.0, 13.0 / 27.0};
  c.x_achieved = c.x_target;
  c.n_shells = 5;
  c.gamma = 1.0;
  c.population = 16;
  c.outputs = 4;
  c.e_tol = 0.0;  // unreachable (Pi denominators never hit 14/27 exactly)
  c.max_generations = 25;
  c.retry_budget = 60;
  c.seed = 7;
  c.log_info = false;
  return c;
}
}  // namespace

TEST_CASE("constructive seeding: exact composition, guaranteed residual symmetry [D4]",
          "[evolution]") {
  const RunConfig cfg = sc27_cfg();
  const RunContext ctx = RunContext::build(cfg);
  REQUIRE_FALSE(ctx.seed_perms.empty());
  int ok_count = 0;
  for (int slot = 0; slot < 12; ++slot) {
    bool ok = false;
    const auto sig = seed_sigma_constructive(cfg, ctx, 0, 0, slot, ok);
    std::vector<int> cnt(2, 0);
    for (int t : sig) cnt[static_cast<size_t>(t)]++;
    REQUIRE(cnt[0] == 14);
    REQUIRE(cnt[1] == 13);
    if (ok) {
      ++ok_count;
      const Structure s = decorate(ctx.geom, sig, cfg.species);
      REQUIRE(get_symmetry(s, cfg.symprec).sg_number > 1);
    }
  }
  REQUIRE(ok_count >= 6);  // rich op pool on sc27; most slots must find a feasible subgroup
}

TEST_CASE("symmetry-preserving mutation: child keeps a cyclic subgroup of the parent or falls "
          "back to a pair swap [D6][A12]",
          "[evolution]") {
  const RunConfig cfg = sc27_cfg();
  const RunContext ctx = RunContext::build(cfg);
  bool ok = false;
  std::vector<int> sig;
  int slot = 0;
  while (!ok && slot < 50) sig = seed_sigma_constructive(cfg, ctx, 0, 0, slot++, ok);
  REQUIRE(ok);
  Individual parent;
  parent.sigma = sig;
  const Structure ps = decorate(ctx.geom, sig, cfg.species);
  const SymmetryInfo pinfo = get_symmetry(ps, cfg.symprec);
  REQUIRE(pinfo.sg_number > 1);
  parent.eq_atoms = pinfo.equivalent_atoms;
  for (const auto& op : pinfo.ops) {
    bool ident = true;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        if (op.R[i][j] != ((i == j) ? 1 : 0)) ident = false;
    if (!ident) parent.stab_ops.push_back(op);
  }
  REQUIRE_FALSE(parent.stab_ops.empty());  // non-P1 => a nontrivial op exists
  CounterRng rng(99, 0, 1, 0, RngPurpose::Mutation);
  int subgroup_children = 0;
  for (int rep = 0; rep < 20; ++rep) {
    const auto child = mutate_sigma(parent, cfg, ctx, rng);
    std::vector<int> cnt(2, 0);
    for (int t : child) cnt[static_cast<size_t>(t)]++;
    REQUIRE(cnt[0] == 14);
    REQUIRE(cnt[1] == 13);
    // the child must stay invariant under AT LEAST ONE nontrivial parent op
    // (the preserved <g>); it MAY drop the rest of the parent's group.
    bool keeps_some_g = false;
    for (const auto& op : parent.stab_ops) {
      const auto perm = permutation_of_op(ctx.geom, op.R, op.t);
      bool inv = true;
      for (size_t i = 0; i < child.size() && inv; ++i) inv = (child[perm[i]] == child[i]);
      if (inv) {
        keeps_some_g = true;
        break;
      }
    }
    if (keeps_some_g) {
      ++subgroup_children;
    } else {
      int nd = 0;  // plain-swap fallback: exactly one unlike pair swapped
      for (size_t i = 0; i < child.size(); ++i) nd += (child[i] != parent.sigma[i]);
      REQUIRE(nd == 2);
    }
  }
  REQUIRE(subgroup_children >= 15);  // orbit swaps dominate on a symmetric parent
}

TEST_CASE("plain mutation (k=1): unlike-pair swap preserves composition [A12]", "[evolution]") {
  RunConfig cfg = sc27_cfg();
  cfg.mut_sympres = false;
  const RunContext ctx = RunContext::build(cfg);
  Individual parent;
  parent.sigma = seed_sigma_rejection(cfg, ctx, 0, 0);
  CounterRng rng(5, 0, 1, 3, RngPurpose::Mutation);
  for (int rep = 0; rep < 25; ++rep) {
    const auto child = mutate_sigma(parent, cfg, ctx, rng);
    std::vector<int> cnt(2, 0);
    for (int t : child) cnt[static_cast<size_t>(t)]++;
    REQUIRE(cnt[0] == 14);
    REQUIRE(cnt[1] == 13);
    int nd = 0;
    for (size_t i = 0; i < child.size(); ++i) nd += (child[i] != parent.sigma[i]);
    REQUIRE(nd == 2);
  }
}

TEST_CASE("evolve_island: monotone elite, unique non-P1 outputs, budget stop [D5][A7][A9][D7]",
          "[evolution]") {
  const RunConfig cfg = sc27_cfg();
  const RunContext ctx = RunContext::build(cfg);
  const IslandResult r = evolve_island(cfg, ctx, 0);
  REQUIRE_FALSE(r.success);
  REQUIRE(r.stop_reason == "max_generations");
  REQUIRE(r.log.size() == static_cast<size_t>(cfg.max_generations) + 1);
  for (size_t g = 1; g < r.log.size(); ++g)
    REQUIRE(r.log[g].best_e_obj <= r.log[g - 1].best_e_obj);  // [D5] E_min monotone
  REQUIRE_FALSE(r.pool.empty());
  REQUIRE(r.pool.size() <= static_cast<size_t>(cfg.outputs));
  std::set<std::vector<uint8_t>> keys;
  for (const auto& I : r.pool) {
    REQUIRE(I.sg > 1);
    keys.insert(I.canonical);
  }
  REQUIRE(keys.size() == r.pool.size());
}

TEST_CASE("T-R1: bit-identical reruns, seed sensitivity, sequential-island determinism",
          "[evolution][T-R1]") {
  RunConfig cfg = sc27_cfg();
  cfg.islands = 2;
  cfg.max_generations = 12;
  const RunContext ctx = RunContext::build(cfg);
  const RunOutput a = run_evolution(cfg, ctx);
  const RunOutput b = run_evolution(cfg, ctx);
  REQUIRE(a.log.size() == b.log.size());
  for (size_t i = 0; i < a.log.size(); ++i) {
    REQUIRE(a.log[i].island == b.log[i].island);
    REQUIRE(a.log[i].gen == b.log[i].gen);
    REQUIRE(a.log[i].best_e_pure == b.log[i].best_e_pure);  // bitwise double equality
    REQUIRE(a.log[i].best_e_obj == b.log[i].best_e_obj);
    REQUIRE(a.log[i].median_e_obj == b.log[i].median_e_obj);
    REQUIRE(a.log[i].best_D == b.log[i].best_D);
    REQUIRE(a.log[i].best_sg == b.log[i].best_sg);
    REQUIRE(a.log[i].killed == b.log[i].killed);
    REQUIRE(a.log[i].evals == b.log[i].evals);
  }
  REQUIRE(a.outputs.size() == b.outputs.size());
  for (size_t i = 0; i < a.outputs.size(); ++i) {
    REQUIRE(a.outputs[i].canonical == b.outputs[i].canonical);
    REQUIRE(a.outputs[i].e_pure == b.outputs[i].e_pure);
    REQUIRE(a.outputs[i].e_obj == b.outputs[i].e_obj);
    REQUIRE(a.outputs[i].D == b.outputs[i].D);
  }
  RunConfig cfg8 = cfg;
  cfg8.seed = 8;
  const RunOutput c = run_evolution(cfg8, ctx);
  bool differs = (c.outputs.size() != a.outputs.size());
  if (!differs)
    for (size_t i = 0; i < a.outputs.size() && !differs; ++i)
      differs = !(a.outputs[i].canonical == c.outputs[i].canonical);
  for (size_t i = 0; i < std::min(a.log.size(), c.log.size()) && !differs; ++i)
    differs = (a.log[i].best_e_obj != c.log[i].best_e_obj) || (a.log[i].killed != c.log[i].killed);
  REQUIRE(differs);
}

TEST_CASE("E_floor: bound computed in RunContext, matches recompute; K=2 full = 2x diag; e_tol auto",
          "[evolution][floor]") {
  const RunConfig cfg = sc27_cfg();
  const RunContext ctx = RunContext::build(cfg);
  REQUIRE(ctx.e_floor > 0.0);
  double ref = 0.0;
  for (int n = 0; n < ctx.zones.n_shells; ++n)
    for (size_t t = 0; t < cfg.counts.size(); ++t) {
      const double P = double(cfg.counts[t]) * ctx.zones.coord_num[static_cast<size_t>(n)];
      const double v = cfg.x_achieved[t] * P;
      ref += ctx.weights[static_cast<size_t>(n)] * std::abs(v - std::round(v)) / P;
    }
  REQUIRE(ctx.e_floor == Approx(ref).epsilon(1e-13));
  const double ff = e_floor_full(ctx.zones, cfg.counts, cfg.x_achieved, ctx.weights);
  REQUIRE(ff == Approx(2.0 * ctx.e_floor).epsilon(1e-12));
  RunConfig c2 = cfg;
  c2.e_tol = -1.0;
  REQUIRE(effective_e_tol(c2, ctx) == Approx(3.0 * ctx.e_floor));
  c2.e_tol = 0.5;
  REQUIRE(effective_e_tol(c2, ctx) == Approx(0.5));
}
