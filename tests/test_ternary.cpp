// v1.4 (step 5): first K>=3 coverage. The spec's machinery (largest-remainder
// composition, full-pairs error [A16], the K-general floor of section 4.1,
// constructive seeding, orbit mutation) was written for K species but only
// ever exercised at K=2. These tests pin the ternary behavior, including the
// commensurate-composition case where E_floor is exactly zero.
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <set>
#include <vector>

#include "exsqs/correlation.hpp"
#include "exsqs/evolution.hpp"
#include "exsqs/lattice.hpp"

using namespace exsqs;
using Catch::Approx;

namespace {

RunConfig tern27_cfg() {
  RunConfig c;
  c.proto = make_sc(3.0);
  c.H = {{{3, 0, 0}, {0, 3, 0}, {0, 0, 3}}};
  c.species = {"W", "Mo", "Cr"};
  c.counts = {10, 9, 8};
  c.x_target = {10.0 / 27.0, 9.0 / 27.0, 8.0 / 27.0};
  c.x_achieved = c.x_target;
  c.n_shells = 5;
  c.full_pairs = true;  // [A16]: auto -> full for K >= 3
  c.gamma = 0.0;
  c.population = 12;
  c.outputs = 4;
  c.e_tol = 0.0;  // unreachable
  c.max_generations = 6;
  c.retry_budget = 60;
  c.seed = 7;
  c.log_info = false;
  return c;
}

RunConfig tern_bcc_cfg() {
  RunConfig c;
  c.proto = make_bcc(3.165);
  c.H = {{{4, 0, 0}, {0, 4, 0}, {0, 0, 4}}};
  c.species = {"W", "Mo", "Cr"};
  c.counts = {64, 32, 32};
  c.x_target = {0.50, 0.25, 0.25};
  c.x_achieved = {64.0 / 128.0, 32.0 / 128.0, 32.0 / 128.0};
  c.n_shells = 7;
  c.full_pairs = true;  // [A16]
  c.gamma = 0.0;
  c.population = 64;
  c.outputs = 5;
  c.e_tol = 1e-12;  // absolute metric: E_floor is exactly 0 here
  c.max_generations = 20;
  c.islands = 2;
  c.metropolis = true;
  c.beta = 3000.0;
  c.seed = 42;
  c.log_info = false;
  return c;
}

std::vector<int> species_counts(const std::vector<int>& sigma, int K) {
  std::vector<int> n(static_cast<size_t>(K), 0);
  for (int t : sigma) n[static_cast<size_t>(t)]++;
  return n;
}

}  // namespace

TEST_CASE("K=3 config: largest remainder, mode auto -> full_pairs [A16]", "[config][ternary]") {
  const std::string path = "/tmp/exsqs_tern_cfg.yaml";
  {
    std::ofstream f(path);
    f << "lattice: {type: bcc, a: 3.165}\nsupercell: {diag: [4,4,4]}\n"
         "composition: {W: 0.50, Mo: 0.25, Cr: 0.25}\n";
  }
  RunConfig c = load_config(path);
  REQUIRE(c.species == std::vector<std::string>{"W", "Mo", "Cr"});
  REQUIRE(c.counts == std::vector<int>{64, 32, 32});
  REQUIRE(c.full_pairs);  // auto at K=3
  c = load_config(path, {"error.mode=diagonal"});
  REQUIRE_FALSE(c.full_pairs);

  {
    std::ofstream f(path);
    f << "lattice: {type: sc, a: 3.0}\nsupercell: {diag: [3,3,3]}\n"
         "composition:\n  A: 0.3333333333333333\n  B: 0.3333333333333333\n"
         "  C: 0.3333333333333333\n";
  }
  c = load_config(path);
  REQUIRE(c.counts == std::vector<int>{9, 9, 9});
}

TEST_CASE("K=3 floor: full-pairs recompute; commensurate composition has floor exactly 0",
          "[evolution][floor][ternary]") {
  const RunConfig cfg = tern27_cfg();
  const RunContext ctx = RunContext::build(cfg);
  REQUIRE(ctx.e_floor > 0.0);
  double full = 0.0, diag = 0.0;
  for (int n = 0; n < ctx.zones.n_shells; ++n)
    for (size_t t = 0; t < cfg.counts.size(); ++t) {
      const double P = double(cfg.counts[t]) * ctx.zones.coord_num[static_cast<size_t>(n)];
      for (size_t u = 0; u < cfg.counts.size(); ++u) {
        const double v = cfg.x_achieved[u] * P;
        const double term = ctx.weights[static_cast<size_t>(n)] * std::abs(v - std::round(v)) / P;
        full += term;
        if (u == t) diag += term;
      }
    }
  REQUIRE(ctx.e_floor == Approx(full).epsilon(1e-13));
  REQUIRE(full > diag);  // K=3: off-diagonal terms contribute (unlike K=2's 2x identity)
  REQUIRE(diag == Approx(e_floor_diagonal(ctx.zones, cfg.counts, cfg.x_achieved, ctx.weights))
                      .epsilon(1e-13));

  // dyadic x on bcc: every x*P integral => the bound vanishes; a perfect SQS
  // is not excluded and floor-relative gates degenerate (absolute e_tol needed)
  const RunConfig big = tern_bcc_cfg();
  const RunContext bctx = RunContext::build(big);
  REQUIRE(bctx.e_floor == 0.0);
}

TEST_CASE("K=3 engine: composition preserved, non-P1, deterministic", "[evolution][ternary]") {
  const RunConfig cfg = tern27_cfg();
  const RunContext ctx = RunContext::build(cfg);
  const RunOutput a = run_evolution(cfg, ctx);
  REQUIRE_FALSE(a.outputs.empty());
  std::set<std::vector<uint8_t>> keys;
  for (const auto& I : a.outputs) {
    REQUIRE(I.sg > 1);
    REQUIRE(species_counts(I.sigma, 3) == cfg.counts);
    keys.insert(I.canonical);
  }
  REQUIRE(keys.size() == a.outputs.size());
  const RunOutput b = run_evolution(cfg, ctx);
  REQUIRE(a.evals == b.evals);
  REQUIRE(a.outputs.size() == b.outputs.size());
  for (size_t i = 0; i < a.outputs.size(); ++i) {
    REQUIRE(a.outputs[i].canonical == b.outputs[i].canonical);
    REQUIRE(a.outputs[i].e_pure == b.outputs[i].e_pure);
  }
}

TEST_CASE("T-E3: ternary bcc W64Mo32Cr32 (full_pairs) reaches E_pure <= 1.2e-1, C-class output",
          "[e2e][ternary]") {
  const RunConfig cfg = tern_bcc_cfg();
  const RunContext ctx = RunContext::build(cfg);
  REQUIRE(ctx.e_floor == 0.0);
  const RunOutput out = run_evolution(cfg, ctx);
  REQUIRE(out.outputs.size() >= static_cast<size_t>(cfg.outputs));
  const Individual& best = out.outputs[0];
  REQUIRE(best.sg > 1);
  REQUIRE(best.e_pure <= 1.2e-1);  // seed-pinned: 1.1125e-1 observed at gen 20
  REQUIRE(best.D < 6 * ctx.geom.natoms());
  REQUIRE(species_counts(best.sigma, 3) == cfg.counts);
}
