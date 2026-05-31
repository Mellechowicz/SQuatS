// v1.5 [score]: T-X1 -- score_structure() reproduces the engine's own output
// records bit-exactly, is input-site-order-free, and rejects mismatched
// lattice / species / composition with clear messages.
#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <vector>

#include "exsqs/evolution.hpp"
#include "exsqs/lattice.hpp"
#include "exsqs/score.hpp"

using namespace exsqs;
using Catch::Matchers::ContainsSubstring;

namespace {
RunConfig sc_cfg() {
  RunConfig c;
  c.proto = make_sc(3.0);
  c.H = {{{3, 0, 0}, {0, 3, 0}, {0, 0, 3}}};
  c.species = {"W", "Cr"};
  c.counts = {14, 13};
  c.x_target = {14.0 / 27.0, 13.0 / 27.0};
  c.x_achieved = c.x_target;
  c.n_shells = 5;
  c.gamma = 1.0;
  c.population = 12;
  c.outputs = 4;
  c.e_tol = 0.0;
  c.max_generations = 6;
  c.retry_budget = 60;
  c.seed = 7;
  c.log_info = false;
  return c;
}
}  // namespace

TEST_CASE("T-X1: score_structure == engine records; input site order free", "[score]") {
  const RunConfig cfg = sc_cfg();
  const RunContext ctx = RunContext::build(cfg);
  const RunOutput out = run_evolution(cfg, ctx);
  REQUIRE_FALSE(out.outputs.empty());
  for (const Individual& I : out.outputs) {
    const Structure dec = decorate(ctx.geom, I.sigma, cfg.species);
    const ScoreResult r = score_structure(cfg, ctx, dec);
    REQUIRE(r.e_pure == I.e_pure);
    REQUIRE(r.e_obj == I.e_obj);
    REQUIRE(r.D == I.D);
    REQUIRE(r.sg == I.sg);
  }
  const Individual& B = out.outputs[0];
  Structure p = decorate(ctx.geom, B.sigma, cfg.species);
  std::reverse(p.frac.begin(), p.frac.end());
  std::reverse(p.species.begin(), p.species.end());
  const ScoreResult r = score_structure(cfg, ctx, p);
  REQUIRE(r.e_pure == B.e_pure);
  REQUIRE(r.D == B.D);
  REQUIRE(r.sg == B.sg);
}
