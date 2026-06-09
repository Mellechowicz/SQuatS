// v1.5 [schedule]: [A11] geometric beta schedule -- parsing/validation, the
// growth = 1.0 == const bitwise identity, and deterministic divergence for
// growth > 1.
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <fstream>

#include "exsqs/evolution.hpp"
#include "exsqs/lattice.hpp"

using namespace exsqs;
using Catch::Approx;
using Catch::Matchers::ContainsSubstring;

namespace {
RunConfig base_cfg() {
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
  c.max_generations = 10;
  c.retry_budget = 60;
  c.metropolis = true;
  c.beta = 50.0;
  c.seed = 7;
  c.log_info = false;
  return c;
}

bool outputs_equal(const RunOutput& a, const RunOutput& b) {
  if (a.evals != b.evals || a.outputs.size() != b.outputs.size()) return false;
  for (size_t i = 0; i < a.outputs.size(); ++i)
    if (a.outputs[i].canonical != b.outputs[i].canonical ||
        a.outputs[i].e_pure != b.outputs[i].e_pure)
      return false;
  return true;
}
}  // namespace

TEST_CASE("[A11] schedule parsing and validation", "[schedule][config]") {
  const std::string path = "/tmp/exsqs_sched.yaml";
  {
    std::ofstream f(path);
    f << "lattice: {type: sc, a: 3.0}\nsupercell: {diag: [3,3,3]}\n"
         "composition: {W: 0.52, Cr: 0.48}\n"
         "evolution:\n"
         "  survival: {mode: metropolis, beta: 50, schedule: geometric, beta_growth: 1.3}\n";
  }
  RunConfig c = load_config(path);
  REQUIRE(c.metropolis);
  REQUIRE(c.beta == Approx(50.0));
  REQUIRE(c.beta_schedule == 1);
  REQUIRE(c.beta_growth == Approx(1.3));
  REQUIRE_THROWS_WITH(load_config(path, {"evolution.survival.mode=ratio"}),
                      ContainsSubstring("metropolis"));
  REQUIRE_THROWS_WITH(load_config(path, {"evolution.survival.beta_growth=0"}),
                      ContainsSubstring("beta_growth"));
  REQUIRE_THROWS_WITH(load_config(path, {"evolution.survival.beta=auto"}),
                      ContainsSubstring("numeric beta"));
}

TEST_CASE("geometric growth 1.0 == const bitwise; growth > 1 diverges deterministically",
          "[schedule]") {
  const RunConfig a = base_cfg();
  const RunContext ctx = RunContext::build(a);
  const RunOutput ra = run_evolution(a, ctx);

  RunConfig b = base_cfg();
  b.beta_schedule = 1;
  b.beta_growth = 1.0;
  const RunOutput rb = run_evolution(b, ctx);
  REQUIRE(outputs_equal(ra, rb));

  RunConfig c = base_cfg();
  c.beta_schedule = 1;
  c.beta_growth = 1.5;
  const RunOutput rc = run_evolution(c, ctx);
  REQUIRE_FALSE(outputs_equal(ra, rc));
  const RunOutput rc2 = run_evolution(c, ctx);
  REQUIRE(outputs_equal(rc, rc2));
}
