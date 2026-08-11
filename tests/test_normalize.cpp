// SPEC 4.3 gates: the analytic E[random] of the pair sector and the
// pair-only normalization (error.normalize).
//   T-N1: the closed formula against an independent in-test recompute, and
//         the K = 2 identity full == 2 x diagonal.
//   T-N2: composition of the normalized objective -- e_pure equals the raw
//         pair error divided by E[random], plus the RAW multiplet sectors.
//   T-N3: the flag enters the trajectory signature (and nothing else does),
//         and the normalized trajectory stays thread-invariant.
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <string>
#include <vector>

#include "exsqs/config.hpp"
#include "exsqs/correlation.hpp"
#include "exsqs/evolution.hpp"
#include "exsqs/lattice.hpp"
#include "exsqs/score.hpp"
#include "exsqs/serialize.hpp"
#include "exsqs/structure.hpp"
#include "exsqs/zones.hpp"

using namespace exsqs;
using Catch::Approx;

namespace {

RunConfig norm_cfg() {
  RunConfig c;
  c.proto = make_bcc(3.165);
  c.H = {{{3, 0, 0}, {0, 3, 0}, {0, 0, 3}}};
  c.species = {"W", "Cr"};
  c.counts = {38, 16};
  c.x_target = {38.0 / 54, 16.0 / 54};
  c.x_achieved = c.x_target;
  c.n_shells = 5;
  c.gamma = 0.0;
  c.population = 12;
  c.outputs = 4;
  c.e_tol = 1e-12;
  c.max_generations = 6;
  c.retry_budget = 60;
  c.seed = 20260718;
  c.log_info = false;
  return c;
}

}  // namespace

TEST_CASE("T-N1 analytic E[random]: formula and K=2 identity", "[normalize]") {
  auto g = make_supercell_diag(make_bcc(3.165), 3, 3, 3, {"W", "Cr"});
  auto zt = build_zones(g, 5);
  const std::vector<int> counts = {38, 16};
  const std::vector<double> x = {38.0 / 54, 16.0 / 54};
  const auto w = make_weights(WeightForm::InvN, zt);

  // independent recompute: half-normal MAD of a binomial proportion,
  // E|p - x| = sqrt(2 x (1 - x) / (pi P)), P = N_t Z_n
  double expect = 0.0;
  for (int n = 0; n < zt.n_shells; ++n)
    for (size_t t = 0; t < counts.size(); ++t) {
      const double P = counts[t] * static_cast<double>(zt.coord_num[static_cast<size_t>(n)]);
      expect += w[static_cast<size_t>(n)] *
                std::sqrt(2.0 * x[t] * (1.0 - x[t]) / (M_PI * P));
    }
  REQUIRE(analytic_e_random_diagonal(zt, counts, x, w) == Approx(expect).epsilon(1e-12));
  // K = 2: every (t, u) term repeats the same x(1-x), so full = 2 x diagonal
  REQUIRE(analytic_e_random_full(zt, counts, x, w) ==
          Approx(2.0 * analytic_e_random_diagonal(zt, counts, x, w)).epsilon(1e-12));
}

TEST_CASE("T-N2 pair-only composition of the normalized objective", "[normalize]") {
  RunConfig off = norm_cfg();
  off.lambda3 = 0.7;      // sectors ON in both configs
  off.mshell3 = 2;
  RunConfig on = off;
  on.normalize_epure = true;

  const RunContext cx_off = RunContext::build(off);
  const RunContext cx_on = RunContext::build(on);
  REQUIRE(cx_on.e_random == Approx(cx_off.e_random).epsilon(1e-15));
  REQUIRE(cx_on.e_random > 0.0);

  // one fixed decoration scored under both configs
  const auto sigma = seed_sigma_rejection(off, cx_off, 0, 0);
  const Structure dec = decorate(cx_off.geom, sigma, off.species);
  const ScoreResult r_off = score_structure(off, cx_off, dec);
  const ScoreResult r_on = score_structure(on, cx_on, dec);

  REQUIRE(r_on.e_pair == Approx(r_off.e_pair).epsilon(1e-15));  // raw pair kept
  REQUIRE(r_on.e3 == Approx(r_off.e3).epsilon(1e-15));          // sectors raw
  const double composed =
      r_off.e_pair / cx_on.e_random + on.lambda3 * r_off.e3;
  REQUIRE(r_on.e_pure == Approx(composed).epsilon(1e-12));
  // the floor composes the same way
  const double floor_composed =
      cx_off.e_floor_pair / cx_on.e_random + on.lambda3 * cx_on.clusters.floor3;
  REQUIRE(cx_on.e_floor == Approx(floor_composed).epsilon(1e-12));
}

TEST_CASE("T-N3 signature coverage and thread invariance", "[normalize]") {
  RunConfig a = norm_cfg();
  RunConfig b = a;
  b.normalize_epure = true;
  REQUIRE(trajectory_signature(a) != trajectory_signature(b));
  REQUIRE(trajectory_signature(b) == trajectory_signature(b));

  const RunContext ctx = RunContext::build(b);
  b.omp_threads = 1;
  const RunOutput r1 = run_evolution(b, ctx);
  b.omp_threads = 4;
  const RunOutput r4 = run_evolution(b, ctx);
  REQUIRE(r1.evals == r4.evals);
  REQUIRE(r1.outputs.size() == r4.outputs.size());
  for (size_t i = 0; i < r1.outputs.size(); ++i) {
    REQUIRE(r1.outputs[i].canonical == r4.outputs[i].canonical);
    REQUIRE(r1.outputs[i].e_pure == r4.outputs[i].e_pure);  // exact
  }
}
