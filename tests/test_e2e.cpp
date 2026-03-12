// Integration tests on the paper's reference system: W70Cr30, bcc 4x4x4
// (-> W90Cr38, 128 sites). Tagged [e2e]; registered as a separate ctest entry.
//
// v1.1: acceptance is floor-relative (SPEC 4.1, 12) -- the paper's absolute
// error scale lies below the provable L1 quantization floor and is therefore
// not comparable. Recipe: 4 sequential islands of high-pressure metropolis
// (fixed beta), the serial (1+lambda) exploitation regime; pooled best [D8].
#include <catch2/catch_test_macros.hpp>
#include <set>
#include <vector>

#include "exsqs/evolution.hpp"
#include "exsqs/lattice.hpp"

using namespace exsqs;

namespace {
RunConfig w_cfg(double gamma, int gens) {
  RunConfig c;
  c.proto = make_bcc(3.165);
  c.H = {{{4, 0, 0}, {0, 4, 0}, {0, 0, 4}}};
  c.species = {"W", "Cr"};
  c.counts = {90, 38};  // largest-remainder of 0.70/0.30 on 128 sites [A5]
  c.x_target = {0.70, 0.30};
  c.x_achieved = {90.0 / 128.0, 38.0 / 128.0};
  c.n_shells = 7;
  c.gamma = gamma;
  c.population = 100;
  c.outputs = 5;
  c.e_tol = -1.0;  // auto = 3.0 * E_floor (v1.1, SPEC 4.1)
  c.max_generations = gens;
  c.islands = 4;
  c.metropolis = true;
  c.beta = 3000.0;  // fixed high pressure: serial (1+lambda) exploitation
  c.seed = 42;
  c.log_info = false;
  return c;
}
}  // namespace

TEST_CASE("T-E1: gamma=0 reaches E_pure <= 3.0*E_floor with >= M unique non-P1 outputs", "[e2e]") {
  const RunConfig cfg = w_cfg(0.0, 25);
  const RunContext ctx = RunContext::build(cfg);
  const RunOutput out = run_evolution(cfg, ctx);
  REQUIRE(out.success);
  REQUIRE(out.outputs.size() >= static_cast<size_t>(cfg.outputs));
  std::set<std::vector<uint8_t>> keys;
  for (const auto& I : out.outputs) {
    REQUIRE(I.sg > 1);
    keys.insert(I.canonical);
  }
  REQUIRE(keys.size() == out.outputs.size());
  REQUIRE(out.outputs[0].e_pure <= 3.0 * ctx.e_floor);
}

TEST_CASE("T-E2: gamma=1 best output reaches D <= D(P1)/3 at E_pure <= 3x the T-E1 bound",
          "[e2e]") {
  // Spec-default dynamics (ratio survival) populate the (D, E_pure) knee of the
  // gamma=1 Pareto front; the high-pressure regime polarizes to its extremes
  // (D=40 R-3m at high E vs D=416 Cm at low E -- see STEP2_REPORT).
  RunConfig cfg = w_cfg(1.0, 300);
  cfg.population = 64;
  cfg.islands = 1;
  cfg.metropolis = false;  // [A11] default: ratio survival
  cfg.beta = -1.0;
  cfg.e_tol = 1e-9;  // disable the success stop; T-E2 asserts on outputs only
  const RunContext ctx = RunContext::build(cfg);
  const RunOutput out = run_evolution(cfg, ctx);
  REQUIRE_FALSE(out.outputs.empty());
  const int D_p1 = 6 * ctx.geom.natoms();  // 768 for 128 sites
  REQUIRE(out.outputs[0].D * 3 <= D_p1);
  REQUIRE(out.outputs[0].e_pure <= 9.0 * ctx.e_floor);
}
