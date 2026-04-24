// v1.3 checkpoint/restart tests: serializer round trip, T-K1 (resume after a
// budget stop is bit-exact vs an uninterrupted run, with a migration boundary
// crossing the split), T-K2 (signature mismatch refused; budget raise allowed).
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>
#include <string>
#include <vector>

#include "exsqs/evolution.hpp"
#include "exsqs/island_engine.hpp"
#include "exsqs/lattice.hpp"
#include "exsqs/serialize.hpp"

using namespace exsqs;

namespace {

RunConfig ck_cfg() {
  RunConfig c;
  c.proto = make_sc(3.0);
  c.H = {{{3, 0, 0}, {0, 3, 0}, {0, 0, 3}}};
  c.species = {"W", "Cr"};
  c.counts = {14, 13};
  c.x_target = {14.0 / 27.0, 13.0 / 27.0};
  c.x_achieved = {14.0 / 27.0, 13.0 / 27.0};
  c.n_shells = 5;
  c.gamma = 1.0;
  c.population = 12;
  c.outputs = 4;
  c.e_tol = 0.0;  // unreachable
  c.retry_budget = 60;
  c.islands = 2;
  c.migration_every = 3;  // migrations at rounds 3, 6, 9, 12 (split crosses 6)
  c.migrants = 1;
  c.seed = 7;
  c.log_info = false;
  return c;
}

bool same_output(const RunOutput& a, const RunOutput& b) {
  if (a.success != b.success || a.evals != b.evals) return false;
  if (a.log.size() != b.log.size() || a.outputs.size() != b.outputs.size()) return false;
  for (size_t i = 0; i < a.log.size(); ++i) {
    const auto& x = a.log[i];
    const auto& y = b.log[i];
    if (x.island != y.island || x.gen != y.gen || x.best_e_pure != y.best_e_pure ||
        x.best_e_obj != y.best_e_obj || x.median_e_obj != y.median_e_obj ||
        x.best_D != y.best_D || x.best_sg != y.best_sg || x.killed != y.killed ||
        x.evals != y.evals || x.dup_rejected != y.dup_rejected ||
        x.p1_rejected != y.p1_rejected || x.fallback_seeds != y.fallback_seeds)
      return false;
  }
  for (size_t i = 0; i < a.outputs.size(); ++i)
    if (a.outputs[i].canonical != b.outputs[i].canonical ||
        a.outputs[i].e_pure != b.outputs[i].e_pure || a.outputs[i].D != b.outputs[i].D ||
        a.outputs[i].sg != b.outputs[i].sg)
      return false;
  if (a.island_generations != b.island_generations) return false;
  if (a.island_migrants_in != b.island_migrants_in) return false;
  return true;
}

}  // namespace

TEST_CASE("v1.3 serializer: IslandResult/Individual bit-exact round trip", "[checkpoint]") {
  RunConfig c = ck_cfg();
  c.islands = 1;
  c.migrants = 0;
  c.max_generations = 4;
  const RunContext ctx = RunContext::build(c);
  const IslandResult r = evolve_island(c, ctx, 0);
  ByteWriter w;
  put_island_result(w, r);
  ByteReader rd(w.data());
  const IslandResult r2 = get_island_result(rd);
  REQUIRE(rd.remaining() == 0);
  REQUIRE(r2.pool.size() == r.pool.size());
  for (size_t i = 0; i < r.pool.size(); ++i) {
    const Individual& a = r.pool[i];
    const Individual& b = r2.pool[i];
    REQUIRE(a.sigma == b.sigma);
    REQUIRE(a.canonical == b.canonical);
    REQUIRE(a.hash == b.hash);
    REQUIRE(a.sg == b.sg);
    REQUIRE(a.sg_symbol == b.sg_symbol);
    REQUIRE(a.eq_atoms == b.eq_atoms);
    REQUIRE(a.stab_ops.size() == b.stab_ops.size());
    for (size_t k = 0; k < a.stab_ops.size(); ++k) {
      REQUIRE(a.stab_ops[k].R == b.stab_ops[k].R);
      for (int j = 0; j < 3; ++j) REQUIRE(a.stab_ops[k].t[j] == b.stab_ops[k].t[j]);
    }
    REQUIRE(a.e_pure == b.e_pure);
    REQUIRE(a.e_obj == b.e_obj);
    REQUIRE(a.D == b.D);
    REQUIRE(a.birth_gen == b.birth_gen);
    REQUIRE(a.origin == b.origin);
  }
  REQUIRE(r2.log.size() == r.log.size());
  REQUIRE(r2.stop_reason == r.stop_reason);
  REQUIRE(r2.evals == r.evals);
}

TEST_CASE("T-K1: resume after a budget stop is bit-exact vs an uninterrupted run",
          "[checkpoint]") {
  const std::string dir = "/tmp/exsqs_ck_tk1";
  std::filesystem::remove_all(dir);
  RunConfig c = ck_cfg();
  c.max_generations = 12;
  const RunContext ctx = RunContext::build(c);
  const RunOutput straight = run_evolution(c, ctx);

  RunConfig c5 = c;
  c5.max_generations = 5;  // budget stop mid-run, before the round-6 migration
  (void)run_evolution(c5, ctx, {}, "", dir);
  const RunOutput resumed = run_evolution(c, ctx, {}, dir + "/state.ckpt", "");
  REQUIRE(same_output(straight, resumed));
}

TEST_CASE("T-K2: resume refuses a trajectory-signature mismatch; budget raise allowed",
          "[checkpoint]") {
  const std::string dir = "/tmp/exsqs_ck_tk2";
  std::filesystem::remove_all(dir);
  RunConfig c = ck_cfg();
  c.max_generations = 4;
  const RunContext ctx = RunContext::build(c);
  (void)run_evolution(c, ctx, {}, "", dir);

  RunConfig bad = c;
  bad.seed = 8;
  REQUIRE_THROWS_WITH(run_evolution(bad, ctx, {}, dir + "/state.ckpt", ""),
                      Catch::Matchers::ContainsSubstring("signature"));

  RunConfig ok = c;
  ok.max_generations = 9;  // raising a budget cap is the intended chain pattern
  const RunOutput more = run_evolution(ok, ctx, {}, dir + "/state.ckpt", "");
  REQUIRE_FALSE(more.log.empty());
  for (int g : more.island_generations) REQUIRE(g == 9);
}
