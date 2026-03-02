// load_config regression: SPEC section-11 parsing + CLI overrides. The engine
// tests construct RunConfig programmatically, so this covers the YAML path.
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdio>
#include <fstream>

#include "exsqs/config.hpp"

using namespace exsqs;
using Catch::Approx;

TEST_CASE("config: section-11 defaults, largest remainder, overrides", "[config]") {
  const std::string path = "/tmp/exsqs_test_cfg.yaml";
  {
    std::ofstream f(path);
    f << "lattice:    {type: bcc, a: 3.165}\n"
         "supercell:  {diag: [4, 4, 4]}\n"
         "composition: {W: 0.70, Cr: 0.30}\n"
         "evolution:\n"
         "  population: 200\n"
         "  survival: {mode: ratio}\n"
         "rng: {seed: 42}\n"
         "output: {dir: ./runs/x, log_level: info}\n";
  }

  SECTION("plain parse: frozen defaults + [A5] composition") {
    const RunConfig c = load_config(path);
    REQUIRE(c.species == std::vector<std::string>{"W", "Cr"});
    REQUIRE(c.counts == std::vector<int>{90, 38});  // paper supercell W90Cr38
    REQUIRE(c.x_achieved[0] == Approx(90.0 / 128.0).epsilon(1e-15));
    REQUIRE(c.n_shells == 7);
    REQUIRE(c.gamma == 1.0);
    REQUIRE(c.population == 200);
    REQUIRE(c.outputs == 10);
    REQUIRE(c.e_tol < 0);  // v1.1 default: auto (resolved to 3.0*E_floor at run start)
    REQUIRE(c.max_generations == 5000);
    REQUIRE_FALSE(c.metropolis);
    REQUIRE(c.elitism_best == 1);
    REQUIRE(c.mut_swaps == 1);
    REQUIRE(c.mut_sympres);
    REQUIRE(c.retry_budget == 100);
    REQUIRE(c.seed_mode == 2);
    REQUIRE_FALSE(c.full_pairs);  // auto, K=2 [A16]
    REQUIRE(c.symprec == Approx(1e-5));
    REQUIRE(c.seed == 42);
    REQUIRE(c.log_info);
  }

  SECTION("overrides rebind, not assign (yaml-cpp reset regression)") {
    const RunConfig c = load_config(
        path, {"error.gamma=0", "evolution.population=64", "evolution.max_generations=40",
               "output.dir=/tmp/exsqs_ovr", "evolution.survival.mode=metropolis"});
    REQUIRE(c.gamma == 0.0);
    REQUIRE(c.population == 64);
    REQUIRE(c.max_generations == 40);
    REQUIRE(c.outdir == "/tmp/exsqs_ovr");
    REQUIRE(c.metropolis);
    // the document must have survived every override
    REQUIRE(c.species.size() == 2);
    REQUIRE(c.counts == std::vector<int>{90, 38});
  }

  SECTION("validation: bad composition rejected") {
    const std::string bad = "/tmp/exsqs_bad_cfg.yaml";
    {
      std::ofstream f(bad);
      f << "lattice: {type: bcc, a: 3.0}\nsupercell: {diag: [1,1,1]}\n"
           "composition: {W: 0.95, Cr: 0.05}\n";  // Cr rounds to 0 on 2 sites
    }
    REQUIRE_THROWS(load_config(bad));
    std::remove(bad.c_str());
  }

  std::remove(path.c_str());
}

TEST_CASE("config v1.1: e_tol auto|number, poisson swaps, stagnation_stop", "[config]") {
  const std::string path = "/tmp/exsqs_test_cfg11.yaml";
  {
    std::ofstream f(path);
    f << "lattice: {type: sc, a: 3.0}\nsupercell: {diag: [3,3,3]}\n"
         "composition: {W: 0.52, Cr: 0.48}\n"
         "evolution:\n  e_tol: auto\n  stagnation_stop: 5\n"
         "  mutation: {swaps: poisson, lambda: 0.25}\n";
  }
  RunConfig c = load_config(path);
  REQUIRE(c.e_tol < 0);
  REQUIRE(c.stagnation_stop == 5);
  REQUIRE(c.mut_poisson_lambda == Approx(0.25));
  c = load_config(path, {"evolution.e_tol=7.5e-3", "evolution.mutation.swaps=2"});
  REQUIRE(c.e_tol == Approx(7.5e-3));
  REQUIRE(c.mut_poisson_lambda < 0);
  REQUIRE(c.mut_swaps == 2);
}
