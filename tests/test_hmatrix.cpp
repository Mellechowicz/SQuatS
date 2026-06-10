// v1.5 [supercell]: non-diagonal integer H (parsed since v0.2, first exercised
// here) -- counts from |det H|, zone/context build, sane empty-cell symmetry.
#include <catch2/catch_test_macros.hpp>
#include <fstream>
#include <vector>

#include "exsqs/evolution.hpp"
#include "exsqs/symmetry.hpp"

using namespace exsqs;

TEST_CASE("supercell.matrix: non-diagonal H builds a consistent run context",
          "[config][supercell]") {
  const std::string path = "/tmp/exsqs_hmat.yaml";
  {
    std::ofstream f(path);
    f << "lattice: {type: sc, a: 3.0}\n"
         "supercell: {matrix: [[2,1,0],[0,2,0],[0,0,2]]}\n"
         "composition: {W: 0.5, Cr: 0.5}\nzones: {n_shells: 3}\n";
  }
  const RunConfig c = load_config(path);
  REQUIRE(c.counts == std::vector<int>{4, 4});  // |det H| = 8 sites
  const RunContext ctx = RunContext::build(c);
  REQUIRE(ctx.geom.natoms() == 8);
  REQUIRE(ctx.zones.n_shells >= 1);
  const SymmetryInfo empty = get_symmetry(ctx.geom, c.symprec);
  REQUIRE(empty.sg_number > 1);  // undecorated sheared cell keeps symmetry
  REQUIRE(ctx.e_floor >= 0.0);
}
