#include <catch2/catch_test_macros.hpp>

#include <fstream>
#include <map>
#include <string>

#include "exsqs/displacements.hpp"
#include "exsqs/lattice.hpp"
#include "exsqs/structure.hpp"
#include "exsqs/testcells.hpp"

using namespace exsqs;

TEST_CASE("displacement sanity on pure cubic cells", "[displacements]") {
  // m-3m site symmetry: one direction, minus related by inversion -> D = 1
  REQUIRE(displacement_count(make_supercell_diag(make_bcc(3.165), 2, 2, 2, {"W"})) == 1);
  REQUIRE(displacement_count(make_supercell_diag(make_fcc(4.05), 2, 2, 2, {"W"})) == 1);
  REQUIRE(displacement_count(make_supercell_diag(make_sc(3.0), 2, 2, 2, {"W"})) == 1);
}

TEST_CASE("T-D1 displacement counts equal phonopy (gate)", "[displacements][gate]") {
  std::ifstream f(std::string(EXSQS_TEST_DATA_DIR) + "/dcounts_ref.txt");
  if (!f) SKIP("dcounts_ref.txt missing: run exsqs_gentests and tools/py/gen_reference.py");
  std::map<std::string, int> ref;
  std::string name;
  int d = 0;
  while (f >> name >> d) ref[name] = d;

  const auto cells = gate_cells();
  REQUIRE(cells.size() >= 12);
  for (const auto& c : cells) {
    INFO(c.name);
    REQUIRE(ref.count(c.name) == 1);
    REQUIRE(displacement_count(c.s) == ref.at(c.name));
  }
}
