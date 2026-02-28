#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <fstream>
#include <map>
#include <random>
#include <string>
#include <vector>

#include "exsqs/correlation.hpp"
#include "exsqs/dedup.hpp"
#include "exsqs/lattice.hpp"
#include "exsqs/structure.hpp"
#include "exsqs/symmetry.hpp"
#include "exsqs/testcells.hpp"
#include "exsqs/zones.hpp"

using namespace exsqs;

TEST_CASE("T-S1 known space groups", "[symmetry]") {
  for (const auto& c : gate_cells()) {
    if (c.expected_sg < 0) continue;
    INFO(c.name);
    REQUIRE(get_symmetry(c.s).sg_number == c.expected_sg);
  }
}

TEST_CASE("T-S1 cross-check vs external reference", "[symmetry]") {
  std::ifstream f(std::string(EXSQS_TEST_DATA_DIR) + "/sg_ref.txt");
  if (!f) SKIP("sg_ref.txt missing: run exsqs_gentests and tools/py/gen_reference.py");
  std::map<std::string, int> ref;
  std::string name;
  int sg = 0;
  while (f >> name >> sg) ref[name] = sg;
  for (const auto& c : gate_cells()) {
    INFO(c.name);
    REQUIRE(ref.count(c.name) == 1);
    REQUIRE(get_symmetry(c.s).sg_number == ref.at(c.name));
  }
}

TEST_CASE("T-S2 canonical form: invariance and discrimination", "[dedup]") {
  auto g = make_supercell_diag(make_sc(3.0), 3, 3, 3, {"W", "Cr"});
  const auto perms = site_permutations(g);
  REQUIRE(perms.size() > 1);
  const int N = g.natoms();
  std::mt19937_64 rng(2024);

  SECTION("invariance under the empty-cell group action") {
    for (int trial = 0; trial < 500; ++trial) {
      std::vector<int> sig(static_cast<size_t>(N), 0);
      std::fill(sig.begin(), sig.begin() + 13, 1);
      std::shuffle(sig.begin(), sig.end(), rng);
      const auto& p = perms[static_cast<size_t>(rng() % perms.size())];
      std::vector<int> sig2(static_cast<size_t>(N));
      for (int i = 0; i < N; ++i)
        sig2[static_cast<size_t>(p[static_cast<size_t>(i)])] = sig[static_cast<size_t>(i)];
      REQUIRE(canonical_labels(sig, perms) == canonical_labels(sig2, perms));
    }
  }
  SECTION("idempotence: canonical of canonical is itself") {
    for (int trial = 0; trial < 100; ++trial) {
      std::vector<int> sig(static_cast<size_t>(N), 0);
      std::fill(sig.begin(), sig.begin() + 13, 1);
      std::shuffle(sig.begin(), sig.end(), rng);
      const auto can = canonical_labels(sig, perms);
      const std::vector<int> sigc(can.begin(), can.end());
      REQUIRE(canonical_labels(sigc, perms) == can);
    }
  }
  SECTION("discrimination: invariant-differing decorations never collide") {
    auto zt = build_zones(g, 3);
    auto invariant = [&](const std::vector<int>& sig) {
      return count_pairs(decorate(g, sig, {"W", "Cr"}), zt).C;
    };
    int checked = 0;
    for (int trial = 0; trial < 2000; ++trial) {
      std::vector<int> a(static_cast<size_t>(N), 0);
      std::fill(a.begin(), a.begin() + 13, 1);
      std::vector<int> b = a;
      std::shuffle(a.begin(), a.end(), rng);
      std::shuffle(b.begin(), b.end(), rng);
      if (invariant(a) == invariant(b)) continue;  // possibly equivalent; skip
      ++checked;
      REQUIRE(canonical_labels(a, perms) != canonical_labels(b, perms));
    }
    REQUIRE(checked > 1000);
  }
}
