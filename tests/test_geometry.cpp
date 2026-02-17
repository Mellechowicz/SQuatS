#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <map>
#include <tuple>
#include <vector>

#include "exsqs/lattice.hpp"
#include "exsqs/structure.hpp"
#include "exsqs/zones.hpp"

using namespace exsqs;
using Catch::Approx;

TEST_CASE("T-G1 shell radii and coordination sequences", "[geometry]") {
  SECTION("bcc") {
    const double a = 3.165;
    auto g = make_supercell_diag(make_bcc(a), 2, 2, 2, {"W"});
    auto zt = build_zones(g, 5);
    REQUIRE(zt.n_shells == 5);
    const std::vector<double> r = {std::sqrt(3.0) / 2 * a, a, std::sqrt(2.0) * a,
                                   std::sqrt(11.0) / 2 * a, std::sqrt(3.0) * a};
    for (int n = 0; n < 5; ++n) REQUIRE(zt.radii[n] == Approx(r[n]).epsilon(1e-9));
    REQUIRE(zt.coord_num == std::vector<int>{8, 6, 12, 24, 8});
  }
  SECTION("fcc") {
    const double a = 4.05;
    auto g = make_supercell_diag(make_fcc(a), 2, 2, 2, {"W"});
    auto zt = build_zones(g, 5);
    const std::vector<double> r = {std::sqrt(0.5) * a, a, std::sqrt(1.5) * a, std::sqrt(2.0) * a,
                                   std::sqrt(2.5) * a};
    for (int n = 0; n < 5; ++n) REQUIRE(zt.radii[n] == Approx(r[n]).epsilon(1e-9));
    REQUIRE(zt.coord_num == std::vector<int>{12, 6, 24, 12, 24});
  }
  SECTION("sc") {
    const double a = 3.0;
    auto g = make_supercell_diag(make_sc(a), 3, 3, 3, {"W"});
    auto zt = build_zones(g, 5);
    const std::vector<double> r = {a, std::sqrt(2.0) * a, std::sqrt(3.0) * a, 2 * a,
                                   std::sqrt(5.0) * a};
    for (int n = 0; n < 5; ++n) REQUIRE(zt.radii[n] == Approx(r[n]).epsilon(1e-9));
    REQUIRE(zt.coord_num == std::vector<int>{6, 12, 8, 6, 24});
  }
  SECTION("hcp nearest neighbours") {
    auto g = make_supercell_diag(make_hcp(3.2), 2, 2, 1, {"Mg"});
    auto zt = build_zones(g, 1);
    REQUIRE(zt.coord_num[0] == 12);
  }
}

TEST_CASE("T-G2 pair table identities", "[geometry]") {
  auto g = make_supercell_diag(make_bcc(3.165), 2, 2, 2, {"W"});
  auto zt = build_zones(g, 5);
  const int N = g.natoms();

  SECTION("sum rule: sum_i n_i(shell) = N * Z_n") {
    std::vector<long long> tot(static_cast<size_t>(zt.n_shells), 0);
    for (int i = 0; i < N; ++i)
      for (const auto& nb : zt.nbrs[static_cast<size_t>(i)]) tot[static_cast<size_t>(nb.shell)]++;
    for (int n = 0; n < zt.n_shells; ++n)
      REQUIRE(tot[static_cast<size_t>(n)] == static_cast<long long>(N) * zt.coord_num[static_cast<size_t>(n)]);
  }
  SECTION("directional symmetry: count(i->j, n) == count(j->i, n)") {
    std::map<std::tuple<int, int, int>, long long> cnt;
    for (int i = 0; i < N; ++i)
      for (const auto& nb : zt.nbrs[static_cast<size_t>(i)]) cnt[{i, nb.j, nb.shell}]++;
    for (const auto& kv : cnt) {
      const auto [i, j, n] = kv.first;
      REQUIRE(cnt[{j, i, n}] == kv.second);
    }
  }
  SECTION("no zero-distance self pairs within the first five shells") {
    // 2x2x2 bcc: nearest self-image sits at 2a = 6.33 A > shell-5 radius 5.48 A
    for (int i = 0; i < N; ++i)
      for (const auto& nb : zt.nbrs[static_cast<size_t>(i)]) REQUIRE(nb.j != i);
  }
  SECTION("nonzero self-images are legitimate neighbours [A1]") {
    auto g1 = make_supercell_diag(make_sc(3.0), 1, 1, 1, {"W"});
    auto z1 = build_zones(g1, 1);
    REQUIRE(z1.coord_num[0] == 6);
    for (const auto& nb : z1.nbrs[0]) REQUIRE(nb.j == 0);
  }
}
