#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <random>
#include <string>
#include <vector>

#include "exsqs/correlation.hpp"
#include "exsqs/lattice.hpp"
#include "exsqs/structure.hpp"
#include "exsqs/zones.hpp"

using namespace exsqs;
using Catch::Approx;

namespace {
const std::vector<std::string> kWCr = {"W", "Cr"};
}

TEST_CASE("T-C1 B2 analytic correlations", "[correlation]") {
  auto g = make_supercell_diag(make_bcc(3.165), 1, 1, 1, kWCr);
  auto s = decorate(g, {0, 1}, kWCr);
  auto zt = build_zones(s, 2);
  auto cd = count_pairs(s, zt);

  REQUIRE(cd.pi(0, 0, 0) == Approx(0.0).margin(1e-15));
  REQUIRE(cd.pi(0, 1, 1) == Approx(0.0).margin(1e-15));
  REQUIRE(cd.pi(1, 0, 0) == Approx(1.0).epsilon(1e-12));
  REQUIRE(cd.pi(1, 1, 1) == Approx(1.0).epsilon(1e-12));

  const auto w = make_weights(WeightForm::InvN, zt);
  const auto x = s.concentrations();
  REQUIRE(e_pure_diagonal(cd, x, w) == Approx(1.5).epsilon(1e-12));
  REQUIRE(e_pure_full(cd, x, w) == Approx(3.0).epsilon(1e-12));
}

TEST_CASE("T-C2 L1_0 and L1_2 analytic correlations", "[correlation]") {
  SECTION("L1_0") {
    auto g = make_supercell_diag(make_fcc(4.05), 1, 1, 1, kWCr);
    auto s = decorate(g, {0, 0, 1, 1}, kWCr);
    auto zt = build_zones(s, 2);
    auto cd = count_pairs(s, zt);
    REQUIRE(cd.pi(0, 0, 0) == Approx(1.0 / 3.0).epsilon(1e-12));
    REQUIRE(cd.pi(0, 1, 1) == Approx(1.0 / 3.0).epsilon(1e-12));
    REQUIRE(cd.pi(1, 0, 0) == Approx(1.0).epsilon(1e-12));
    REQUIRE(cd.pi(1, 1, 1) == Approx(1.0).epsilon(1e-12));
    const auto w = make_weights(WeightForm::InvN, zt);
    REQUIRE(e_pure_diagonal(cd, s.concentrations(), w) == Approx(5.0 / 6.0).epsilon(1e-12));
  }
  SECTION("L1_2") {
    auto g = make_supercell_diag(make_fcc(4.05), 1, 1, 1, kWCr);
    auto s = decorate(g, {0, 1, 1, 1}, kWCr);
    auto zt = build_zones(s, 2);
    auto cd = count_pairs(s, zt);
    REQUIRE(cd.pi(0, 0, 0) == Approx(0.0).margin(1e-15));
    REQUIRE(cd.pi(0, 1, 1) == Approx(2.0 / 3.0).epsilon(1e-12));
    REQUIRE(cd.pi(1, 0, 0) == Approx(1.0).epsilon(1e-12));
    REQUIRE(cd.pi(1, 1, 1) == Approx(1.0).epsilon(1e-12));
    const auto w = make_weights(WeightForm::InvN, zt);
    REQUIRE(e_pure_diagonal(cd, s.concentrations(), w) == Approx(5.0 / 6.0).epsilon(1e-12));
  }
}

TEST_CASE("T-C3 hypergeometric mean of random canonical decorations", "[correlation]") {
  auto g = make_supercell_diag(make_bcc(3.165), 2, 2, 2, kWCr);
  auto zt = build_zones(g, 5);
  const int N = g.natoms();
  const int nB = 5, nA = N - nB;
  const int M = 10000;

  std::vector<int> sig(static_cast<size_t>(N), 0);
  std::fill(sig.begin(), sig.begin() + nB, 1);
  std::mt19937_64 rng(12345);

  std::vector<double> sumA(static_cast<size_t>(zt.n_shells), 0.0),
      sumsqA(static_cast<size_t>(zt.n_shells), 0.0);
  std::vector<double> sumB(static_cast<size_t>(zt.n_shells), 0.0),
      sumsqB(static_cast<size_t>(zt.n_shells), 0.0);
  for (int m = 0; m < M; ++m) {
    std::shuffle(sig.begin(), sig.end(), rng);
    auto s = decorate(g, sig, kWCr);
    auto cd = count_pairs(s, zt);
    for (int n = 0; n < zt.n_shells; ++n) {
      const double pa = cd.pi(n, 0, 0), pb = cd.pi(n, 1, 1);
      sumA[static_cast<size_t>(n)] += pa;
      sumsqA[static_cast<size_t>(n)] += pa * pa;
      sumB[static_cast<size_t>(n)] += pb;
      sumsqB[static_cast<size_t>(n)] += pb * pb;
    }
  }
  const double expA = static_cast<double>(nA - 1) / (N - 1);
  const double expB = static_cast<double>(nB - 1) / (N - 1);
  for (int n = 0; n < zt.n_shells; ++n) {
    INFO("shell " << n);
    const double mA = sumA[static_cast<size_t>(n)] / M;
    const double vA = sumsqA[static_cast<size_t>(n)] / M - mA * mA;
    const double seA = std::sqrt(std::max(vA, 0.0) / M);
    REQUIRE(std::abs(mA - expA) <= 3.0 * seA + 1e-12);
    const double mB = sumB[static_cast<size_t>(n)] / M;
    const double vB = sumsqB[static_cast<size_t>(n)] / M - mB * mB;
    const double seB = std::sqrt(std::max(vB, 0.0) / M);
    REQUIRE(std::abs(mB - expB) <= 3.0 * seB + 1e-12);
  }
}

TEST_CASE("T-C4 K=2 identity: E_full == 2 * E_diag", "[correlation]") {
  auto g = make_supercell_diag(make_bcc(3.165), 2, 2, 2, kWCr);
  auto zt = build_zones(g, 5);
  const auto w = make_weights(WeightForm::InvN, zt);
  std::mt19937_64 rng(777);
  for (int trial = 0; trial < 200; ++trial) {
    const int nB = 1 + static_cast<int>(rng() % 15);
    std::vector<int> sig(static_cast<size_t>(g.natoms()), 0);
    std::fill(sig.begin(), sig.begin() + nB, 1);
    std::shuffle(sig.begin(), sig.end(), rng);
    auto s = decorate(g, sig, kWCr);
    auto cd = count_pairs(s, zt);
    for (int n = 0; n < zt.n_shells; ++n) REQUIRE(cd.getC(n, 0, 1) == cd.getC(n, 1, 0));
    const auto x = s.concentrations();
    const double ed = e_pure_diagonal(cd, x, w);
    const double ef = e_pure_full(cd, x, w);
    REQUIRE(ef == Approx(2.0 * ed).epsilon(1e-12).margin(1e-15));
  }
}
