#include <catch2/catch_test_macros.hpp>
#include <numeric>
#include <vector>

#include "exsqs/rng.hpp"

using namespace exsqs;

TEST_CASE("CounterRng: identical keys give identical streams [A14]", "[rng]") {
  CounterRng a(42, 1, 3, 7, RngPurpose::Mutation);
  CounterRng b(42, 1, 3, 7, RngPurpose::Mutation);
  for (int i = 0; i < 16; ++i) REQUIRE(a() == b());
}

TEST_CASE("CounterRng: every key component separates streams [A14]", "[rng]") {
  const uint64_t ref = CounterRng(42, 1, 3, 7, RngPurpose::Mutation)();
  REQUIRE(CounterRng(43, 1, 3, 7, RngPurpose::Mutation)() != ref);
  REQUIRE(CounterRng(42, 2, 3, 7, RngPurpose::Mutation)() != ref);
  REQUIRE(CounterRng(42, 1, 4, 7, RngPurpose::Mutation)() != ref);
  REQUIRE(CounterRng(42, 1, 3, 8, RngPurpose::Mutation)() != ref);
  REQUIRE(CounterRng(42, 1, 3, 7, RngPurpose::ParentPick)() != ref);
}

TEST_CASE("CounterRng: uniform() in [0,1) with sane mean; below(n) bounded", "[rng]") {
  CounterRng r(123, 0, 0, 0, RngPurpose::Generic);
  double sum = 0;
  for (int i = 0; i < 10000; ++i) {
    const double u = r.uniform();
    REQUIRE(u >= 0.0);
    REQUIRE(u < 1.0);
    sum += u;
  }
  const double mean = sum / 10000.0;
  REQUIRE(mean > 0.48);
  REQUIRE(mean < 0.52);
  REQUIRE(r.below(1) == 0);
  for (uint64_t n : {2ULL, 3ULL, 10ULL, 1000ULL})
    for (int i = 0; i < 200; ++i) REQUIRE(r.below(n) < n);
}

TEST_CASE("CounterRng: shuffle is reproducible and slot-sensitive", "[rng]") {
  std::vector<int> v0(27), v1(27), v2(27);
  std::iota(v0.begin(), v0.end(), 0);
  v1 = v0;
  v2 = v0;
  CounterRng a(9, 0, 0, 5, RngPurpose::SeedInit);
  CounterRng b(9, 0, 0, 5, RngPurpose::SeedInit);
  CounterRng c(9, 0, 0, 6, RngPurpose::SeedInit);
  a.shuffle(v0.begin(), v0.end());
  b.shuffle(v1.begin(), v1.end());
  c.shuffle(v2.begin(), v2.end());
  REQUIRE(v0 == v1);
  REQUIRE(v0 != v2);
}
