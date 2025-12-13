#pragma once
// Counter-based RNG [A14]: stateless streams keyed by
// (master_seed, island, generation, slot, purpose). Trajectories are therefore
// independent of thread count and scheduling, and reruns are bit-identical.
//
// Design: SplitMix64-class key derivation + counter mixing (same finalizer
// family as PCG's output permutations). Statistically strong for evolutionary
// search; not cryptographic. Pure 64-bit integer arithmetic => portable and
// bit-reproducible across platforms/compilers (GCC/Clang: uses __int128 for
// the unbiased-enough bounded draw).

#include <cstdint>
#include <utility>

namespace exsqs {

enum class RngPurpose : uint64_t {
  SeedInit = 1,
  ExtinctionDraw = 2,
  ParentPick = 3,
  Mutation = 4,
  ConstructiveSeed = 5,
  Generic = 6,
};

namespace detail {
inline uint64_t mix64(uint64_t x) {
  x ^= x >> 30;
  x *= 0xbf58476d1ce4e5b9ULL;
  x ^= x >> 27;
  x *= 0x94d049bb133111ebULL;
  x ^= x >> 31;
  return x;
}
}  // namespace detail

class CounterRng {
 public:
  using result_type = uint64_t;

  CounterRng(uint64_t master_seed, uint64_t island, uint64_t generation, uint64_t slot,
             RngPurpose purpose) {
    uint64_t k = detail::mix64(master_seed + 0x9e3779b97f4a7c15ULL);
    k = detail::mix64(k ^ (island * 0xd1342543de82ef95ULL + 0x2545f4914f6cdd1dULL));
    k = detail::mix64(k ^ (generation * 0xda942042e4dd58b5ULL + 0x9e6c63d0876a9a67ULL));
    k = detail::mix64(k ^ (slot * 0xc2b2ae3d27d4eb4fULL + 0x165667b19e3779f9ULL));
    key_ = detail::mix64(k ^ (static_cast<uint64_t>(purpose) * 0x27d4eb2f165667c5ULL));
  }

  static constexpr result_type min() { return 0; }
  static constexpr result_type max() { return ~0ULL; }

  result_type operator()() { return detail::mix64(key_ + 0x9e3779b97f4a7c15ULL * ++ctr_); }

  // uniform double in [0, 1), 53 significant bits
  double uniform() { return static_cast<double>(operator()() >> 11) * 0x1.0p-53; }

  // deterministic bounded draw (Lemire multiply-shift; bias negligible for n << 2^64)
  uint64_t below(uint64_t n) {
    return static_cast<uint64_t>((static_cast<unsigned __int128>(operator()()) * n) >> 64);
  }

  // deterministic Fisher-Yates
  template <typename RandomIt>
  void shuffle(RandomIt first, RandomIt last) {
    const auto n = static_cast<uint64_t>(last - first);
    for (uint64_t i = n; i > 1; --i) {
      const uint64_t j = below(i);
      std::swap(first[i - 1], first[j]);
    }
  }

 private:
  uint64_t key_ = 0;
  uint64_t ctr_ = 0;
};

}  // namespace exsqs
