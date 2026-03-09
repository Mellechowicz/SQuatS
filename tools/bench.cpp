// T-B1 [SPEC section 12]: evaluations/s vs N, recorded, no gate. Serial, one core.
#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "exsqs/correlation.hpp"
#include "exsqs/dedup.hpp"
#include "exsqs/displacements.hpp"
#include "exsqs/lattice.hpp"
#include "exsqs/rng.hpp"
#include "exsqs/structure.hpp"
#include "exsqs/symmetry.hpp"
#include "exsqs/zones.hpp"

using namespace exsqs;
using Clock = std::chrono::steady_clock;

static double ms(Clock::time_point a, Clock::time_point b) {
  return std::chrono::duration<double, std::milli>(b - a).count();
}

int main() {
  std::printf("T-B1 benchmark: bcc n x n x n, W70Cr30-like decorations, symprec 1e-5\n");
  std::printf("%6s %8s %10s %11s %9s %9s %10s\n", "N", "|Pi|", "canon_ms", "spglib_ms", "D_ms",
              "corr_ms", "evals/s");
  for (int n : {3, 4, 5, 6, 8}) {
    const Prototype proto = make_bcc(3.165);
    const std::vector<std::string> names = {"W", "Cr"};
    const Mat3i H = {{{n, 0, 0}, {0, n, 0}, {0, 0, n}}};
    const Structure geom = make_supercell(proto, H, names, 0);
    const int N = geom.natoms();
    const ZoneTable zt = build_zones(geom, 7, 1e-3);
    const auto w = make_weights(WeightForm::InvN, zt);
    const auto perms = site_permutations(geom, 1e-5);
    const int nW = static_cast<int>(std::lround(0.7 * N));
    std::vector<int> base(static_cast<size_t>(N), 0);
    for (int i = nW; i < N; ++i) base[static_cast<size_t>(i)] = 1;
    const int reps = (N <= 128) ? 20 : (N <= 250 ? 10 : (N <= 432 ? 6 : 3));
    double t_can = 0, t_spg = 0, t_D = 0, t_corr = 0, t_all = 0;
    for (int r = 0; r < reps; ++r) {
      auto sig = base;
      CounterRng rng(1234, static_cast<uint64_t>(n), 0, static_cast<uint64_t>(r),
                     RngPurpose::Generic);
      rng.shuffle(sig.begin(), sig.end());
      const auto a0 = Clock::now();
      const auto can = canonical_labels(sig, perms);
      const auto a1 = Clock::now();
      const Structure s = decorate(geom, sig, names);
      const SymmetryInfo info = get_symmetry(s, 1e-5);
      const auto a2 = Clock::now();
      const int D = displacement_count(s, info, 1e-5);
      const auto a3 = Clock::now();
      const CorrData cd = count_pairs(s, zt);
      const double ep = e_pure_diagonal(cd, {0.7, 0.3}, w);
      const auto a4 = Clock::now();
      (void)can;
      (void)D;
      (void)ep;
      t_can += ms(a0, a1);
      t_spg += ms(a1, a2);
      t_D += ms(a2, a3);
      t_corr += ms(a3, a4);
      t_all += ms(a0, a4);
    }
    std::printf("%6d %8zu %10.2f %11.2f %9.2f %9.3f %10.1f\n", N, perms.size(), t_can / reps,
                t_spg / reps, t_D / reps, t_corr / reps, 1000.0 * reps / t_all);
  }
  return 0;
}
