// v1.9 multiplet-sector gates (SPEC 4.2).
//   T-M1: anchored instance counts against the lattice inventories (bcc has a
//         single isosceles triangle class and NO compact quad; fcc has the
//         equilateral NN triangle and the NN tetrahedron), plus a golden
//         recount of the decorated multiset histogram by an independent
//         brute-force enumerator.
//   T-M2: the per-sector floor bound against exhaustive decoration
//         enumeration on a tiny cell (bound holds; sectors need not attain it).
//   T-M3: Pi-invariance of the sector errors (the invariance lemma that keeps
//         the archive and [D6] valid), and the lambda = 0 evaluation identity.
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <random>
#include <string>
#include <vector>

#include "exsqs/cluster.hpp"
#include "exsqs/config.hpp"
#include "exsqs/dedup.hpp"
#include "exsqs/evolution.hpp"
#include "exsqs/lattice.hpp"
#include "exsqs/linalg.hpp"
#include "exsqs/structure.hpp"
#include "exsqs/zones.hpp"

using namespace exsqs;
using Catch::Approx;

namespace {

const std::vector<std::string> kWCr = {"W", "Cr"};
const std::vector<std::string> kWCrTa = {"W", "Cr", "Ta"};

std::vector<int> random_counts_sigma(int n, const std::vector<int>& counts, uint64_t seed) {
  std::vector<int> s;
  for (size_t t = 0; t < counts.size(); ++t)
    for (int k = 0; k < counts[t]; ++k) s.push_back(static_cast<int>(t));
  REQUIRE(static_cast<int>(s.size()) == n);
  std::mt19937_64 rng(seed);
  std::shuffle(s.begin(), s.end(), rng);
  return s;
}

// Independent brute-force recount: anchored triples over the 27-image
// background, histogram over the species multiset (classes aggregated).
// Deliberately reimplemented from scratch -- shares no helper with cluster.cpp.
std::map<std::array<int, 3>, long long> brute_triple_histogram(const Structure& g,
                                                               const std::vector<int>& sigma,
                                                               double rcut) {
  std::vector<Vec3> cart(static_cast<size_t>(g.natoms()));
  for (int i = 0; i < g.natoms(); ++i) {
    Vec3 r{0, 0, 0};
    for (int c = 0; c < 3; ++c) r = r + g.frac[static_cast<size_t>(i)][c] * g.cell[static_cast<size_t>(c)];
    cart[static_cast<size_t>(i)] = r;
  }
  std::map<std::array<int, 3>, long long> h;
  for (int i = 0; i < g.natoms(); ++i) {
    std::vector<std::pair<int, Vec3>> nb;
    for (int j = 0; j < g.natoms(); ++j)
      for (int a = -1; a <= 1; ++a)
        for (int b = -1; b <= 1; ++b)
          for (int c = -1; c <= 1; ++c) {
            const Vec3 p = cart[static_cast<size_t>(j)] + static_cast<double>(a) * g.cell[0] +
                           static_cast<double>(b) * g.cell[1] + static_cast<double>(c) * g.cell[2];
            const double d = norm(p - cart[static_cast<size_t>(i)]);
            if (d > 1e-9 && d <= rcut) nb.emplace_back(j, p);
          }
    for (size_t u = 0; u < nb.size(); ++u)
      for (size_t v = u + 1; v < nb.size(); ++v) {
        const double duv = norm(nb[u].second - nb[v].second);
        if (duv <= 1e-9 || duv > rcut) continue;
        std::array<int, 3> tau = {sigma[static_cast<size_t>(i)],
                                  sigma[static_cast<size_t>(nb[u].first)],
                                  sigma[static_cast<size_t>(nb[v].first)]};
        std::sort(tau.begin(), tau.end());
        h[tau]++;
      }
  }
  return h;
}

}  // namespace

TEST_CASE("T-M1 bcc inventory: one isosceles class; disphenoid quad", "[multiplets][gate]") {
  auto g = make_supercell_diag(make_bcc(3.165), 3, 3, 3, kWCr);
  const int N = g.natoms();
  REQUIRE(N == 54);
  auto zt = build_zones(g, 3);
  auto ct = build_clusters(g, zt, {0.5, 0.5}, /*shell3=*/2, /*shell4=*/2, /*quads=*/true);

  // the (r1, r1, r2) isosceles is the ONLY triangle within two zones; 12
  // geometric per site -> 36 N anchored. The two-zone 4-clique is the
  // r1^4 r2^2 disphenoid (6 per site -> 24 N anchored); bcc admits no
  // EQUILATERAL compact cluster, and the NN-only cutoff admits nothing.
  REQUIRE(ct.c3.size() == 1);
  REQUIRE(ct.inst3 == 36LL * N);
  REQUIRE(ct.c3[0].edges[0] == Approx(zt.radii[0]).epsilon(1e-9));
  REQUIRE(ct.c3[0].edges[1] == Approx(zt.radii[0]).epsilon(1e-9));
  REQUIRE(ct.c3[0].edges[2] == Approx(zt.radii[1]).epsilon(1e-9));
  REQUIRE(ct.c4.size() == 1);
  REQUIRE(ct.inst4 == 24LL * N);
  REQUIRE(ct.c4[0].edges[3] == Approx(zt.radii[0]).epsilon(1e-9));  // 4 x r1 ...
  REQUIRE(ct.c4[0].edges[4] == Approx(zt.radii[1]).epsilon(1e-9));  // ... + 2 x r2

  auto nn = build_clusters(g, zt, {0.5, 0.5}, /*shell3=*/1, /*shell4=*/1, /*quads=*/true);
  REQUIRE(nn.c3.empty());  // two bcc NN are never NN of each other
  REQUIRE(nn.c4.empty());
}

TEST_CASE("T-M1 fcc inventory: equilateral triangle + NN tetrahedron", "[multiplets][gate]") {
  auto g = make_supercell_diag(make_fcc(3.6), 2, 2, 2, kWCr);
  const int N = g.natoms();
  REQUIRE(N == 32);
  auto zt = build_zones(g, 2);
  auto ct = build_clusters(g, zt, {0.5, 0.5}, /*shell3=*/1, /*shell4=*/1, /*quads=*/true);

  // 8 equilateral NN triangles and 2 NN tetrahedra per site (textbook fcc):
  // anchored 3 * 8 N and 4 * 2 N.
  REQUIRE(ct.c3.size() == 1);
  REQUIRE(ct.inst3 == 24LL * N);
  REQUIRE(ct.c3[0].edges[2] == Approx(zt.radii[0]).epsilon(1e-9));
  REQUIRE(ct.c4.size() == 1);
  REQUIRE(ct.inst4 == 8LL * N);
}

TEST_CASE("T-M1 golden recount on sc, binary and ternary", "[multiplets][gate]") {
  for (int K = 2; K <= 3; ++K) {
    auto g = make_supercell_diag(make_sc(3.0), 3, 3, 3,
                                 K == 2 ? kWCr : kWCrTa);
    const int N = g.natoms();
    auto zt = build_zones(g, 3);
    const std::vector<double> x = K == 2 ? std::vector<double>{14.0 / 27, 13.0 / 27}
                                         : std::vector<double>{10.0 / 27, 9.0 / 27, 8.0 / 27};
    auto ct = build_clusters(g, zt, x, /*shell3=*/2, /*shell4=*/1, /*quads=*/false);
    const std::vector<int> counts =
        K == 2 ? std::vector<int>{14, 13} : std::vector<int>{10, 9, 8};

    for (uint64_t seed = 1; seed <= 5; ++seed) {
      const auto sigma = random_counts_sigma(N, counts, seed);
      // kernel counts, aggregated over classes
      std::map<std::array<int, 3>, long long> agg;
      for (const auto& cc : ct.c3)
        for (long long n = 0; n < cc.instances; ++n) {
          std::array<int, 3> tau = {sigma[static_cast<size_t>(cc.sites[static_cast<size_t>(3 * n)])],
                                    sigma[static_cast<size_t>(cc.sites[static_cast<size_t>(3 * n + 1)])],
                                    sigma[static_cast<size_t>(cc.sites[static_cast<size_t>(3 * n + 2)])]};
          std::sort(tau.begin(), tau.end());
          agg[tau]++;
        }
      const auto brute = brute_triple_histogram(g, sigma, ct.r3);
      REQUIRE(agg == brute);
    }
  }
}

TEST_CASE("T-M2 floor bound holds under exhaustive enumeration", "[multiplets][gate]") {
  auto g = make_supercell_diag(make_sc(3.0), 2, 2, 2, kWCr);
  const int N = g.natoms();
  REQUIRE(N == 8);
  auto zt = build_zones(g, 2);

  // equiatomic x = 1/2: every p(tau) I_c is integer, the bound is 0 and the
  // 70 decorations must contain one attaining it exactly
  auto ceq = build_clusters(g, zt, {0.5, 0.5}, /*shell3=*/2, /*shell4=*/1, /*quads=*/false);
  REQUIRE(ceq.floor3 == 0.0);

  // x = (3/8, 5/8): the targets quantize, the bound is positive
  auto ct = build_clusters(g, zt, {3.0 / 8, 5.0 / 8}, /*shell3=*/2, /*shell4=*/1, false);
  REQUIRE(ct.floor3 > 0.0);

  double emin = 1e300, emin_eq = 1e300;
  for (int mask = 0; mask < 256; ++mask) {
    const int pc = __builtin_popcount(static_cast<unsigned>(mask));
    std::vector<int> sigma(static_cast<size_t>(N), 0);
    for (int i = 0; i < N; ++i) sigma[static_cast<size_t>(i)] = (mask >> i) & 1;
    double e3 = 0.0, e4 = 0.0;
    if (pc == 5) {  // 3 W + 5 Cr
      e_multiplets(sigma, ct, e3, e4);
      emin = std::min(emin, e3);
    } else if (pc == 4) {  // 4 W + 4 Cr
      e_multiplets(sigma, ceq, e3, e4);
      emin_eq = std::min(emin_eq, e3);
    }
  }
  REQUIRE(emin >= ct.floor3 - 1e-12);  // bound: never above the true minimum
  REQUIRE(emin > 0.0);
  REQUIRE(emin_eq == Approx(0.0).margin(1e-12));  // 0-bound attained at x = 1/2
}

TEST_CASE("T-M3 sector errors are Pi-invariant; lambda 0 is inert", "[multiplets][gate]") {
  auto g = make_supercell_diag(make_bcc(3.165), 2, 2, 2, kWCrTa);
  const int N = g.natoms();
  auto zt = build_zones(g, 3);
  auto ct = build_clusters(g, zt, {6.0 / 16, 5.0 / 16, 5.0 / 16}, 2, 1, true);
  const auto perms = site_permutations(g, 1e-5);
  REQUIRE(perms.size() >= 2);

  const auto sigma = random_counts_sigma(N, {6, 5, 5}, 20260717);
  double e3 = 0.0, e4 = 0.0;
  e_multiplets(sigma, ct, e3, e4);
  REQUIRE(e3 > 0.0);

  // relabelling by any op of the dedup group leaves every sector unchanged
  int checked = 0;
  for (size_t p = 0; p < perms.size(); p += std::max<size_t>(1, perms.size() / 7)) {
    std::vector<int> moved(static_cast<size_t>(N));
    for (int i = 0; i < N; ++i) moved[static_cast<size_t>(perms[p][static_cast<size_t>(i)])] =
        sigma[static_cast<size_t>(i)];
    double m3 = 0.0, m4 = 0.0;
    e_multiplets(moved, ct, m3, m4);
    REQUIRE(m3 == Approx(e3).margin(1e-12));
    REQUIRE(m4 == Approx(e4).margin(1e-12));
    ++checked;
  }
  REQUIRE(checked >= 2);
}

TEST_CASE("T-M4 trajectory with active sectors is thread-invariant", "[multiplets][gate]") {
  RunConfig c;
  c.proto = make_sc(3.0);
  c.H = {{{3, 0, 0}, {0, 3, 0}, {0, 0, 3}}};
  c.species = {"W", "Cr"};
  c.counts = {14, 13};
  c.x_target = {14.0 / 27.0, 13.0 / 27.0};
  c.x_achieved = c.x_target;
  c.n_shells = 5;
  c.lambda3 = 0.6;   // sectors ON: the point of the gate
  c.lambda4 = 0.2;
  c.mshell3 = 2;
  c.mshell4 = 2;
  c.population = 12;
  c.outputs = 4;
  c.e_tol = 1e-12;  // unreachable
  c.max_generations = 8;
  c.islands = 2;
  c.migration_every = 3;
  c.retry_budget = 60;
  c.seed = 20260717;
  c.log_info = false;
  const RunContext ctx = RunContext::build(c);
  REQUIRE(ctx.multiplets);
  REQUIRE(ctx.e_floor > ctx.e_floor_pair);

  c.omp_threads = 1;
  const RunOutput a = run_evolution(c, ctx);
  c.omp_threads = 4;
  const RunOutput b = run_evolution(c, ctx);
  REQUIRE(a.evals == b.evals);
  REQUIRE(a.outputs.size() == b.outputs.size());
  for (size_t i = 0; i < a.outputs.size(); ++i) {
    REQUIRE(a.outputs[i].canonical == b.outputs[i].canonical);  // bit-identical
    REQUIRE(a.outputs[i].e_pure == b.outputs[i].e_pure);        // exact, no eps
    REQUIRE(a.outputs[i].e_obj == b.outputs[i].e_obj);
    REQUIRE(a.outputs[i].D == b.outputs[i].D);
  }
  REQUIRE(a.log.size() == b.log.size());
  for (size_t i = 0; i < a.log.size(); ++i) {
    REQUIRE(a.log[i].best_e_obj == b.log[i].best_e_obj);
    REQUIRE(a.log[i].killed == b.log[i].killed);
    REQUIRE(a.log[i].dup_rejected == b.log[i].dup_rejected);
  }
}
