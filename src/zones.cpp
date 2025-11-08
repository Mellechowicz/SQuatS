#include "exsqs/zones.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace exsqs {

ZoneTable build_zones(const Structure& geom, int n_shells, double shell_rel_tol) {
  const int N = geom.natoms();
  if (N <= 0) throw std::runtime_error("build_zones: empty structure");
  if (n_shells <= 0) throw std::runtime_error("build_zones: n_shells must be positive");

  struct Pair {
    int i, j;
    double d;
  };
  std::vector<Pair> pairs;
  pairs.reserve(static_cast<size_t>(N) * N * 27);

  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j)
      for (int a = -1; a <= 1; ++a)
        for (int b = -1; b <= 1; ++b)
          for (int c = -1; c <= 1; ++c) {
            const Vec3 df = {geom.frac[j][0] + a - geom.frac[i][0],
                             geom.frac[j][1] + b - geom.frac[i][1],
                             geom.frac[j][2] + c - geom.frac[i][2]};
            const double d = norm(frac_to_cart(df, geom.cell));
            if (d < 1e-9) continue;  // exclude only the exact zero self pair [A1]
            pairs.push_back({i, j, d});
          }

  std::vector<double> ds;
  ds.reserve(pairs.size());
  for (const auto& p : pairs) ds.push_back(p.d);
  std::sort(ds.begin(), ds.end());

  ZoneTable zt;
  (void)n_shells;
  zt.radii.push_back(ds.front());  // first coordination shell only
  zt.n_shells = 1;

  const double rmax = zt.radii.back() * (1.0 + shell_rel_tol);

  zt.nbrs.assign(N, {});
  std::vector<std::vector<int>> per_site(N, std::vector<int>(zt.n_shells, 0));
  for (const auto& p : pairs) {
    if (p.d > rmax) continue;
    int sh = -1;
    for (int n = 0; n < zt.n_shells; ++n)
      if (p.d <= zt.radii[n] * (1.0 + shell_rel_tol)) {
        sh = n;
        break;
      }
    if (sh < 0) continue;
    zt.nbrs[p.i].push_back({p.j, sh});
    per_site[p.i][sh]++;
  }

  zt.coord_num = per_site[0];
  for (int i = 1; i < N; ++i)
    if (per_site[i] != zt.coord_num)
      throw std::runtime_error(
          "build_zones: sites do not form a single orbit (unequal coordination); "
          "v1 requires one Wyckoff orbit");

  zt.half_min_width = 0.0;  // minimum-image guard comes later
  zt.shells_exceed_half_width = false;
  return zt;
}

}  // namespace exsqs
