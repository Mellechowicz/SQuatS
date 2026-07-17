#include "exsqs/cluster.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>

#include "exsqs/linalg.hpp"

namespace exsqs {
namespace {

// One background neighbour of an anchor site: real site index j, cartesian
// position (image-resolved) and distance to the anchor. The zone table cannot
// serve here -- its entries drop the image vector, and two images of the same
// site can close two different triangles.
struct BgNbr {
  int j;
  Vec3 pos;
  double d;
};

Vec3 cart_of(const Structure& g, int i) {
  Vec3 r{0, 0, 0};
  for (int c = 0; c < 3; ++c)
    r = r + g.frac[static_cast<size_t>(i)][c] * g.cell[static_cast<size_t>(c)];
  return r;
}

// Neighbours of site i within rcut over the 27-image background [A1]; the
// zero-distance self excluded, nonzero self-images included.
std::vector<BgNbr> bg_neighbours(const Structure& g, const std::vector<Vec3>& cart, int i,
                                 double rcut) {
  std::vector<BgNbr> out;
  for (int j = 0; j < g.natoms(); ++j)
    for (int a = -1; a <= 1; ++a)
      for (int b = -1; b <= 1; ++b)
        for (int c = -1; c <= 1; ++c) {
          const Vec3 p = cart[static_cast<size_t>(j)] + static_cast<double>(a) * g.cell[0] +
                         static_cast<double>(b) * g.cell[1] + static_cast<double>(c) * g.cell[2];
          const double d = norm(p - cart[static_cast<size_t>(i)]);
          if (d < 1e-9 || d > rcut) continue;
          out.push_back({j, p, d});
        }
  // deterministic order: by (j, image) via position lexicographic
  std::sort(out.begin(), out.end(), [](const BgNbr& u, const BgNbr& v) {
    if (u.j != v.j) return u.j < v.j;
    for (int c = 0; c < 3; ++c)
      if (u.pos[c] != v.pos[c]) return u.pos[c] < v.pos[c];
    return false;
  });
  return out;
}

// Quantize an edge length into a bin index; bins merge within rel_tol exactly
// like the zone radii do. Representatives are APPEND-ONLY: indices issued
// earlier must stay valid (they are already baked into class keys), so the
// list keeps discovery order -- deterministic, and the classes are re-sorted
// by their edge lengths at the end of the build anyway.
int edge_bin(std::vector<double>& bins, double d, double rel_tol) {
  for (size_t b = 0; b < bins.size(); ++b)
    if (std::abs(d - bins[b]) <= rel_tol * bins[b]) return static_cast<int>(b);
  bins.push_back(d);
  return static_cast<int>(bins.size()) - 1;
}

double factorial(int n) {
  double f = 1.0;
  for (int i = 2; i <= n; ++i) f *= i;
  return f;
}

// Multiset slots of k species out of K, lexicographic in the sorted tuple;
// fills the dense sorted-tuple lookup and the multinomial target p(tau).
void multiset_tables(int K, int k, std::vector<int>& rank, std::vector<double>& target,
                     const std::vector<double>& x, int& nmul) {
  int size = 1;
  for (int i = 0; i < k; ++i) size *= K;
  rank.assign(static_cast<size_t>(size), -1);
  target.clear();
  nmul = 0;
  std::vector<int> t(static_cast<size_t>(k), 0);
  while (true) {
    // dense index of the sorted tuple + its multinomial probability
    int idx = 0;
    for (int i = 0; i < k; ++i) idx = idx * K + t[static_cast<size_t>(i)];
    rank[static_cast<size_t>(idx)] = nmul;
    std::vector<int> cnt(static_cast<size_t>(K), 0);
    double p = factorial(k);
    for (int i = 0; i < k; ++i) {
      p *= x[static_cast<size_t>(t[static_cast<size_t>(i)])];
      cnt[static_cast<size_t>(t[static_cast<size_t>(i)])]++;
    }
    for (int s = 0; s < K; ++s) p /= factorial(cnt[static_cast<size_t>(s)]);
    target.push_back(p);
    ++nmul;
    // next non-decreasing tuple
    int pos = k - 1;
    while (pos >= 0 && t[static_cast<size_t>(pos)] == K - 1) --pos;
    if (pos < 0) break;
    const int v = t[static_cast<size_t>(pos)] + 1;
    for (int i = pos; i < k; ++i) t[static_cast<size_t>(i)] = v;
  }
}

double class_floor(const ClusterClass& cc, const std::vector<double>& target) {
  double f = 0.0;
  const double I = static_cast<double>(cc.instances);
  for (double p : target) {
    const double v = p * I;
    f += std::abs(v - std::round(v)) / I;
  }
  return f;
}

}  // namespace

ClusterTable build_clusters(const Structure& geom, const ZoneTable& zt,
                            const std::vector<double>& x, int shell3, int shell4, bool quads) {
  ClusterTable ct;
  ct.K = static_cast<int>(x.size());
  if (shell3 < 1 || shell3 > zt.n_shells || (quads && (shell4 < 1 || shell4 > zt.n_shells)))
    throw std::runtime_error("cluster: multiplet shell out of the zone range");
  ct.shell3 = shell3;
  ct.shell4 = quads ? shell4 : 0;
  const double rel_tol = 1e-3;  // matches build_zones' default merge tolerance
  ct.r3 = zt.radii[static_cast<size_t>(shell3 - 1)] * (1.0 + rel_tol);
  ct.r4 = quads ? zt.radii[static_cast<size_t>(shell4 - 1)] * (1.0 + rel_tol) : 0.0;

  multiset_tables(ct.K, 3, ct.rank3, ct.target3, x, ct.nmul3);
  if (quads) multiset_tables(ct.K, 4, ct.rank4, ct.target4, x, ct.nmul4);

  std::vector<Vec3> cart(static_cast<size_t>(geom.natoms()));
  for (int i = 0; i < geom.natoms(); ++i) cart[static_cast<size_t>(i)] = cart_of(geom, i);

  // class key: sorted edge-bin tuple (isometry invariant => union of Pi-orbits)
  std::vector<double> bins;
  std::map<std::vector<int>, size_t> key3, key4;

  for (int i = 0; i < geom.natoms(); ++i) {
    const auto nb = bg_neighbours(geom, cart, i, std::max(ct.r3, ct.r4));
    for (size_t u = 0; u < nb.size(); ++u)
      for (size_t v = u + 1; v < nb.size(); ++v) {
        const double duv = norm(nb[u].pos - nb[v].pos);
        // ---- triplets: all three edges inside r3 ----
        if (nb[u].d <= ct.r3 && nb[v].d <= ct.r3 && duv <= ct.r3 && duv > 1e-9) {
          std::vector<int> kk = {edge_bin(bins, nb[u].d, rel_tol),
                                 edge_bin(bins, nb[v].d, rel_tol),
                                 edge_bin(bins, duv, rel_tol)};
          std::sort(kk.begin(), kk.end());
          auto it = key3.find(kk);
          if (it == key3.end()) {
            it = key3.emplace(kk, ct.c3.size()).first;
            ct.c3.push_back({3, {nb[u].d, nb[v].d, duv}, 0, {}, 0.0});
            std::sort(ct.c3.back().edges.begin(), ct.c3.back().edges.end());
          }
          ClusterClass& cc = ct.c3[it->second];
          cc.sites.push_back(i);
          cc.sites.push_back(nb[u].j);
          cc.sites.push_back(nb[v].j);
          cc.instances++;
        }
        // ---- quadruplets: pick a third neighbour, all six edges inside r4 ----
        if (!quads || nb[u].d > ct.r4 || nb[v].d > ct.r4 || duv > ct.r4 || duv <= 1e-9)
          continue;
        for (size_t w = v + 1; w < nb.size(); ++w) {
          if (nb[w].d > ct.r4) continue;
          const double duw = norm(nb[u].pos - nb[w].pos);
          const double dvw = norm(nb[v].pos - nb[w].pos);
          if (duw > ct.r4 || duw <= 1e-9 || dvw > ct.r4 || dvw <= 1e-9) continue;
          std::vector<int> kk = {edge_bin(bins, nb[u].d, rel_tol), edge_bin(bins, nb[v].d, rel_tol),
                                 edge_bin(bins, nb[w].d, rel_tol), edge_bin(bins, duv, rel_tol),
                                 edge_bin(bins, duw, rel_tol),     edge_bin(bins, dvw, rel_tol)};
          std::sort(kk.begin(), kk.end());
          auto it = key4.find(kk);
          if (it == key4.end()) {
            it = key4.emplace(kk, ct.c4.size()).first;
            ct.c4.push_back({4, {nb[u].d, nb[v].d, nb[w].d, duv, duw, dvw}, 0, {}, 0.0});
            std::sort(ct.c4.back().edges.begin(), ct.c4.back().edges.end());
          }
          ClusterClass& cc = ct.c4[it->second];
          cc.sites.push_back(i);
          cc.sites.push_back(nb[u].j);
          cc.sites.push_back(nb[v].j);
          cc.sites.push_back(nb[w].j);
          cc.instances++;
        }
      }
  }

  // deterministic class order: by sorted edge lengths (map order is by bins,
  // which depend on discovery order for equal-tol merges -- re-sort explicitly)
  const auto by_edges = [](const ClusterClass& a, const ClusterClass& b) {
    return a.edges < b.edges;
  };
  std::sort(ct.c3.begin(), ct.c3.end(), by_edges);
  std::sort(ct.c4.begin(), ct.c4.end(), by_edges);

  for (auto& cc : ct.c3) {
    cc.floor_bound = class_floor(cc, ct.target3);
    ct.floor3 += cc.floor_bound;
    ct.inst3 += cc.instances;
  }
  for (auto& cc : ct.c4) {
    cc.floor_bound = class_floor(cc, ct.target4);
    ct.floor4 += cc.floor_bound;
    ct.inst4 += cc.instances;
  }
  return ct;
}

void e_multiplets(const std::vector<int>& sigma, const ClusterTable& ct, double& e3, double& e4) {
  e3 = e4 = 0.0;
  const int K = ct.K;
  std::vector<long long> cnt(static_cast<size_t>(std::max(ct.nmul3, ct.nmul4)));
  for (const ClusterClass& cc : ct.c3) {
    std::fill(cnt.begin(), cnt.begin() + ct.nmul3, 0);
    for (long long n = 0; n < cc.instances; ++n) {
      int a = sigma[static_cast<size_t>(cc.sites[static_cast<size_t>(3 * n)])];
      int b = sigma[static_cast<size_t>(cc.sites[static_cast<size_t>(3 * n + 1)])];
      int c = sigma[static_cast<size_t>(cc.sites[static_cast<size_t>(3 * n + 2)])];
      if (a > b) std::swap(a, b);  // 3-element sorting network
      if (b > c) std::swap(b, c);
      if (a > b) std::swap(a, b);
      cnt[static_cast<size_t>(ct.rank3[static_cast<size_t>((a * K + b) * K + c)])]++;
    }
    const double I = static_cast<double>(cc.instances);
    for (int m = 0; m < ct.nmul3; ++m)
      e3 += std::abs(static_cast<double>(cnt[static_cast<size_t>(m)]) / I -
                     ct.target3[static_cast<size_t>(m)]);
  }
  for (const ClusterClass& cc : ct.c4) {
    std::fill(cnt.begin(), cnt.begin() + ct.nmul4, 0);
    for (long long n = 0; n < cc.instances; ++n) {
      int q[4];
      for (int j = 0; j < 4; ++j)
        q[j] = sigma[static_cast<size_t>(cc.sites[static_cast<size_t>(4 * n + j)])];
      std::sort(q, q + 4);
      cnt[static_cast<size_t>(
          ct.rank4[static_cast<size_t>(((q[0] * K + q[1]) * K + q[2]) * K + q[3])])]++;
    }
    const double I = static_cast<double>(cc.instances);
    for (int m = 0; m < ct.nmul4; ++m)
      e4 += std::abs(static_cast<double>(cnt[static_cast<size_t>(m)]) / I -
                     ct.target4[static_cast<size_t>(m)]);
  }
}

}  // namespace exsqs
