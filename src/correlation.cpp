#include "exsqs/correlation.hpp"

#include <cmath>
#include <stdexcept>

namespace exsqs {

std::vector<double> make_weights(WeightForm form, const ZoneTable& zt, double p,
                                 const std::vector<double>& custom) {
  std::vector<double> w(zt.n_shells, 0.0);
  switch (form) {
    case WeightForm::InvN:
      for (int n = 0; n < zt.n_shells; ++n) w[n] = 1.0 / (n + 1);
      break;
    case WeightForm::InvNPow:
      for (int n = 0; n < zt.n_shells; ++n) w[n] = 1.0 / std::pow(n + 1.0, p);
      break;
    case WeightForm::InvR:
      for (int n = 0; n < zt.n_shells; ++n) w[n] = 1.0 / zt.radii[n];
      break;
    case WeightForm::Custom:
      if (static_cast<int>(custom.size()) != zt.n_shells)
        throw std::runtime_error("make_weights: custom weight size mismatch");
      w = custom;
      break;
  }
  return w;
}

double CorrData::pi(int n, int t, int t2) const {
  const long long p = getP(n, t);
  if (p == 0) throw std::runtime_error("CorrData::pi: empty population");
  return static_cast<double>(getC(n, t, t2)) / static_cast<double>(p);
}

CorrData count_pairs(const Structure& s, const ZoneTable& zt) {
  CorrData d;
  d.n_shells = zt.n_shells;
  d.K = s.nspecies();
  d.C.assign(static_cast<size_t>(d.n_shells) * d.K * d.K, 0);
  d.P.assign(static_cast<size_t>(d.n_shells) * d.K, 0);
  for (int i = 0; i < s.natoms(); ++i) {
    const int ti = s.species[i];
    for (const auto& nb : zt.nbrs[i]) {
      const int tj = s.species[nb.j];
      d.C[(static_cast<size_t>(nb.shell) * d.K + ti) * d.K + tj]++;
      d.P[static_cast<size_t>(nb.shell) * d.K + ti]++;
    }
  }
  return d;
}

double e_pure_diagonal(const CorrData& cd, const std::vector<double>& x,
                       const std::vector<double>& w) {
  double E = 0.0;
  for (int n = 0; n < cd.n_shells; ++n)
    for (int t = 0; t < cd.K; ++t) {
      const long long p = cd.getP(n, t);
      if (p == 0) continue;  // absent species contributes nothing
      const double pi = static_cast<double>(cd.getC(n, t, t)) / p;
      E += w[n] * std::abs(pi - x[t]);
    }
  return E;
}

double e_pure_full(const CorrData& cd, const std::vector<double>& x,
                   const std::vector<double>& w) {
  double E = 0.0;
  for (int n = 0; n < cd.n_shells; ++n)
    for (int t = 0; t < cd.K; ++t) {
      const long long p = cd.getP(n, t);
      if (p == 0) continue;
      for (int t2 = 0; t2 < cd.K; ++t2) {
        const double pi = static_cast<double>(cd.getC(n, t, t2)) / p;
        E += w[n] * std::abs(pi - x[t2]);
      }
    }
  return E;
}


double e_floor_diagonal(const ZoneTable& zt, const std::vector<int>& counts,
                        const std::vector<double>& x, const std::vector<double>& w) {
  double f = 0.0;
  for (int n = 0; n < zt.n_shells; ++n)
    for (size_t t = 0; t < counts.size(); ++t) {
      if (counts[t] <= 0) continue;
      const double P = static_cast<double>(counts[t]) * zt.coord_num[static_cast<size_t>(n)];
      const double v = x[t] * P;
      f += w[static_cast<size_t>(n)] * std::abs(v - std::round(v)) / P;
    }
  return f;
}

double e_floor_full(const ZoneTable& zt, const std::vector<int>& counts,
                    const std::vector<double>& x, const std::vector<double>& w) {
  double f = 0.0;
  for (int n = 0; n < zt.n_shells; ++n)
    for (size_t t = 0; t < counts.size(); ++t) {
      if (counts[t] <= 0) continue;
      const double P = static_cast<double>(counts[t]) * zt.coord_num[static_cast<size_t>(n)];
      for (size_t u = 0; u < counts.size(); ++u) {
        const double v = x[u] * P;
        f += w[static_cast<size_t>(n)] * std::abs(v - std::round(v)) / P;
      }
    }
  return f;
}

}  // namespace exsqs
