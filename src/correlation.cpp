#include "exsqs/correlation.hpp"

#include <cmath>
#include <stdexcept>

namespace exsqs {

std::vector<double> make_weights(WeightForm form, const ZoneTable& zt, double p,
                                 const std::vector<double>& custom) {
  std::vector<double> w(zt.n_shells, 0.0);
  (void)form;
  (void)p;
  (void)custom;
  for (int n = 0; n < zt.n_shells; ++n) w[n] = 1.0 / (n + 1);  // inv_n only for now
  return w;
}

double CorrData::pi(int n, int t, int t2) const {
  const long long p = getP(n, t);
  if (p == 0) throw std::runtime_error("CorrData::pi: empty population");
  return static_cast<double>(getC(n, t, t2)) / static_cast<double>(p);
}

// C indexed as (shell, t, t2); P as (shell, t)
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

}  // namespace exsqs
