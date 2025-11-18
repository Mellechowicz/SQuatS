#pragma once
#include <vector>

#include "exsqs/structure.hpp"
#include "exsqs/zones.hpp"

namespace exsqs {

enum class WeightForm { InvN, InvNPow, InvR, Custom };

std::vector<double> make_weights(WeightForm form, const ZoneTable& zt, double p = 1.0,
                                 const std::vector<double>& custom = {});

// Pair-correlation counters over ordered pairs [SPEC section 4].
//   C_{t,t'}(n) : ordered pairs with centre type t, neighbour type t'
//   P_t(n)      : pairs whose central atom is type t  [D1]
struct CorrData {
  int n_shells = 0;
  int K = 0;
  std::vector<long long> C;  // [(n*K + t)*K + t2]
  std::vector<long long> P;  // [n*K + t]
  long long getC(int n, int t, int t2) const {
    return C[(static_cast<size_t>(n) * K + t) * K + t2];
  }
  long long getP(int n, int t) const { return P[static_cast<size_t>(n) * K + t]; }
  double pi(int n, int t, int t2) const;  // conditional probability Pi_{t,t'}(n)
};

CorrData count_pairs(const Structure& s, const ZoneTable& zt);

// E_pure in the two modes [A16]; weights w indexed by shell (0-based).
// K = 2 identity: e_pure_full == 2 * e_pure_diagonal (test T-C4).
double e_pure_diagonal(const CorrData& cd, const std::vector<double>& x,
                       const std::vector<double>& w);
double e_pure_full(const CorrData& cd, const std::vector<double>& x, const std::vector<double>& w);


// v1.1 (SPEC 4.1): exact L1 lower bound from count quantization. C counters are
// integers while P_t(n) = N_t*Z_n is fixed, so each term obeys
// |Pi - x| >= dist(x*P, Z)/P. Independent per-term minima => computable bound.
// For K = 2, e_floor_full == 2 * e_floor_diagonal exactly.
double e_floor_diagonal(const ZoneTable& zt, const std::vector<int>& counts,
                        const std::vector<double>& x, const std::vector<double>& w);
double e_floor_full(const ZoneTable& zt, const std::vector<int>& counts,
                    const std::vector<double>& x, const std::vector<double>& w);

}  // namespace exsqs
