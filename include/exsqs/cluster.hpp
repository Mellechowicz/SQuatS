#pragma once
#include <cstdint>
#include <vector>

#include "exsqs/structure.hpp"
#include "exsqs/zones.hpp"

namespace exsqs {

// v1.9 multiplet sectors (SPEC section 4.2). Cluster classes of the
// UNDECORATED geometry, k = 3 (triplets) and k = 4 (quadruplets), built from
// the same 27-image background as the zones [A1]. Classes partition the
// instances by the sorted edge-length multiset -- an isometry invariant, so
// every class is a union of orbits of the dedup group Pi [A9] and every
// class-resolved count is invariant under the relabelling group: the
// canonical-label archive, the [D6] mutation and the bit-reproducibility
// contract carry over with no further argument.
//
// Instances are ANCHORED: every geometric k-cluster whose edges all fit the
// cutoff appears exactly k times (once per vertex able to anchor it), so
// I_c = k x the geometric count and every occupancy Pi_c(tau) = C_c(tau)/I_c
// is unaffected by the uniform multiplicity.
struct ClusterClass {
  int k = 3;
  std::vector<double> edges;   // sorted representative edge lengths (Angstrom)
  long long instances = 0;     // I_c (anchored)
  std::vector<int32_t> sites;  // flat [instances * k] site indices
  double floor_bound = 0.0;    // sum_tau dist(p(tau) I_c, Z)/I_c  (bound)
};

struct ClusterTable {
  int K = 0;
  int shell3 = 0, shell4 = 0;  // zones included in the cutoffs (1-based)
  double r3 = 0.0, r4 = 0.0;   // resolved cutoff radii
  std::vector<ClusterClass> c3, c4;
  // dense sorted-tuple -> multiset-slot lookups (valid for t1 <= t2 <= ...)
  std::vector<int> rank3;  // [K*K*K]
  std::vector<int> rank4;  // [K*K*K*K]
  int nmul3 = 0, nmul4 = 0;              // C(K+2,3), C(K+3,4)
  std::vector<double> target3, target4;  // p(tau) per multiset slot
  double floor3 = 0.0, floor4 = 0.0;     // per-sector floor bounds (A_c = 1)
  long long inst3 = 0, inst4 = 0;        // total anchored instances
};

// Build the table on the undecorated supercell. shell3/shell4 select how many
// coordination zones the triplet/quadruplet cutoff radius spans (>= 1); x is
// the achieved composition [A5] (the multinomial targets use it, matching the
// pair sectors). quads = false skips the k = 4 enumeration entirely.
ClusterTable build_clusters(const Structure& geom, const ZoneTable& zt,
                            const std::vector<double>& x, int shell3, int shell4,
                            bool quads);

// L1 sector errors of a decoration sigma (site -> species index):
//   e_k = sum_c sum_tau | C_c(tau)/I_c - p(tau) |.
// Serial by design: the engine parallelizes ACROSS candidates, one evaluation
// stays on one thread (integer counting, deterministic in any order).
void e_multiplets(const std::vector<int>& sigma, const ClusterTable& ct, double& e3, double& e4);

}  // namespace exsqs
