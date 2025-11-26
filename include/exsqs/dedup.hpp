#pragma once
#include <cstdint>
#include <vector>

#include "exsqs/structure.hpp"

namespace exsqs {

// Site permutation induced by one symmetry op of the geometry: perm[i] = image of i.
// Throws if the op does not map the site set onto itself.
std::vector<int> permutation_of_op(const Structure& geom, const Mat3i& R, const Vec3& t);

// Full site-permutation group of the UNDECORATED supercell [A9]: species are
// ignored, every spglib op of the empty geometry is converted to a permutation,
// duplicates removed. Includes the identity and all internal translations.
std::vector<std::vector<int>> site_permutations(const Structure& geom, double symprec = 1e-5);

// Canonical decoration label [A10]: for each perm build sigma'(perm[i]) = sigma(i)
// and keep the lexicographic minimum. Two decorations of the same geometry are
// symmetry-duplicates iff their canonical labels are equal.
std::vector<uint8_t> canonical_labels(const std::vector<int>& sigma,
                                      const std::vector<std::vector<int>>& perms);

// FNV-1a 64-bit over the canonical label (fast prefilter; exact compare is the
// canonical vector itself).
uint64_t hash_labels(const std::vector<uint8_t>& labels);

}  // namespace exsqs
