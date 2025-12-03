#include <cstdio>
#include <string>
#include <vector>

#include "exsqs/lattice.hpp"
#include "exsqs/structure.hpp"
#include "exsqs/zones.hpp"
#include "exsqs/symmetry.hpp"
#include "exsqs/dedup.hpp"

using namespace exsqs;

namespace {
Structure demo_cell() {
  const Prototype p = make_sc(3.0);
  return make_supercell(p, {{{3, 0, 0}, {0, 3, 0}, {0, 0, 3}}}, {"W", "Cr"}, 0);
}
}  // namespace

int main(int argc, char** argv) {
  const std::string a = argc > 1 ? argv[1] : "";
  if (a == "--demo-geom") {
    const Structure s = demo_cell();
    std::printf("demo-geom: sc 3x3x3 supercell -> %d sites\n", s.natoms());
    return 0;
  }
  if (a == "--demo-zones") {
    const Structure s = demo_cell();
    const ZoneTable zt = build_zones(s, 5, 1e-3);
    std::printf("demo-zones: %d shells:", zt.n_shells);
    for (double r : zt.radii) std::printf(" %.4f", r);
    std::printf("  (z1=%d)\n", zt.coord_num.empty() ? 0 : zt.coord_num[0]);
    return 0;
  }
  if (a == "--demo-sym") {
    const Structure s = demo_cell();
    const SymmetryInfo info = get_symmetry(s, 1e-5);
    std::printf("demo-sym: empty cell SG=%d (%s), %zu ops\n", info.sg_number,
                info.sg_symbol.c_str(), info.ops.size());
    return 0;
  }
  if (a == "--demo-group") {
    const Structure g = demo_cell();
    std::printf("demo-group: |Pi| =%4zu supercell site permutations\n", site_permutations(g, 1e-5).size());  // species-blind
    return 0;
  }
  std::printf("exsqs 0.2.0 (development)\n");
  return 0;
}
