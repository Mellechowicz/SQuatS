#include <cstdio>
#include <string>
#include <vector>

#include "exsqs/lattice.hpp"
#include "exsqs/structure.hpp"
#include <exception>
#include "exsqs/zones.hpp"
#include "exsqs/symmetry.hpp"
#include "exsqs/correlation.hpp"
#include "exsqs/config.hpp"

using namespace exsqs;

namespace {
Structure demo_cell() {
  const Prototype p = make_sc(3.0);
  return make_supercell(p, {{{3, 0, 0}, {0, 3, 0}, {0, 0, 3}}}, {"W", "Cr"}, 0);
}
}  // namespace

int main(int argc, char** argv) {
  const std::string a = argc > 1 ? argv[1] : "";
  if (!a.empty() && a[0] != '-') {
    std::vector<std::string> ovr;
    for (int i = 2; i + 1 < argc; ++i)
      if (std::string(argv[i]) == "--set") ovr.push_back(argv[++i]);
    try {
      const RunConfig cfg = load_config(a, ovr);
      std::printf("config: %zu species |", cfg.species.size());
      for (size_t t = 0; t < cfg.species.size(); ++t)
        std::printf(" %s=%d", cfg.species[t].c_str(), cfg.counts[t]);
      std::printf("\n");
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "exsqs error: %s\n", e.what());
      return 1;
    }
  }
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
  if (a == "--demo-corr") {
    const Structure g = demo_cell();
    const ZoneTable zt = build_zones(g, 3, 1e-3);
    std::vector<int> sigma(static_cast<size_t>(g.natoms()));
    for (size_t i = 0; i < sigma.size(); ++i) sigma[i] = static_cast<int>(i % 2);
    const Structure dec = decorate(g, sigma, {"W", "Cr"});
    const CorrData cd = count_pairs(dec, zt);
    const std::vector<double> w = make_weights(WeightForm::InvN, zt);
    const std::vector<double> x = {0.5, 0.5};
    std::printf("demo-corr: alternating decoration E_pure=%.6e\n", e_pure_diagonal(cd, x, w));
    return 0;
  }
  std::printf("exsqs 0.3.0 (development)\n");
  return 0;
}
