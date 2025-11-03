#include <cstdio>
#include <string>
#include <vector>

#include "exsqs/lattice.hpp"
#include "exsqs/structure.hpp"

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
  std::printf("exsqs 0.0.1 (development)\n");
  return 0;
}
