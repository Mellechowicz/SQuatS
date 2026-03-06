#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>

#include "exsqs/displacements.hpp"
#include "exsqs/structure.hpp"
#include "exsqs/symmetry.hpp"
#include "exsqs/testcells.hpp"

int main(int argc, char** argv) {
  using namespace exsqs;
  const std::string outdir = (argc > 1) ? argv[1] : "tests/data";
  std::filesystem::create_directories(outdir + "/cells");

  std::ofstream man(outdir + "/manifest.txt");
  std::ofstream dc(outdir + "/dcounts_cpp.txt");
  std::printf("%-16s %6s %5s %-10s %5s\n", "name", "natoms", "sg", "symbol", "D");
  for (const auto& c : gate_cells()) {
    write_poscar(c.s, outdir + "/cells/" + c.name + ".vasp", c.name);
    const auto info = get_symmetry(c.s);
    const int D = displacement_count(c.s);
    man << c.name << " " << c.s.natoms() << " " << info.sg_number << " " << info.sg_symbol << " "
        << D << "\n";
    dc << c.name << " " << D << "\n";
    std::printf("%-16s %6d %5d %-10s %5d\n", c.name.c_str(), c.s.natoms(), info.sg_number,
                info.sg_symbol.c_str(), D);
  }
  std::printf("wrote %s/cells/*.vasp, manifest.txt, dcounts_cpp.txt\n", outdir.c_str());
  return 0;
}
