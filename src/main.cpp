#include <cstdio>
#include <exception>
#include <string>
#include <vector>

#include "exsqs/config.hpp"
#include "exsqs/evolution.hpp"

int main(int argc, char** argv) {
  std::string path;
  std::vector<std::string> ovr;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "--set" && i + 1 < argc) {
      ovr.push_back(argv[++i]);
    } else if (a == "--out" && i + 1 < argc) {
      ovr.push_back(std::string("output.dir=") + argv[++i]);
    } else if (path.empty()) {
      path = a;
    }
  }
  if (path.empty()) {
    std::fprintf(stderr, "usage: exsqs <config.yaml> [--set key=value ...] [--out DIR]\nruns the extinction loop and prints the archived structures\n");
    return 1;
  }
  try {
    const exsqs::RunConfig cfg = exsqs::load_config(path, ovr);
    const exsqs::RunContext ctx = exsqs::RunContext::build(cfg);
    std::printf("exsqs 0.5.0 | %d sites | run | E_floor=%.6e\n", ctx.geom.natoms(), ctx.e_floor);
    const exsqs::RunOutput out = exsqs::run_evolution(cfg, ctx);
    for (const exsqs::Individual& I : out.outputs)
      std::printf("  E_pure=%.6e D=%4d E_obj=%.6e SG=%d (%s)\n", I.e_pure, I.D, I.e_obj, I.sg,
                  I.sg_symbol.c_str());
    return 0;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "exsqs error: %s\n", e.what());
    return 1;
  }
}
