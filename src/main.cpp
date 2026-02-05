#include <cstdio>
#include <exception>
#include <string>
#include <vector>

#include "exsqs/config.hpp"
#include "exsqs/evolution.hpp"

static void usage() {
  std::fprintf(stderr,
               "exsqs-mini: extinction-based evolutionary SQS generator\n"
               "usage: exsqs <config.yaml> [--set key=value ...] [--out DIR]\n"
               "exit codes: 0 converged to e_tol, 3 budget exhausted, 1 error\n");
}

int main(int argc, char** argv) {
  std::string path;
  std::vector<std::string> ovr;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "--set") {
      if (i + 1 >= argc) {
        usage();
        return 1;
      }
      ovr.push_back(argv[++i]);
    } else if (a == "--out") {
      if (i + 1 >= argc) {
        usage();
        return 1;
      }
      ovr.push_back(std::string("output.dir=") + argv[++i]);
    } else if (path.empty()) {
      path = a;
    } else {
      usage();
      return 1;
    }
  }
  if (path.empty()) {
    usage();
    return 1;
  }
  try {
    const exsqs::RunConfig cfg = exsqs::load_config(path, ovr);
    const exsqs::RunContext ctx = exsqs::RunContext::build(cfg);
    std::printf("exsqs-mini 1.0.0 | %d sites | E_floor=%.6e\n", ctx.geom.natoms(), ctx.e_floor);
    const exsqs::RunOutput out = exsqs::run_evolution(cfg, ctx);
    exsqs::write_outputs(cfg, ctx, out);
    return out.success ? 0 : 3;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "exsqs error: %s\n", e.what());
    return 1;
  }
}
