#include <cstdio>
#include <string>
#include <vector>

#include "exsqs/config.hpp"
#include "exsqs/evolution.hpp"
#include "exsqs/score.hpp"

static void usage() {
  std::printf(
      "usage: exsqs <config.yaml> [--set key.path=value ...] [--out DIR] [--resume DIR]\n"
      "  Extinction-based evolutionary SQS search (arXiv:2602.10872 workflow; SPEC v1.0).\n"
      "  --set   override any config key, e.g. --set evolution.population=64\n"
      "  --out    shorthand for --set output.dir=DIR\n"
      "  --resume continue from DIR/state.ckpt (v1.3; bit-exact, budget caps may be raised)\n"
      "subcommand (v1.5): exsqs score <config> [--set ...] <POSCAR>... [--json PATH]\n"
      "  scores external structures on the config geometry: E_pure, E/E_floor, D, SG, E_obj\n"
      "exit codes: 0 success (min E_pure <= e_tol) | 3 budget exhausted | 1 error\n");
}

int main(int argc, char** argv) {
  if (argc >= 2 && std::string(argv[1]) == "score") {
    std::string cfgpath, json;
    std::vector<std::string> ovr2, files;
    for (int i = 2; i < argc; ++i) {
      const std::string a = argv[i];
      if (a == "--set" && i + 1 < argc) {
        ovr2.push_back(argv[++i]);
      } else if (a == "--json" && i + 1 < argc) {
        json = argv[++i];
      } else if (cfgpath.empty()) {
        cfgpath = a;
      } else {
        files.push_back(a);
      }
    }
    if (cfgpath.empty() || files.empty()) {
      usage();
      return 1;
    }
    try {
      const exsqs::RunConfig cfg = exsqs::load_config(cfgpath, ovr2);
      return exsqs::run_score_cli(cfg, files, json);
    } catch (const std::exception& e) {
      std::fprintf(stderr, "exsqs error: %s\n", e.what());
      return 1;
    }
  }
  std::string path;
  std::vector<std::string> ovr;
  std::string resume_dir;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "-h" || a == "--help") {
      usage();
      return 0;
    }
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
    } else if (a == "--resume") {
      if (i + 1 >= argc) {
        usage();
        return 1;
      }
      resume_dir = argv[++i];
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
    return exsqs::run_from_config(
        cfg, resume_dir.empty() ? std::string() : resume_dir + "/state.ckpt");
  } catch (const std::exception& e) {
    std::fprintf(stderr, "exsqs error: %s\n", e.what());
    return 1;
  }
}
