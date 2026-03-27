// T-B1 scaling section (v1.2): wall time of a fixed workload vs omp_threads
// and islands. Recorded, no gate. Thread counts never change results [A14] --
// the best_e_pure column doubles as an invariance check.
#include <chrono>
#include <cstdio>

#include "exsqs/evolution.hpp"
#include "exsqs/lattice.hpp"

using namespace exsqs;
using Clock = std::chrono::steady_clock;

static RunConfig base_cfg(int islands) {
  RunConfig c;
  c.proto = make_bcc(3.165);
  c.H = {{{4, 0, 0}, {0, 4, 0}, {0, 0, 4}}};
  c.species = {"W", "Cr"};
  c.counts = {90, 38};
  c.x_target = {0.70, 0.30};
  c.x_achieved = {90.0 / 128.0, 38.0 / 128.0};
  c.n_shells = 7;
  c.gamma = 0.0;
  c.population = 48;
  c.outputs = 4;
  c.e_tol = 1e-12;  // unreachable: fixed workload
  c.max_generations = 6;
  c.islands = islands;
  c.metropolis = true;
  c.beta = 3000.0;
  c.migration_every = 0;
  c.seed = 42;
  c.log_info = false;
  return c;
}

int main() {
  std::printf("T-B1 scaling (fixed workload: W90Cr38 bcc 4x4x4, pop 48, 6 generations)\n");
  std::printf("%-8s %-12s %-9s %-14s\n", "islands", "omp_threads", "wall_s", "best_e_pure");
  for (int isl : {1, 4}) {
    RunConfig c = base_cfg(isl);
    const RunContext ctx = RunContext::build(c);
    for (int t : {1, 2, 4}) {
      c.omp_threads = t;
      const auto t0 = Clock::now();
      const RunOutput o = run_evolution(c, ctx);
      const double w = std::chrono::duration<double>(Clock::now() - t0).count();
      std::printf("%-8d %-12d %-9.2f %-14.6e\n", isl, t, w,
                  o.outputs.empty() ? 0.0 : o.outputs[0].e_pure);
    }
  }
  std::printf("(single-core containers show ~flat wall times; rerun on a multicore node)\n");
  return 0;
}
