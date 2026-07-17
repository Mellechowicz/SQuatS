#pragma once
#include <cstdint>
#include <string>
#include <vector>

#include "exsqs/correlation.hpp"
#include "exsqs/lattice.hpp"
#include "exsqs/linalg.hpp"

namespace exsqs {

// Parsed + validated form of the SPEC section-11 YAML schema. All defaults are
// the frozen ones; composition is resolved to counts by largest remainder [A5].
struct RunConfig {
  // system
  Prototype proto;
  Mat3i H{};
  std::vector<std::string> species;
  std::vector<double> x_target;
  std::vector<double> x_achieved;
  std::vector<int> counts;

  // zones / error
  int n_shells = 7;
  double shell_tol = 1e-3;
  WeightForm wform = WeightForm::InvN;
  double wpow = 1.0;
  std::vector<double> wcustom;
  bool full_pairs = false;  // resolved from error.mode auto [A16]
  double gamma = 1.0;       // E_obj = E_pure * D^gamma [D3]
  // v1.9 multiplet sectors (SPEC 4.2): E = E_2 + l3*E_3 + l4*E_4; both
  // lambdas 0 (the default) reproduces the pair-only engine bitwise.
  double lambda3 = 0.0;
  double lambda4 = 0.0;
  int mshell3 = 2;  // zones spanned by the triplet cutoff radius
  int mshell4 = 1;  // zones spanned by the quadruplet cutoff radius

  // evolution
  int population = 200;
  int outputs = 10;  // M
  double e_tol = -1.0;  // < 0 => auto = 3.0 * E_floor (v1.1, SPEC 4.1)
  int max_generations = 5000;
  bool metropolis = false;  // survival mode [A11]
  double beta = -1.0;       // <= 0 => auto = ln2 / median(E_i - E_min)
  int beta_schedule = 0;    // [A11] v1.5: 0 = const, 1 = geometric (beta_g = beta*growth^g)
  double beta_growth = 1.0;  // [A11] v1.5: geometric factor per generation (> 0)
  int elitism_best = 1;     // [D5]
  int p1_elite_quota = 0;   // [A8]
  int mut_swaps = 1;
  double mut_poisson_lambda = -1.0;  // >=0: k ~ 1+Poisson(lambda) [A12]
  bool mut_sympres = true;  // [D6]
  int retry_budget = 100;   // [A12]
  int stagnation_stop = 3;  // v1.1 [A13]: stop after S all-fallback repopulations (0 = off)
  int seed_mode = 2;        // 0 rejection | 1 constructive | 2 mixed [D4]

  // symmetry / displacements
  double symprec = 1e-5;  // [A6]; filter policy fixed to reject_p1 [A7]

  // parallel (v1.2: lockstep islands, synchronous ring migration, OpenMP)
  int islands = 1;
  int omp_threads = -1;     // -1 = auto (OpenMP default)
  int migration_every = 50; // rounds between synchronous ring migrations (0 = off)
  int migrants = 2;         // pool-best copies sent per island per migration (0 = off)

  // rng / output
  uint64_t seed = 42;
  std::string outdir = "./run1";
  int checkpoint_every = 100;
  double max_wall_s = -1.0;  // < 0 => unlimited [A13]
  bool log_info = true;      // per-generation stdout rows

  std::string config_echo;  // raw YAML text for summary.json
};

// Load + validate a YAML config; overrides are "a.b.c=value" strings applied
// to the document before parsing (CLI-overridable per section 11).
RunConfig load_config(const std::string& path, const std::vector<std::string>& overrides = {});

}  // namespace exsqs
