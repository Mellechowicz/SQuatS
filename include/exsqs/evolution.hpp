#pragma once
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "exsqs/config.hpp"
#include "exsqs/rng.hpp"
#include "exsqs/structure.hpp"
#include "exsqs/symmetry.hpp"
#include "exsqs/zones.hpp"

namespace exsqs {

// One population member. `canonical` is the [A10] lex-min label (exact dedup
// key; strictly safer than the 128-bit-hash variant the spec allows), and
// eq_atoms is the spglib orbit map reused by symmetry-preserving mutation [D6].
struct Individual {
  std::vector<int> sigma;
  std::vector<uint8_t> canonical;
  uint64_t hash = 0;
  int sg = 0;
  std::string sg_symbol;
  std::vector<int> eq_atoms;
  std::vector<SymOp> stab_ops;  // nontrivial ops (R != I) of the decorated cell;
                                // mutation picks one and preserves its cyclic subgroup [D6]
  double e_pure = 0.0;
  double e_obj = 0.0;  // E_pure * D^gamma [D3]
  int D = 0;
  int birth_gen = 0;
  int birth_island = 0;
  char origin = '?';  // r rejection seed | c constructive seed | m mutant | f fallback seed
};

struct GenStats {
  int island = 0, gen = 0;
  double best_e_pure = 0.0, best_e_obj = 0.0, median_e_obj = 0.0;
  int best_D = 0, best_sg = 0;
  int killed = 0, evals = 0, dup_rejected = 0, p1_rejected = 0, fallback_seeds = 0;
  double elapsed_s = 0.0;
};

struct IslandResult {
  std::vector<Individual> pool;  // best <= M pairwise-inequivalent, ranked by E_obj
  std::vector<GenStats> log;
  bool success = false;  // min E_pure <= e_tol reached [D7]
  int generations = 0;
  long long evals = 0;
  std::string stop_reason;  // e_tol | max_generations | wall_time
};

struct RunOutput {
  std::vector<Individual> outputs;  // pooled + deduped across islands [D8]
  std::vector<GenStats> log;
  bool success = false;
  long long evals = 0;
  double wall_s = 0.0;
  int islands = 1;
  std::vector<int> island_generations;
  std::vector<int> island_success;
  std::vector<std::string> island_stop;
};

// Immutable per-run context shared by all islands (geometry, zones, weights,
// the [A9] dedup permutation group, and the [D4] constructive-seeding op pool).
struct RunContext {
  Structure geom;  // undecorated supercell
  ZoneTable zones;
  std::vector<double> weights;
  double e_floor = 0.0;  // v1.1 (SPEC 4.1): L1 quantization lower bound
  std::vector<std::vector<int>> perms;       // empty-cell site permutations [A9]
  SymmetryInfo empty_info;                   // ops of the empty supercell
  std::vector<std::vector<int>> seed_perms;  // perms of ops with R != I (non-P1 guarantee)
  static RunContext build(const RunConfig& cfg);
};

// v1.1: resolves e_tol 'auto' (cfg.e_tol < 0) to 3.0 * ctx.e_floor.
double effective_e_tol(const RunConfig& cfg, const RunContext& ctx);

// Building blocks (exposed for unit tests):
std::vector<int> composition_sigma(const RunConfig& cfg);
std::vector<int> seed_sigma_rejection(const RunConfig& cfg, const RunContext& ctx, int island,
                                      int slot);
std::vector<int> seed_sigma_constructive(const RunConfig& cfg, const RunContext& ctx, int island,
                                         int gen, int slot, bool& constructive_ok);
std::vector<int> mutate_sigma(const Individual& parent, const RunConfig& cfg,
                              const RunContext& ctx, CounterRng& rng);

using CheckpointFn = std::function<void(int island, int gen, const std::vector<Individual>& pool)>;

IslandResult evolve_island(const RunConfig& cfg, const RunContext& ctx, int island,
                           const CheckpointFn& cb = {});
RunOutput run_evolution(const RunConfig& cfg, const RunContext& ctx, const CheckpointFn& cb = {});

// Merge per-island results into the pooled, deduplicated top-M output [D8].
RunOutput merge_island_results(const RunConfig& cfg, std::vector<IslandResult>&& rs);
void write_outputs(const RunConfig& cfg, const RunContext& ctx, const RunOutput& out,
                   bool checkpoint);
// Full serial driver (banner, run, outputs). Exit code: 0 success, 3 budget.
int run_from_config(const RunConfig& cfg);

}  // namespace exsqs
