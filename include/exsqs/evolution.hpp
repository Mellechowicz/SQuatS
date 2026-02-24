#pragma once
// exsqs-mini: sequential extinction evolutionary SQS engine (minimal build).
// One island, random seeding, plain unlike-pair swap mutation, ratio or
// metropolis extinction, canonical-label deduplication,
// objective E_obj = E_pure * D^gamma and the quantization floor of SPEC 4.1.

#include <cstdint>
#include <string>
#include <vector>

#include "exsqs/config.hpp"
#include "exsqs/structure.hpp"
#include "exsqs/symmetry.hpp"
#include "exsqs/zones.hpp"

namespace exsqs {

struct RunContext {
  Structure geom;                        // undecorated supercell
  ZoneTable zones;                       // coordination shells
  std::vector<double> weights;           // A_n
  double e_floor = 0.0;                  // provable L1 quantization bound
  std::vector<std::vector<int>> perms;   // full site-permutation group [A9]
  SymmetryInfo empty_info;                   // ops of the empty supercell
  std::vector<std::vector<int>> seed_perms;  // perms of ops with R != I (non-P1 guarantee)
  static RunContext build(const RunConfig& cfg);
};

struct Individual {
  std::vector<int> sigma;
  std::vector<uint8_t> canonical;
  double e_pure = 0.0;
  double e_obj = 0.0;
  int D = 0;
  int sg = 0;
  std::string sg_symbol;
};

struct GenStat {
  int gen = 0;
  double best = 0.0;  // best E_pure in the population after survival/refill
  long evals = 0;     // cumulative evaluations
  int killed = 0;
};

struct RunOutput {
  std::vector<Individual> outputs;  // top-M by E_obj, symmetry-unique
  std::vector<GenStat> log;
  long evals = 0;
  int generations = 0;
  bool success = false;  // reached e_tol
};

double effective_e_tol(const RunConfig& cfg, const RunContext& ctx);
std::vector<int> composition_sigma(const RunConfig& cfg);
std::vector<int> seed_sigma_rejection(const RunConfig& cfg, const RunContext& ctx, int island,
                                      int slot);
std::vector<int> seed_sigma_constructive(const RunConfig& cfg, const RunContext& ctx, int island,
                                         int gen, int slot, bool& constructive_ok);
Individual evaluate(const RunConfig& cfg, const RunContext& ctx, std::vector<int> sigma);
RunOutput run_evolution(const RunConfig& cfg, const RunContext& ctx);
void write_outputs(const RunConfig& cfg, const RunContext& ctx, const RunOutput& out);

}  // namespace exsqs
