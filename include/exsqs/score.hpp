#pragma once
// v1.5 (step 6): scoring of externally supplied structures under a run config
// -- the executable half of the T-V1 doctrine: absolute error scalars are not
// comparable across codes (SPEC 4.1), raw structures are. `exsqs score`
// evaluates any POSCAR on the config's geometry with the config's zones,
// weights, error mode and gamma, reporting E_pure, E/E_floor, D, SG and E_obj.
// Input site ORDER is free; sites are matched by fractional coordinate.

#include <string>
#include <vector>

#include "exsqs/config.hpp"
#include "exsqs/evolution.hpp"
#include "exsqs/structure.hpp"

namespace exsqs {

struct ScoreResult {
  std::string file;
  double e_pure = 0.0;  // total: pair + lambda-weighted sectors (v1.9)
  double e_obj = 0.0;
  double e_pair = 0.0;  // pair sector alone
  double e3 = 0.0, e4 = 0.0;  // multiplet sectors (0 when off)
  int D = 0;
  int sg = 0;
  std::string sg_symbol;
  int pg_order = 0;
};

// Maps an external structure onto ctx.geom: the lattice must match elementwise
// (tolerance 1e-6 relative + 1e-9 absolute) and every geometry site must match
// exactly one input site by wrapped fractional distance (<= 1e-4). Species
// names must belong to cfg.species and the composition must equal cfg.counts.
// Returns sigma in geometry site order; throws with a clear message otherwise.
std::vector<int> sigma_on_geometry(const RunConfig& cfg, const RunContext& ctx,
                                   const Structure& s);

ScoreResult score_structure(const RunConfig& cfg, const RunContext& ctx, const Structure& s);

// CLI backend: exsqs score <config.yaml> [--set k=v ...] <POSCAR> [...]
//              [--json PATH]
// Prints a table (and optional JSON); returns 0 iff every file scored.
int run_score_cli(const RunConfig& cfg, const std::vector<std::string>& files,
                  const std::string& json_path);

}  // namespace exsqs
