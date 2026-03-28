// v1.2 parallel-layer tests: spglib/pipeline thread safety (T-P0), thread-count
// invariance incl. migration (T-P1, the "1 vs 16 threads" clause of T-R1), and
// ring-migration semantics (T-M1).
#include <catch2/catch_test_macros.hpp>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "exsqs/correlation.hpp"
#include "exsqs/displacements.hpp"
#include "exsqs/evolution.hpp"
#include "exsqs/lattice.hpp"
#include "exsqs/rng.hpp"
#include "exsqs/structure.hpp"
#include "exsqs/symmetry.hpp"

using namespace exsqs;

namespace {

RunConfig par_cfg() {
  RunConfig c;
  c.proto = make_sc(3.0);
  c.H = {{{3, 0, 0}, {0, 3, 0}, {0, 0, 3}}};
  c.species = {"W", "Cr"};
  c.counts = {14, 13};
  c.x_target = {14.0 / 27.0, 13.0 / 27.0};
  c.x_achieved = {14.0 / 27.0, 13.0 / 27.0};
  c.n_shells = 5;
  c.gamma = 1.0;
  c.population = 12;
  c.outputs = 4;
  c.e_tol = 0.0;  // unreachable
  c.max_generations = 10;
  c.retry_budget = 60;
  c.seed = 7;
  c.log_info = false;
  return c;
}

struct EvalRef {
  int sg = 0;
  size_t nops = 0;
  int D = 0;
  std::vector<int> eq;
  double e = 0;
};

EvalRef eval_one(const Structure& s, const ZoneTable& zt, const std::vector<double>& x,
                 const std::vector<double>& w) {
  EvalRef r;
  const SymmetryInfo info = get_symmetry(s, 1e-5);
  r.sg = info.sg_number;
  r.nops = info.ops.size();
  r.eq = info.equivalent_atoms;
  r.D = displacement_count(s, info, 1e-5);
  r.e = e_pure_diagonal(count_pairs(s, zt), x, w);
  return r;
}

bool same_output(const RunOutput& a, const RunOutput& b) {
  if (a.success != b.success || a.evals != b.evals) return false;
  if (a.log.size() != b.log.size() || a.outputs.size() != b.outputs.size()) return false;
  for (size_t i = 0; i < a.log.size(); ++i) {
    const auto& x = a.log[i];
    const auto& y = b.log[i];
    if (x.island != y.island || x.gen != y.gen || x.best_e_pure != y.best_e_pure ||
        x.best_e_obj != y.best_e_obj || x.median_e_obj != y.median_e_obj ||
        x.best_D != y.best_D || x.best_sg != y.best_sg || x.killed != y.killed ||
        x.evals != y.evals || x.dup_rejected != y.dup_rejected ||
        x.p1_rejected != y.p1_rejected || x.fallback_seeds != y.fallback_seeds)
      return false;
  }
  for (size_t i = 0; i < a.outputs.size(); ++i)
    if (a.outputs[i].canonical != b.outputs[i].canonical ||
        a.outputs[i].e_pure != b.outputs[i].e_pure || a.outputs[i].D != b.outputs[i].D ||
        a.outputs[i].sg != b.outputs[i].sg)
      return false;
  if (a.island_migrants_in != b.island_migrants_in) return false;
  return true;
}

}  // namespace

TEST_CASE("T-P0: symmetry/displacement/correlation pipeline is thread-safe", "[parallel]") {
  const RunConfig c = par_cfg();
  const Structure geom = make_supercell(c.proto, c.H, c.species, 0);
  const ZoneTable zt = build_zones(geom, c.n_shells, c.shell_tol);
  const auto w = make_weights(WeightForm::InvN, zt, 1.0, {});
  const int M = 48;
  std::vector<Structure> cases;
  cases.reserve(M);
  for (int r = 0; r < M; ++r) {
    std::vector<int> sig;
    for (size_t t = 0; t < c.counts.size(); ++t) sig.insert(sig.end(), c.counts[t], (int)t);
    CounterRng rng(999, 0, 0, (uint64_t)r, RngPurpose::Generic);
    rng.shuffle(sig.begin(), sig.end());
    cases.push_back(decorate(geom, sig, c.species));
  }
  std::vector<EvalRef> serial(M);
  for (int i = 0; i < M; ++i) serial[i] = eval_one(cases[i], zt, c.x_achieved, w);
  std::vector<EvalRef> par(M);
#ifdef _OPENMP
#pragma omp parallel for num_threads(4) schedule(dynamic)
#endif
  for (int i = 0; i < M; ++i) par[i] = eval_one(cases[i], zt, c.x_achieved, w);
  for (int i = 0; i < M; ++i) {
    REQUIRE(par[i].sg == serial[i].sg);
    REQUIRE(par[i].nops == serial[i].nops);
    REQUIRE(par[i].eq == serial[i].eq);
    REQUIRE(par[i].D == serial[i].D);
    REQUIRE(par[i].e == serial[i].e);
  }
}

TEST_CASE("T-P1: results are independent of omp_threads (islands + migration)", "[parallel]") {
  RunConfig c = par_cfg();
  c.islands = 3;
  c.migration_every = 3;
  c.migrants = 1;
  const RunContext ctx = RunContext::build(c);
  c.omp_threads = 1;
  const RunOutput a = run_evolution(c, ctx);
  c.omp_threads = 4;
  const RunOutput b = run_evolution(c, ctx);
  REQUIRE(same_output(a, b));

  // single-island path too
  RunConfig s = par_cfg();
  s.islands = 1;
  const RunContext ctx1 = RunContext::build(s);
  s.omp_threads = 1;
  const RunOutput a1 = run_evolution(s, ctx1);
  s.omp_threads = 4;
  const RunOutput b1 = run_evolution(s, ctx1);
  REQUIRE(same_output(a1, b1));
}

TEST_CASE("T-M1: synchronous ring migration exchanges pool-best individuals", "[parallel]") {
  RunConfig c = par_cfg();
  c.islands = 3;
  c.migration_every = 3;
  c.migrants = 1;
  const RunContext ctx = RunContext::build(c);
  const RunOutput out = run_evolution(c, ctx);
  REQUIRE(out.island_migrants_in.size() == 3);
  int tin = 0, tout = 0;
  for (int v : out.island_migrants_in) tin += v;
  for (int v : out.island_migrants_out) tout += v;
  REQUIRE(tout >= 2);        // exchanges happened at rounds 3, 6, 9
  REQUIRE(tin >= 1);         // at least one migrant was novel to its receiver
  REQUIRE(tin <= tout);      // dedup can only drop, never invent
  // determinism of the migrating system (rerun)
  const RunOutput out2 = run_evolution(c, ctx);
  REQUIRE(same_output(out, out2));
}
