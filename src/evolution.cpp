#include "exsqs/evolution.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <set>
#include <stdexcept>

#include "exsqs/correlation.hpp"
#include "exsqs/dedup.hpp"
#include "exsqs/displacements.hpp"
#include "exsqs/lattice.hpp"
#include "exsqs/rng.hpp"

namespace exsqs {

RunContext RunContext::build(const RunConfig& cfg) {
  RunContext ctx;
  ctx.geom = make_supercell(cfg.proto, cfg.H, cfg.species, 0);
  const int N = ctx.geom.natoms();
  int csum = 0;
  for (int c : cfg.counts) csum += c;
  if (csum != N) throw std::runtime_error("evolution: composition counts do not sum to N");
  ctx.zones = build_zones(ctx.geom, cfg.n_shells, cfg.shell_tol);
  ctx.weights = make_weights(cfg.wform, ctx.zones, cfg.wpow, cfg.wcustom);
  ctx.e_floor = cfg.full_pairs
                    ? e_floor_full(ctx.zones, cfg.counts, cfg.x_achieved, ctx.weights)
                    : e_floor_diagonal(ctx.zones, cfg.counts, cfg.x_achieved, ctx.weights);
  ctx.perms = site_permutations(ctx.geom, cfg.symprec);
  return ctx;
}

double effective_e_tol(const RunConfig& cfg, const RunContext& ctx) {
  return cfg.e_tol < 0 ? 3.0 * ctx.e_floor : cfg.e_tol;
}

std::vector<int> composition_sigma(const RunConfig& cfg) {
  std::vector<int> s;
  for (size_t t = 0; t < cfg.counts.size(); ++t)
    for (int k = 0; k < cfg.counts[t]; ++k) s.push_back(static_cast<int>(t));
  return s;
}

Individual evaluate(const RunConfig& cfg, const RunContext& ctx, std::vector<int> sigma) {
  Individual I;
  const Structure dec = decorate(ctx.geom, sigma, cfg.species);
  const CorrData cd = count_pairs(dec, ctx.zones);
  I.e_pure = cfg.full_pairs ? e_pure_full(cd, cfg.x_achieved, ctx.weights)
                            : e_pure_diagonal(cd, cfg.x_achieved, ctx.weights);
  const SymmetryInfo info = get_symmetry(dec, cfg.symprec);
  I.sg = info.sg_number;
  I.sg_symbol = info.sg_symbol;
  I.D = displacement_count(dec, info, cfg.symprec);
  I.e_obj = I.e_pure * std::pow(static_cast<double>(I.D), cfg.gamma);
  I.sigma = std::move(sigma);
  return I;
}

namespace {

void consider_output(std::vector<Individual>& top, const Individual& I, int M) {
  top.push_back(I);
  std::sort(top.begin(), top.end(),
            [](const Individual& a, const Individual& b) { return a.e_obj < b.e_obj; });
  if (static_cast<int>(top.size()) > M) top.resize(static_cast<size_t>(M));
}

}  // namespace

RunOutput run_evolution(const RunConfig& cfg, const RunContext& ctx) {
  RunOutput out;
  std::set<std::vector<uint8_t>> archive;
  std::vector<Individual> pop;
  const double tol = effective_e_tol(cfg, ctx);

  // admit: canonical dedup [A10] -> evaluate. The minimal build admits P1
  // decorations (random seeds at large N are always P1); the D^gamma term of
  // the objective still rewards symmetry whenever gamma > 0. The strict P1
  // rejection filter and constructive space-group seeding are extensions.
  auto admit = [&](std::vector<int> sigma, Individual& I) -> bool {
    std::vector<uint8_t> labels = canonical_labels(sigma, ctx.perms);
    if (archive.count(labels)) return false;
    archive.insert(labels);
    out.evals++;
    I = evaluate(cfg, ctx, std::move(sigma));
    I.canonical = std::move(labels);
    consider_output(out.outputs, I, cfg.outputs);
    return true;
  };

  // seeding: random shuffles of the exact composition [A5]
  {
    uint64_t slot = 0;
    long budget = static_cast<long>(cfg.retry_budget) * cfg.population;
    while (static_cast<int>(pop.size()) < cfg.population && budget-- > 0) {
      std::vector<int> s = composition_sigma(cfg);
      CounterRng rng(cfg.seed, 0, 0, slot++, RngPurpose::SeedInit);
      rng.shuffle(s.begin(), s.end());
      Individual I;
      if (admit(std::move(s), I)) pop.push_back(std::move(I));
    }
  }
  if (pop.empty()) throw std::runtime_error("evolution: seeding produced no admissible individual");

  int stagnant = 0;
  int g = 0;
  for (g = 1; g <= cfg.max_generations; ++g) {
    std::sort(pop.begin(), pop.end(),
              [](const Individual& a, const Individual& b) { return a.e_obj < b.e_obj; });
    if (pop.front().e_pure <= tol) {
      out.success = true;
      break;
    }
    // extinction [A11]: ratio keeps the better half; metropolis draws survival
    const size_t before = pop.size();
    std::vector<Individual> next;
    if (cfg.metropolis) {
      double beta = cfg.beta;
      if (beta <= 0) {
        const double emin = pop.front().e_obj;
        std::vector<double> gaps;
        for (const Individual& I : pop) gaps.push_back(I.e_obj - emin);
        std::nth_element(gaps.begin(), gaps.begin() + gaps.size() / 2, gaps.end());
        const double med = std::max(1e-300, gaps[gaps.size() / 2]);
        beta = std::log(2.0) / med;
      }
      for (size_t i = 0; i < pop.size(); ++i) {
        CounterRng rng(cfg.seed, 0, static_cast<uint64_t>(g), i, RngPurpose::ExtinctionDraw);
        const double p = std::exp(-beta * (pop[i].e_obj - pop.front().e_obj));
        if (static_cast<int>(i) < cfg.elitism_best || rng.uniform() < p)
          next.push_back(std::move(pop[i]));
      }
    } else {
      const size_t keep = std::max<size_t>(static_cast<size_t>(cfg.elitism_best),
                                           pop.size() / 2);
      for (size_t i = 0; i < keep && i < pop.size(); ++i) next.push_back(std::move(pop[i]));
    }
    const int killed = static_cast<int>(before - next.size());
    pop.swap(next);

    // repopulation: unlike-pair swap mutation of uniformly picked parents [A12]
    bool any_admitted = false;
    uint64_t slot = 0;
    int budget = cfg.retry_budget;
    while (static_cast<int>(pop.size()) < cfg.population && budget > 0) {
      CounterRng pick(cfg.seed, 0, static_cast<uint64_t>(g), slot, RngPurpose::ParentPick);
      const Individual& parent =
          pop[static_cast<size_t>(pick.uniform() * static_cast<double>(pop.size())) % pop.size()];
      std::vector<int> child = parent.sigma;
      CounterRng mut(cfg.seed, 0, static_cast<uint64_t>(g), slot, RngPurpose::Mutation);
      for (int sw = 0; sw < std::max(1, cfg.mut_swaps); ++sw) {
        for (int attempt = 0; attempt < 64; ++attempt) {
          const size_t i = static_cast<size_t>(mut.uniform() * child.size()) % child.size();
          const size_t j = static_cast<size_t>(mut.uniform() * child.size()) % child.size();
          if (child[i] != child[j]) {
            std::swap(child[i], child[j]);
            break;
          }
        }
      }
      ++slot;
      Individual I;
      if (admit(std::move(child), I)) {
        pop.push_back(std::move(I));
        any_admitted = true;
      } else {
        --budget;
      }
    }
    std::sort(pop.begin(), pop.end(),
              [](const Individual& a, const Individual& b) { return a.e_obj < b.e_obj; });
    out.log.push_back({g, pop.front().e_pure, out.evals, killed});
    if (cfg.log_info)
      std::printf("[gen %4d] E_pure*=%.6e E_obj*=%.6e D*=%4d SG*=%3d kill=%3d ev=%ld\n", g,
                  pop.front().e_pure, pop.front().e_obj, pop.front().D, pop.front().sg, killed,
                  out.evals);
    if (!any_admitted) {
      if (cfg.stagnation_stop > 0 && ++stagnant >= cfg.stagnation_stop) break;  // [A13]
    } else {
      stagnant = 0;
    }
  }
  out.generations = std::min(g, cfg.max_generations);
  return out;
}

}  // namespace exsqs
