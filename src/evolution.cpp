#include "exsqs/evolution.hpp"

#include "exsqs/island_engine.hpp"
#include "exsqs/serialize.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <memory>
#include <mutex>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <tuple>
#include <unordered_set>

#include "exsqs/correlation.hpp"
#include "exsqs/dedup.hpp"
#include "exsqs/displacements.hpp"

namespace exsqs {
namespace {

namespace fs = std::filesystem;
using Clock = std::chrono::steady_clock;

// v1.2: OpenMP is a performance layer only -- results are identical for any
// thread count (counter-keyed RNG [A14] + deterministic commit ordering).
int total_threads(const RunConfig& cfg) {
#ifdef _OPENMP
  return cfg.omp_threads > 0 ? cfg.omp_threads : omp_get_max_threads();
#else
  (void)cfg;
  return 1;
#endif
}

std::mutex& checkpoint_mutex() {
  static std::mutex m;
  return m;
}

double seconds_since(const Clock::time_point& t0) {
  return std::chrono::duration<double>(Clock::now() - t0).count();
}

bool is_identity(const Mat3i& R) {
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      if (R[i][j] != (i == j ? 1 : 0)) return false;
  return true;
}

std::string key_of(const std::vector<uint8_t>& c) { return std::string(c.begin(), c.end()); }

double median_of(std::vector<double> v) {
  if (v.empty()) return 0.0;
  std::sort(v.begin(), v.end());
  const size_t n = v.size();
  return (n % 2 == 1) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

std::string json_escape(const std::string& s) {
  std::string o;
  o.reserve(s.size() + 16);
  for (unsigned char c : s) {
    switch (c) {
      case '"': o += "\\\""; break;
      case '\\': o += "\\\\"; break;
      case '\n': o += "\\n"; break;
      case '\r': o += "\\r"; break;
      case '\t': o += "\\t"; break;
      default:
        if (c < 0x20) {
          char buf[8];
          std::snprintf(buf, sizeof(buf), "\\u%04x", c);
          o += buf;
        } else {
          o += static_cast<char>(c);
        }
    }
  }
  return o;
}

// DP subset-sum with randomized backtrack: choose a subset of `sizes` summing
// to `target`; when both include/exclude are feasible the branch is a coin
// flip, giving a randomized (deterministic under CounterRng) solution.
std::vector<int> subset_pick(const std::vector<int>& sizes, int target, CounterRng& rng,
                             bool& ok) {
  ok = false;
  if (target < 0) return {};
  const int m = static_cast<int>(sizes.size());
  std::vector<std::vector<char>> reach(m + 1, std::vector<char>(target + 1, 0));
  reach[0][0] = 1;
  for (int j = 1; j <= m; ++j) {
    const int sj = sizes[j - 1];
    for (int t = 0; t <= target; ++t)
      reach[j][t] = reach[j - 1][t] || (t >= sj && reach[j - 1][t - sj]);
  }
  if (!reach[m][target]) return {};
  std::vector<int> chosen;
  int t = target;
  for (int j = m; j >= 1; --j) {
    const int sj = sizes[j - 1];
    const bool can_exclude = reach[j - 1][t] != 0;
    const bool can_include = (t >= sj) && reach[j - 1][t - sj] != 0;
    const bool take = can_include && (!can_exclude || ((rng() & 1ULL) != 0));
    if (take) {
      chosen.push_back(j - 1);
      t -= sj;
    }
  }
  ok = true;
  return chosen;
}

// Constructive decoration [D4]: pick an op with R != I from the pool, use the
// cycles of its site permutation as orbits, and subset-sum whole orbits to the
// exact composition. The decoration is invariant under that op, so the child
// carries a point-group element beyond identity => guaranteed non-P1.
std::vector<int> constructive_from_rng(const RunConfig& cfg, const RunContext& ctx,
                                       CounterRng& rng, bool& ok) {
  ok = false;
  const int N = ctx.geom.natoms();
  if (!ctx.seed_perms.empty()) {
    const size_t P = ctx.seed_perms.size();
    const size_t start = rng.below(P);
    const size_t tries = std::min<size_t>(32, P);
    for (size_t a = 0; a < tries; ++a) {
      const auto& perm = ctx.seed_perms[(start + a) % P];
      std::vector<std::vector<int>> orbs;
      std::vector<char> seen(N, 0);
      for (int i = 0; i < N; ++i)
        if (!seen[i]) {
          std::vector<int> cyc;
          int j = i;
          while (!seen[j]) {
            seen[j] = 1;
            cyc.push_back(j);
            j = perm[j];
          }
          orbs.push_back(std::move(cyc));
        }
      std::vector<int> sigma(N, 0);
      std::vector<int> avail(orbs.size());
      std::iota(avail.begin(), avail.end(), 0);
      bool good = true;
      for (int t = static_cast<int>(cfg.counts.size()) - 1; t >= 1 && good; --t) {
        std::vector<int> sizes;
        sizes.reserve(avail.size());
        for (int oi : avail) sizes.push_back(static_cast<int>(orbs[oi].size()));
        bool okk = false;
        const auto chosen = subset_pick(sizes, cfg.counts[t], rng, okk);
        if (!okk) {
          good = false;
          break;
        }
        std::vector<char> used(avail.size(), 0);
        for (int ci : chosen) {
          used[ci] = 1;
          for (int site : orbs[avail[ci]]) sigma[site] = t;
        }
        std::vector<int> navail;
        navail.reserve(avail.size());
        for (size_t q = 0; q < avail.size(); ++q)
          if (!used[q]) navail.push_back(avail[q]);
        avail.swap(navail);
      }
      if (good) {
        ok = true;
        return sigma;
      }
    }
  }
  // no feasible op found: plain shuffle from the same stream (caller filters)
  auto s = composition_sigma(cfg);
  rng.shuffle(s.begin(), s.end());
  return s;
}

struct Tally {
  int ev = 0, dup = 0, p1 = 0;
};

}  // namespace

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
  ctx.empty_info = get_symmetry(ctx.geom, cfg.symprec);
  std::set<std::vector<int>> sp;
  for (const auto& op : ctx.empty_info.ops)
    if (!is_identity(op.R)) sp.insert(permutation_of_op(ctx.geom, op.R, op.t));
  ctx.seed_perms.assign(sp.begin(), sp.end());
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

std::vector<int> seed_sigma_rejection(const RunConfig& cfg, const RunContext& ctx, int island,
                                      int slot) {
  (void)ctx;
  auto s = composition_sigma(cfg);
  CounterRng rng(cfg.seed, static_cast<uint64_t>(island), 0, static_cast<uint64_t>(slot),
                 RngPurpose::SeedInit);
  rng.shuffle(s.begin(), s.end());
  return s;
}

std::vector<int> seed_sigma_constructive(const RunConfig& cfg, const RunContext& ctx, int island,
                                         int gen, int slot, bool& constructive_ok) {
  CounterRng rng(cfg.seed, static_cast<uint64_t>(island), static_cast<uint64_t>(gen),
                 static_cast<uint64_t>(slot), RngPurpose::ConstructiveSeed);
  return constructive_from_rng(cfg, ctx, rng, constructive_ok);
}

// [D6] symmetry-preserving mutation: swap species assignments between two
// equal-size orbits (different species) of H = the parent's full stabilizer
// (spglib equivalent_atoms). Decorations constant on H-orbits stay constant,
// so the child keeps H by construction. When no such orbit pair exists the
// plain unlike-pair swap of [A12] is the fallback.
namespace {
// deterministic inverse-CDF Poisson draw (small lambda regime)
int poisson_draw(CounterRng& rng, double lam) {
  const double u = rng.uniform();
  double p = std::exp(-lam), F = p;
  int k = 0;
  while (u >= F && k < 64) {
    ++k;
    p *= lam / k;
    F += p;
  }
  return k;
}
}  // namespace

std::vector<int> mutate_sigma(const Individual& parent, const RunConfig& cfg,
                              const RunContext& ctx, CounterRng& rng) {
  const int N = ctx.geom.natoms();
  const int kswaps = (cfg.mut_poisson_lambda >= 0)
                         ? 1 + poisson_draw(rng, cfg.mut_poisson_lambda)
                         : cfg.mut_swaps;
  std::vector<int> s = parent.sigma;
  // [D6] symmetry-preserving move, revised: choose ONE nontrivial op g of the
  // parent and swap species between equal-size cycles of <g>. The child keeps
  // <g> by construction -- so it stays non-P1 and never enters the
  // mutate->P1->reject loop -- but MAY drop the rest of the parent's group.
  // (H = full stabilizer proved too rigid: symmetry could never decrease and
  // the search was trapped in high-symmetry strata; see STEP2_REPORT.)
  std::vector<std::vector<int>> orbs;
  if (cfg.mut_sympres && !parent.stab_ops.empty()) {
    for (int a = 0; a < 4 && orbs.empty(); ++a) {
      const SymOp& g = parent.stab_ops[rng.below(parent.stab_ops.size())];
      const auto perm = permutation_of_op(ctx.geom, g.R, g.t);
      std::vector<char> seen(N, 0);
      std::vector<std::vector<int>> cyc;
      for (int i = 0; i < N; ++i) {
        if (seen[i]) continue;
        std::vector<int> c;
        int j = i;
        while (!seen[j]) {
          seen[j] = 1;
          c.push_back(j);
          j = perm[j];
        }
        cyc.push_back(std::move(c));
      }
      std::map<int, std::map<int, int>> probe;  // cycle size -> species -> count
      for (const auto& c : cyc) probe[static_cast<int>(c.size())][s[c[0]]]++;
      for (const auto& kv : probe)
        if (kv.second.size() >= 2) {  // an equal-size unlike-species pair exists
          orbs = std::move(cyc);
          break;
        }
    }
  }
  for (int k = 0; k < kswaps; ++k) {
    bool done = false;
    if (!orbs.empty()) {
      // size -> species -> cycle ids (std::map keeps the ordering deterministic)
      std::map<int, std::map<int, std::vector<int>>> by;
      for (size_t oi = 0; oi < orbs.size(); ++oi)
        by[static_cast<int>(orbs[oi].size())][s[orbs[oi][0]]].push_back(static_cast<int>(oi));
      std::vector<std::tuple<int, int, int>> cand;  // (size, species a, species b)
      for (const auto& sz_sp : by)
        if (sz_sp.second.size() >= 2) {
          std::vector<int> ts;
          for (const auto& t_ids : sz_sp.second) ts.push_back(t_ids.first);
          for (size_t i = 0; i < ts.size(); ++i)
            for (size_t j = i + 1; j < ts.size(); ++j)
              cand.emplace_back(sz_sp.first, ts[i], ts[j]);
        }
      if (!cand.empty()) {
        const auto pick = cand[rng.below(cand.size())];
        const int sz = std::get<0>(pick), ta = std::get<1>(pick), tb = std::get<2>(pick);
        const auto& la = by[sz][ta];
        const auto& lb = by[sz][tb];
        const int ia = la[rng.below(la.size())];
        const int ib = lb[rng.below(lb.size())];
        for (int i : orbs[static_cast<size_t>(ia)]) s[i] = tb;
        for (int i : orbs[static_cast<size_t>(ib)]) s[i] = ta;
        done = true;
      }
    }
    if (!done) {
      const int i = static_cast<int>(rng.below(static_cast<uint64_t>(N)));
      std::vector<int> others;
      for (int j = 0; j < N; ++j)
        if (s[j] != s[i]) others.push_back(j);
      if (others.empty()) break;  // degenerate single-species case
      const int j = others[rng.below(others.size())];
      std::swap(s[i], s[j]);
    }
  }
  return s;
}

// ---------------------------------------------------------------------------
// IslandEngine (v1.2): resumable per-island state machine. The counter-keyed
// RNG [A14] carries no state, so (cfg, island, gen, population, archive, pool,
// stall counter) fully determine the trajectory -- enabling lockstep
// scheduling, synchronous migration, OpenMP execution and, later,
// checkpoint/restart. seed_generation0()/advance() are the v1.1 sequential
// algorithm transposed verbatim; a direct evolve_island() call is bit-exact
// with v1.1.
// ---------------------------------------------------------------------------
namespace {

class IslandEngine {
 public:
  IslandEngine(const RunConfig& cfg, const RunContext& ctx, int island, Clock::time_point t0,
               CheckpointFn cb)
      : cfg_(cfg), ctx_(ctx), island_(island), t0_(t0), cb_(std::move(cb)) {
    pop_.reserve(static_cast<size_t>(cfg.population));
    etol_ = effective_e_tol(cfg, ctx);
  }

  int generation() const { return gen_; }
  int island() const { return island_; }
  bool done() const { return done_; }
  IslandResult&& take_result() { return std::move(res_); }

  // ---- generation 0: mixed seeding [D4] (odd slots constructive) ----
  // v1.2 round-based execution, bit-equivalent to the v1.1 sequential loop:
  // per round exactly `remaining` slots are consumed; candidate generation and
  // canonicalization run in parallel (pure per-slot functions of the counter
  // RNG [A14]), dup-check + provisional insert and the quota commit run
  // serially in slot order.
  void seed_generation0() {
    Tally tl;
    long long attempts = 0;
    const long long cap =
        static_cast<long long>(cfg_.population) * std::max(100, cfg_.retry_budget);
    int slot = 0;
    while (static_cast<int>(pop_.size()) < cfg_.population) {
      const int B = cfg_.population - static_cast<int>(pop_.size());
      std::vector<Individual> cand(static_cast<size_t>(B));
      std::vector<char> need(static_cast<size_t>(B), 0);
#ifdef _OPENMP
#pragma omp parallel for num_threads(inner_threads_) schedule(dynamic) if (inner_threads_ > 1)
#endif
      for (int b = 0; b < B; ++b) {
        const int ts = slot + b;
        Individual& I = cand[static_cast<size_t>(b)];
        const bool constructive = cfg_.seed_mode == 1 || (cfg_.seed_mode == 2 && ts % 2 == 1);
        if (constructive) {
          bool cok = false;
          I.sigma = seed_sigma_constructive(cfg_, ctx_, island_, 0, ts, cok);
          I.origin = cok ? 'c' : 'r';
        } else {
          I.sigma = seed_sigma_rejection(cfg_, ctx_, island_, ts);
          I.origin = 'r';
        }
        I.birth_gen = 0;
        I.birth_island = island_;
        fill_canonical(I);
      }
      slot += B;
      for (int b = 0; b < B; ++b) {
        if (++attempts > cap)
          throw std::runtime_error(
              "evolution: seeding exhausted -- population too large for the unique non-P1 "
              "space?");
        if (dup_or_insert(cand[static_cast<size_t>(b)], tl)) need[static_cast<size_t>(b)] = 1;
      }
#ifdef _OPENMP
#pragma omp parallel for num_threads(inner_threads_) schedule(dynamic) if (inner_threads_ > 1)
#endif
      for (int b = 0; b < B; ++b)
        if (need[static_cast<size_t>(b)]) evaluate_individual(cand[static_cast<size_t>(b)]);
      int p1_children = 0;
      for (int b = 0; b < B; ++b) {
        if (!need[static_cast<size_t>(b)]) continue;
        Individual& I = cand[static_cast<size_t>(b)];
        if (!commit_p1(I, tl, p1_children)) continue;
        if (I.sg == 1) ++p1_children;
        pop_.push_back(std::move(I));
        pool_insert(pop_.back());
      }
    }
    record(0, 0, tl, 0);
  }

  // One termination check + one generation (extinction [A11] + repopulation
  // [A12] + record) -- the v1.1 loop body. No-op once done().
  void advance() {
    if (done_) return;

    double min_ep = 1e300;
    for (const auto& I : pop_) min_ep = std::min(min_ep, I.e_pure);
    if (min_ep <= etol_) {
      res_.success = true;
      finish("e_tol");
      return;
    }
    if (gen_ >= cfg_.max_generations) {
      finish("max_generations");
      return;
    }
    if (cfg_.max_wall_s > 0 && seconds_since(t0_) > cfg_.max_wall_s) {
      finish("wall_time");
      return;
    }

    // extinction [A11] with elitism [D5]
    const int P = static_cast<int>(pop_.size());
    std::vector<int> order(static_cast<size_t>(P));
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
      return std::tie(pop_[static_cast<size_t>(a)].e_obj, pop_[static_cast<size_t>(a)].e_pure,
                      pop_[static_cast<size_t>(a)].hash, a) <
             std::tie(pop_[static_cast<size_t>(b)].e_obj, pop_[static_cast<size_t>(b)].e_pure,
                      pop_[static_cast<size_t>(b)].hash, b);
    });
    std::vector<char> elite(static_cast<size_t>(P), 0);
    for (int e = 0; e < std::min(cfg_.elitism_best, P); ++e)
      elite[static_cast<size_t>(order[static_cast<size_t>(e)])] = 1;
    const double e_min = pop_[static_cast<size_t>(order[0])].e_obj;
    double beta = cfg_.beta;
    if (cfg_.metropolis && beta > 0 && cfg_.beta_schedule == 1)
      beta *= std::pow(cfg_.beta_growth, static_cast<double>(gen_));  // [A11] geometric, v1.5
    if (cfg_.metropolis && beta <= 0) {  // auto: median structure survives 50% [A11]
      std::vector<double> diffs;
      diffs.reserve(pop_.size());
      for (const auto& I : pop_) diffs.push_back(I.e_obj - e_min);
      const double med = median_of(std::move(diffs));
      beta = (med > 1e-300) ? (0.6931471805599453 / med) : 1e300;
    }
    std::vector<Individual> survivors;
    survivors.reserve(pop_.size());
    int killed = 0;
    for (int i = 0; i < P; ++i) {
      bool keep;
      if (elite[static_cast<size_t>(i)]) {
        keep = true;
      } else {
        const double ei = pop_[static_cast<size_t>(i)].e_obj;
        double p;
        if (!cfg_.metropolis)
          p = (ei <= 0) ? 1.0 : std::min(1.0, e_min / ei);
        else
          p = std::min(1.0, std::exp(-beta * (ei - e_min)));
        CounterRng r(cfg_.seed, static_cast<uint64_t>(island_), static_cast<uint64_t>(gen_),
                     static_cast<uint64_t>(i), RngPurpose::ExtinctionDraw);
        keep = r.uniform() < p;
      }
      if (keep)
        survivors.push_back(std::move(pop_[static_cast<size_t>(i)]));
      else
        ++killed;
    }
    pop_.swap(survivors);
    const int kept = static_cast<int>(pop_.size());

    // repopulation [A12]; RNG streams keyed by the children's birth generation.
    // v1.2 round-based execution: each unresolved slot spins (serially, in
    // slot order) to its next non-duplicate candidate (provisional archive
    // insert), the round's provisionals are evaluated in parallel, then
    // committed serially in slot order; a P1-rejected slot resumes its own
    // attempt sequence next round. With the default quota 0 and no post-seed
    // P1 events (the empirical case) this is bit-equivalent to the v1.1
    // sequential loop.
    Tally tl;
    int fb = 0;
    const int nslots = cfg_.population - kept;
    struct SlotState {
      CounterRng prng, mrng;
      std::unique_ptr<CounterRng> crng;
      int stage = 0;  // 0 mutate, 1 constructive fallback, 2 shuffle fallback
      int at_m = 0, at_c = 0, at_s = 0;
      bool resolved = false;
      Individual child;
      SlotState(CounterRng p, CounterRng m) : prng(p), mrng(m) {}
    };
    std::vector<SlotState> ss;
    ss.reserve(static_cast<size_t>(std::max(0, nslots)));
    for (int v = kept; v < cfg_.population; ++v)
      ss.emplace_back(
          CounterRng(cfg_.seed, static_cast<uint64_t>(island_), static_cast<uint64_t>(gen_) + 1,
                     static_cast<uint64_t>(v), RngPurpose::ParentPick),
          CounterRng(cfg_.seed, static_cast<uint64_t>(island_), static_cast<uint64_t>(gen_) + 1,
                     static_cast<uint64_t>(v), RngPurpose::Mutation));
    int unresolved = nslots;
    int p1_children = 0;
    std::vector<int> pend;
    while (unresolved > 0) {
      pend.clear();
      for (int q = 0; q < nslots; ++q) {
        SlotState& S = ss[static_cast<size_t>(q)];
        if (S.resolved) continue;
        bool prov = false;
        while (!prov) {
          std::vector<int> sig;
          char org = 'r';
          if (S.stage == 0) {
            if (S.at_m >= cfg_.retry_budget) {  // constructive fallback [A12]
              ++fb;
              S.stage = 1;
              S.crng.reset(new CounterRng(
                  cfg_.seed, static_cast<uint64_t>(island_), static_cast<uint64_t>(gen_) + 1,
                  static_cast<uint64_t>(kept + q), RngPurpose::ConstructiveSeed));
              continue;
            }
            const Individual& parent =
                pop_[static_cast<size_t>(S.prng.below(static_cast<uint64_t>(kept)))];
            sig = mutate_sigma(parent, cfg_, ctx_, S.mrng);
            org = 'm';
            ++S.at_m;
          } else if (S.stage == 1) {
            if (S.at_c >= cfg_.retry_budget) {
              S.stage = 2;
              continue;
            }
            bool cok = false;
            sig = constructive_from_rng(cfg_, ctx_, *S.crng, cok);
            org = cok ? 'f' : 'r';
            ++S.at_c;
          } else {
            if (S.at_s >= cfg_.retry_budget * 10)
              throw std::runtime_error(
                  "evolution: repopulation exhausted -- unique non-P1 space likely saturated; "
                  "reduce population or generations");
            sig = composition_sigma(cfg_);
            S.crng->shuffle(sig.begin(), sig.end());
            org = 'r';
            ++S.at_s;
          }
          S.child = Individual{};
          S.child.sigma = std::move(sig);
          S.child.origin = org;
          S.child.birth_gen = gen_ + 1;
          S.child.birth_island = island_;
          fill_canonical(S.child);
          if (!dup_or_insert(S.child, tl)) continue;
          prov = true;
        }
        pend.push_back(q);
      }
#ifdef _OPENMP
#pragma omp parallel for num_threads(inner_threads_) schedule(dynamic) if (inner_threads_ > 1)
#endif
      for (int j = 0; j < static_cast<int>(pend.size()); ++j)
        evaluate_individual(ss[static_cast<size_t>(pend[static_cast<size_t>(j)])].child);
      for (int j = 0; j < static_cast<int>(pend.size()); ++j) {
        SlotState& S = ss[static_cast<size_t>(pend[static_cast<size_t>(j)])];
        if (!commit_p1(S.child, tl, p1_children)) continue;  // P1 slot resumes next round
        if (S.child.sg == 1) ++p1_children;
        S.resolved = true;
        --unresolved;
      }
    }
    for (int q = 0; q < nslots; ++q) {
      pop_.push_back(std::move(ss[static_cast<size_t>(q)].child));
      pool_insert(pop_.back());
    }

    ++gen_;
    record(gen_, killed, tl, fb);
    if (cb_ && cfg_.checkpoint_every > 0 && gen_ % cfg_.checkpoint_every == 0) {
      std::lock_guard<std::mutex> lk(checkpoint_mutex());  // islands may run concurrently (v1.2)
      cb_(island_, gen_, res_.pool);
    }
    // v1.1 [A13]: stagnation stop -- repopulation degenerated to all-fallback seeding
    if (cfg_.stagnation_stop > 0 && killed > 0 && fb >= killed) {
      if (++stall_gens_ >= cfg_.stagnation_stop) {
        finish("stagnation");
        return;
      }
    } else {
      stall_gens_ = 0;
    }
  }

  // ---- synchronous ring migration (v1.2): best-k pool copies, dedup-aware ----
  std::vector<Individual> emigrants(int k) const {
    const int n = std::min<int>(k, static_cast<int>(res_.pool.size()));
    return std::vector<Individual>(res_.pool.begin(), res_.pool.begin() + n);
  }
  void note_emigrants(int n) { res_.migrants_out += n; }
  void receive_migrants(const std::vector<Individual>& in) {
    if (done_ || pop_.empty()) return;
    int acc = 0;
    for (const Individual& m : in) {
      const std::string k = key_of(m.canonical);
      if (archive_.count(k)) continue;  // receiver already explored it
      archive_.insert(k);
      // replace the current worst member; the elite [D5] is never the worst
      // while P >= 2, so E_min stays monotone
      size_t wi = 0;
      for (size_t i = 1; i < pop_.size(); ++i)
        if (std::tie(pop_[i].e_obj, pop_[i].e_pure, pop_[i].hash) >
            std::tie(pop_[wi].e_obj, pop_[wi].e_pure, pop_[wi].hash))
          wi = i;
      pop_[wi] = m;
      pool_insert(pop_[wi]);
      ++acc;
    }
    res_.migrants_in += acc;
    if (acc > 0 && cfg_.log_info)
      std::printf("[isl %d gen %4d] migration: accepted %d/%zu migrant(s)\n", island_, gen_, acc,
                  in.size());
  }


  // ---- v1.3 checkpoint/restart ----
  void serialize(ByteWriter& w) const {
    w.i32(island_);
    w.i32(gen_);
    w.i32(stall_gens_);
    w.u8(done_ ? 1 : 0);
    put_island_result(w, res_);
    put_individuals(w, pop_);
    std::vector<std::string> keys(archive_.begin(), archive_.end());
    std::sort(keys.begin(), keys.end());  // deterministic state files
    w.u64(static_cast<uint64_t>(keys.size()));
    for (const auto& k : keys) w.str(k);
  }

  void deserialize(ByteReader& r) {
    const int isl = r.i32();
    if (isl != island_) throw std::runtime_error("state: island id mismatch in blob");
    gen_ = r.i32();
    stall_gens_ = r.i32();
    done_ = r.u8() != 0;
    res_ = get_island_result(r);
    pop_ = get_individuals(r);
    archive_.clear();
    const uint64_t na = r.u64();
    archive_.reserve(static_cast<size_t>(na) * 2);
    for (uint64_t i = 0; i < na; ++i) archive_.insert(r.str());
  }

  bool reopen_for_resume() {
    if (done_ &&
        (res_.stop_reason == "max_generations" || res_.stop_reason == "wall_time")) {
      done_ = false;  // raisable caps re-arm; e_tol/stagnation stay terminal
      res_.stop_reason.clear();
    }
    return !done_;
  }

 private:
  void finish(const char* why) {
    res_.stop_reason = why;
    res_.generations = gen_;
    done_ = true;
  }

  int p1_in_pop() const {
    int c = 0;
    for (const auto& I : pop_)
      if (I.sg == 1) ++c;
    return c;
  }

  // Filtration + evaluation pipeline, split for the v1.2 round-based parallel
  // execution: fill_canonical (parallel-safe, pure), dup_or_insert (serial;
  // provisional archive insert -- correct for BOTH later outcomes, since v1.1
  // archived admitted AND P1-rejected candidates alike), evaluate_individual
  // (parallel-safe, one spglib call), commit_p1 (serial quota check [A7][A8]).
  void fill_canonical(Individual& out) const {
    out.canonical = canonical_labels(out.sigma, ctx_.perms);
    out.hash = hash_labels(out.canonical);
  }

  bool dup_or_insert(const Individual& out, Tally& tl) {
    const std::string k = key_of(out.canonical);
    if (archive_.count(k)) {
      ++tl.dup;
      return false;
    }
    archive_.insert(k);
    return true;
  }

  void evaluate_individual(Individual& out) const {
    const Structure s = decorate(ctx_.geom, out.sigma, cfg_.species);
    const SymmetryInfo info = get_symmetry(s, cfg_.symprec);
    out.sg = info.sg_number;
    out.sg_symbol = info.sg_symbol;
    out.eq_atoms = info.equivalent_atoms;
    out.stab_ops.clear();
    for (const auto& op : info.ops)
      if (!is_identity(op.R)) out.stab_ops.push_back(op);
    out.D = displacement_count(s, info, cfg_.symprec);
    const CorrData cd = count_pairs(s, ctx_.zones);
    out.e_pure = cfg_.full_pairs ? e_pure_full(cd, cfg_.x_achieved, ctx_.weights)
                                 : e_pure_diagonal(cd, cfg_.x_achieved, ctx_.weights);
    out.e_obj = out.e_pure * std::pow(static_cast<double>(out.D), cfg_.gamma);
  }

  bool commit_p1(const Individual& out, Tally& tl, int p1_children) {
    ++tl.ev;
    if (out.sg == 1 && p1_in_pop() + p1_children >= cfg_.p1_elite_quota) {
      ++tl.p1;  // canonical already archived at dup_or_insert time (= v1.1)
      return false;
    }
    return true;
  }

  void pool_insert(const Individual& I) {
    res_.pool.push_back(I);
    std::sort(res_.pool.begin(), res_.pool.end(), [](const Individual& a, const Individual& b) {
      return std::tie(a.e_obj, a.e_pure, a.hash) < std::tie(b.e_obj, b.e_pure, b.hash);
    });
    if (static_cast<int>(res_.pool.size()) > cfg_.outputs)
      res_.pool.resize(static_cast<size_t>(cfg_.outputs));
  }

  void record(int gen, int killed, const Tally& tl, int fb) {
    GenStats st;
    st.island = island_;
    st.gen = gen;
    double bp = 1e300;
    size_t bi = 0;
    for (size_t i = 0; i < pop_.size(); ++i) {
      bp = std::min(bp, pop_[i].e_pure);
      if (std::tie(pop_[i].e_obj, pop_[i].e_pure, pop_[i].hash) <
          std::tie(pop_[bi].e_obj, pop_[bi].e_pure, pop_[bi].hash))
        bi = i;
    }
    st.best_e_pure = bp;
    st.best_e_obj = pop_[bi].e_obj;
    st.best_D = pop_[bi].D;
    st.best_sg = pop_[bi].sg;
    std::vector<double> eo;
    eo.reserve(pop_.size());
    for (const auto& I : pop_) eo.push_back(I.e_obj);
    st.median_e_obj = median_of(std::move(eo));
    st.killed = killed;
    st.evals = tl.ev;
    st.dup_rejected = tl.dup;
    st.p1_rejected = tl.p1;
    st.fallback_seeds = fb;
    st.elapsed_s = seconds_since(t0_);
    res_.log.push_back(st);
    res_.evals += tl.ev;
    if (cfg_.log_info)
      std::printf(
          "[isl %d gen %4d] E_pure*=%.6e E_obj*=%.6e D*=%4d SG*=%3d med=%.3e kill=%3d ev=%3d "
          "dup=%d p1=%d fb=%d t=%.1fs\n",
          island_, gen, st.best_e_pure, st.best_e_obj, st.best_D, st.best_sg, st.median_e_obj,
          st.killed, st.evals, st.dup_rejected, st.p1_rejected, st.fallback_seeds, st.elapsed_s);
  }

  const RunConfig& cfg_;
  const RunContext& ctx_;
  int island_ = 0;
  Clock::time_point t0_;
  CheckpointFn cb_;
  double etol_ = 0;
  std::unordered_set<std::string> archive_;  // [A9][A10] canonical labels ever seen
  std::vector<Individual> pop_;
  IslandResult res_;
  int gen_ = 0;
  int stall_gens_ = 0;
  bool done_ = false;
  int inner_threads_ = 1;

 public:
  void set_inner_threads(int n) { inner_threads_ = std::max(1, n); }
};

}  // namespace

// ---------------------------------------------------------------------------
// EngineHandle (v1.3): thin public pimpl over IslandEngine for external
// drivers (exsqs_mpi, tests) and checkpoint/restart.
// ---------------------------------------------------------------------------
struct EngineHandle::Impl {
  IslandEngine e;
  Impl(const RunConfig& cfg, const RunContext& ctx, int island, CheckpointFn cb)
      : e(cfg, ctx, island, Clock::now(), std::move(cb)) {}
};

EngineHandle::EngineHandle(const RunConfig& cfg, const RunContext& ctx, int island,
                           CheckpointFn cb)
    : impl_(new Impl(cfg, ctx, island, std::move(cb))) {}
EngineHandle::EngineHandle(EngineHandle&&) noexcept = default;
EngineHandle& EngineHandle::operator=(EngineHandle&&) noexcept = default;
EngineHandle::~EngineHandle() = default;
void EngineHandle::set_inner_threads(int n) { impl_->e.set_inner_threads(n); }
void EngineHandle::seed_generation0() { impl_->e.seed_generation0(); }
void EngineHandle::advance() { impl_->e.advance(); }
bool EngineHandle::done() const { return impl_->e.done(); }
int EngineHandle::generation() const { return impl_->e.generation(); }
int EngineHandle::island() const { return impl_->e.island(); }
std::vector<Individual> EngineHandle::emigrants(int k) const { return impl_->e.emigrants(k); }
void EngineHandle::receive_migrants(const std::vector<Individual>& in) {
  impl_->e.receive_migrants(in);
}
void EngineHandle::note_emigrants(int n) { impl_->e.note_emigrants(n); }
IslandResult EngineHandle::take_result() { return impl_->e.take_result(); }
void EngineHandle::serialize(ByteWriter& w) const { impl_->e.serialize(w); }
void EngineHandle::deserialize(ByteReader& r) { impl_->e.deserialize(r); }
bool EngineHandle::reopen_for_resume() { return impl_->e.reopen_for_resume(); }

void save_run_state_blobs(const std::string& path, const RunConfig& cfg,
                          const std::vector<std::string>& blobs_by_island) {
  ensure_little_endian();
  ByteWriter w;
  w.raw("EXSQSTAT", 8);
  w.u32(1);  // format version
  w.str(trajectory_signature(cfg));
  w.u32(static_cast<uint32_t>(blobs_by_island.size()));
  for (size_t i = 0; i < blobs_by_island.size(); ++i) {
    w.u32(static_cast<uint32_t>(i));
    w.str(blobs_by_island[i]);
  }
  namespace fs = std::filesystem;
  const fs::path fp(path);
  if (fp.has_parent_path()) fs::create_directories(fp.parent_path());
  const std::string tmp = path + ".tmp";
  {
    std::ofstream f(tmp, std::ios::binary | std::ios::trunc);
    if (!f) throw std::runtime_error("state: cannot write " + tmp);
    f.write(w.data().data(), static_cast<std::streamsize>(w.data().size()));
  }
  fs::rename(tmp, fp);
}

void save_run_state(const std::string& path, const RunConfig& cfg,
                    const std::vector<EngineHandle>& engines) {
  std::vector<std::string> blobs(engines.size());
  for (const auto& e : engines) {
    ByteWriter b;
    e.serialize(b);
    blobs[static_cast<size_t>(e.island())] = b.data();
  }
  save_run_state_blobs(path, cfg, blobs);
}

int load_run_state(const std::string& path, const RunConfig& cfg,
                   std::vector<EngineHandle>& engines) {
  ensure_little_endian();
  std::ifstream f(path, std::ios::binary);
  if (!f) throw std::runtime_error("state: cannot open " + path);
  const std::string buf((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
  ByteReader r(buf);
  char magic[8];
  r.raw_read(magic, 8);
  if (std::memcmp(magic, "EXSQSTAT", 8) != 0) throw std::runtime_error("state: bad magic");
  if (r.u32() != 1) throw std::runtime_error("state: unsupported format version");
  if (r.str() != trajectory_signature(cfg))
    throw std::runtime_error(
        "state: trajectory signature mismatch -- the config differs from the checkpointed run "
        "in a trajectory-relevant field (only budget caps may change on resume)");
  const uint32_t nisl = r.u32();
  if (nisl != engines.size() || static_cast<int>(nisl) != cfg.islands)
    throw std::runtime_error("state: island count mismatch");
  std::vector<char> seen(nisl, 0);
  for (uint32_t k = 0; k < nisl; ++k) {
    const uint32_t id = r.u32();
    const std::string blob = r.str();
    if (id >= nisl || seen[id]) throw std::runtime_error("state: duplicate/invalid island record");
    ByteReader br(blob);
    engines[id].deserialize(br);
    seen[id] = 1;
  }
  int round = 0;
  for (auto& e : engines)
    if (e.reopen_for_resume()) round = std::max(round, e.generation());
  return round;
}

IslandResult evolve_island(const RunConfig& cfg, const RunContext& ctx, int island,
                           const CheckpointFn& cb) {
  EngineHandle e(cfg, ctx, island, cb);
  e.set_inner_threads(total_threads(cfg));
  e.seed_generation0();
  while (!e.done()) e.advance();
  return e.take_result();
}

RunOutput merge_island_results(const RunConfig& cfg, std::vector<IslandResult>&& rs) {
  RunOutput out;
  out.islands = cfg.islands;
  std::vector<Individual> pooled;
  for (auto& r : rs) {
    out.success = out.success || r.success;
    out.evals += r.evals;
    out.island_generations.push_back(r.generations);
    out.island_success.push_back(r.success ? 1 : 0);
    out.island_stop.push_back(r.stop_reason);
    out.island_migrants_in.push_back(r.migrants_in);
    out.island_migrants_out.push_back(r.migrants_out);
    for (auto& g : r.log) out.log.push_back(g);
    for (auto& I : r.pool) pooled.push_back(std::move(I));
  }
  std::sort(pooled.begin(), pooled.end(), [](const Individual& a, const Individual& b) {
    return std::tie(a.e_obj, a.e_pure, a.hash) < std::tie(b.e_obj, b.e_pure, b.hash);
  });
  std::unordered_set<std::string> seen;
  for (auto& I : pooled) {
    const std::string k = key_of(I.canonical);
    if (seen.count(k)) continue;
    seen.insert(k);
    out.outputs.push_back(std::move(I));
    if (static_cast<int>(out.outputs.size()) >= cfg.outputs) break;
  }
  return out;
}

RunOutput run_evolution(const RunConfig& cfg, const RunContext& ctx, const CheckpointFn& cb,
                        const std::string& resume_from, const std::string& state_dir) {
  const auto t0 = Clock::now();
  std::vector<EngineHandle> eng;
  eng.reserve(static_cast<size_t>(cfg.islands));
  for (int i = 0; i < cfg.islands; ++i) eng.emplace_back(cfg, ctx, i, cb);

  // Thread budget (v1.2): outer over islands, remainder inside each island's
  // evaluation rounds. Thread counts affect wall time only, never results.
  const int tt = total_threads(cfg);
  const int outer = std::min<int>(cfg.islands, std::max(1, tt));
  const int inner = std::max(1, tt / std::max(1, outer));
#ifdef _OPENMP
  if (outer > 1 && inner > 1) omp_set_max_active_levels(2);
#endif
  for (auto& e : eng) e.set_inner_threads(inner);

  int round = 0;
  if (!resume_from.empty()) {
    round = load_run_state(resume_from, cfg, eng);  // v1.3: bit-exact resume
  } else {
    std::vector<std::exception_ptr> errs(eng.size());
#ifdef _OPENMP
#pragma omp parallel for num_threads(outer) schedule(static, 1)
#endif
    for (int i = 0; i < static_cast<int>(eng.size()); ++i) {
      try {
        eng[static_cast<size_t>(i)].seed_generation0();
      } catch (...) {
        errs[static_cast<size_t>(i)] = std::current_exception();
      }
    }
    for (auto& ep : errs)
      if (ep) std::rethrow_exception(ep);
  }

  const auto save_state = [&]() {
    if (!state_dir.empty()) save_run_state(state_dir + "/state.ckpt", cfg, eng);
  };

  // Lockstep generation rounds (v1.2); synchronous ring migration (SPEC 8.1).
  while (true) {
    bool any = false;
    for (auto& e : eng)
      if (!e.done()) any = true;
    if (!any) break;
    std::vector<std::exception_ptr> errs(eng.size());
    {
#ifdef _OPENMP
#pragma omp parallel for num_threads(outer) schedule(static, 1)
#endif
      for (int i = 0; i < static_cast<int>(eng.size()); ++i) {
        try {
          eng[static_cast<size_t>(i)].advance();
        } catch (...) {
          errs[static_cast<size_t>(i)] = std::current_exception();
        }
      }
    }
    for (auto& ep : errs)
      if (ep) std::rethrow_exception(ep);
    ++round;
    if (cfg.islands > 1 && cfg.migrants > 0 && cfg.migration_every > 0 &&
        round % cfg.migration_every == 0) {
      std::vector<int> act;
      for (int i = 0; i < cfg.islands; ++i)
        if (!eng[static_cast<size_t>(i)].done()) act.push_back(i);
      if (act.size() >= 2) {
        std::vector<std::vector<Individual>> sends;
        sends.reserve(act.size());
        for (int i : act) sends.push_back(eng[static_cast<size_t>(i)].emigrants(cfg.migrants));
        for (size_t j = 0; j < act.size(); ++j) {
          eng[static_cast<size_t>(act[(j + 1) % act.size()])].receive_migrants(sends[j]);
          eng[static_cast<size_t>(act[j])].note_emigrants(static_cast<int>(sends[j].size()));
        }
      }
    }
    if (!state_dir.empty() && cfg.checkpoint_every > 0 && round % cfg.checkpoint_every == 0)
      save_state();
  }
  save_state();  // final state enables budget-raising resume chains [A13]

  std::vector<IslandResult> rs;
  rs.reserve(eng.size());
  for (auto& e : eng) rs.push_back(e.take_result());
  RunOutput out = merge_island_results(cfg, std::move(rs));
  out.wall_s = seconds_since(t0);
  return out;
}

void write_outputs(const RunConfig& cfg, const RunContext& ctx, const RunOutput& out,
                   bool checkpoint) {
  fs::create_directories(cfg.outdir);

  {
    std::ofstream csv(cfg.outdir + "/generations.csv");
    csv << "island,gen,best_e_pure,best_e_obj,median_e_obj,best_D,best_sg,killed,evals,"
           "dup_rejected,p1_rejected,fallback_seeds,elapsed_s\n";
    char buf[256];
    for (const auto& g : out.log) {
      std::snprintf(buf, sizeof(buf), "%d,%d,%.10e,%.10e,%.10e,%d,%d,%d,%d,%d,%d,%d,%.3f\n",
                    g.island, g.gen, g.best_e_pure, g.best_e_obj, g.median_e_obj, g.best_D,
                    g.best_sg, g.killed, g.evals, g.dup_rejected, g.p1_rejected,
                    g.fallback_seeds, g.elapsed_s);
      csv << buf;
    }
  }

  // POSCARs + section-5 per-structure symmetry report (SG, pg order, ineq sites)
  struct Extra {
    int pg = 0;
    int ineq = 0;
  };
  std::vector<Extra> extra(out.outputs.size());
  for (size_t i = 0; i < out.outputs.size(); ++i) {
    const auto& I = out.outputs[i];
    const Structure s = decorate(ctx.geom, I.sigma, cfg.species);
    const SymmetryInfo info = get_symmetry(s, cfg.symprec);
    extra[i].pg = pointgroup_order(info);
    extra[i].ineq = static_cast<int>(info.independent_atoms().size());
    char name[64];
    std::snprintf(name, sizeof(name), "best_%02zu.vasp", i);
    char comment[256];
    std::snprintf(comment, sizeof(comment),
                  "exsqs E_pure=%.6e D=%d E_obj=%.6e SG=%d (%s) pg_order=%d ineq_sites=%d",
                  I.e_pure, I.D, I.e_obj, I.sg, I.sg_symbol.c_str(), extra[i].pg, extra[i].ineq);
    write_poscar(grouped_by_species(s), cfg.outdir + "/" + name, comment);
  }

  std::ofstream j(cfg.outdir + (checkpoint ? "/checkpoint.json" : "/summary.json"));
  char nb[64];
  auto num = [&nb](double v) {
    std::snprintf(nb, sizeof(nb), "%.17g", v);
    return std::string(nb);
  };
  j << "{\n";
  j << "  \"exsqs_version\": \"1.6.0\",\n";
  j << "  \"checkpoint\": " << (checkpoint ? "true" : "false") << ",\n";
  j << "  \"success\": " << (out.success ? "true" : "false") << ",\n";
  j << "  \"islands\": " << out.islands << ",\n";
  j << "  \"island_generations\": [";
  for (size_t i = 0; i < out.island_generations.size(); ++i)
    j << (i ? "," : "") << out.island_generations[i];
  j << "],\n  \"island_stop\": [";
  for (size_t i = 0; i < out.island_stop.size(); ++i)
    j << (i ? "," : "") << '"' << json_escape(out.island_stop[i]) << '"';
  j << "],\n";
  j << "  \"island_migrants_in\": [";
  for (size_t i = 0; i < out.island_migrants_in.size(); ++i)
    j << (i ? "," : "") << out.island_migrants_in[i];
  j << "],\n  \"island_migrants_out\": [";
  for (size_t i = 0; i < out.island_migrants_out.size(); ++i)
    j << (i ? "," : "") << out.island_migrants_out[i];
  j << "],\n";
  j << "  \"total_evaluations\": " << out.evals << ",\n";
  j << "  \"wall_s\": " << num(out.wall_s) << ",\n";
  j << "  \"species\": [";
  for (size_t t = 0; t < cfg.species.size(); ++t)
    j << (t ? "," : "") << '"' << json_escape(cfg.species[t]) << '"';
  j << "],\n  \"counts\": [";
  for (size_t t = 0; t < cfg.counts.size(); ++t) j << (t ? "," : "") << cfg.counts[t];
  j << "],\n  \"x_achieved\": [";
  for (size_t t = 0; t < cfg.x_achieved.size(); ++t) j << (t ? "," : "") << num(cfg.x_achieved[t]);
  j << "],\n";
  j << "  \"e_tol\": " << num(cfg.e_tol) << ",\n";
  j << "  \"e_tol_effective\": " << num(effective_e_tol(cfg, ctx)) << ",\n";
  j << "  \"e_floor\": " << num(ctx.e_floor) << ",\n";
  j << "  \"gamma\": " << num(cfg.gamma) << ",\n";
  j << "  \"error_mode\": \"" << (cfg.full_pairs ? "full_pairs" : "diagonal") << "\",\n";
  j << "  \"weights\": [";
  for (size_t n = 0; n < ctx.weights.size(); ++n) j << (n ? "," : "") << num(ctx.weights[n]);
  j << "],\n";
  j << "  \"n_shells\": " << ctx.zones.n_shells << ",\n";
  j << "  \"shell_tol\": " << num(cfg.shell_tol) << ",\n";
  j << "  \"shell_radii\": [";
  for (size_t n = 0; n < ctx.zones.radii.size(); ++n)
    j << (n ? "," : "") << num(ctx.zones.radii[n]);
  j << "],\n";
  j << "  \"shells_exceed_half_width\": "
    << (ctx.zones.shells_exceed_half_width ? "true" : "false") << ",\n";
  j << "  \"symprec\": " << num(cfg.symprec) << ",\n";
  j << "  \"seed\": " << cfg.seed << ",\n";
  j << "  \"spglib_version\": \"" << json_escape(spglib_version()) << "\",\n";
  j << "  \"displacement_convention\": \"phonopy_default (ported from phonopy 4.3.1) [A15]\",\n";
  j << "  \"outputs\": [\n";
  for (size_t i = 0; i < out.outputs.size(); ++i) {
    const auto& I = out.outputs[i];
    char name[64];
    std::snprintf(name, sizeof(name), "best_%02zu.vasp", i);
    j << "    {\"file\": \"" << name << "\", \"e_pure\": " << num(I.e_pure) << ", \"D\": " << I.D
      << ", \"e_obj\": " << num(I.e_obj) << ", \"sg\": " << I.sg << ", \"sg_symbol\": \""
      << json_escape(I.sg_symbol) << "\", \"pg_order\": " << extra[i].pg
      << ", \"ineq_sites\": " << extra[i].ineq << ", \"origin\": \"" << I.origin
      << "\", \"birth_island\": " << I.birth_island << ", \"birth_gen\": " << I.birth_gen << "}"
      << (i + 1 < out.outputs.size() ? "," : "") << "\n";
  }
  j << "  ],\n";
  j << "  \"config_echo\": \"" << json_escape(cfg.config_echo) << "\"\n";
  j << "}\n";
}

int run_from_config(const RunConfig& cfg, const std::string& resume_from) {
  std::printf("exsqs 1.6.0 | spglib %s\n", spglib_version().c_str());
  const RunContext ctx = RunContext::build(cfg);
  const int N = ctx.geom.natoms();
  std::printf("system: %d sites |", N);
  for (size_t t = 0; t < cfg.species.size(); ++t)
    std::printf(" %s=%d (x=%.6f)", cfg.species[t].c_str(), cfg.counts[t], cfg.x_achieved[t]);
  std::printf("\nzones: %d shells, radii(A):", ctx.zones.n_shells);
  for (double r : ctx.zones.radii) std::printf(" %.4f", r);
  std::printf("%s\n", ctx.zones.shells_exceed_half_width
                          ? "  [WARNING: outer shells exceed half cell width, A3]"
                          : "");
  std::printf("dedup group |Pi| = %zu perms; constructive pool = %zu ops (R != I)\n",
              ctx.perms.size(), ctx.seed_perms.size());
  std::printf("mode=%s gamma=%.3g survival=%s seeding=%s sympres=%s e_tol=%.3g pop=%d islands=%d\n",
              cfg.full_pairs ? "full_pairs" : "diagonal", cfg.gamma,
              cfg.metropolis ? "metropolis" : "ratio",
              cfg.seed_mode == 0 ? "rejection" : (cfg.seed_mode == 1 ? "constructive" : "mixed"),
              cfg.mut_sympres ? "on" : "off", effective_e_tol(cfg, ctx), cfg.population, cfg.islands);
  std::printf("E_floor = %.6e (L1 quantization bound, SPEC 4.1)%s\n", ctx.e_floor,
              cfg.e_tol < 0 ? "  [e_tol auto -> 3.0 x E_floor]" : "");
  if (effective_e_tol(cfg, ctx) < ctx.e_floor)
    std::fprintf(stderr,
                 "exsqs warning: e_tol=%.3g is below E_floor=%.3g (infeasible); "
                 "the run will terminate on generation/wall caps only\n",
                 effective_e_tol(cfg, ctx), ctx.e_floor);
  std::printf("reference: D(P1, %d sites) = %d\n", N, 6 * N);

  const CheckpointFn cb = [&cfg, &ctx](int island, int gen, const std::vector<Individual>& pool) {
    RunOutput partial;
    partial.islands = cfg.islands;
    partial.outputs = pool;
    partial.island_generations = {gen};
    partial.island_stop = {"checkpoint@island" + std::to_string(island)};
    write_outputs(cfg, ctx, partial, true);
  };

  if (!resume_from.empty())
    std::printf("resuming from %s\n", resume_from.c_str());
  RunOutput out = run_evolution(cfg, ctx, cb, resume_from, cfg.outdir);
  write_outputs(cfg, ctx, out, false);

  std::printf("\n== %s | evals=%lld | wall=%.1fs ==\n",
              out.success ? "SUCCESS (min E_pure <= e_tol)"
                          : "budget exhausted (best-effort output) [A13]",
              out.evals, out.wall_s);
  for (size_t i = 0; i < out.outputs.size(); ++i) {
    const auto& I = out.outputs[i];
    std::printf("  best_%02zu: E_pure=%.6e D=%4d (P1:%d) E_obj=%.6e SG=%d (%s)\n", i, I.e_pure,
                I.D, 6 * N, I.e_obj, I.sg, I.sg_symbol.c_str());
  }
  if (!out.outputs.empty() && ctx.e_floor > 0)
    std::printf("E_floor=%.6e | best E_pure/E_floor = %.2f\n", ctx.e_floor,
                out.outputs[0].e_pure / ctx.e_floor);
  std::printf("outputs written to %s\n", cfg.outdir.c_str());
  return out.success ? 0 : 3;
}

}  // namespace exsqs
