#include "exsqs/serialize.hpp"

#include <cstring>
#include <stdexcept>

namespace exsqs {

void ensure_little_endian() {
  const uint32_t probe = 0x01020304u;
  uint8_t b[4];
  std::memcpy(b, &probe, 4);
  if (b[0] != 0x04)
    throw std::runtime_error("serialize: big-endian host unsupported for state/wire format v1");
}

void ByteReader::raw_read(void* dst, size_t n) {
  const char* p = take(n);
  std::memcpy(dst, p, n);
}

const char* ByteReader::take(size_t n) {
  if (static_cast<size_t>(end_ - p_) < n)
    throw std::runtime_error("serialize: truncated or corrupt payload");
  const char* p = p_;
  p_ += n;
  return p;
}

void put_individual(ByteWriter& w, const Individual& I) {
  w.ints(I.sigma);
  w.bytes(I.canonical);
  w.u64(I.hash);
  w.i32(I.sg);
  w.str(I.sg_symbol);
  w.ints(I.eq_atoms);
  w.u32(static_cast<uint32_t>(I.stab_ops.size()));
  for (const SymOp& op : I.stab_ops) {
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) w.i32(op.R[i][j]);
    for (int j = 0; j < 3; ++j) w.f64(op.t[j]);
  }
  w.f64(I.e_pure);
  w.f64(I.e_obj);
  w.i32(I.D);
  w.i32(I.birth_gen);
  w.i32(I.birth_island);
  w.u8(static_cast<uint8_t>(I.origin));
}

Individual get_individual(ByteReader& r) {
  Individual I;
  I.sigma = r.ints();
  I.canonical = r.bytes();
  I.hash = r.u64();
  I.sg = r.i32();
  I.sg_symbol = r.str();
  I.eq_atoms = r.ints();
  const uint32_t nops = r.u32();
  I.stab_ops.resize(nops);
  for (uint32_t k = 0; k < nops; ++k) {
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) I.stab_ops[k].R[i][j] = r.i32();
    for (int j = 0; j < 3; ++j) I.stab_ops[k].t[j] = r.f64();
  }
  I.e_pure = r.f64();
  I.e_obj = r.f64();
  I.D = r.i32();
  I.birth_gen = r.i32();
  I.birth_island = r.i32();
  I.origin = static_cast<char>(r.u8());
  return I;
}

void put_individuals(ByteWriter& w, const std::vector<Individual>& v) {
  w.u32(static_cast<uint32_t>(v.size()));
  for (const auto& I : v) put_individual(w, I);
}

std::vector<Individual> get_individuals(ByteReader& r) {
  const uint32_t n = r.u32();
  std::vector<Individual> v;
  v.reserve(n);
  for (uint32_t i = 0; i < n; ++i) v.push_back(get_individual(r));
  return v;
}

void put_genstats(ByteWriter& w, const GenStats& g) {
  w.i32(g.island);
  w.i32(g.gen);
  w.f64(g.best_e_pure);
  w.f64(g.best_e_obj);
  w.f64(g.median_e_obj);
  w.i32(g.best_D);
  w.i32(g.best_sg);
  w.i32(g.killed);
  w.i32(g.evals);
  w.i32(g.dup_rejected);
  w.i32(g.p1_rejected);
  w.i32(g.fallback_seeds);
  w.f64(g.elapsed_s);
}

GenStats get_genstats(ByteReader& r) {
  GenStats g;
  g.island = r.i32();
  g.gen = r.i32();
  g.best_e_pure = r.f64();
  g.best_e_obj = r.f64();
  g.median_e_obj = r.f64();
  g.best_D = r.i32();
  g.best_sg = r.i32();
  g.killed = r.i32();
  g.evals = r.i32();
  g.dup_rejected = r.i32();
  g.p1_rejected = r.i32();
  g.fallback_seeds = r.i32();
  g.elapsed_s = r.f64();
  return g;
}

void put_island_result(ByteWriter& w, const IslandResult& r) {
  put_individuals(w, r.pool);
  w.u32(static_cast<uint32_t>(r.log.size()));
  for (const auto& g : r.log) put_genstats(w, g);
  w.u8(r.success ? 1 : 0);
  w.i32(r.generations);
  w.i64(r.evals);
  w.i32(r.migrants_in);
  w.i32(r.migrants_out);
  w.str(r.stop_reason);
}

IslandResult get_island_result(ByteReader& rd) {
  IslandResult r;
  r.pool = get_individuals(rd);
  const uint32_t nl = rd.u32();
  r.log.reserve(nl);
  for (uint32_t i = 0; i < nl; ++i) r.log.push_back(get_genstats(rd));
  r.success = rd.u8() != 0;
  r.generations = rd.i32();
  r.evals = rd.i64();
  r.migrants_in = rd.i32();
  r.migrants_out = rd.i32();
  r.stop_reason = rd.str();
  return r;
}

std::string trajectory_signature(const RunConfig& c) {
  ByteWriter w;
  w.raw("EXSQSIG1", 8);
  // system geometry + composition
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) w.f64(c.proto.cell[i][j]);
  w.u32(static_cast<uint32_t>(c.proto.sites.size()));
  for (const auto& s : c.proto.sites)
    for (int j = 0; j < 3; ++j) w.f64(s[j]);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) w.i32(c.H[i][j]);
  w.u32(static_cast<uint32_t>(c.species.size()));
  for (const auto& sp : c.species) w.str(sp);
  w.ints(c.counts);
  for (double x : c.x_achieved) w.f64(x);
  // zones / error model
  w.i32(c.n_shells);
  w.f64(c.shell_tol);
  w.i32(static_cast<int32_t>(c.wform));
  w.f64(c.wpow);
  w.u32(static_cast<uint32_t>(c.wcustom.size()));
  for (double x : c.wcustom) w.f64(x);
  w.u8(c.full_pairs ? 1 : 0);
  w.f64(c.gamma);
  // evolution dynamics
  w.i32(c.population);
  w.i32(c.outputs);
  w.f64(c.e_tol);
  w.u8(c.metropolis ? 1 : 0);
  w.f64(c.beta);
  w.i32(c.beta_schedule);  // v1.5 [A11]
  w.f64(c.beta_growth);
  w.i32(c.elitism_best);
  w.i32(c.p1_elite_quota);
  w.i32(c.mut_swaps);
  w.f64(c.mut_poisson_lambda);
  w.u8(c.mut_sympres ? 1 : 0);
  w.i32(c.retry_budget);
  w.i32(c.seed_mode);
  w.i32(c.stagnation_stop);
  // symmetry / rng / topology
  w.f64(c.symprec);
  w.u64(c.seed);
  w.i32(c.islands);
  w.i32(c.migration_every);
  w.i32(c.migrants);
  return w.data();
}

}  // namespace exsqs
