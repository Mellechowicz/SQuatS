#pragma once
// v1.3: public handle over the internal island engine (defined in
// evolution.cpp). One island's complete evolutionary state; because the
// counter-keyed RNG [A14] is stateless, (population, archive, pool,
// generation, stall counter) serialize to a bit-exact resume. Consumed by
// run_evolution (serial lockstep driver), exsqs_mpi (rank-distributed lockstep
// driver, SPEC 8.2) and the checkpoint tests.

#include <memory>
#include <string>
#include <vector>

#include "exsqs/config.hpp"
#include "exsqs/evolution.hpp"
#include "exsqs/serialize.hpp"

namespace exsqs {

class EngineHandle {
 public:
  EngineHandle(const RunConfig& cfg, const RunContext& ctx, int island, CheckpointFn cb = {});
  EngineHandle(EngineHandle&&) noexcept;
  EngineHandle& operator=(EngineHandle&&) noexcept;
  ~EngineHandle();

  void set_inner_threads(int n);
  void seed_generation0();
  void advance();
  bool done() const;
  int generation() const;
  int island() const;

  // synchronous ring migration (SPEC 8.1)
  std::vector<Individual> emigrants(int k) const;
  void receive_migrants(const std::vector<Individual>& in);
  void note_emigrants(int n);

  IslandResult take_result();

  // ---- checkpoint/restart (v1.3) ----
  void serialize(ByteWriter& w) const;
  void deserialize(ByteReader& r);  // full state; throws on corruption/mismatch
  // Re-arms an island stopped on a raisable cap (max_generations, wall_time);
  // e_tol success and stagnation stay terminal. Returns true if now active.
  bool reopen_for_resume();

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace exsqs
