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

// Whole-run state file (atomic tmp+rename): magic "EXSQSTAT", format version,
// trajectory signature, island count, then (island id, engine blob) records.
void save_run_state(const std::string& path, const RunConfig& cfg,
                    const std::vector<EngineHandle>& engines);

// Same file format from pre-serialized per-island blobs (index = island id);
// used by the MPI driver's rank-0 writer.
void save_run_state_blobs(const std::string& path, const RunConfig& cfg,
                          const std::vector<std::string>& blobs_by_island);

// Loads into pre-constructed engines (size must equal cfg.islands); verifies
// the trajectory signature (throws on mismatch -- refuses to resume a
// different run) and re-arms raisable stops. Returns the lockstep round to
// resume from (= max generation over the still-active islands, 0 if none).
int load_run_state(const std::string& path, const RunConfig& cfg,
                   std::vector<EngineHandle>& engines);

}  // namespace exsqs
