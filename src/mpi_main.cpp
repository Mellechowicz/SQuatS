// exsqs_mpi (v1.3, SPEC 8.2): rank-distributed lockstep driver.
//
// Islands are assigned round-robin (island i -> rank i % R). Trajectories are
// island-keyed [A14] and migration is deterministic (SPEC 8.1), so results are
// bit-identical to the single-process `exsqs` for ANY rank count (T-MPI1).
// The wire and state formats are the v1.3 little-endian encoding; a
// homogeneous cluster is assumed (ensure_little_endian guards every boundary).
// Checkpoint/restart: rank 0 writes the same state.ckpt as the serial driver;
// on --resume every rank reads it from the shared filesystem and keeps its
// owned islands.

#include <mpi.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <exception>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "exsqs/config.hpp"
#include "exsqs/evolution.hpp"
#include "exsqs/island_engine.hpp"
#include "exsqs/serialize.hpp"

using Clock = std::chrono::steady_clock;

namespace {

void usage() {
  std::printf(
      "usage: mpirun -n R exsqs_mpi <config.yaml> [--set k.p=v ...] [--out DIR]\n"
      "  Rank-distributed islands (SPEC 8.2); results are bit-identical to `exsqs`\n"
      "  for any rank count. exit codes: 0 success | 3 budget exhausted | 1 error\n");
}

int local_threads(const exsqs::RunConfig& cfg) {
#ifdef _OPENMP
  return cfg.omp_threads > 0 ? cfg.omp_threads : omp_get_max_threads();
#else
  (void)cfg;
  return 1;
#endif
}

// Allgatherv of variable-length local byte buffers -> concatenation on all ranks.
std::string allgather_bytes(const std::string& local) {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  const int n = static_cast<int>(local.size());
  std::vector<int> counts(static_cast<size_t>(size)), displs(static_cast<size_t>(size));
  MPI_Allgather(&n, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
  int total = 0;
  for (int i = 0; i < size; ++i) {
    displs[static_cast<size_t>(i)] = total;
    total += counts[static_cast<size_t>(i)];
  }
  std::string out(static_cast<size_t>(total), '\0');
  MPI_Allgatherv(local.data(), n, MPI_BYTE, out.empty() ? nullptr : &out[0], counts.data(),
                 displs.data(), MPI_BYTE, MPI_COMM_WORLD);
  return out;
}

// Gatherv to rank 0; other ranks receive "".
std::string gather_bytes_root(const std::string& local, int rank) {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  const int n = static_cast<int>(local.size());
  std::vector<int> counts(static_cast<size_t>(size)), displs(static_cast<size_t>(size));
  MPI_Gather(&n, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
  int total = 0;
  if (rank == 0)
    for (int i = 0; i < size; ++i) {
      displs[static_cast<size_t>(i)] = total;
      total += counts[static_cast<size_t>(i)];
    }
  std::string out(static_cast<size_t>(total), '\0');
  MPI_Gatherv(local.data(), n, MPI_BYTE, out.empty() ? nullptr : &out[0], counts.data(),
              displs.data(), MPI_BYTE, 0, MPI_COMM_WORLD);
  return out;
}

// parse a concatenated (u32 island id, str blob) record stream
std::vector<std::string> parse_records(const std::string& all, int islands) {
  std::vector<std::string> blobs(static_cast<size_t>(islands));
  exsqs::ByteReader r(all);
  while (r.remaining() > 0) {
    const uint32_t id = r.u32();
    if (static_cast<int>(id) >= islands) throw std::runtime_error("mpi: bad island record id");
    blobs[id] = r.str();
  }
  return blobs;
}

}  // namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int rank = 0, size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  std::string path;
  std::vector<std::string> ovr;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "-h" || a == "--help") {
      if (rank == 0) usage();
      MPI_Finalize();
      return 0;
    }
    if (a == "--set" && i + 1 < argc) {
      ovr.push_back(argv[++i]);
    } else if (a == "--out" && i + 1 < argc) {
      ovr.push_back(std::string("output.dir=") + argv[++i]);
    } else if (path.empty()) {
      path = a;
    } else {
      if (rank == 0) usage();
      MPI_Finalize();
      return 1;
    }
  }
  if (path.empty()) {
    if (rank == 0) usage();
    MPI_Finalize();
    return 1;
  }

  try {
    exsqs::ensure_little_endian();
    exsqs::RunConfig cfg = exsqs::load_config(path, ovr);
    if (rank != 0) cfg.log_info = false;  // per-generation rows from rank 0's islands only
    const exsqs::RunContext ctx = exsqs::RunContext::build(cfg);
    const auto t0 = Clock::now();

    std::vector<int> owned;
    for (int i = 0; i < cfg.islands; ++i)
      if (i % size == rank) owned.push_back(i);

    std::vector<exsqs::EngineHandle> eng;
    eng.reserve(static_cast<size_t>(cfg.islands));
    for (int i = 0; i < cfg.islands; ++i) eng.emplace_back(cfg, ctx, i);

    const int T = local_threads(cfg);
    const int outer = std::max(1, std::min<int>(static_cast<int>(owned.size()), T));
    const int inner = std::max(1, T / outer);
#ifdef _OPENMP
    if (outer > 1 && inner > 1) omp_set_max_active_levels(2);
#endif
    for (int i : owned) eng[static_cast<size_t>(i)].set_inner_threads(inner);

    if (rank == 0) {
      std::printf("exsqs_mpi 1.3.0 | ranks=%d islands=%d (round-robin)\n", size, cfg.islands);
      std::printf("E_floor = %.6e | e_tol(effective) = %.4g%s\n", ctx.e_floor,
                  exsqs::effective_e_tol(cfg, ctx), cfg.e_tol < 0 ? " [auto]" : "");
    }

    int round = 0;
    {
      std::vector<std::exception_ptr> errs(owned.size());
#ifdef _OPENMP
#pragma omp parallel for num_threads(outer) schedule(static, 1)
#endif
      for (int k = 0; k < static_cast<int>(owned.size()); ++k) {
        try {
          eng[static_cast<size_t>(owned[static_cast<size_t>(k)])].seed_generation0();
        } catch (...) {
          errs[static_cast<size_t>(k)] = std::current_exception();
        }
      }
      for (auto& ep : errs)
        if (ep) std::rethrow_exception(ep);
    }

    const auto save_state_root = [&]() {
      exsqs::ByteWriter lw;
      for (int i : owned) {
        exsqs::ByteWriter b;
        eng[static_cast<size_t>(i)].serialize(b);
        lw.u32(static_cast<uint32_t>(i));
        lw.str(b.data());
      }
      const std::string all = gather_bytes_root(lw.data(), rank);
      if (rank == 0)
        exsqs::save_run_state_blobs(cfg.outdir + "/state.ckpt", cfg,
                                    parse_records(all, cfg.islands));
      MPI_Barrier(MPI_COMM_WORLD);
    };

    // ---- lockstep rounds (SPEC 8.1/8.2) ----
    while (true) {
      std::vector<int> flags(static_cast<size_t>(cfg.islands), 0);
      for (int i : owned)
        flags[static_cast<size_t>(i)] = eng[static_cast<size_t>(i)].done() ? 0 : 1;
      MPI_Allreduce(MPI_IN_PLACE, flags.data(), cfg.islands, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      bool any = false;
      for (int f : flags) any = any || (f != 0);
      if (!any) break;

      std::vector<std::exception_ptr> errs(owned.size());
#ifdef _OPENMP
#pragma omp parallel for num_threads(outer) schedule(static, 1)
#endif
      for (int k = 0; k < static_cast<int>(owned.size()); ++k) {
        try {
          eng[static_cast<size_t>(owned[static_cast<size_t>(k)])].advance();
        } catch (...) {
          errs[static_cast<size_t>(k)] = std::current_exception();
        }
      }
      for (auto& ep : errs)
        if (ep) std::rethrow_exception(ep);
      ++round;

      if (cfg.islands > 1 && cfg.migrants > 0 && cfg.migration_every > 0 &&
          round % cfg.migration_every == 0) {
        std::vector<int> af(static_cast<size_t>(cfg.islands), 0);
        for (int i : owned)
          af[static_cast<size_t>(i)] = eng[static_cast<size_t>(i)].done() ? 0 : 1;
        MPI_Allreduce(MPI_IN_PLACE, af.data(), cfg.islands, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        std::vector<int> act;
        for (int i = 0; i < cfg.islands; ++i)
          if (af[static_cast<size_t>(i)]) act.push_back(i);
        if (act.size() >= 2) {
          exsqs::ByteWriter lw;
          std::vector<std::pair<int, int>> sent;
          for (int i : owned)
            if (af[static_cast<size_t>(i)]) {
              const auto m = eng[static_cast<size_t>(i)].emigrants(cfg.migrants);
              exsqs::ByteWriter b;
              exsqs::put_individuals(b, m);
              lw.u32(static_cast<uint32_t>(i));
              lw.str(b.data());
              sent.emplace_back(i, static_cast<int>(m.size()));
            }
          const std::vector<std::string> blob =
              parse_records(allgather_bytes(lw.data()), cfg.islands);
          for (size_t j = 0; j < act.size(); ++j) {
            const int dst = act[(j + 1) % act.size()];
            if (dst % size == rank) {
              exsqs::ByteReader br(blob[static_cast<size_t>(act[j])]);
              eng[static_cast<size_t>(dst)].receive_migrants(exsqs::get_individuals(br));
            }
          }
          for (const auto& pr : sent)
            eng[static_cast<size_t>(pr.first)].note_emigrants(pr.second);
        }
      }
      if (cfg.checkpoint_every > 0 && round % cfg.checkpoint_every == 0) save_state_root();
    }
    save_state_root();  // final state: resubmission chains [A13]

    // ---- gather island results to rank 0; write pooled outputs [D8] ----
    exsqs::ByteWriter lw;
    for (int i : owned) {
      const exsqs::IslandResult r = eng[static_cast<size_t>(i)].take_result();
      exsqs::ByteWriter b;
      exsqs::put_island_result(b, r);
      lw.u32(static_cast<uint32_t>(i));
      lw.str(b.data());
    }
    const std::string all = gather_bytes_root(lw.data(), rank);
    int code = 3;
    if (rank == 0) {
      const std::vector<std::string> blobs = parse_records(all, cfg.islands);
      std::vector<exsqs::IslandResult> rs;
      rs.reserve(blobs.size());
      for (const auto& b : blobs) {
        exsqs::ByteReader br(b);
        rs.push_back(exsqs::get_island_result(br));
      }
      exsqs::RunOutput out = exsqs::merge_island_results(cfg, std::move(rs));
      out.wall_s = std::chrono::duration<double>(Clock::now() - t0).count();
      exsqs::write_outputs(cfg, ctx, out, false);
      if (!out.outputs.empty())
        std::printf("%s | best E_pure=%.6e (%.2fx floor) D=%d SG=%d | outputs in %s\n",
                    out.success ? "SUCCESS" : "BUDGET EXHAUSTED", out.outputs[0].e_pure,
                    ctx.e_floor > 0 ? out.outputs[0].e_pure / ctx.e_floor : 0.0,
                    out.outputs[0].D, out.outputs[0].sg, cfg.outdir.c_str());
      code = out.success ? 0 : 3;
    }
    MPI_Bcast(&code, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Finalize();
    return code;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "exsqs_mpi[rank %d] error: %s\n", rank, e.what());
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  return 1;
}
