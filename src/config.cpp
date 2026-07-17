#include "exsqs/config.hpp"

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <numeric>
#include <sstream>
#include <stdexcept>

#include "exsqs/structure.hpp"

namespace exsqs {
namespace {

[[noreturn]] void fail(const std::string& msg) { throw std::runtime_error("config: " + msg); }

std::vector<std::string> split(const std::string& s, char sep) {
  std::vector<std::string> out;
  std::string cur;
  for (char c : s) {
    if (c == sep) {
      out.push_back(cur);
      cur.clear();
    } else {
      cur += c;
    }
  }
  out.push_back(cur);
  return out;
}

void apply_override(YAML::Node root, const std::string& spec) {
  const auto eq = spec.find('=');
  if (eq == std::string::npos) fail("override must be key.path=value: " + spec);
  const auto keys = split(spec.substr(0, eq), '.');
  const std::string value = spec.substr(eq + 1);
  YAML::Node cur;
  cur.reset(root);  // reset() rebinds the wrapper; operator= would ASSIGN INTO root
  for (size_t i = 0; i + 1 < keys.size(); ++i) {
    YAML::Node next = cur[keys[i]];
    cur.reset(next);
  }
  cur[keys.back()] = YAML::Load(value);  // parse scalar/typed value
}

// largest-remainder rounding [A5]
std::vector<int> largest_remainder(const std::vector<double>& x, int N) {
  const int K = static_cast<int>(x.size());
  std::vector<int> cnt(K);
  std::vector<std::pair<double, int>> rem(K);
  int assigned = 0;
  for (int t = 0; t < K; ++t) {
    const double v = x[t] * N;
    cnt[t] = static_cast<int>(std::floor(v));
    rem[t] = {v - cnt[t], t};
    assigned += cnt[t];
  }
  std::sort(rem.begin(), rem.end(), [](const auto& a, const auto& b) {
    if (a.first != b.first) return a.first > b.first;
    return a.second < b.second;  // deterministic tie break
  });
  for (int k = 0; k < N - assigned; ++k) cnt[rem[static_cast<size_t>(k)].second]++;
  return cnt;
}

}  // namespace

RunConfig load_config(const std::string& path, const std::vector<std::string>& overrides) {
  RunConfig c;
  {
    std::ifstream f(path);
    if (!f) fail("cannot open " + path);
    std::stringstream ss;
    ss << f.rdbuf();
    c.config_echo = ss.str();
  }
  YAML::Node root;
  try {
    root = YAML::Load(c.config_echo);
    for (const auto& ov : overrides) apply_override(root, ov);
  } catch (const YAML::Exception& e) {
    fail(std::string("YAML error: ") + e.what());
  }

  // ---- lattice ----
  const auto lat = root["lattice"];
  if (!lat) fail("missing lattice section");
  const std::string type = lat["type"] ? lat["type"].as<std::string>() : "";
  const bool has_file = lat["file"] && !lat["file"].IsNull();
  if (has_file) {
    const Structure ps = read_poscar(lat["file"].as<std::string>());
    c.proto.name = "file";
    c.proto.cell = ps.cell;
    c.proto.sites = ps.frac;  // prototype = geometry; decorations replace species
  } else {
    const double a = lat["a"] ? lat["a"].as<double>() : -1.0;
    if (a <= 0) fail("lattice.a must be positive");
    if (type == "sc")
      c.proto = make_sc(a);
    else if (type == "bcc")
      c.proto = make_bcc(a);
    else if (type == "fcc")
      c.proto = make_fcc(a);
    else if (type == "hcp")
      c.proto = make_hcp(a, lat["c"] ? lat["c"].as<double>() : -1.0);
    else
      fail("lattice.type must be sc|bcc|fcc|hcp or lattice.file given");
  }

  // ---- supercell ----
  const auto sc = root["supercell"];
  if (!sc) fail("missing supercell section");
  if (sc["diag"]) {
    const auto d = sc["diag"].as<std::vector<int>>();
    if (d.size() != 3) fail("supercell.diag needs 3 integers");
    c.H = {{{d[0], 0, 0}, {0, d[1], 0}, {0, 0, d[2]}}};
  } else if (sc["matrix"]) {
    const auto m = sc["matrix"].as<std::vector<std::vector<int>>>();
    if (m.size() != 3 || m[0].size() != 3) fail("supercell.matrix must be 3x3");
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) c.H[i][j] = m[static_cast<size_t>(i)][static_cast<size_t>(j)];
  } else {
    fail("supercell needs diag or matrix");
  }
  if (det(c.H) <= 0) fail("supercell matrix must have positive determinant");

  // ---- composition (document order preserved by yaml-cpp) ----
  const auto comp = root["composition"];
  if (!comp || !comp.IsMap() || comp.size() < 2) fail("composition needs >= 2 species");
  double xsum = 0;
  for (auto it = comp.begin(); it != comp.end(); ++it) {
    c.species.push_back(it->first.as<std::string>());
    c.x_target.push_back(it->second.as<double>());
    xsum += c.x_target.back();
  }
  if (std::abs(xsum - 1.0) > 1e-6) fail("composition must sum to 1");
  const long long Nll = det(c.H) * static_cast<long long>(c.proto.sites.size());
  const int N = static_cast<int>(Nll);
  c.counts = largest_remainder(c.x_target, N);
  c.x_achieved.resize(c.species.size());
  const double comp_tol = 0.02;  // section 9 default
  for (size_t t = 0; t < c.species.size(); ++t) {
    if (c.counts[t] <= 0) fail("species " + c.species[t] + " rounds to zero atoms");
    c.x_achieved[t] = static_cast<double>(c.counts[t]) / N;
    if (std::abs(c.x_achieved[t] - c.x_target[t]) > comp_tol)
      fail("achieved composition deviates > comp_tol for " + c.species[t]);
  }

  // ---- zones ----
  if (const auto z = root["zones"]) {
    if (z["n_shells"]) c.n_shells = z["n_shells"].as<int>();
    if (z["shell_tol"]) c.shell_tol = z["shell_tol"].as<double>();
  }

  // ---- error ----
  if (const auto e = root["error"]) {
    if (const auto w = e["weights"]) {
      const std::string form = w["form"] ? w["form"].as<std::string>() : "inv_n";
      if (form == "inv_n")
        c.wform = WeightForm::InvN;
      else if (form == "inv_n_pow") {
        c.wform = WeightForm::InvNPow;
        if (!w["p"]) fail("weights.p required for inv_n_pow");
        c.wpow = w["p"].as<double>();
      } else if (form == "inv_R")
        c.wform = WeightForm::InvR;
      else if (form == "custom") {
        c.wform = WeightForm::Custom;
        if (!w["values"]) fail("weights.values required for custom");
        c.wcustom = w["values"].as<std::vector<double>>();
      } else
        fail("unknown weights.form " + form);
    }
    const std::string mode = e["mode"] ? e["mode"].as<std::string>() : "auto";
    if (mode == "auto")
      c.full_pairs = c.species.size() >= 3;  // [A16]
    else if (mode == "diagonal")
      c.full_pairs = false;
    else if (mode == "full_pairs")
      c.full_pairs = true;
    else
      fail("error.mode must be auto|diagonal|full_pairs");
    if (e["gamma"]) c.gamma = e["gamma"].as<double>();
    // v1.9 multiplet sectors (SPEC 4.2); keys parsed and validated even when
    // the lambdas stay 0, so misuse fails loudly (project config policy)
    if (const auto mu = e["multiplets"]) {
      if (mu["lambda3"]) c.lambda3 = mu["lambda3"].as<double>();
      if (mu["lambda4"]) c.lambda4 = mu["lambda4"].as<double>();
      if (c.lambda3 < 0 || c.lambda4 < 0) fail("multiplets.lambda3/4 must be >= 0");
      if (mu["shell3"]) c.mshell3 = mu["shell3"].as<int>();
      if (mu["shell4"]) c.mshell4 = mu["shell4"].as<int>();
      if (c.mshell3 < 1 || c.mshell3 > c.n_shells)
        fail("multiplets.shell3 must be in [1, zones.n_shells]");
      if (c.mshell4 < 1 || c.mshell4 > c.n_shells)
        fail("multiplets.shell4 must be in [1, zones.n_shells]");
    }
  } else {
    c.full_pairs = c.species.size() >= 3;
  }

  // ---- evolution ----
  if (const auto ev = root["evolution"]) {
    if (ev["population"]) c.population = ev["population"].as<int>();
    if (ev["outputs"]) c.outputs = ev["outputs"].as<int>();
    if (ev["e_tol"]) {
      const std::string et = ev["e_tol"].as<std::string>();
      if (et == "auto") {
        c.e_tol = -1.0;  // resolved to 3.0 * E_floor at run start (v1.1, SPEC 4.1)
      } else {
        c.e_tol = ev["e_tol"].as<double>();
        if (c.e_tol <= 0) fail("e_tol must be > 0 or 'auto'");
      }
    }
    if (ev["max_generations"]) c.max_generations = ev["max_generations"].as<int>();
    if (ev["max_wall_s"]) c.max_wall_s = ev["max_wall_s"].as<double>();
    if (ev["stagnation_stop"]) {
      c.stagnation_stop = ev["stagnation_stop"].as<int>();
      if (c.stagnation_stop < 0) fail("stagnation_stop must be >= 0");
    }
    if (const auto sv = ev["survival"]) {
      const std::string m = sv["mode"] ? sv["mode"].as<std::string>() : "ratio";
      if (m == "ratio")
        c.metropolis = false;
      else if (m == "metropolis") {
        c.metropolis = true;
        if (sv["beta"] && !sv["beta"].IsNull()) {
          const std::string b = sv["beta"].as<std::string>();
          c.beta = (b == "auto") ? -1.0 : sv["beta"].as<double>();
        }
      } else
        fail("survival.mode must be ratio|metropolis");
      // [A11] v1.5: schedule keys apply to metropolis only, but are parsed and
      // validated for every mode so misuse fails loudly (project config policy)
      if (sv["schedule"]) {
        const std::string sch = sv["schedule"].as<std::string>();
        if (sch == "const")
          c.beta_schedule = 0;
        else if (sch == "geometric")
          c.beta_schedule = 1;
        else
          fail("survival.schedule must be const|geometric");
      }
      if (sv["beta_growth"]) c.beta_growth = sv["beta_growth"].as<double>();
      if (c.beta_schedule == 1) {
        if (!c.metropolis) fail("survival.schedule geometric requires mode metropolis");
        if (c.beta <= 0) fail("survival.schedule geometric requires a numeric beta");
        if (c.beta_growth <= 0) fail("survival.beta_growth must be > 0");
      }
    }
    if (ev["elitism_best"]) c.elitism_best = ev["elitism_best"].as<int>();
    if (ev["p1_elite_quota"]) c.p1_elite_quota = ev["p1_elite_quota"].as<int>();
    if (const auto mu = ev["mutation"]) {
      if (mu["swaps"]) {
        const std::string sw = mu["swaps"].as<std::string>();
        if (sw == "poisson") {  // k ~ 1 + Poisson(lambda) [A12], v1.1
          c.mut_poisson_lambda = mu["lambda"] ? mu["lambda"].as<double>() : 1.0;
          if (c.mut_poisson_lambda < 0) fail("mutation.lambda must be >= 0");
        } else {
          c.mut_swaps = mu["swaps"].as<int>();
          if (c.mut_swaps < 1) fail("mutation.swaps must be >= 1 or 'poisson'");
        }
      }
      if (mu["symmetry_preserving"]) c.mut_sympres = mu["symmetry_preserving"].as<bool>();
      if (mu["retry_budget"]) c.retry_budget = mu["retry_budget"].as<int>();
    }
    if (const auto sd = ev["seeding"]) {
      const std::string m = sd["mode"] ? sd["mode"].as<std::string>() : "mixed";
      if (m == "rejection")
        c.seed_mode = 0;
      else if (m == "constructive")
        c.seed_mode = 1;
      else if (m == "mixed")
        c.seed_mode = 2;
      else
        fail("seeding.mode must be rejection|constructive|mixed");
      if (sd["sg_pool"] && sd["sg_pool"].as<std::string>() != "auto")
        fail("seeding.sg_pool: only 'auto' is defined (the spec names no alternatives)");
    }
  }
  if (c.population < 2) fail("population must be >= 2");
  if (c.outputs < 1) fail("outputs must be >= 1");
  if (c.elitism_best < 1 || c.elitism_best >= c.population) fail("elitism_best out of range");

  // ---- symmetry ----
  if (const auto sy = root["symmetry"]) {
    if (sy["symprec"]) c.symprec = sy["symprec"].as<double>();
    if (sy["filter"] && sy["filter"]["policy"] &&
        sy["filter"]["policy"].as<std::string>() != "reject_p1")
      fail("filter.policy: only reject_p1 is defined [A7]");
  }

  // ---- displacements ----
  if (root["displacements"] && root["displacements"]["convention"] &&
      root["displacements"]["convention"].as<std::string>() != "phonopy_default")
    fail("displacements.convention must be phonopy_default [A15]");

  // ---- parallel ----
  if (const auto pa = root["parallel"]) {
    if (pa["islands"]) c.islands = pa["islands"].as<int>();
    if (c.islands < 1) fail("islands must be >= 1");
    if (pa["omp_threads"] && !pa["omp_threads"].IsNull()) {
      const std::string ot = pa["omp_threads"].as<std::string>();
      if (ot == "auto") {
        c.omp_threads = -1;
      } else {
        c.omp_threads = pa["omp_threads"].as<int>();
        if (c.omp_threads < 1) fail("omp_threads must be >= 1 or 'auto'");
      }
    }
    if (pa["migration_every"]) {
      c.migration_every = pa["migration_every"].as<int>();
      if (c.migration_every < 0) fail("migration_every must be >= 0");
    }
    if (pa["migrants"]) {
      c.migrants = pa["migrants"].as<int>();
      if (c.migrants < 0) fail("migrants must be >= 0");
    }
  }

  // ---- rng / output ----
  if (root["rng"] && root["rng"]["seed"]) c.seed = root["rng"]["seed"].as<uint64_t>();
  if (const auto o = root["output"]) {
    if (o["dir"]) c.outdir = o["dir"].as<std::string>();
    if (o["checkpoint_every"]) c.checkpoint_every = o["checkpoint_every"].as<int>();
    if (o["log_level"]) c.log_info = (o["log_level"].as<std::string>() == "info" ||
                                      o["log_level"].as<std::string>() == "debug");
    if (o["formats"])
      for (const auto& fm : o["formats"])
        if (fm.as<std::string>() != "poscar")
          std::fprintf(stderr, "config: output format '%s' ignored (poscar is the only writer)\n",
                       fm.as<std::string>().c_str());
  }

  if (c.species.size() >= 4)
    std::fprintf(stderr,
                 "config: K >= 4 detected -- exercised but not campaign-validated "
                 "(K = 3 validated per SPEC section 12; cross-check runs with "
                 "tools/py/validate.py, T-V1)\n");
  return c;
}

}  // namespace exsqs
