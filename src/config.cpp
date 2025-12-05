#include "exsqs/config.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include <yaml-cpp/yaml.h>

#include "exsqs/lattice.hpp"

namespace exsqs {
namespace {

[[noreturn]] void fail(const std::string& m) { throw std::runtime_error("config: " + m); }

void apply_override(YAML::Node root, const std::string& kv) {
  const size_t eq = kv.find('=');
  if (eq == std::string::npos) fail("--set expects key=value, got '" + kv + "'");
  const std::string path = kv.substr(0, eq), val = kv.substr(eq + 1);
  YAML::Node cur = root;
  size_t pos = 0;
  while (true) {
    const size_t dot = path.find('.', pos);
    const std::string key = path.substr(pos, dot == std::string::npos ? dot : dot - pos);
    if (dot == std::string::npos) {
      cur[key] = val;
      break;
    }
    YAML::Node next = cur[key];
    if (!next) {
      cur[key] = YAML::Node(YAML::NodeType::Map);
      next = cur[key];
    }
    cur.reset(next);
    pos = dot + 1;
  }
}

}  // namespace

RunConfig load_config(const std::string& path, const std::vector<std::string>& overrides) {
  YAML::Node root = YAML::LoadFile(path);
  for (const std::string& kv : overrides) apply_override(root, kv);
  RunConfig c;

  const auto lat = root["lattice"];
  if (!lat) fail("missing lattice");
  const std::string type = lat["type"] ? lat["type"].as<std::string>() : "";
  const double a = lat["a"] ? lat["a"].as<double>() : -1.0;
  if (a <= 0) fail("lattice.a must be positive");
  if (type == "sc")
    c.proto = make_sc(a);
  else if (type == "bcc")
    c.proto = make_bcc(a);
  else if (type == "fcc")
    c.proto = make_fcc(a);
  else
    fail("lattice.type must be sc|bcc|fcc");

  const auto sc = root["supercell"];
  if (!sc || !sc["diag"]) fail("supercell.diag required");
  const auto d = sc["diag"].as<std::vector<int>>();
  if (d.size() != 3 || d[0] <= 0 || d[1] <= 0 || d[2] <= 0)
    fail("supercell.diag needs three positive integers");
  c.H = {{{d[0], 0, 0}, {0, d[1], 0}, {0, 0, d[2]}}};

  const auto comp = root["composition"];
  if (!comp || !comp.IsMap()) fail("missing composition map");
  double xsum = 0.0;
  for (const auto& kv : comp) {
    c.species.push_back(kv.first.as<std::string>());
    c.x_target.push_back(kv.second.as<double>());
    xsum += c.x_target.back();
  }
  if (std::abs(xsum - 1.0) > 1e-6) fail("composition fractions must sum to 1");

  const int cells = c.H[0][0] * c.H[1][1] * c.H[2][2];
  const int N = cells * static_cast<int>(c.proto.sites.size());
  std::vector<double> ideal;
  int assigned = 0;
  for (double x : c.x_target) {
    ideal.push_back(x * N);
    c.counts.push_back(static_cast<int>(std::floor(x * N)));
    assigned += c.counts.back();
  }
  std::vector<size_t> order(c.counts.size());
  for (size_t i = 0; i < order.size(); ++i) order[i] = i;
  std::sort(order.begin(), order.end(), [&](size_t p, size_t q) {
    return ideal[p] - std::floor(ideal[p]) > ideal[q] - std::floor(ideal[q]);
  });
  for (int r = 0; r < N - assigned; ++r) c.counts[order[static_cast<size_t>(r) % order.size()]]++;
  for (int k : c.counts) c.x_achieved.push_back(static_cast<double>(k) / N);

  if (root["output"] && root["output"]["dir"]) c.outdir = root["output"]["dir"].as<std::string>();
  return c;
}

}  // namespace exsqs
