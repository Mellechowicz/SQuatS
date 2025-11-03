#include "exsqs/structure.hpp"

#include <algorithm>
#include <climits>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>

namespace exsqs {

std::vector<int> Structure::counts() const {
  std::vector<int> c(names.size(), 0);
  for (int t : species) c.at(static_cast<size_t>(t))++;
  return c;
}

std::vector<double> Structure::concentrations() const {
  const auto c = counts();
  std::vector<double> x(c.size());
  const double n = static_cast<double>(natoms());
  for (size_t i = 0; i < c.size(); ++i) x[i] = c[i] / n;
  return x;
}

Structure make_supercell(const Prototype& p, const Mat3i& H,
                         const std::vector<std::string>& names, int species_fill) {
  const long long dH = det(H);
  if (dH <= 0) throw std::runtime_error("make_supercell: det(H) must be positive");

  Structure s;
  s.names = names;
  for (int i = 0; i < 3; ++i)
    for (int k = 0; k < 3; ++k) {
      double v = 0.0;
      for (int j = 0; j < 3; ++j) v += H[i][j] * p.cell[j][k];
      s.cell[i][k] = v;
    }

  Mat3 Hd;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) Hd[i][j] = static_cast<double>(H[i][j]);
  const Mat3 Hinv = inverse(Hd);

  // Bounding box of the supercell corners in old-lattice integer coordinates.
  int lo[3] = {INT_MAX, INT_MAX, INT_MAX}, hi[3] = {INT_MIN, INT_MIN, INT_MIN};
  for (int u0 = 0; u0 <= 1; ++u0)
    for (int u1 = 0; u1 <= 1; ++u1)
      for (int u2 = 0; u2 <= 1; ++u2)
        for (int j = 0; j < 3; ++j) {
          const int x = u0 * H[0][j] + u1 * H[1][j] + u2 * H[2][j];
          lo[j] = std::min(lo[j], x);
          hi[j] = std::max(hi[j], x);
        }
  for (int j = 0; j < 3; ++j) {
    lo[j] -= 1;
    hi[j] += 1;
  }

  const double eps = 1e-9;
  for (int n0 = lo[0]; n0 <= hi[0]; ++n0)
    for (int n1 = lo[1]; n1 <= hi[1]; ++n1)
      for (int n2 = lo[2]; n2 <= hi[2]; ++n2)
        for (const Vec3& f0 : p.sites) {
          const Vec3 f = {n0 + f0[0], n1 + f0[1], n2 + f0[2]};
          Vec3 fr{0, 0, 0};  // row vector times matrix: fr = f * Hinv
          for (int j = 0; j < 3; ++j)
            fr[j] = f[0] * Hinv[0][j] + f[1] * Hinv[1][j] + f[2] * Hinv[2][j];
          bool inside = true;
          for (int j = 0; j < 3; ++j)
            if (!(fr[j] >= -eps && fr[j] < 1.0 - eps)) {
              inside = false;
              break;
            }
          if (!inside) continue;
          for (int j = 0; j < 3; ++j)
            if (fr[j] < 0.0) fr[j] = 0.0;
          s.frac.push_back(fr);
        }

  const long long expected = dH * static_cast<long long>(p.sites.size());
  if (static_cast<long long>(s.frac.size()) != expected)
    throw std::runtime_error("make_supercell: site enumeration mismatch");
  s.species.assign(s.frac.size(), species_fill);
  return s;
}

Structure make_supercell_diag(const Prototype& p, int n1, int n2, int n3,
                              const std::vector<std::string>& names, int species_fill) {
  const Mat3i H = {{{n1, 0, 0}, {0, n2, 0}, {0, 0, n3}}};
  return make_supercell(p, H, names, species_fill);
}

Structure decorate(const Structure& geom, const std::vector<int>& sigma,
                   const std::vector<std::string>& names) {
  if (sigma.size() != geom.frac.size()) throw std::runtime_error("decorate: sigma size mismatch");
  for (int t : sigma)
    if (t < 0 || t >= static_cast<int>(names.size()))
      throw std::runtime_error("decorate: species index out of range");
  Structure s = geom;
  s.species = sigma;
  s.names = names;
  return s;
}

Structure grouped_by_species(const Structure& s) {
  Structure o = s;
  o.frac.clear();
  o.species.clear();
  for (int t = 0; t < s.nspecies(); ++t)
    for (int i = 0; i < s.natoms(); ++i)
      if (s.species[i] == t) {
        o.frac.push_back(s.frac[i]);
        o.species.push_back(t);
      }
  return o;
}

void write_poscar(const Structure& s, const std::string& path, const std::string& comment) {
  std::ofstream f(path);
  if (!f) throw std::runtime_error("write_poscar: cannot open " + path);
  f << (comment.empty() ? "exsqs" : comment) << "\n";
  f << "1.0\n";
  f << std::setprecision(16) << std::fixed;
  for (int i = 0; i < 3; ++i)
    f << "  " << s.cell[i][0] << " " << s.cell[i][1] << " " << s.cell[i][2] << "\n";
  const auto cnt = s.counts();
  std::string line1, line2;
  for (size_t t = 0; t < s.names.size(); ++t) {
    if (cnt[t] == 0) continue;
    if (!line1.empty()) {
      line1 += " ";
      line2 += " ";
    }
    line1 += s.names[t];
    line2 += std::to_string(cnt[t]);
  }
  f << line1 << "\n" << line2 << "\n" << "Direct\n";
  for (size_t t = 0; t < s.names.size(); ++t) {
    if (cnt[t] == 0) continue;
    for (int i = 0; i < s.natoms(); ++i)
      if (s.species[i] == static_cast<int>(t))
        f << "  " << s.frac[i][0] << " " << s.frac[i][1] << " " << s.frac[i][2] << "\n";
  }
}


Structure read_poscar(const std::string& path) {
  std::ifstream f(path);
  if (!f) throw std::runtime_error("read_poscar: cannot open " + path);
  auto next_line = [&]() {
    std::string line;
    if (!std::getline(f, line)) throw std::runtime_error("read_poscar: truncated " + path);
    return line;
  };
  next_line();  // comment
  const double scale = std::stod(next_line());
  if (scale <= 0.0) throw std::runtime_error("read_poscar: only positive scale supported");
  Structure s;
  for (int i = 0; i < 3; ++i) {
    std::istringstream ls(next_line());
    for (int k = 0; k < 3; ++k) {
      ls >> s.cell[i][k];
      s.cell[i][k] *= scale;
    }
    if (!ls) throw std::runtime_error("read_poscar: bad lattice row");
  }
  std::string species_line = next_line();
  {
    std::istringstream probe(species_line);
    double x;
    if (probe >> x) throw std::runtime_error("read_poscar: VASP4 (no species line) unsupported");
  }
  std::istringstream sl(species_line);
  for (std::string w; sl >> w;) s.names.push_back(w);
  std::vector<int> cnt;
  {
    std::istringstream cl(next_line());
    for (int c; cl >> c;) cnt.push_back(c);
  }
  if (cnt.size() != s.names.size() || cnt.empty())
    throw std::runtime_error("read_poscar: species/count mismatch");
  std::string mode = next_line();
  if (!mode.empty() && (mode[0] == 'S' || mode[0] == 's')) mode = next_line();  // selective dyn.
  const bool cart = !mode.empty() && (mode[0] == 'C' || mode[0] == 'c' || mode[0] == 'K' || mode[0] == 'k');
  const Mat3 inv = inverse(s.cell);
  int total = 0;
  for (int c : cnt) total += c;
  for (size_t t = 0; t < cnt.size(); ++t)
    for (int i = 0; i < cnt[t]; ++i) {
      std::istringstream ls(next_line());
      Vec3 v{};
      ls >> v[0] >> v[1] >> v[2];
      if (!ls) throw std::runtime_error("read_poscar: bad coordinate line");
      if (cart) {
        const Vec3 r = {scale * v[0], scale * v[1], scale * v[2]};
        Vec3 fr{};
        for (int j = 0; j < 3; ++j) fr[j] = r[0] * inv[0][j] + r[1] * inv[1][j] + r[2] * inv[2][j];
        v = fr;
      }
      for (int j = 0; j < 3; ++j) {
        v[j] -= std::floor(v[j]);
        if (v[j] > 1.0 - 1e-9) v[j] = 0.0;
      }
      s.frac.push_back(v);
      s.species.push_back(static_cast<int>(t));
    }
  (void)total;
  return s;
}

}  // namespace exsqs

