#include "exsqs/lattice.hpp"

#include <cmath>

namespace exsqs {

Prototype make_sc(double a) {
  Prototype p;
  p.name = "sc";
  p.cell = {{{a, 0, 0}, {0, a, 0}, {0, 0, a}}};
  p.sites = {{0.0, 0.0, 0.0}};
  return p;
}

Prototype make_bcc(double a) {
  Prototype p;
  p.name = "bcc";
  p.cell = {{{a, 0, 0}, {0, a, 0}, {0, 0, a}}};
  p.sites = {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}};
  return p;
}

Prototype make_fcc(double a) {
  Prototype p;
  p.name = "fcc";
  p.cell = {{{a, 0, 0}, {0, a, 0}, {0, 0, a}}};
  p.sites = {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.0}, {0.5, 0.0, 0.5}, {0.0, 0.5, 0.5}};
  return p;
}

Prototype make_hcp(double a, double c) {
  if (c <= 0.0) c = a * std::sqrt(8.0 / 3.0);
  Prototype p;
  p.name = "hcp";
  const double s3 = std::sqrt(3.0);
  p.cell = {{{a, 0, 0}, {-0.5 * a, 0.5 * s3 * a, 0}, {0, 0, c}}};
  p.sites = {{1.0 / 3.0, 2.0 / 3.0, 0.25}, {2.0 / 3.0, 1.0 / 3.0, 0.75}};
  return p;
}

}  // namespace exsqs
