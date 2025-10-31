#pragma once
// Minimal 3D linear algebra for exsqs.
//
// Conventions (match phonopy / pymatgen):
//   * Mat3 rows are the basis vectors a, b, c  =>  cartesian = fractional * L.
//   * Integer rotation matrices act on fractional/lattice vectors as v' = R v,
//     which equals numpy's np.dot(v, R.T) used throughout phonopy.
//
// NOTE (deviation from the initial plan): Eigen3 is hosted on gitlab.com, which
// is unreachable from the sandboxed build network; a 3x3 in-house
// implementation is all the core needs and removes an external dependency.

#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>

namespace exsqs {

using Vec3 = std::array<double, 3>;
using Mat3 = std::array<Vec3, 3>;                 // rows = basis vectors
using Mat3i = std::array<std::array<int, 3>, 3>;  // integer rotation (fractional)

inline Vec3 operator+(const Vec3& a, const Vec3& b) {
  return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}
inline Vec3 operator-(const Vec3& a, const Vec3& b) {
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}
inline Vec3 operator*(double s, const Vec3& a) { return {s * a[0], s * a[1], s * a[2]}; }

inline double dot(const Vec3& a, const Vec3& b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
inline Vec3 cross(const Vec3& a, const Vec3& b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
inline double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }

// cartesian = fractional * L  (L rows are basis vectors)
inline Vec3 frac_to_cart(const Vec3& f, const Mat3& L) {
  Vec3 c{0.0, 0.0, 0.0};
  for (int k = 0; k < 3; ++k) c[k] = f[0] * L[0][k] + f[1] * L[1][k] + f[2] * L[2][k];
  return c;
}

inline double det(const Mat3& m) {
  return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
         m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
         m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

inline long long det(const Mat3i& m) {
  auto M = [&](int i, int j) { return static_cast<long long>(m[i][j]); };
  return M(0, 0) * (M(1, 1) * M(2, 2) - M(1, 2) * M(2, 1)) -
         M(0, 1) * (M(1, 0) * M(2, 2) - M(1, 2) * M(2, 0)) +
         M(0, 2) * (M(1, 0) * M(2, 1) - M(1, 1) * M(2, 0));
}

// determinant of three integer row vectors (phonopy's determinant([d, r1, r2]))
inline long long det_rows(const std::array<int, 3>& a, const std::array<int, 3>& b,
                          const std::array<int, 3>& c) {
  auto A = [&](int i) { return static_cast<long long>(a[i]); };
  auto B = [&](int i) { return static_cast<long long>(b[i]); };
  auto C = [&](int i) { return static_cast<long long>(c[i]); };
  return A(0) * (B(1) * C(2) - B(2) * C(1)) - A(1) * (B(0) * C(2) - B(2) * C(0)) +
         A(2) * (B(0) * C(1) - B(1) * C(0));
}

inline Mat3 inverse(const Mat3& m) {
  const double d = det(m);
  if (std::abs(d) < 1e-14) throw std::runtime_error("linalg: singular 3x3 matrix");
  Mat3 inv;
  inv[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) / d;
  inv[0][1] = -(m[0][1] * m[2][2] - m[0][2] * m[2][1]) / d;
  inv[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / d;
  inv[1][0] = -(m[1][0] * m[2][2] - m[1][2] * m[2][0]) / d;
  inv[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / d;
  inv[1][2] = -(m[0][0] * m[1][2] - m[0][2] * m[1][0]) / d;
  inv[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) / d;
  inv[2][1] = -(m[0][0] * m[2][1] - m[0][1] * m[2][0]) / d;
  inv[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / d;
  return inv;
}

// v' = R v (equals numpy's np.dot(v, R.T))
inline Vec3 apply_rot(const Mat3i& R, const Vec3& v) {
  Vec3 o{0.0, 0.0, 0.0};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) o[i] += R[i][j] * v[j];
  return o;
}
inline std::array<int, 3> apply_rot(const Mat3i& R, const std::array<int, 3>& v) {
  std::array<int, 3> o{0, 0, 0};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) o[i] += R[i][j] * v[j];
  return o;
}

}  // namespace exsqs
