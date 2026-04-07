#pragma once
// v1.3 binary encoding for checkpoint/restart state and MPI wire transfer.
// Fixed-width little-endian primitives; doubles travel as IEEE-754 bit
// patterns, so round trips are bit-exact. Assumes a homogeneous little-endian
// cluster (the practical HPC case) -- ensure_little_endian() guards at every
// file/wire boundary. The state file carries a format version for evolution.

#include <cstdint>
#include <string>
#include <vector>

#include "exsqs/config.hpp"
#include "exsqs/evolution.hpp"

namespace exsqs {

void ensure_little_endian();  // throws on big-endian hosts

class ByteWriter {
 public:
  void u8(uint8_t v) { buf_.push_back(static_cast<char>(v)); }
  void u32(uint32_t v) { raw(&v, sizeof v); }
  void u64(uint64_t v) { raw(&v, sizeof v); }
  void i32(int32_t v) { raw(&v, sizeof v); }
  void i64(int64_t v) { raw(&v, sizeof v); }
  void f64(double v) { raw(&v, sizeof v); }
  void str(const std::string& s) {
    u32(static_cast<uint32_t>(s.size()));
    raw(s.data(), s.size());
  }
  void ints(const std::vector<int>& v) {
    u32(static_cast<uint32_t>(v.size()));
    for (int x : v) i32(x);
  }
  void bytes(const std::vector<uint8_t>& v) {
    u32(static_cast<uint32_t>(v.size()));
    if (!v.empty()) raw(v.data(), v.size());
  }
  void raw(const void* p, size_t n) { buf_.append(static_cast<const char*>(p), n); }
  const std::string& data() const { return buf_; }

 private:
  std::string buf_;
};

class ByteReader {
 public:
  ByteReader(const char* p, size_t n) : p_(p), end_(p + n) {}
  explicit ByteReader(const std::string& s) : ByteReader(s.data(), s.size()) {}

  uint8_t u8() { return static_cast<uint8_t>(*take(1)); }
  uint32_t u32() { return get<uint32_t>(); }
  uint64_t u64() { return get<uint64_t>(); }
  int32_t i32() { return get<int32_t>(); }
  int64_t i64() { return get<int64_t>(); }
  double f64() { return get<double>(); }
  std::string str() {
    const uint32_t n = u32();
    const char* p = take(n);
    return std::string(p, n);
  }
  std::vector<int> ints() {
    const uint32_t n = u32();
    std::vector<int> v(n);
    for (uint32_t i = 0; i < n; ++i) v[i] = i32();
    return v;
  }
  std::vector<uint8_t> bytes() {
    const uint32_t n = u32();
    const char* p = take(n);
    return std::vector<uint8_t>(p, p + n);
  }
  void raw_read(void* dst, size_t n);
  size_t remaining() const { return static_cast<size_t>(end_ - p_); }

 private:
  const char* take(size_t n);  // bounds-checked; throws on truncation
  template <typename T>
  T get() {
    T v;
    const char* p = take(sizeof(T));
    __builtin_memcpy(&v, p, sizeof(T));
    return v;
  }
  const char* p_;
  const char* end_;
};

}  // namespace exsqs
