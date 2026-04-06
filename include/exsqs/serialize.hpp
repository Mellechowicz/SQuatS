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

}  // namespace exsqs
