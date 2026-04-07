#include "exsqs/serialize.hpp"

#include <cstring>
#include <stdexcept>

namespace exsqs {

void ensure_little_endian() {
  const uint32_t probe = 0x01020304u;
  uint8_t b[4];
  std::memcpy(b, &probe, 4);
  if (b[0] != 0x04)
    throw std::runtime_error("serialize: big-endian host unsupported for state/wire format v1");
}

void ByteReader::raw_read(void* dst, size_t n) {
  const char* p = take(n);
  std::memcpy(dst, p, n);
}

const char* ByteReader::take(size_t n) {
  if (static_cast<size_t>(end_ - p_) < n)
    throw std::runtime_error("serialize: truncated or corrupt payload");
  const char* p = p_;
  p_ += n;
  return p;
}

}  // namespace exsqs
