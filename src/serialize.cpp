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

}  // namespace exsqs
