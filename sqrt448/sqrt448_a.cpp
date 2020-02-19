#include <cstdint>
#include <cstring>
#include "sqrt448.h"


extern "C" void asm_sqrt448(uint64_t* __restrict dst, const uint64_t* __restrict src, int exp1);
cpp_bin_float_132 my_sqrt(const cpp_bin_float_132& x) {
  if (x <= 0) {
    return 0; // do whatever you feel correct for negatives
  }
  int e = x.backend().exponent();
  cpp_bin_float_132 res(x);
  res.backend().exponent() = int(e & -2) / 2;
  asm_sqrt448((uint64_t*)&res.backend().bits(), (uint64_t*)&x.backend().bits(), e & 1);

  return res;
}
