#include <fenv.h>

static const int rm[4] = { FE_DOWNWARD,  FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD };
void mulq_all_rm(__float128 res[4], __float128 values[2])
{
  for (int i = 0; i < 4; ++i) {
    fesetround(rm[i]);
    res[i] = values[0] * values[1];
  }
  fesetround(FE_TONEAREST); // restore default
}
