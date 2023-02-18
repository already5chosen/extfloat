#include <fenv.h>

static const int rm[4] = { FE_DOWNWARD,  FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD };
void addq_all_rm(__float128 res[4], __float128 values[2], int exceptions[4])
{
  for (int i = 0; i < 4; ++i) {
    fesetround(rm[i]);
    feclearexcept(FE_ALL_EXCEPT);
    res[i] = values[0] + values[1];
    exceptions[i] = fetestexcept(FE_ALL_EXCEPT);
  }
  fesetround(FE_TONEAREST); // restore default
}


void subq_all_rm(__float128 res[4], __float128 values[2], int exceptions[4])
{
  for (int i = 0; i < 4; ++i) {
    fesetround(rm[i]);
    feclearexcept(FE_ALL_EXCEPT);
    res[i] = values[0] - values[1];
    exceptions[i] = fetestexcept(FE_ALL_EXCEPT);
  }
  fesetround(FE_TONEAREST); // restore default
}

