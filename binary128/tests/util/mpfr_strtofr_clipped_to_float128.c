#include "mpfr_strtofr_clipped_to_float128.h"

int mpfr_strtofr_clipped_to_float128(mpfr_t rop, const char *nptr, char **endptr)
{
  int ret = mpfr_strtofr(rop, nptr, endptr, 0, GMP_RNDN);
  if (mpfr_regular_p(rop)) {
    MPFR_DECL_INIT(flt128_nrm_min, 113);
    mpfr_set_ui_2exp(flt128_nrm_min, 1, -16382, GMP_RNDN);
    if (mpfr_cmpabs(rop, flt128_nrm_min) < 0) { // subnormal or zero
      int sign = mpfr_signbit(rop);
      mpfr_setsign(rop, rop, 0, GMP_RNDN);
      mpfr_add(rop, rop, flt128_nrm_min, GMP_RNDN);
      mpfr_sub(rop, rop, flt128_nrm_min, GMP_RNDN);
      mpfr_setsign(rop, rop, sign, GMP_RNDN);
    } else {
      MPFR_DECL_INIT(flt128_max, 113);
      mpfr_set_ui_2exp(flt128_max, 1, 16384, GMP_RNDN);
      if (mpfr_cmpabs(rop, flt128_max) >= 0) { // very big number
        int sign = mpfr_signbit(rop);
        mpfr_set(rop, flt128_max, GMP_RNDN);
        mpfr_nextbelow(rop);
        mpfr_setsign(rop, rop, sign, GMP_RNDN);
      }
    }
  }
  return ret;
}
