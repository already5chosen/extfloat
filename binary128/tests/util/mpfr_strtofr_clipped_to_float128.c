#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <quadmath.h>

#include "mpfr_strtofr_clipped_to_float128.h"
#include "qmx2mpfr.h"

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

void manual_input_parser(mpfr_t rop, _Float128* rop128, const char *nptr)
{
  mpfr_init2(rop, 113);
  char* endp;
  mpfr_strtofr_clipped_to_float128(rop, nptr, &endp);
  mpfr_to_float128(rop128, rop);
  if (isnanq(*rop128) && endp[0] == '.') {
    bool signaling = false;
    if (endp[1] == 's') {
      ++endp;
      signaling = true;
    }
    char tmp[128];
    snprintf(tmp, sizeof(tmp), "NAN(%s)", endp+1);
    *rop128 = strtoflt128(tmp, NULL);
    if (signaling) { // convert xy[i] to signaling NaN
       uint32_t dw[4];
       memcpy(dw, rop128, sizeof(dw));
       #if (__FLOAT_WORD_ORDER__ == __ORDER_LITTLE_ENDIAN__)
       dw[3] &= ~((uint32_t)1 << 15);
       #else
       dw[0] &= ~((uint32_t)1 << 15);
       #endif
       memcpy(rop128, dw, sizeof(*rop128));
    }
  }
}
