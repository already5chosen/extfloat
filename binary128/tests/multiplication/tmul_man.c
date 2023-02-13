#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <fenv.h>
#include <mpfr.h>
#include <quadmath.h>

#include "qmx2mpfr.h"
#include "mpfr_strtofr_clipped_to_float128.h"

int main(int argz, char** argv)
{
  mpfr_t xa[2];
  mpfr_init2(xa[0], 113);
  mpfr_init2(xa[1], 113);

  for (int i = 0; i < 2 && i < argz-1; ++i)
    mpfr_strtofr_clipped_to_float128(xa[i], argv[i+1], NULL);

  __float128 x, y, res;
  mpfr_to_float128(&x, xa[0]);
  mpfr_to_float128(&y, xa[1]);
  uint64_t xu[2], yu[2];
  memcpy(xu, &x, sizeof(xu));
  memcpy(yu, &y, sizeof(yu));
  mpfr_printf(
   "x    %-+45.28Ra %016llx:%016llx %+-54.40Re\n"
   "y    %-+45.28Ra %016llx:%016llx %+-54.40Re\n"
   , xa[0], xu[1], xu[0], xa[0]
   , xa[1], yu[1], yu[0], xa[1]
  );

  MPFR_DECL_INIT(res_max,         113);
  MPFR_DECL_INIT(res_nrm_min,     113);
  MPFR_DECL_INIT(res_subnormal_h, 113);

  __float128 flt128_mx = FLT128_MAX;
  float128_to_mpfr(res_max, &flt128_mx);

  __float128 flt128_mn = FLT128_MIN;
  float128_to_mpfr(res_nrm_min, &flt128_mn);

  __float128 flt128_subnormal_min = nextafterq((__float128)(0.0), FLT128_MAX);
  float128_to_mpfr(res_subnormal_h, &flt128_subnormal_min);
  mpfr_mul_d(res_subnormal_h, res_subnormal_h, 0.5, GMP_RNDN);

  MPFR_DECL_INIT(ref, 113);
  int ref_inexact = mpfr_mul(ref, xa[0], xa[1], GMP_RNDN);
  int ref_ex = ref_inexact ? FE_INEXACT : 0;
  if (mpfr_regular_p(ref)) {
    int r_sign = mpfr_signbit(ref);
    if (mpfr_cmpabs(ref, res_max) > 0) {
      mpfr_set_inf(ref, r_sign ? -1 : 1);
      ref_ex |= FE_OVERFLOW | FE_INEXACT;
    } else if (mpfr_cmpabs(ref, res_nrm_min) <= 0) {
      ref_ex |= FE_UNDERFLOW;
      MPFR_DECL_INIT(ref_1st, 113);
      mpfr_set(ref_1st, ref, GMP_RNDN);
      if (mpfr_cmpabs(ref, res_subnormal_h) <= 0) {
        mpfr_set_d(ref, r_sign ? -0.0 : 0.0, GMP_RNDN);
      } else { // subnormal
        if (r_sign) {
          mpfr_fms(ref, xa[0], xa[1], res_nrm_min, GMP_RNDN);
          mpfr_add(ref, ref, res_nrm_min, GMP_RNDN);
        } else {
          mpfr_fma(ref, xa[0], xa[1], res_nrm_min, GMP_RNDN);
          mpfr_sub(ref, ref, res_nrm_min, GMP_RNDN);
        }
      }
      if (!mpfr_equal_p(ref, ref_1st))
        ref_ex |= FE_INEXACT;
    }
  }

  memset(&res, 42, sizeof(res));
  feclearexcept(FE_ALL_EXCEPT);
  res = x * y;
  int res_ex = fetestexcept(FE_ALL_EXCEPT);

  MPFR_DECL_INIT(resx, 113);
  float128_to_mpfr(resx, &res);
  int succ = mpfr_equal_p(ref, resx);
  if (!succ) {
    if (mpfr_nan_p(ref) && mpfr_nan_p(resx))
      succ = true;
  }

  char resbuf[256];//, ybuf[256], rbuf[256];
  // quadmath_snprintf(xbuf, sizeof(xbuf), "%+-.28Qa", x);
  // quadmath_snprintf(ybuf, sizeof(ybuf), "%+-.28Qa", y);
  quadmath_snprintf(resbuf, sizeof(resbuf), "%+-45.28Qa", res);
  uint64_t resu[2];
  memcpy(resu, &res, sizeof(resu));
  mpfr_printf(
   "res  %s %016llx:%016llx\n"
   "resx %-+45.28Ra %+-54.40Re%s%s%s%s\n"
   "ref  %-+45.28Ra %+-54.40Re%s%s%s%s\n"
   "%s\n"
   , resbuf, resu[1], resu[0]
   , resx, resx
   , res_ex & FE_INVALID   ? " Invalid"   : ""
   , res_ex & FE_INEXACT   ? " Inexact"   : ""
   , res_ex & FE_OVERFLOW  ? " Overflow"  : ""
   , res_ex & FE_UNDERFLOW ? " Underflow" : ""
   , ref,  ref
   , ref_ex & FE_INVALID   ? " Invalid"   : ""
   , ref_ex & FE_INEXACT   ? " Inexact"   : ""
   , ref_ex & FE_OVERFLOW  ? " Overflow"  : ""
   , ref_ex & FE_UNDERFLOW ? " Underflow" : ""
   , succ ? "o.k." : "fail"
  );

  if (!succ) {
    MPFR_DECL_INIT(xx, 256); float128_to_mpfr(xx, &x);
    MPFR_DECL_INIT(yx, 113); float128_to_mpfr(yx, &y);
    mpfr_printf(
     "xx   %-+45.28Ra %d\n"
     "yx   %-+45.28Ra %d\n"
     , xx, mpfr_cmp(xa[0], xx)
     , yx, mpfr_cmp(xa[1], yx)
    );
  }

  return !succ;
}
