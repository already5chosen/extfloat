#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpfr.h>
#include <quadmath.h>

#include "qmx2mpfr.h"

int main(int argz, char** argv)
{
  mpfr_t xa[2];
  mpfr_init2(xa[0], 113);
  mpfr_init2(xa[1], 113);
  MPFR_DECL_INIT(flt128_nrm_min, 113);

  mpfr_set_ui_2exp(flt128_nrm_min, 1, -16382, GMP_RNDN);
  for (int i = 0; i < 2 && i < argz-1; ++i) {
    mpfr_strtofr(xa[i], argv[i+1], NULL, 0, GMP_RNDN);
    if (mpfr_cmpabs(xa[i], flt128_nrm_min) < 0) {
      int sign = mpfr_signbit(xa[i]);
      mpfr_setsign(xa[i], xa[i], 0, GMP_RNDN);
      mpfr_add(xa[i], xa[i], flt128_nrm_min, GMP_RNDN);
      mpfr_sub(xa[i], xa[i], flt128_nrm_min, GMP_RNDN);
      mpfr_setsign(xa[i], xa[i], sign, GMP_RNDN);
    }
  }


  MPFR_DECL_INIT(ref, 113);
  mpfr_add(ref, xa[0], xa[1], GMP_RNDN);

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
  fflush(stdout);

  res = x + y;

  MPFR_DECL_INIT(resx, 113);
  float128_to_mpfr(resx, &res);
  int succ = mpfr_equal_p(ref, resx);
  if (!succ && mpfr_nan_p(ref) && mpfr_nan_p(resx))
    succ = true;

  uint64_t resu[2];
  memcpy(resu, &res, sizeof(resu));
  char resbuf[256];
  quadmath_snprintf(resbuf, sizeof(resbuf), "%+-.28Qa", res);
  mpfr_printf(
   "res  %-45s %016llx:%016llx\n"
   "resx %-+45.28Ra %+-54.40Re\n"
   "ref  %-+45.28Ra %+-54.40Re\n"
   "%s\n"
   , resbuf
   , resu[1], resu[0]
   , resx, resx
   , ref,  ref
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
