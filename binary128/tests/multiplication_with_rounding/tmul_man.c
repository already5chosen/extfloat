#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpfr.h>
#include <quadmath.h>

#include "qmx2mpfr.h"

void mulq_all_rm(__float128 res[4], __float128 values[2]);

int main(int argz, char** argv)
{
  mpfr_t xa[2];
  mpfr_init2(xa[0], 113);
  mpfr_init2(xa[1], 113);

  for (int i = 0; i < 2 && i < argz-1; ++i)
    mpfr_strtofr(xa[i], argv[i+1], NULL, 0, GMP_RNDN);

  __float128 xy[2], results[4];
  mpfr_to_float128(&xy[0], xa[0]);
  mpfr_to_float128(&xy[1], xa[1]);
  uint64_t xu[2], yu[2];
  memcpy(xu, &xy[0], sizeof(xu));
  memcpy(yu, &xy[1], sizeof(yu));
  mpfr_printf(
   "x    %-+45.28Ra %016llx:%016llx %+-54.40Re\n"
   "y    %-+45.28Ra %016llx:%016llx %+-54.40Re\n"
   , xa[0], xu[1], xu[0], xa[0]
   , xa[1], yu[1], yu[0], xa[1]
  );

  MPFR_DECL_INIT(res_max,         113);
  MPFR_DECL_INIT(res_nrm_min,     113);
  MPFR_DECL_INIT(subnormal_adjust,113);

  __float128 flt128_mx = FLT128_MAX;
  float128_to_mpfr(res_max, &flt128_mx);

  __float128 flt128_mn = FLT128_MIN;
  float128_to_mpfr(res_nrm_min, &flt128_mn);

  MPFR_DECL_INIT(ref, 113);
  MPFR_DECL_INIT(resx, 113);

  mulq_all_rm(results, xy);
  int err = 0;
  for (int mode_i = 0; mode_i < 4; ++mode_i) {
    static const struct {
      const char* str;
      mpfr_rnd_t  mpfr_mode;
    } rm_db[4] = {
      {"rounding towards negative infinity",           MPFR_RNDD},
      {"rounding towards nearest representable value", MPFR_RNDN},
      {"rounding towards zero",                        MPFR_RNDZ},
      {"rounding towards positive infinity",           MPFR_RNDU},
    };

    mpfr_rnd_t mpfr_mode = rm_db[mode_i].mpfr_mode;
    mpfr_mul(ref, xa[0], xa[1], mpfr_mode);
    if (mpfr_regular_p(ref)) {
      int r_sign = mpfr_signbit(ref);
      bool round_toward_zero =
          mpfr_mode == MPFR_RNDZ ||
         (mpfr_mode == MPFR_RNDD && r_sign==0) ||
         (mpfr_mode == MPFR_RNDU && r_sign==1);
      if (mpfr_cmpabs(ref, res_max) > 0) {
        if (round_toward_zero)
          mpfr_setsign(ref, res_max, r_sign, MPFR_RNDN);
        else
          mpfr_set_inf(ref, r_sign ? -1 : 1);
      } else if (mpfr_cmpabs(ref, res_nrm_min) < 0) {
        mpfr_setsign(subnormal_adjust, res_nrm_min, r_sign, MPFR_RNDN);
        mpfr_fma(ref, xa[0], xa[1], subnormal_adjust, mpfr_mode);
        mpfr_sub(ref, ref,          subnormal_adjust, mpfr_mode);
        mpfr_setsign(ref, ref, r_sign, MPFR_RNDN);
      }
    }

    float128_to_mpfr(resx, &results[mode_i]);
    int succ = mpfr_equal_p(ref, resx);
    if (!succ) {
      if (mpfr_nan_p(ref) && mpfr_nan_p(resx))
        succ = true;
    }

    char resbuf[256];
    quadmath_snprintf(resbuf, sizeof(resbuf), "%+-45.28Qa", results[mode_i]);
    uint64_t resu[2];
    memcpy(resu, &results[mode_i], sizeof(resu));
    mpfr_printf(
     "%s\n"
     "res  %s %016llx:%016llx\n"
     "resx %-+45.28Ra %+-54.40Re\n"
     "ref  %-+45.28Ra %+-54.40Re\n"
     "%s\n"
     , rm_db[mode_i].str
     , resbuf, resu[1], resu[0]
     , resx, resx
     , ref,  ref
     , succ ? "o.k." : "fail"
    );

    if (!succ)
      err = 1;
  }

  if (err) {
    MPFR_DECL_INIT(xx, 256); float128_to_mpfr(xx, &xy[0]);
    MPFR_DECL_INIT(yx, 113); float128_to_mpfr(yx, &xy[1]);
    mpfr_printf(
     "\n"
     "xx   %-+45.28Ra %d\n"
     "yx   %-+45.28Ra %d\n"
     "fail!!!\n"
     , xx, mpfr_cmp(xa[0], xx)
     , yx, mpfr_cmp(xa[1], yx)
    );
  }
  return err;
}
