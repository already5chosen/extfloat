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

void addq_all_rm(__float128 res[4], __float128 values[2], int exceptions[4]);

int main(int argz, char** argv)
{
  if (argz < 3)
    return 1;
    
  mpfr_t xa[2];
  __float128 xy[2];
  for (int i = 0; i < 2 && i < argz-1; ++i)
    manual_input_parser(xa[i], &xy[i], argv[i+1]);

  MPFR_DECL_INIT(ref,  113);
  MPFR_DECL_INIT(resx, 113);
  MPFR_DECL_INIT(res_max, 113);
  int err = 0;

  __float128 results[4];
  int results_ex[4];
  addq_all_rm(results, xy, results_ex);

  __float128 flt128_mx = FLT128_MAX;
  float128_to_mpfr(res_max, &flt128_mx);

  uint64_t xu[2], yu[2];
  memcpy(xu, &xy[0], sizeof(xu));
  memcpy(yu, &xy[1], sizeof(yu));
  mpfr_printf(
   "x    %-+45.28Ra %016llx:%016llx %+-54.40Re\n"
   "y    %-+45.28Ra %016llx:%016llx %+-54.40Re\n"
   , xa[0], xu[1], xu[0], xa[0]
   , xa[1], yu[1], yu[0], xa[1]
  );
  fflush(stdout);

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
    int ref_inexact = mpfr_add(ref, xa[0], xa[1], mpfr_mode);
    int ref_ex = ref_inexact ? FE_INEXACT : 0;
    if (mpfr_regular_p(ref) && mpfr_cmpabs(ref, res_max) > 0) {
      int r_sign = mpfr_signbit(ref);
      if (mpfr_mode == MPFR_RNDZ ||
         (mpfr_mode == MPFR_RNDD && r_sign==0) ||
         (mpfr_mode == MPFR_RNDU && r_sign==1)) {
        mpfr_setsign(ref, res_max, r_sign, MPFR_RNDN);
      } else {
        mpfr_set_inf(ref, r_sign ? -1 : 1);
      }
      ref_ex |= FE_OVERFLOW | FE_INEXACT;
    }

    float128_to_mpfr(resx, &results[mode_i]);
    int succ = mpfr_total_order_p(ref, resx) && mpfr_total_order_p(resx, ref);
    if (!succ && mpfr_nan_p(ref) && mpfr_nan_p(resx))
      succ = true;

    uint64_t resu[2];
    memcpy(resu, &results[mode_i], sizeof(resu));
    char resbuf[256];
    quadmath_snprintf(resbuf, sizeof(resbuf), "%+-.28Qa", results[mode_i]);
    int res_ex = results_ex[mode_i];
    mpfr_printf(
     "%s\n"
     "res  %-45s %016llx:%016llx\n"
     "resx %-+45.28Ra %+-54.40Re%s%s%s%s\n"
     "ref  %-+45.28Ra %+-54.40Re%s%s%s%s\n"
     "%s\n"
     , rm_db[mode_i].str
     , resbuf
     , resu[1], resu[0]
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
