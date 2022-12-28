#include <stdint.h>
#include <quadmath.h>
#include <gmp.h>
#include <mpfr.h>

#include "qmx2mpfr.h"

void float128_to_mpfr(mpfr_t dst, const __float128* src)
{
  if (isnanq(*src)) {
    mpfr_set_nan(dst);
    return;
  }

  if (isinfq(*src)) {
    mpfr_set_inf(dst, signbitq(*src) ? -1 : 1);
    return;
  }

  if (*src==0) {
    mpfr_set_zero(dst, signbitq(*src) ? -1 : 1);
    return;
  }

  int exp;
  __float128 srcnrm = frexpq(*src, &exp);
  __float128 srcnrm60 = ldexpq(srcnrm, 60);
  __float128 srcHi = truncq(srcnrm60);
  __float128 srcLo = ldexpq(srcnrm60 - srcHi, 60);
  long long iSrcHi = llroundq(srcHi);
  long long iSrcLo = llroundq(srcLo);
  MPFR_DECL_INIT(mpHi, 128);
  MPFR_DECL_INIT(mpLo, 128);
  mpfr_set_sj_2exp(mpHi, iSrcHi, exp-60,  GMP_RNDN);
  mpfr_set_sj_2exp(mpLo, iSrcLo, exp-120, GMP_RNDN);
  mpfr_add(mpHi, mpHi, mpLo, GMP_RNDN);
  mpfr_set(dst, mpHi, GMP_RNDN);
}

void mpfr_to_float128(__float128* dst, mpfr_t src)
{
  if (mpfr_regular_p(src)) {
    MPFR_DECL_INIT(tmp, 113);
    mpfr_exp_t exp;
    mpfr_frexp(&exp, tmp, src, GMP_RNDN);
    if (exp < INT_MIN) {
      *dst = mpfr_signbit(src) ? -0.0q : +0.0q;
    } else if (exp > INT_MAX) {
      *dst = mpfr_signbit(src) ? -HUGE_VALQ : HUGE_VALQ;
    } else {
      double d1 = mpfr_get_d(tmp, GMP_RNDN);
      mpfr_sub_d(tmp, tmp, d1, GMP_RNDN);
      double d2 = mpfr_get_d(tmp, GMP_RNDN);
      mpfr_sub_d(tmp, tmp, d2, GMP_RNDN);
      *dst = mpfr_get_d(tmp, GMP_RNDN);
      *dst += d2;
      *dst += d1;
      *dst = ldexpq(*dst, (int)exp);
    }
  } else {
    if        (mpfr_zero_p(src)) {
      *dst = mpfr_signbit(src) ? -0.0q : +0.0q;
    } else if (mpfr_nan_p (src)) {
      *dst = nanq("");
    } else {
      *dst = mpfr_signbit(src) ? -HUGE_VALQ : HUGE_VALQ;
    }
  }
}
