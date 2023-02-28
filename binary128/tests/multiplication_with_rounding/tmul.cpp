#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <cfenv>
#include <random>
#include <thread>
#include <mutex>
#include <algorithm>

#include <quadmath.h>
#include <mpfr.h>

#include "qmx2mpfr.h"
#include "mulq_test_generator.h"

extern "C" {
void mulq_all_rm(__float128 res[4], __float128 values[2], int exceptions[4]);
}

struct test_context {
  std::mutex sema;
  long long  nIter;
  long long  nMism[4];
  long long  nMismTotal;
  long long  nSamples;
  long long  nZeros;
  long long  nUnderflow;
  long long  nSubnormal;
  long long  nNormal;
  long long  nOverflow;
  long long  nInf;
  long long  nNan;
};

void test(long long it0, long long it1, test_context* context)
{
  // unsigned long long nIter = context->nIter;
  std::mt19937_64 rndGen(it0+1);

  MPFR_DECL_INIT(xx,   113);
  MPFR_DECL_INIT(yx,   113);
  MPFR_DECL_INIT(ref,  113);
  MPFR_DECL_INIT(resx, 113);
  MPFR_DECL_INIT(res_max, 113);
  MPFR_DECL_INIT(res_nrm_min, 113);
  MPFR_DECL_INIT(subnormal_adj, 113);
  MPFR_DECL_INIT(ref_1st, 113);

  __float128 flt128_mx = FLT128_MAX;
  float128_to_mpfr(res_max, &flt128_mx);

  __float128 flt128_mn = FLT128_MIN;
  float128_to_mpfr(res_nrm_min, &flt128_mn);

  long long nSamples   = 0;
  long long nZeros     = 0;
  long long nUnderflow = 0;
  long long nSubnormal = 0;
  long long nNormal    = 0;
  long long nOverflow  = 0;
  long long nInf       = 0;
  long long nNan       = 0;

  for (long long i = it0; i < it1; ++i) {
    const int RA_LEN = 1+2+2;
    uint64_t ra[RA_LEN];
    for (int k = 0; k < RA_LEN; ++k)
      ra[k] = rndGen();

    __float128 xy[2];
    make_test_values_for_mulq(xy, ra);

    float128_to_mpfr(xx, &xy[0]);
    float128_to_mpfr(yx, &xy[1]);

    int results_ex[4];
    __float128 results[4];
    mulq_all_rm(results, xy, results_ex);

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
      int ref_inexact = mpfr_mul(ref, xx, yx, mpfr_mode);
      int ref_ex = ref_inexact ? FE_INEXACT : 0;
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
          ref_ex |= FE_OVERFLOW | FE_INEXACT;
          ++nOverflow;
        } else if (mpfr_cmpabs(ref, res_nrm_min) < 0) {
          // subnormal or zero
          mpfr_set(ref_1st, ref, GMP_RNDN);
          mpfr_setsign(subnormal_adj, res_nrm_min, r_sign, MPFR_RNDN);
          mpfr_fma(ref, xx, yx, subnormal_adj, mpfr_mode);
          mpfr_sub(ref, ref,    subnormal_adj, mpfr_mode);
          if (mpfr_zero_p(ref)) {
            mpfr_setsign(ref, ref, r_sign, MPFR_RNDN);
            ++nUnderflow;
          } else {
            ++nSubnormal;
          }
          if (!mpfr_equal_p(ref, ref_1st))
            ref_ex |= FE_INEXACT;
          if (ref_ex & FE_INEXACT)
            ref_ex |= FE_UNDERFLOW;
          // Underflow detected after rounding as permitted by
          // option (a) of IEEE Std 754-2008 paragraph 7.5
        } else {
          ++nNormal;
        }
      } else {
        if (mpfr_nan_p(ref))
          ++nNan;
        else if (mpfr_zero_p(ref))
          ++nZeros;
        else
          ++nInf;
      }

      float128_to_mpfr(resx, &results[mode_i]);
      int equal = mpfr_total_order_p(ref, resx) && mpfr_total_order_p(resx, ref);
      if ((ref_ex ^ results_ex[mode_i]) & FE_OVERFLOW)
        equal = 0;
      if (!equal) {
        if (!mpfr_nan_p(ref) || !mpfr_nan_p(resx)) {
          context->sema.lock();
          if (context->nMismTotal < 10) {
            char xbuf[256], ybuf[256], rbuf[256];
            quadmath_snprintf(xbuf, sizeof(xbuf), "%+-.28Qa", xy[0]);
            quadmath_snprintf(ybuf, sizeof(ybuf), "%+-.28Qa", xy[1]);
            quadmath_snprintf(rbuf, sizeof(rbuf), "%+-.28Qa", results[mode_i]);
            mpfr_printf(
             "%s\n"
             "x    %-+45.28Ra %-50.36Re %s\n"
             "y    %-+45.28Ra %-50.36Re %s\n"
             "res  %-+45.28Ra %-50.36Re %s%s\n"
             "ref  %-+45.28Ra %-50.36Re%s\n"
             ,rm_db[mode_i].str
             ,xx,   xx, xbuf
             ,yx,   yx, ybuf
             ,resx, resx, rbuf
             ,results_ex[mode_i] & FE_OVERFLOW  ? " Overflow"  : ""
             ,ref,  ref
             ,ref_ex & FE_OVERFLOW  ? " Overflow"  : ""
            );
            fflush(stdout);
          }
          ++context->nMism[mode_i];
          ++context->nMismTotal;
          context->sema.unlock();
        }
      }
    }
    ++nSamples;
  }
  context->sema.lock();
  context->nSamples   += nSamples;
  context->nZeros     += nZeros;
  context->nUnderflow += nUnderflow;
  context->nSubnormal += nSubnormal;
  context->nNormal    += nNormal;
  context->nOverflow  += nOverflow;
  context->nInf       += nInf;
  context->nNan       += nNan;
  context->sema.unlock();
}


int main(int argz, char** argv)
{
  test_context tc;
  tc.nIter = 100000;
  if (argz > 1) {
    long long val = strtoll(argv[1], 0, 0);
    if (val > 0 && val < 1000ll*1000ll*1000ll*100)
      tc.nIter = val;
  }

  tc.nMismTotal = 0;
  tc.nSamples   = 0;
  tc.nZeros     = 0;
  tc.nUnderflow = 0;
  tc.nSubnormal = 0;
  tc.nNormal    = 0;
  tc.nOverflow  = 0;
  tc.nInf       = 0;
  tc.nNan       = 0;
  memset(tc.nMism, 0, sizeof(tc.nMism));
  const int N_THREADS = 64;
  long long itPerThread = (tc.nIter-1)/N_THREADS + 1;
  std::thread tid[N_THREADS];
  for (int i = 0; i < N_THREADS; ++i)
    tid[i] = std::thread(test, itPerThread*i, std::min(itPerThread*(i+1), tc.nIter), &tc);
  for (int i = 0; i < N_THREADS; ++i)
    tid[i].join();

  printf("%lld iterations, %lld (%lld %lld/%lld/%lld/%lld %lld %lld) samples, %lld (%lld+%lld+%lld+%lld) mismatches.\n"
    , tc.nIter
    , tc.nSamples
    , tc.nZeros
    , tc.nUnderflow
    , tc.nSubnormal
    , tc.nNormal
    , tc.nOverflow
    , tc.nInf
    , tc.nNan
    , tc.nMismTotal
    , tc.nMism[0], tc.nMism[1], tc.nMism[2], tc.nMism[3]);

  return tc.nMismTotal != 0;
}
