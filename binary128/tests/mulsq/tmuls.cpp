#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <random>
#include <thread>
#include <mutex>
#include <algorithm>

#include <quadmath.h>
#include <mpfr.h>

#include "qmx2mpfr.h"
#include "mulsq.h"


struct test_context {
  std::mutex sema;
  long long  nIter;
  long long  nMism;
  long long  nSamples;
  long long  nZeros;
  long long  nUnderflow;
  long long  nSubnormal;
  long long  nNormal;
  long long  nOverflow;
  long long  nInf;
  long long  nNan;
};

static bool make_float128(__float128* dst, const uint64_t ra[2], int fp_class)
{
  const uint64_t BIT_63 = (uint64_t)1 << 63;
  const uint64_t BIT_48 = (uint64_t)1 << 48;
  const uint64_t MSK_48 = BIT_48 - 1;
  const uint64_t MSK_15 = (uint64_t)(-1) >> (64-15);
  uint64_t w[2]={0};
  bool mundane = false;
  if (fp_class >= 4) { // 252 out of 256 are fully random, so, mostly normal
    w[0] = ra[0];
    w[1] = ra[1];
    uint64_t mh = ra[1] & MSK_48;
    uint64_t s  = ra[1] & BIT_63;
    uint64_t exp = (ra[1] >> 48) & MSK_15;
    if (fp_class < 224) { // 220 out of 252 has reduced range to prevent over/under flow
      exp = exp/2 + MSK_15/4;
      mundane = true;
    }
    w[1] = s | (exp << 48) | mh;
  } else {
    if (fp_class == 0)        { // zero
    } else if (fp_class == 1) { // subnormal
      uint64_t exp = (ra[1] >> 48) & MSK_15;
      int nbits = ((exp*112)>>15)+1;
      if (nbits <= 64) {
        w[0] = (ra[0] >> (64-nbits)) | (uint64_t(1) << (nbits-1));
        w[1] = 0;
      } else {
        w[0] = ra[0];
        w[1] = (ra[1] >> (128-nbits)) | (uint64_t(1) << (nbits-65));
      }
    } else if (fp_class == 3) { // infinity
      w[1] = uint64_t(MSK_15) << 48;
    } else                    { // NaN
      w[0] = ra[0];
      w[1] = (uint64_t(MSK_15) << 48) | (ra[1] >> 48);
    }
    w[1] |= (ra[1] & BIT_63); // copy sign
  }
  memcpy(dst, w, sizeof(*dst));
  return mundane;
}

static void force_y_to_subnormal(__float128* y, int x_exp)
{
  if (x_exp <= 0) {
    const int r_exp_max = -0x4000 + 10;
    int y_exp;
    frexpq(*y, &y_exp);
    int de = ((r_exp_max - x_exp - y_exp)/256 - 1)*256;
    *y = ldexpq(*y, de);
  }
}

void test(long long it0, long long it1, test_context* context)
{
  std::mt19937_64 rndGen(it0+1);

  MPFR_DECL_INIT(xx,   113);
  MPFR_DECL_INIT(yx,   113);
  MPFR_DECL_INIT(ref,  113);
  MPFR_DECL_INIT(resx, 113);
  MPFR_DECL_INIT(res_max, 113);
  MPFR_DECL_INIT(res_nrm_min, 113);
  MPFR_DECL_INIT(res_subnormal_h, 113);

  __float128 flt128_mx = FLT128_MAX;
  float128_to_mpfr(res_max, &flt128_mx);

  __float128 flt128_mn = FLT128_MIN;
  float128_to_mpfr(res_nrm_min, &flt128_mn);

  __float128 flt128_subnormal_min = nextafterq(__float128(0.0), FLT128_MAX);
  float128_to_mpfr(res_subnormal_h, &flt128_subnormal_min);
  mpfr_mul_d(res_subnormal_h, res_subnormal_h, 0.5, GMP_RNDN);

  long long nSamples   = 0;
  long long nZeros     = 0;
  long long nUnderflow = 0;
  long long nSubnormal = 0;
  long long nNormal    = 0;
  long long nOverflow  = 0;
  long long nInf       = 0;
  long long nNan       = 0;

  const int MAX_VLEN = 10;
  for (long long i = it0; i < it1; ++i) {
    int vlen = rndGen() % MAX_VLEN + 1;
    uint64_t ra[(MAX_VLEN+1)*3];
    for (int k = 0; k < (vlen+1)*3; ++k)
      ra[k] = rndGen();

    __float128 x, y[MAX_VLEN];
    uint8_t  x_class = ra[0];
    bool x_mundane = make_float128(&x, &ra[1], x_class);
    for (int k = 0; k < vlen; ++k) {
      uint8_t y_class = ra[(k+1)*3+0] >> (8*0);
      uint8_t r_class = ra[(k+1)*3+0] >> (8*1);
      bool y_mundane = make_float128(&y[k], &ra[(k+1)*3+1], y_class);
      if (x_mundane && y_mundane && r_class == 0) { // force result to be in or near subnormal range
        int x_exp;
        frexpq(x, &x_exp);
        force_y_to_subnormal(&y[k], x_exp);
      }
    }

    __float128 res[MAX_VLEN];
    float128_to_mpfr(xx, &x);
    mulSq(res, y, x, vlen);
    for (int k = 0; k < vlen; ++k) {
      float128_to_mpfr(yx,   &y[k]);
      float128_to_mpfr(resx, &res[k]);

      mpfr_mul(ref, xx, yx, GMP_RNDN);
      if (mpfr_regular_p(ref)) {
        int r_sign = mpfr_signbit(ref);
        if (mpfr_cmpabs(ref, res_max) > 0) {
          mpfr_set_inf(ref, r_sign ? -1 : 1);
          ++nOverflow;
        } else if (mpfr_cmpabs(ref, res_nrm_min) < 0) {
          if (mpfr_cmpabs(ref, res_subnormal_h) <= 0) {
            mpfr_set_d(ref, r_sign ? -0.0 : 0.0, GMP_RNDN);
            ++nUnderflow;
          } else { // subnormal
            if (r_sign) {
              mpfr_fms(ref, xx, yx, res_nrm_min, GMP_RNDN);
              mpfr_add(ref, ref, res_nrm_min, GMP_RNDN);
            } else {
              mpfr_fma(ref, xx, yx, res_nrm_min, GMP_RNDN);
              mpfr_sub(ref, ref, res_nrm_min, GMP_RNDN);
            }
            ++nSubnormal;
          }
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

      int equal = mpfr_total_order_p(ref, resx) && mpfr_total_order_p(resx, ref);
      if (!equal) {
        if (!mpfr_nan_p(ref) || !mpfr_nan_p(resx)) {
          context->sema.lock();
          if (context->nMism < 10) {
            char xbuf[256], ybuf[256], rbuf[256];
            quadmath_snprintf(xbuf, sizeof(xbuf), "%+-.28Qa", x);
            quadmath_snprintf(ybuf, sizeof(ybuf), "%+-.28Qa", y[k]);
            quadmath_snprintf(rbuf, sizeof(rbuf), "%+-.28Qa", res[k]);
            mpfr_printf(
             "x    %-+45.28Ra %-50.36Re %s\n"
             "y    %-+45.28Ra %-50.36Re %s\n"
             "res  %-+45.28Ra %-50.36Re %s\n"
             "ref  %-+45.28Ra %-50.36Re\n"
             ,xx,   xx, xbuf
             ,yx,   yx, ybuf
             ,resx, resx, rbuf
             ,ref,  ref
            );
            fflush(stdout);
          }
          ++context->nMism;
          context->sema.unlock();
        }
      }
      ++nSamples;
    }
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

  tc.nMism      = 0;
  tc.nSamples   = 0;
  tc.nZeros     = 0;
  tc.nUnderflow = 0;
  tc.nSubnormal = 0;
  tc.nNormal    = 0;
  tc.nOverflow  = 0;
  tc.nInf       = 0;
  tc.nNan       = 0;
  const int N_THREADS = 64;
  long long itPerThread = (tc.nIter-1)/N_THREADS + 1;
  std::thread tid[N_THREADS];
  for (int i = 0; i < N_THREADS; ++i)
    tid[i] = std::thread(test, itPerThread*i, std::min(itPerThread*(i+1), tc.nIter), &tc);
  for (int i = 0; i < N_THREADS; ++i)
    tid[i].join();

  printf("%lld iterations, %lld (%lld %lld/%lld/%lld/%lld %lld %lld) samples, %lld mismatches.\n"
    , tc.nIter
    , tc.nSamples
    , tc.nZeros
    , tc.nUnderflow
    , tc.nSubnormal
    , tc.nNormal
    , tc.nOverflow
    , tc.nInf
    , tc.nNan
    , tc.nMism);

  return tc.nMism != 0;
}
