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

struct test_context {
  std::mutex sema;
  long long  nIter;
  long long  nMism;
  long long  nSamples;
};

static void make_subnormal_mantissa(uint64_t w[2], const uint64_t ra[2], int nbits)
{
  if (nbits <= 64) {
    w[0] = (ra[0] >> (64-nbits)) | (uint64_t(1) << (nbits-1));
    w[1] = 0;
  } else {
    w[0] = ra[0];
    w[1] = (ra[1] >> (128-nbits)) | (uint64_t(1) << (nbits-65));
  }
}

static int make_float128(__float128* dst, const uint64_t ra[2], int fp_class, int min_exp, int max_exp)
{
  const uint64_t BIT_63 = (uint64_t)1 << 63;
  const uint64_t BIT_48 = (uint64_t)1 << 48;
  const uint64_t MSK_48 = BIT_48 - 1;
  const uint64_t MSK_15 = (uint64_t)(-1) >> (64-15);
  uint64_t w[2]={0};
  int retexp = 0x7FFF;
  if (fp_class >= 4) { // 252 out of 256 are fully random, so, mostly normal
    w[0] = ra[0];
    w[1] = ra[1];
    uint64_t mh = ra[1] & MSK_48;
    uint64_t s  = ra[1] & BIT_63;
    uint64_t exp = (ra[1] >> 48) & MSK_15;
    retexp = exp;
    if (fp_class < 224) { // 220 out of 252 has exponent in reduced range to make arithmetic more involved
      int exp1 = int((exp * (max_exp-min_exp+1)) >> 15) + min_exp;
      retexp = exp1;
      exp = exp1;
      if (exp1 <= 0) { // subnormal
        int nbits = exp1+112;
        if (nbits < 1) nbits=1;
        make_subnormal_mantissa(w, ra, nbits);
      }
    }
    w[1] = s | (exp << 48) | mh;
  } else {
    if (fp_class == 0)        { // zero
    } else if (fp_class == 1) { // subnormal
      uint64_t exp = (ra[1] >> 48) & MSK_15;
      int nbits = ((exp*112)>>15)+1;
      make_subnormal_mantissa(w, ra, nbits);
      retexp = nbits - 112;
    } else if (fp_class == 3) { // infinity
      w[1] = uint64_t(MSK_15) << 48;
    } else                    { // NaN
      w[0] = ra[0];
      w[1] = (uint64_t(MSK_15) << 48) | (ra[1] >> 48);
    }
    w[1] |= (ra[1] & BIT_63); // copy sign
  }
  memcpy(dst, w, sizeof(*dst));
  return retexp;
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

  long long nSamples = 0;
  for (long long i = it0; i < it1; ++i) {
    const int RA_LEN = 1+2+2;
    uint64_t ra[RA_LEN];
    for (int k = 0; k < RA_LEN; ++k)
      ra[k] = rndGen();

    uint8_t  x_class = ra[0] >> (8*0);
    uint8_t  y_class = ra[0] >> (8*1);
    uint8_t  r_class = ra[0] >> (8*2);

    __float128 x, y;
    int x_exp = make_float128(&x, &ra[1], x_class, 0, 0x7FFE);
    int y_exp_min = 0;
    int y_exp_max = 0x7FFE;
    if (x_exp < 0x7FFF) {
      y_exp_min = x_exp - 128;
      if (y_exp_min < -111)
        y_exp_min = -111;
      y_exp_max = x_exp + 128;
      if (y_exp_max > 0x7FFE)
        y_exp_min = 0x7FFE;
    }
    int y_exp = make_float128(&y, &ra[3], y_class, y_exp_min, y_exp_max);
    if (x_exp < 0x7FFF && y_exp < 0x7FFF) {
      if (r_class == 0) { // force result to be near subnormal range
        int max_exp = x_exp > y_exp ? x_exp : y_exp;
        int min_exp = x_exp + y_exp - max_exp;
        int de = ((min_exp - 10)/0x0800 - 1)*0x800;
        x = ldexpq(x, de);
        y = ldexpq(y, de);
      }
    }

    float128_to_mpfr(xx, &x);
    float128_to_mpfr(yx, &y);

    __float128 res;
#ifdef SUBTRACT
  res = x - y;
  mpfr_sub(ref, xx, yx, GMP_RNDN);
#else
  res = x + y;
  mpfr_add(ref, xx, yx, GMP_RNDN);
#endif

    float128_to_mpfr(resx, &res);
    if (mpfr_regular_p(ref)) {
      int r_sign = mpfr_signbit(ref);
      if (mpfr_cmpabs(ref, res_max) > 0) {
        mpfr_set_inf(ref, r_sign ? -1 : 1);
      }
    }

    int equal = mpfr_total_order_p(ref, resx) && mpfr_total_order_p(resx, ref);
    if (!equal) {
      if (!mpfr_nan_p(ref) || !mpfr_nan_p(resx)) {
        context->sema.lock();
        if (context->nMism < 10) {
          char xbuf[256], ybuf[256], rbuf[256];
          quadmath_snprintf(xbuf, sizeof(xbuf), "%+-.28Qa", x);
          quadmath_snprintf(ybuf, sizeof(ybuf), "%+-.28Qa", y);
          quadmath_snprintf(rbuf, sizeof(rbuf), "%+-.28Qa", res);
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
  context->sema.lock();
  context->nSamples += nSamples;
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

  tc.nMism = 0;
  tc.nSamples = 0;
  const int N_THREADS = 64;
  long long itPerThread = (tc.nIter-1)/N_THREADS + 1;
  std::thread tid[N_THREADS];
  for (int i = 0; i < N_THREADS; ++i)
    tid[i] = std::thread(test, itPerThread*i, std::min(itPerThread*(i+1), tc.nIter), &tc);
  for (int i = 0; i < N_THREADS; ++i)
    tid[i].join();

  printf("%lld iterations, %lld samples, %lld mismatches.\n", tc.nIter, tc.nSamples, tc.nMism);

  return tc.nMism != 0;
}
