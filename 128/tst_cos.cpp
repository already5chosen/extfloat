#include <stdio.h>
#include <stdint.h>
#include <cmath>
#include <random>
#include <functional>           // for std::bind
#include <iostream>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/math/constants/constants.hpp>

#include "extfloat128.h"
#include "extfloat128_to_bobinf.h"
#include "bobinf_to_extfloat128.h"

typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<128, boost::multiprecision::backends::digit_base_2, void, boost::int64_t, INT_MIN, INT_MAX> > boost_quadfloat_t;
typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<256, boost::multiprecision::backends::digit_base_2, void, boost::int64_t, INT_MIN, INT_MAX> > boost_octafloat_t;

static const double MIN_N_ITER  = 1;
static const double MAX_N_ITER  = 1e10;

static void make_random_quadfloat(extfloat128_t* dst, std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr, bool saneExponent = false, bool uniformal = false);
static void speed_test(std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr);

static void make_quadfloat(extfloat128_t* dst, int sign, unsigned exp, uint64_t msw, uint64_t lsw) {
  dst->m_sign           = sign;
  dst->m_exponent      = exp;
  dst->m_significand[0] = lsw;
  dst->m_significand[1] = msw;
}

static void print(const boost_quadfloat_t& a) {
  printf("%d %11I64d", a.backend().sign(), a.backend().exponent());
  const size_t alen = a.backend().bits().size()*sizeof(a.backend().bits().limbs()[0]);
  if (alen == 16) {
    const uint64_t* bits = reinterpret_cast<const uint64_t*>(a.backend().bits().limbs());
    printf(" %016I64x:%016I64x", bits[1], bits[0]);
  }
  else
    printf(" bits.size = %u*%u", a.backend().bits().size(), unsigned(sizeof(a.backend().bits().limbs()[0])));
  std::cout << " " << a << "\n";
}

static void print(const boost_octafloat_t& a) {
  printf("%d %11I64d", a.backend().sign(), a.backend().exponent());
  const size_t alen = a.backend().bits().size()*sizeof(a.backend().bits().limbs()[0]);
  if (alen == 32) {
    const uint64_t* bits = reinterpret_cast<const uint64_t*>(a.backend().bits().limbs());
    printf(" %016I64x:%016I64x:%016I64x:%016I64x", bits[3], bits[2], bits[1], bits[0]);
  }
  else
    printf(" bits.size = %u", a.backend().bits().size());
  std::cout << " " << a << "\n";
}


static void print(const extfloat128_t& a) {
  printf("%d %11u %016I64x:%016I64x {%d}\n", a.m_sign, a.m_exponent, a.m_significand[1], a.m_significand[0], a._get_exponent());
}

static bool report_mismatch(extfloat128_t a, int64_t n);
int main(int argz, char** argv)
{
  if (argz < 2)
  {
    fprintf(stderr,
      "cos_test - validate implementation of cosPi().\n"
      "Usage:\n"
      "cos_test nIter [u]\n"
      "where\n"
      " u - test with arguments uniformalyy distributed in range (-0.5..0.5)\n"
      );
    return 1;
  }

  char* endp;
  double dnIter = strtod(argv[1], &endp);
  if (endp == argv[1] || dnIter < MIN_N_ITER || dnIter >= MAX_N_ITER+1)
  {
    fprintf(stderr, "Illegal nIter argument '%s'. Please specify number in range [%.0f..%.0f]\n",
      argv[1], MIN_N_ITER, MAX_N_ITER);
    return 1;
  }
  int64_t nIter = int64_t(dnIter);
  bool uniformal = false;
  if (argz > 2 && argv[2][0] == 'u')
    uniformal = true;

  std::mt19937_64 rndGen;
  std::uniform_int_distribution<uint64_t> rndDistr(0, uint64_t(-1));

  boost_octafloat_t x_o, yRef_o, yRes_o, yDiff_o, yUlp_o;
  extfloat128_t     x;
  int64_t n = 1;
  uint64_t mCnt[4] = {0};
  static const boost_octafloat_t oPi = boost::math::constants::pi<boost_octafloat_t>();
  double maxErr_v = 0;
  extfloat128_t maxErr_x = 0;
  while (nIter > 0) {
    make_random_quadfloat(&x, rndGen, rndDistr, false, uniformal);
    convert_to_boost_bin_float(&x_o, x);
    yRef_o = cos(x_o*oPi);
    extfloat128_t yRef; convert_from_boost_bin_float(&yRef, yRef_o);
    extfloat128_t yRes = cosPi(x);
    if (yRef != yRes) {
      if (!isnan(yRes) || !isnan(yRef)) {
        ++mCnt[0];
        // if (!isfinite(yRes) || !isfinite(yRef)) {
          // report_mismatch(x, n);
          // return 1;
        // }
        extfloat128_t yUlp = yRef.ulp();
        convert_to_boost_bin_float(&yUlp_o, yUlp);
        convert_to_boost_bin_float(&yRes_o, yRes);
        yDiff_o = fabs(yRes_o - yRef_o);
        if (yDiff_o >= yUlp_o)
          if (report_mismatch(x, n))
            return 1;
        double err = (yDiff_o / yUlp_o).convert_to<double>();
        if (err >= 9./16) {
          ++mCnt[1];
          if (err >= 5./8) {
            ++mCnt[2];
            if (err >= 3./4) {
              ++mCnt[3];
            }
          }
        }
        if (err > maxErr_v) {
          maxErr_v = err;
          maxErr_x = x;
        }
      }
    }
    ++n;
    nIter -= 1;
  }

  printf("o.k. %I64u/%I64u/%I64u/%I64u mismatches. Max mismatch = %.6f ULP at x = %.17e\n",
    mCnt[0],
    mCnt[1],
    mCnt[2],
    mCnt[3],
    maxErr_v, maxErr_x.convert_to_double());
  report_mismatch(maxErr_x, -1);
  speed_test(rndGen, rndDistr);

  return 0;
}

static bool report_mismatch(extfloat128_t a, int64_t n)
{
  boost_quadfloat_t b;
  convert_to_boost_bin_float(&b, a);
  extfloat128_t     ra = cosPi(a);
  boost_quadfloat_t rb = cos(b*boost::math::constants::pi<boost_quadfloat_t>());
  boost_octafloat_t o(b);
  boost_octafloat_t ro = cos(o*boost::math::constants::pi<boost_octafloat_t>());

  if (n >= 0) printf("fail at iteration %I64d.\n", n);
  print(b);
  print(a);
  print(rb);
  print(ra);
  print(ro);
  return true;
}

static void make_random_quadfloat(
  extfloat128_t* dst,
  std::mt19937_64& rndGen,
  std::uniform_int_distribution<uint64_t>& rndDistr,
  bool saneExponent,
  bool uniformal)
{
  auto rndFunc = std::bind ( rndDistr, std::ref(rndGen) );
  uint64_t msw = rndFunc();
  uint64_t lsw = rndFunc();
  uint64_t exw = rndFunc();
  if (uniformal) {
    *dst  = extfloat128_t::ldexp(extfloat128_t(msw), -(64*1+1));
    *dst += extfloat128_t::ldexp(extfloat128_t(lsw), -(64*2+1));
    *dst += extfloat128_t::ldexp(extfloat128_t(exw & uint64_t(-2)), -(64*3+1));
    dst->m_sign = exw & 1;
    // printf("%016I64x:%016I64x:%016I64x\n", msw, lsw, exw);
    // print(*dst);
    // boost_quadfloat_t b;
    // convert_to_boost_bin_float(&b, *dst);
    // print(b);
    return;
  }

  msw |= (uint64_t(1) << 63);
  uint32_t expLsw = uint32_t(exw);
  uint32_t expMsw = uint32_t(exw >> 32);
  int sign = expMsw & 1;
  uint32_t biased_exp = std::min(std::max(expLsw, extfloat128_t::min_biased_exponent), extfloat128_t::max_biased_exponent);
  if (saneExponent) {
    msw |= (uint64_t(1) << 63);
    make_quadfloat(dst, sign, extfloat128_t::exponent_bias, msw, lsw);
    *dst *= double(int(expLsw>>1)) * (1.0/(1<<(31-3)));
  } else {
    unsigned expSel = (expMsw >> 1) & 63;
    if (expSel != 0) {
      if ((expSel & 1) != 0) {
        biased_exp = extfloat128_t::exponent_bias + 5 - (expLsw % 16);
      } else {
        // full negative range
        if (biased_exp > extfloat128_t::exponent_bias + 5)
          biased_exp -= extfloat128_t::exponent_bias;
      }
      msw |= (uint64_t(1) << 63);
    } else {
      // generate special value of exponent
      msw = lsw = 0;
      switch (expLsw % 3) {
        case 0:
          biased_exp = dst->zero_biased_exponent; // zero
          break;
        case 1:
          biased_exp = dst->inf_nan_biased_exponent; // inf
          break;
        default:
          biased_exp = dst->inf_nan_biased_exponent; // NaN
          msw  = dst->qnan_bit;
          sign = 0;
          break;
      }
    }
    make_quadfloat(dst, sign, biased_exp, msw, lsw);
  }
}

static void speed_test(std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr)
{
  const int VECLEN = 111;
  const int NITER  = 777;
  extfloat128_t      inpvec_a[VECLEN];
  boost_quadfloat_t  inpvec_b[VECLEN];
  extfloat128_t      outvec_a[VECLEN];
  boost_quadfloat_t  outvec_b[VECLEN];
  uint64_t dummy[3] = {0};
  uint64_t dt_a[NITER];
  uint64_t dt_b[NITER];
  static const boost_quadfloat_t bqfPi = boost::math::constants::pi<boost_quadfloat_t>();
  for (int it = 0; it < NITER; ++it) {
    for (int i = 0; i < VECLEN; ++i)
      make_random_quadfloat(&inpvec_a[i], std::ref(rndGen), rndDistr, true);

    uint64_t t0 = __rdtsc();
    for (int i = 0; i < VECLEN; ++i)
      outvec_a[i] = cosPi(inpvec_a[i]);
    uint64_t t1 = __rdtsc();
    dt_a[it] = t1 - t0;
    for (int i = 0; i < VECLEN; ++i) {
      dummy[0] ^= outvec_a[i].m_significand[0];
      dummy[1] ^= outvec_a[i].m_significand[1];
      dummy[2] ^= (uint64_t(outvec_a[i].m_sign) << 32) | outvec_a[i].m_exponent;
    }

    for (int i = 0; i < VECLEN; ++i) {
      convert_to_boost_bin_float(&inpvec_b[i], inpvec_a[i]);
      inpvec_b[i] *= bqfPi;
    }

    t0 = __rdtsc();
    for (int i = 0; i < VECLEN; ++i)
      outvec_b[i] = cos(inpvec_b[i]);
    t1 = __rdtsc();
    dt_b[it] = t1 - t0;
    for (int i = 0; i < VECLEN; ++i) {
      dummy[0] ^= reinterpret_cast<uint64_t*>(outvec_b[i].backend().bits().limbs())[0];
      dummy[1] ^= reinterpret_cast<uint64_t*>(outvec_b[i].backend().bits().limbs())[1];
      dummy[2] ^= (uint64_t(outvec_b[i].backend().sign()) << 32) | outvec_b[i].backend().exponent();
    }
  }
  std::nth_element(&dt_a[0], &dt_a[NITER/2], &dt_a[NITER]);
  std::nth_element(&dt_b[0], &dt_b[NITER/2], &dt_b[NITER]);
  printf("my %I64u. boost %I64u\n", dt_a[NITER/2]/VECLEN, dt_b[NITER/2]/VECLEN);

  printf("%I64u %I64u %I64u\n", dummy[0], dummy[1], dummy[2]);
}