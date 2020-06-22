#include <stdio.h>
#include <stdint.h>
#include <cmath>
#include <random>
#include <functional>           // for std::bind
#include <iostream>
#include <intrin.h>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/math/constants/constants.hpp>

#include "extfloat128.h"
#include "extfloat128_to_bobinf.h"
#include "bobinf_to_extfloat128.h"

typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<128, boost::multiprecision::backends::digit_base_2, void, boost::int64_t, INT_MIN, INT_MAX> > boost_quadfloat_t;
typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<256, boost::multiprecision::backends::digit_base_2, void, boost::int64_t, INT_MIN, INT_MAX> > boost_octafloat_t;
typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<320, boost::multiprecision::backends::digit_base_2, void, boost::int64_t, INT_MIN, INT_MAX> > boost_decafloat_t;

static const double MIN_N_ITER  = 1;
static const double MAX_N_ITER  = 1e10;

static extfloat128_t make_random_quadfloat(std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr, int distrib = 0);
static void make_nonrandom_quadfloats(std::vector<extfloat128_t>* dst);
static void speed_test(std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr);

static bool isEqual(const extfloat128_t& a, const extfloat128_t& b) {
  if (isnan(a))             return isnan(b);
  if (a.m_sign != b.m_sign) return false;
  return a == b;
}

static extfloat128_t make_quadfloat(int sign, unsigned exp, uint64_t msw, uint64_t lsw) {
  extfloat128_t dst;
  dst.m_sign           = sign;
  dst.m_exponent       = exp;
  dst.m_significand[0] = lsw;
  dst.m_significand[1] = msw;
  return dst;
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
//static void report_mismatch_ex(extfloat128_t a);
int main(int argz, char** argv)
{
  if (argz < 2)
  {
    fprintf(stderr,
      "tst_atan - validate implementation of atanPi().\n"
      "Usage:\n"
      "tst_atan nIter [u | r | b]\n"
      "where\n"
      " u - test with arguments uniformly distributed in range (-1..1)\n"
      " r - test with results uniformly distributed in range (-1/2..1/2)\n"
      " b - test boost implementation with results uniformly distributed in range (-4..4)\n"
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
  int distrib = 0;
  bool testBoost = false;
  if (argz > 2) {
    distrib = argv[2][0];
    if (distrib == 'b') {
      testBoost = true;
      distrib = 'r';
    }
  }

  std::vector<extfloat128_t> specialPoints;
  make_nonrandom_quadfloats(&specialPoints);

  std::mt19937_64 rndGen;
  std::uniform_int_distribution<uint64_t> rndDistr(0, uint64_t(-1));

  boost_octafloat_t x_o, yRef_o, yRes_o, yDiff_o, yUlp_o;
  uint64_t mCnt[6] = {0};
  static const boost_octafloat_t oPi = boost::math::constants::pi<boost_octafloat_t>();
  static const boost_quadfloat_t qPi = boost::math::constants::pi<boost_quadfloat_t>();
  double maxErr_v = 0;
  extfloat128_t maxErr_x = 0;

  int64_t nSpecial = specialPoints.size();
  nIter += nSpecial;
  for (int64_t n = 0;n < nIter; ++n) {
    extfloat128_t x = (n < nSpecial) ? specialPoints[n] : make_random_quadfloat(rndGen, rndDistr, distrib);
    convert_to_boost_bin_float(&x_o, x);
    yRef_o = atan(x_o)/oPi;
    // std::cout << x_o << "=>" << yRef_o << "\n";
    extfloat128_t yRes;
    if (!testBoost) {
      yRes = atanPi(x);
    } else {
      boost_quadfloat_t yRes_q = atan(boost_quadfloat_t(x_o))/qPi;
      convert_from_boost_bin_float(&yRes, yRes_q);
    }
    extfloat128_t yRef; convert_from_boost_bin_float(&yRef, yRef_o);
    if (!isEqual(yRef, yRes)) {
      if (isnan(yRes) || isnan(yRef)) {
        if (report_mismatch(x, n))
          return 1;
      } else {
        ++mCnt[0];
        extfloat128_t yUlp = yRef.ulp();
        convert_to_boost_bin_float(&yUlp_o, yUlp);
        convert_to_boost_bin_float(&yRes_o, yRes);
        yDiff_o = fabs(yRes_o - yRef_o);
        if (!testBoost) {
          if (yDiff_o >= yUlp_o)
            if (report_mismatch(x, n))
              return 1;
        } else {
          if (yDiff_o >= yUlp_o*16)
            if (report_mismatch(x, n))
              return 1;
        }
        double err = (yDiff_o / yUlp_o).convert_to<double>();
        if (err >= 33./64) {
          ++mCnt[1];
          if (err >= 17./32) {
            ++mCnt[2];
            if (err >= 9./16) {
              ++mCnt[3];
              if (err >= 5./8) {
                ++mCnt[4];
                if (err >= 3./4) {
                  ++mCnt[5];
                }
              }
            }
          }
        }
        if (err > maxErr_v) {
          maxErr_v = err;
          maxErr_x = x;
        }
      }
    }
  }

  printf("o.k. %I64u/%I64u/%I64u/%I64u/%I64u/%I64u mismatches. Max mismatch = %.6f ULP at x = ",
    mCnt[0],
    mCnt[1],
    mCnt[2],
    mCnt[3],
    mCnt[4],
    mCnt[5],
    maxErr_v);
  boost_quadfloat_t b;
  convert_to_boost_bin_float(&b, maxErr_x);
  std::cout << b << "\n";
  report_mismatch(maxErr_x, -1);
  //report_mismatch_ex(maxErr_x);
  speed_test(rndGen, rndDistr);

  return 0;
}

static bool report_mismatch(extfloat128_t a, int64_t n)
{
  boost_quadfloat_t b;
  convert_to_boost_bin_float(&b, a);
  boost_octafloat_t o(b);
  extfloat128_t     ra = atanPi(a);
  boost_quadfloat_t rb = atan(b)/boost::math::constants::pi<boost_quadfloat_t>();
  boost_octafloat_t ro = atan(o)/boost::math::constants::pi<boost_octafloat_t>();

  if (n >= 0) printf("fail at iteration %I64d.\n", n);
  print(b);
  print(a);
  print(rb);
  print(ra);
  print(ro);
  convert_to_boost_bin_float(&rb, ra);
  print(rb);
  printf("ix=%.8f\n", (ro*64).convert_to<double>());
  return true;
}

static void make_nonrandom_quadfloats(std::vector<extfloat128_t>* dst)
{
  // positive
  dst->push_back(extfloat128_t::zero());

  extfloat128_t x = extfloat128_t::zero();
  for (int i = 0; i < 256; ++i) {
    x = nextinout(x, true);
    dst->push_back(x);
  }

  x = extfloat128_t::one();
  for (int i = 0; i < 256; ++i) {
    x = nextinout(x, false);
    dst->push_back(x);
  }

  for (int i = -150; i < 5; ++i) {
    x = x.pow2(i);
    dst->push_back(x);
    dst->push_back(nextinout(x,false));
    dst->push_back(nextinout(x,true));
  }

  dst->push_back(extfloat128_t::inf());

  // negative
  size_t n = dst->size();
  dst->resize(n*2);
  for (size_t i = 0; i < n; ++i)
    (*dst)[n+i] = -(*dst)[i];

  dst->push_back(extfloat128_t::nan());
}

static extfloat128_t make_random_quadfloat(
  std::mt19937_64& rndGen,
  std::uniform_int_distribution<uint64_t>& rndDistr,
  int  distrib)
{
  auto rndFunc = std::bind ( rndDistr, std::ref(rndGen) );
  uint64_t msw = rndFunc();
  uint64_t lsw = rndFunc();
  uint64_t exw = rndFunc();
  if (distrib == 'u') {
    extfloat128_t dst = extfloat128_t::ldexp(extfloat128_t(msw),   -(64*1));
    dst += extfloat128_t::ldexp(extfloat128_t(lsw),                -(64*2));
    dst += extfloat128_t::ldexp(extfloat128_t(exw & uint64_t(-2)), -(64*3));
    dst.m_sign = exw & 1;
    return dst;
  } else if (distrib == 'r') {
    static const double SCALE48 = 5.5805895968138108216e-015;
    static const double SCALE63 = 5.5805895968138108216e-015/(1 << 15);
    extfloat128_t dst;
    if (exw > (1u << 20))
      dst = tan(int64_t(exw >> 16)*SCALE48); // argument of tan uniformly distributed in range [0..pi/2)
    else
      dst = extfloat128_t(exw>>1)*SCALE63;
    dst.m_significand[0] = lsw;
    dst.m_significand[1] = (dst.m_significand[1] & (uint64_t(-1) << 32)) | (msw >> 32);
    dst.m_sign = exw & 1;
    return dst;
  }

  msw |= (uint64_t(1) << 63);
  uint32_t expLsw = uint32_t(exw);
  uint32_t expMsw = uint32_t(exw >> 32);
  int sign = expMsw & 1;
  uint32_t biased_exp = std::min(std::max(expLsw, extfloat128_t::min_biased_exponent), extfloat128_t::max_biased_exponent);
  unsigned expSel = (expMsw >> 1) & 63;
  if ((expSel & 1) != 0) {
    biased_exp = extfloat128_t::exponent_bias + 16 - (expLsw % 32);
  }
  else {
    // full range of exponent
    ;
  }
  msw |= (uint64_t(1) << 63);

  return make_quadfloat(sign, biased_exp, msw, lsw);
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
  for (int it = 0; it < NITER; ++it) {
    for (int i = 0; i < VECLEN; ++i)
      inpvec_a[i] = make_random_quadfloat(std::ref(rndGen), rndDistr, 'r');

    uint64_t t0 = __rdtsc();
    for (int i = 0; i < VECLEN; ++i)
      outvec_a[i] = atanPi(inpvec_a[i]);
    uint64_t t1 = __rdtsc();
    dt_a[it] = t1 - t0;
    for (int i = 0; i < VECLEN; ++i) {
      dummy[0] ^= outvec_a[i].m_significand[0];
      dummy[1] ^= outvec_a[i].m_significand[1];
      dummy[2] ^= (uint64_t(outvec_a[i].m_sign) << 32) | outvec_a[i].m_exponent;
    }

    for (int i = 0; i < VECLEN; ++i)
      convert_to_boost_bin_float(&inpvec_b[i], inpvec_a[i]);

    t0 = __rdtsc();
    for (int i = 0; i < VECLEN; ++i)
      outvec_b[i] = atan(inpvec_b[i]);
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
