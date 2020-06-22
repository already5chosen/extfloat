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

static void make_random_quadfloat(extfloat128_t* dst, std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr, int distrib = 0);
static void speed_test(std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr);

static bool isEqual(const extfloat128_t& a, const extfloat128_t& b) {
  if (isnan(a))             return isnan(b);
  if (a.m_sign != b.m_sign) return false;
  return a == b;
}

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

static bool report_mismatch(const extfloat128_t x[2], int64_t n);
static bool mine_is_better(const extfloat128_t a[2], const extfloat128_t& ra);
//static void report_mismatch_ex(extfloat128_t a);
int main(int argz, char** argv)
{
  if (argz < 2)
  {
    fprintf(stderr,
      "tst_atan2 - validate implementation of atan2Pi().\n"
      "Usage:\n"
      "tst_atan2 nIter [u | r | b]\n"
      "where\n"
      " u - test with arguments uniformly distributed in range (-1..1)\n"
      " r - test with results uniformly distributed in range (-1..1)\n"
      " b - test boost implementation with results uniformly distributed in range (-1..1)\n"
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

  std::mt19937_64 rndGen;
  std::uniform_int_distribution<uint64_t> rndDistr(0, uint64_t(-1));

  boost_octafloat_t x_o[2], yRef_o, yRes_o, yDiff_o, yUlp_o;
  extfloat128_t     x[2];
  int64_t n = 1;
  uint64_t mCnt[6] = {0};
  static const boost_octafloat_t oPi = boost::math::constants::pi<boost_octafloat_t>();
  static const boost_quadfloat_t qPi = boost::math::constants::pi<boost_quadfloat_t>();
  double maxErr_v = 0;
  extfloat128_t maxErr_x[2];
  maxErr_x[0]=maxErr_x[1]=0;
  while (nIter > 0) {
    make_random_quadfloat(x, rndGen, rndDistr, distrib);
    for (int i = 0; i < 2;++i) convert_to_boost_bin_float(&x_o[i], x[i]);
    yRef_o = atan2(x_o[0], x_o[1])/oPi;
    // std::cout << x_o << "=>" << yRef_o << "\n";
    extfloat128_t yRef; convert_from_boost_bin_float(&yRef, yRef_o);
    extfloat128_t yRes;
    if (!testBoost) {
      yRes = atan2Pi(x[0], x[1]);
    } else {
      boost_quadfloat_t yRes_q = atan2(boost_quadfloat_t(x_o[0]),boost_quadfloat_t(x_o[1]))/qPi;
      convert_from_boost_bin_float(&yRes, yRes_q);
    }
    if (!isEqual(yRef, yRes) && !mine_is_better(x, yRes)) {
      if (!isfinite(yRes) || !isfinite(yRef)) {
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
          maxErr_x[0] = x[0];
          maxErr_x[1] = x[1];
        }
      }
    }
    ++n;
    nIter -= 1;
  }

  printf("o.k. %I64u/%I64u/%I64u/%I64u/%I64u/%I64u mismatches. Max mismatch = %.6f ULP at y/x = ",
    mCnt[0],
    mCnt[1],
    mCnt[2],
    mCnt[3],
    mCnt[4],
    mCnt[5],
    maxErr_v);
  boost_quadfloat_t b[2];
  for (int i = 0; i < 2; ++i) convert_to_boost_bin_float(&b[i], maxErr_x[i]);
  std::cout << b[0] << "/" << b[1];
  if (maxErr_x[1] != 0) std::cout << " = " << b[0]/b[1];
  std::cout << "\n";
  report_mismatch(maxErr_x, -1);

  speed_test(rndGen, rndDistr);

  return 0;
}

static bool mine_is_better(const extfloat128_t a[2], const extfloat128_t& ra)
{
  if (a[0]==0 && a[1] < 0 && a[0].m_sign != 0 && ra==-1)
    return true; // patch, boost returns +1, my answer is -1. Mine is better
  if (isnan(a[1]) && a[0]==0 && isnan(ra))
    return true; // patch, boost returns 0, my answer is nan. Mine is better
  if (isnan(a[1]) && isinf(a[0]) && isnan(ra))
    return true; // patch, boost returns 0, my answer is nan. Mine is better
  if (a[0]==0 && a[1]==0 && isnan(ra))
    return true; // patch, boost returns 0, my answer is nan. Mine is better
  return false;
}

static bool report_mismatch(const extfloat128_t a[2], int64_t n)
{
  boost_quadfloat_t b[2];
  convert_to_boost_bin_float(&b[0], a[0]);
  convert_to_boost_bin_float(&b[1], a[1]);
  boost_octafloat_t ox(b[0]);
  boost_octafloat_t oy(b[1]);
  extfloat128_t     ra = atan2Pi(a[0], a[1]);
  boost_quadfloat_t rb = atan2(b[0],b[1])/boost::math::constants::pi<boost_quadfloat_t>();
  boost_octafloat_t ro = atan2(ox, oy)/boost::math::constants::pi<boost_octafloat_t>();

  // if (a[0]==0 && a[1] < 0 && a[0].m_sign != 0 && ra==-1)
    // return false; // patch, boost returns +1, my answer is -1. Mine is better
  // if (isnan(a[1]) && isnan(ra))
    // return false; // patch, boost returns 0, my answer is nan. Mine is better
  // if (a[0]==0 && a[1]==0 && isnan(ra))
    // return false; // patch, boost returns 0, my answer is nan. Mine is better

  if (n >= 0) printf("fail at iteration %I64d.\n", n);
  print(b[0]);
  print(b[1]);
  // print(a[0]);
  // print(a[1]);
  print(rb);
  print(ra);
  print(ro);
  convert_to_boost_bin_float(&rb, ra);
  print(rb);
  //printf("ix=%.8f\n", (ro*64).convert_to<double>());
  return true;
}

static void make_random_quadfloat(
  extfloat128_t* dst,
  std::mt19937_64& rndGen,
  std::uniform_int_distribution<uint64_t>& rndDistr,
  int  distrib)
{
  auto rndFunc = std::bind ( rndDistr, std::ref(rndGen) );
  uint64_t rndw[6];
  for (int i = 0; i < 6; ++i)
    rndw[i] = rndFunc();
  if (distrib == 'u') {
    for (int k = 0; k < 2; ++k) {
      dst[k]  = extfloat128_t::ldexp(extfloat128_t(rndw[k*3+0]),                -(64*1));
      dst[k] += extfloat128_t::ldexp(extfloat128_t(rndw[k*3+1]),                -(64*2));
      dst[k] += extfloat128_t::ldexp(extfloat128_t(rndw[k*3+2] & uint64_t(-2)), -(64*3));
      dst[k].m_sign = rndw[k*3+2] & 1;
    }
    return;
  } else if (distrib == 'r') {
    static const double SCALE48 = 1.1161179193627622078333570894698e-14; // pi/2^48
    static const double SCALE63 = SCALE48/(1 << 15);                     // pi/2^63
    int phis = (rndw[2] & 1) == 0 ? 1 : -1;
    if (rndw[2] > (1u << 20)) {
      double phi = int64_t(rndw[2] >> 16)*phis*SCALE48;
      // argument of sin/cos uniformly distributed in range [-pi..pi/2)
      dst[0] = sin(phi);
      dst[1] = cos(phi);
    } else {
      dst[0] = extfloat128_t(rndw[2]>>1)*(phis*SCALE63);
      dst[1] = (-fma(dst[0],dst[0], -extfloat128_t::one())).sqrt();
    }
    for (int k = 0; k < 2; ++k) {
      dst[k].m_significand[1] = (dst[k].m_significand[1] & (uint64_t(-1) << 32)) | (rndw[k] >> 32);
      dst[k].m_significand[0] = rndw[k] << 32;
    }
    extfloat128_t magn;
    magn.m_significand[0] = rndw[3];
    magn.m_significand[1] = rndw[4] | (uint64_t(1) << 63);
    magn.m_exponent = std::min(std::max(uint32_t(rndw[5]>>1), extfloat128_t::min_biased_exponent), extfloat128_t::max_biased_exponent-1);
    magn.m_sign = rndw[5] & 1;
    for (int k = 0; k < 2; ++k)
      dst[k] *= magn;
    return;
  }

  uint32_t sel_ex = uint32_t(rndw[2] >> 33);
  bool bothSpecial  = ((sel_ex & 255) == 0);
  bool firstSpecial = ((sel_ex & 63)  == 0);
  for (int k = 0; k < 2; ++k) {
    uint64_t msw = 0;
    uint64_t lsw = 0;
    uint64_t exw = rndw[k*2+2];
    uint32_t expLsw = uint32_t(exw);
    uint32_t expMsw = uint32_t(exw >> 32);
    int sign = expMsw & 1;
    uint32_t biased_exp = expLsw;
    unsigned expSel = (expMsw >> 1) & 63;
    if (expSel == 0 || bothSpecial) {
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
    } else {
      if (k == 1 && !firstSpecial && (expSel & 1) != 0) {
        // exponent close to buddy
        int32_t dExp = ((expSel & 2) == 0) ?
          (expLsw % 16) - 8 :   // very close
          (expLsw % 512) - 256; // not very close
        biased_exp = dst[0].m_exponent + dExp;
      }
      else {
        // full range of exponent
        ;
      }
      biased_exp = std::min(std::max(biased_exp, extfloat128_t::min_biased_exponent), extfloat128_t::max_biased_exponent);
      msw = rndw[k*2+0] | (uint64_t(1) << 63);
      lsw = rndw[k*2+1];
    }
    make_quadfloat(&dst[k], sign, biased_exp, msw, lsw);
  }
}

static void speed_test(std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr)
{
  const int VECLEN = 111;
  const int NITER  = 777;
  extfloat128_t      inpvec_a[VECLEN*2];
  boost_quadfloat_t  inpvec_b[VECLEN*2];
  extfloat128_t      outvec_a[VECLEN];
  boost_quadfloat_t  outvec_b[VECLEN];
  uint64_t dummy[3] = {0};
  uint64_t dt_a[NITER];
  uint64_t dt_b[NITER];
  static const boost_quadfloat_t bqfPi = boost::math::constants::pi<boost_quadfloat_t>();
  for (int it = 0; it < NITER; ++it) {
    for (int i = 0; i < VECLEN; ++i)
      make_random_quadfloat(&inpvec_a[i*2], std::ref(rndGen), rndDistr, 'r');

    uint64_t t0 = __rdtsc();
    for (int i = 0; i < VECLEN; ++i)
      outvec_a[i] = atan2Pi(inpvec_a[i*2+0],inpvec_a[i*2+1]);
    uint64_t t1 = __rdtsc();
    dt_a[it] = t1 - t0;
    for (int i = 0; i < VECLEN; ++i) {
      dummy[0] ^= outvec_a[i].m_significand[0];
      dummy[1] ^= outvec_a[i].m_significand[1];
      dummy[2] ^= (uint64_t(outvec_a[i].m_sign) << 32) | outvec_a[i].m_exponent;
    }

    for (int i = 0; i < VECLEN*2; ++i) {
      convert_to_boost_bin_float(&inpvec_b[i], inpvec_a[i]);
      inpvec_b[i] *= bqfPi;
    }

    t0 = __rdtsc();
    for (int i = 0; i < VECLEN; ++i)
      outvec_b[i] = atan2(inpvec_b[i*2+0],inpvec_b[i*2+1]);
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
