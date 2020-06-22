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
//static void report_mismatch_ex(extfloat128_t a);
int main(int argz, char** argv)
{
  if (argz < 2)
  {
    fprintf(stderr,
      "tst_sinPi - validate implementation of sinPi().\n"
      "Usage:\n"
      "tst_sinPi nIter [u | b]\n"
      "where\n"
      " u - test with arguments uniformly distributed in range (-0.5..0.5)\n"
      " b - test boost implementation with arguments uniformly distributed in range (-0.5..0.5)\n"
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
  bool testBoost = false;
  if (argz > 2) {
    if (argv[2][0] == 'u')
      uniformal = true;
    else if (argv[2][0] == 'b') {
      testBoost = true;
      uniformal = true;
    }
  }

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
    x_o *= oPi;
    yRef_o = sin(x_o);
    extfloat128_t yRef; convert_from_boost_bin_float(&yRef, yRef_o);
    extfloat128_t yRes;
    if (!testBoost) {
      yRes = sinPi(x);
    } else {
      boost_quadfloat_t yRes_q = sin(boost_quadfloat_t(x_o));
      convert_from_boost_bin_float(&yRes, yRes_q);
    }
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
        if (err >= 9./16) {
          ++mCnt[1];
          if (err >= 5./8) {
            ++mCnt[2];
            if (err >= 3./4) {
              ++mCnt[3];
              // double d = x.convert_to_double();
              // double d2 = fmod(fabs(d), 1.0);
              // if (d2 > 0.5) d2 = 1.0 - d2;
              // printf("%.6f : %.6e %.6f %.6f\n", err, d, d2, d2*64);
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

  printf("o.k. %I64u/%I64u/%I64u/%I64u mismatches. Max mismatch = %.6f ULP at x = ",
    mCnt[0],
    mCnt[1],
    mCnt[2],
    mCnt[3],
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
  extfloat128_t     ra = sinPi(a);
  boost_octafloat_t o(b);
  o *= boost::math::constants::pi<boost_octafloat_t>();
  boost_quadfloat_t rb = sin(boost_quadfloat_t(o));
  boost_octafloat_t ro = sin(o);

  if (n >= 0) printf("fail at iteration %I64d.\n", n);
  print(b);
  print(a);
  print(rb);
  print(ra);
  print(ro);
  convert_to_boost_bin_float(&rb, ra);
  print(rb);
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
      outvec_a[i] = sinPi(inpvec_a[i]);
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
      outvec_b[i] = sin(inpvec_b[i]);
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

#if 0
extern extfloat128_t dbg_xxx[16];
//dbg_xxx[1] = dx * (1.0/64);
//dbg_xxx[2] = x;
//dbg_xxx[3] = ySin;
//dbg_xxx[4] = yCosn;
//dbg_xxx[5] = tabSinA;
//dbg_xxx[6] = tabCosA;
//dbg_xxx[7] = tabSinB;
//dbg_xxx[8] = tabCosB;
//        extfloat128_t acc = (tabSinA + tabSinB)*yCosn;
//dbg_xxx[9] = acc;
//        acc = fma(tabCosB, ySin, acc);
//dbg_xxx[10] = acc;
//        acc = fma(tabCosA, ySin, acc);
//dbg_xxx[11] = acc;
//        acc += tabSinB;
//dbg_xxx[12] = acc;
//        ySin = acc + tabSinA;
//dbg_xxx[13] = ySin;
static boost_octafloat_t ulp(const boost_octafloat_t& x) {
  extfloat128_t y;
  convert_from_boost_bin_float(&y, x);
  boost_octafloat_t z;
  convert_to_boost_bin_float(&z, y.ulp());
  return z;
}
static void report_mismatch_ex(extfloat128_t a)
{
  extfloat128_t x = mod_pow2(a, 1); // reduce to range (-2..+2)
  x.m_sign = 0;
  if (x._get_exponent() >= 0)  x = 2.0 - x;
  if (x._get_exponent() >= -1) x = 1.0 - x;
  sinPi(x);

  boost_quadfloat_t b;
  convert_to_boost_bin_float(&b, dbg_xxx[0]); std::cout << "x0 = "; print(b);
  convert_to_boost_bin_float(&b, dbg_xxx[1]); std::cout << "xa = "; print(b);
  convert_to_boost_bin_float(&b, dbg_xxx[2]); std::cout << "xb = "; print(b);

  boost_octafloat_t x0_o;convert_to_boost_bin_float(&x0_o, dbg_xxx[0]);
  boost_octafloat_t xa_o;convert_to_boost_bin_float(&xa_o, dbg_xxx[1]);
  boost_octafloat_t xb_o;convert_to_boost_bin_float(&xb_o, dbg_xxx[2]);
  boost_octafloat_t y_r = sin(x0_o*boost::math::constants::pi<boost_octafloat_t>());
  boost_octafloat_t xb_r = x0_o - xa_o;
  std::cout << "xb err = " << xb_o/xb_r - 1 <<  "\n";

  convert_to_boost_bin_float(&b, dbg_xxx[3]); std::cout << "sin(xb) =\n"; print(b);
  boost_octafloat_t ySin_o;convert_to_boost_bin_float(&ySin_o, dbg_xxx[3]);
  boost_octafloat_t ySin_r = sin(xb_r*boost::math::constants::pi<boost_octafloat_t>());
  print(ySin_r);
  std::cout << "ySin err = " << (ySin_o-ySin_r)/ulp(ySin_r) << " : " << (ySin_o-ySin_r)/ulp(y_r) << "\n";

  convert_to_boost_bin_float(&b, dbg_xxx[4]); std::cout << "cosn(xb) =\n"; print(b);
  boost_octafloat_t yCosn_o;convert_to_boost_bin_float(&yCosn_o, dbg_xxx[4]);
  boost_decafloat_t yCosn_d = cos(boost_decafloat_t(xb_r)*boost::math::constants::pi<boost_decafloat_t>()) - 1;
  boost_octafloat_t yCosn_r = boost_octafloat_t(yCosn_d);
  print(yCosn_r);
  std::cout << "yCosn err = " << (yCosn_o-yCosn_r)/ulp(yCosn_r) << " : " << (yCosn_o-yCosn_r)/ulp(y_r) << "\n";

  convert_to_boost_bin_float(&b, dbg_xxx[5]); std::cout << "tabSinA = "; print(b);
  convert_to_boost_bin_float(&b, dbg_xxx[6]); std::cout << "tabCosA = "; print(b);
  boost_octafloat_t tabSinA_o;convert_to_boost_bin_float(&tabSinA_o, dbg_xxx[5]);
  boost_octafloat_t tabCosA_o;convert_to_boost_bin_float(&tabCosA_o, dbg_xxx[6]);

  convert_to_boost_bin_float(&b, dbg_xxx[7]); std::cout << "tabSinB =\n"; print(b);
  boost_octafloat_t tabSinB_o;convert_to_boost_bin_float(&tabSinB_o, dbg_xxx[7]);
  boost_octafloat_t tabSinB_r = sin(xa_o*boost::math::constants::pi<boost_octafloat_t>()) - tabSinA_o;
  print(tabSinB_r);
  std::cout << "tabSinB err = " << (tabSinB_o-tabSinB_r)/ulp(tabSinB_r) << " : " << (tabSinB_o-tabSinB_r)/ulp(y_r) << "\n";

  convert_to_boost_bin_float(&b, dbg_xxx[8]); std::cout << "tabCosB =\n"; print(b);
  boost_octafloat_t tabCosB_o;convert_to_boost_bin_float(&tabCosB_o, dbg_xxx[8]);
  boost_octafloat_t tabCosB_r = cos(xa_o*boost::math::constants::pi<boost_octafloat_t>()) - tabCosA_o;
  print(tabCosB_r);
  std::cout << "tabCosB err = " << (tabCosB_o-tabCosB_r)/ulp(tabCosB_r) << " : " << (tabCosB_o-tabCosB_r)/ulp(y_r) << "\n";

  boost_octafloat_t acc, accv[5];
  accv[0] = acc = (tabSinA_o + tabSinB_r)*yCosn_r;
  accv[1] = acc += tabCosB_r * ySin_r;
  accv[2] = acc += tabSinB_r;
  accv[3] = acc += tabCosA_o * ySin_r;
  accv[4] = acc += tabSinA_o;
  for (int i = 0; i < 5; ++i) {
    convert_to_boost_bin_float(&b, dbg_xxx[i+9]);
    std::cout << "acc" << i << " =\n";
    print(b);
    boost_octafloat_t acc_o;
    convert_to_boost_bin_float(&acc_o, dbg_xxx[i+9]);
    boost_octafloat_t acc_r = accv[i];
    print(acc_r);
    std::cout << "acc" << i << " err = " << (acc_o-acc_r)/ulp(acc_r) << " : " << (acc_o-acc_r)/ulp(y_r) << "\n";
  }
}
#endif

