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
typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<320, boost::multiprecision::backends::digit_base_2, void, boost::int64_t, INT_MIN, INT_MAX> > boost_decafloat_t;
typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<2048, boost::multiprecision::backends::digit_base_2, void, boost::int64_t, INT_MIN, INT_MAX> > boost_float2048b_t;

static const double MIN_N_ITER  = 1;
static const double MAX_N_ITER  = 1e10;

static void make_random_quadfloat(extfloat128_t* dst, std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr, int maxExp, bool uniformal = false);
static void speed_test(std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr, int maxExp);
static boost_float2048b_t find_nearest_above_NxPi(int exp, double nPlus);
static boost_float2048b_t find_nearest_below_NxPi(int exp, double nPlus);

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

static boost_octafloat_t to_octafloat(const extfloat128_t& x) {
  boost_octafloat_t b;
  convert_to_boost_bin_float(&b, x);
  return b;
}

static boost_octafloat_t ref_sin(extfloat128_t x) {
  static const boost_float2048b_t Pi    = boost::math::constants::pi<boost_float2048b_t>();
  static const boost_float2048b_t TwoPi = boost::math::constants::pi<boost_float2048b_t>()*2;
  int sign = x.m_sign;
  x.m_sign = 0;
  boost_float2048b_t xx;
  convert_to_boost_bin_float(&xx, x);
  xx = fmod(xx, TwoPi);
  if (xx >= Pi) {
    xx -= Pi;
    sign = !sign;
  }
  boost_decafloat_t xo(xx);
  xo = sin(xo);
  if (sign)
    xo = -xo;
  return static_cast<boost_octafloat_t>(xo);
}

static bool report_mismatch(extfloat128_t a, int64_t n, bool testBoost);
static double NormalizeDiff(const extfloat128_t& ref, const extfloat128_t& res, const boost_octafloat_t& oRef);

int main(int argz, char** argv)
{
  if (argz < 2)
  {
    fprintf(stderr,
      "tst_sin - validate implementation of sin().\n"
      "Usage:\n"
      "tst_sin nIter [u | b]  [number]\n"
      "where\n"
      " number - test arguments in range (-2.0**number..2.0**number). Default number = 1\n"
      " u - test with arguments uniformly distributed in range (-2.0**number..2.0**number)\n"
      " b - test boost implementation with arguments uniformly distributed in range (-2.0**number..2.0**number)\n"
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
  int  maxExp = 1;
  for (int arg_i = 2; arg_i < argz; ++arg_i) {
    char* arg =argv[arg_i];
    switch (arg[0]) {
      case 'u':
      case 'U':
        uniformal = true;
        break;
      case 'b':
      case 'B':
        uniformal = true;
        testBoost = true;
        break;
      default:
      {
        char* endp;
        int val = strtol(arg, &endp, 0);
        if (endp != arg && val <= 1024)
          maxExp = val;
      }
        break;
    }
  }
  printf("maxExp = %d.%s\n", maxExp, uniformal ? " Uniformal distribution." : "");

  std::vector<extfloat128_t> specialPoints;
  if (!uniformal) {
    for (int exp = 2; exp < maxExp; ++exp) {
      extfloat128_t x;
      convert_from_boost_bin_float(&x, find_nearest_above_NxPi(exp, 0.0));
      specialPoints.push_back(x);
      convert_from_boost_bin_float(&x, find_nearest_below_NxPi(exp, 0.0));
      specialPoints.push_back(x);
    }
  }

  std::mt19937_64 rndGen;
  std::uniform_int_distribution<uint64_t> rndDistr(0, uint64_t(-1));

  const double errThr = testBoost ? 64: 2;
  boost_octafloat_t yRef_o, yRes_o, yDiff_o, yUlp_o;
  extfloat128_t     x;
  int64_t n = 0;
  uint64_t mCnt[5] = {0};
  double maxErr_v = 0;
  extfloat128_t maxErr_x = 0;
  int64_t nSpecialPoints = specialPoints.size();
  nIter += nSpecialPoints;
  while (nIter > 0) {
    if (n < nSpecialPoints)
      x = specialPoints[n];
    else
      make_random_quadfloat(&x, rndGen, rndDistr, maxExp, uniformal);
    yRef_o = ref_sin(x);
    extfloat128_t yRef; convert_from_boost_bin_float(&yRef, yRef_o);
    extfloat128_t yRes;
    if (!testBoost) {
      yRes = sin(x);
    } else {
      boost_quadfloat_t x_q;
      convert_to_boost_bin_float(&x_q, x);
      boost_quadfloat_t yRes_q = sin(x_q);
      convert_from_boost_bin_float(&yRes, yRes_q);
    }
    if (yRef != yRes) {
      if (!isnan(yRes) || !isnan(yRef)) {
        ++mCnt[0];
        // if (!isfinite(yRes) || !isfinite(yRef)) {
          // report_mismatch(x, n);
          // return 1;
        // }
        double err = NormalizeDiff(yRef, yRes, yRef_o);
        if (err > errThr) {
          report_mismatch(x, n, testBoost);
          return 1;
        }
        if (err >= 9./16) {
          ++mCnt[1];
          if (err >= 5./8) {
            ++mCnt[2];
            if (err >= 3./4) {
              ++mCnt[3];
              if (err >= 1.) {
                ++mCnt[4];
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
    ++n;
    nIter -= 1;
  }

  printf("o.k. %I64u/%I64u/%I64u/%I64u/%I64u mismatches. Max mismatch = %.6f ULP at x = ",
    mCnt[0],
    mCnt[1],
    mCnt[2],
    mCnt[3],
    mCnt[4],
    maxErr_v);
  boost_quadfloat_t b;
  convert_to_boost_bin_float(&b, maxErr_x);
  std::cout << b << "\n";
  report_mismatch(maxErr_x, -1, testBoost);
  speed_test(rndGen, rndDistr, maxExp);

  return 0;
}

static double NormalizeDiff(const extfloat128_t& ref, const extfloat128_t& res, const boost_octafloat_t& oRef)
{
  return ((to_octafloat(res)-oRef)/to_octafloat(nextafter(ref, res) - ref)).convert_to<double>();
}
static bool report_mismatch(extfloat128_t a, int64_t n, bool testBoost)
{
  boost_quadfloat_t b;
  convert_to_boost_bin_float(&b, a);
  extfloat128_t     ra = sin(a);
  boost_octafloat_t o(b);
  boost_quadfloat_t rb = testBoost ? sin(boost_quadfloat_t(o)) : 0;
  boost_octafloat_t ro = ref_sin(a);

  if (n >= 0) printf("fail at iteration %I64d.\n", n);
  print(a);
  print(b);
  if (testBoost)
    print(rb);
  else
    print(ra);
  convert_to_boost_bin_float(&rb, ra);
  print(rb);
  print(ro);

  // boost_float2048b_t xx = abs(b);
  // xx = fmod(xx, boost::math::constants::pi<boost_float2048b_t>()*2);
  // xx /= boost::math::constants::pi<boost_float2048b_t>();
  // boost_octafloat_t xo(xx-1);
  // print(xo);

  // boost_float2048b_t xx = b;
  // xx = fmod(xx, boost::math::constants::pi<boost_float2048b_t>()*2);
  // boost_octafloat_t xo(xx);
  // print(xo);
  // xo -= boost::math::constants::pi<boost_octafloat_t>();
  // print(xo);
  // xx -= boost::math::constants::pi<boost_float2048b_t>();
  // xo = static_cast<boost_octafloat_t>(xx);
  // print(xo);

  // printf("---\n");
  // print(o/boost::math::constants::pi<boost_octafloat_t>());
  // extfloat128_t aa = fma(a, extfloat128_t::invPi(), a*extfloat128_t::invPi_lo());
  // print(aa);
  // boost_octafloat_t aa_o;
  // convert_to_boost_bin_float(&aa_o, aa);
  // print(sin(aa_o*boost::math::constants::pi<boost_octafloat_t>()));
  return true;
}

static void make_random_quadfloat(
  extfloat128_t*                           dst,
  std::mt19937_64&                         rndGen,
  std::uniform_int_distribution<uint64_t>& rndDistr,
  int                                      maxExp,
  bool                                     uniformal)
{
  auto rndFunc = std::bind ( rndDistr, std::ref(rndGen) );
  uint64_t msw = rndFunc();
  uint64_t lsw = rndFunc();
  uint64_t exw = rndFunc();
  if (uniformal) {
    *dst  = extfloat128_t::ldexp(extfloat128_t(msw),                maxExp-64*1);
    *dst += extfloat128_t::ldexp(extfloat128_t(lsw),                maxExp-64*2);
    *dst += extfloat128_t::ldexp(extfloat128_t(exw & uint64_t(-2)), maxExp-64*3);
    dst->m_sign = exw & 1;
    return;
  }

  msw |= (uint64_t(1) << 63);
  uint32_t expLsw = uint32_t(exw);
  uint32_t expMsw = uint32_t(exw >> 32);
  int sign = expMsw & 1;
  uint32_t biased_exp = 0;
  unsigned expSel = (expMsw >> 1) & 63;
  if (expSel != 0) {
    int32_t minExp = ((expSel & 1) != 0) ?
      -4 :                   // range [-4..maxExp)
      dst->min_exponent_val; // range [min_exponent_val..maxExp)
    biased_exp = extfloat128_t::exponent_bias + (((int64_t(maxExp)-minExp)*expLsw) >> 32) + minExp;
    biased_exp = std::min(std::max(biased_exp, extfloat128_t::min_biased_exponent), extfloat128_t::max_biased_exponent);
    msw |= (uint64_t(1) << 63);
  } else {
    // generate special value of exponent
    msw = lsw = 0;
    switch (expLsw % 3) {
      case 0:
        biased_exp = dst->zero_biased_exponent;    // zero
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
  if ((expSel & 3)==1 && biased_exp > dst->exponent_bias + 1) {
    // generate number in close proximity of N*pi/2
    *dst *= (dst->invPi()*2.0);
    static const boost_decafloat_t halfPi = boost::math::constants::pi<boost_decafloat_t>() * 0.5;
    boost_decafloat_t b;
    convert_to_boost_bin_float(&b, round(*dst));
    convert_from_boost_bin_float(dst, b*halfPi);
  }
}

static void speed_test(std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr, int maxExp)
{
  const int VECLEN = 1111;
  const int NITER  = 77;
  extfloat128_t      inpvec_a[VECLEN];
  boost_quadfloat_t  inpvec_b[VECLEN];
  extfloat128_t      outvec_a[VECLEN];
  boost_quadfloat_t  outvec_b[VECLEN];
  uint64_t dummy[3] = {0};
  uint64_t dt_a[NITER];
  uint64_t dt_b[NITER];
  for (int it = 0; it < NITER; ++it) {
    for (int i = 0; i < VECLEN; ++i)
      make_random_quadfloat(&inpvec_a[i], std::ref(rndGen), rndDistr, maxExp, true);

    uint64_t t0 = __rdtsc();
    for (int i = 0; i < VECLEN; ++i)
      outvec_a[i] = sin(inpvec_a[i]);
    uint64_t t1 = __rdtsc();
    dt_a[it] = t1 - t0;
    for (int i = 0; i < VECLEN; ++i) {
      dummy[0] ^= outvec_a[i].m_significand[0];
      dummy[1] ^= outvec_a[i].m_significand[1];
      dummy[2] ^= (uint64_t(outvec_a[i].m_sign) << 32) | outvec_a[i].m_exponent;
    }

    for (int i = 0; i < VECLEN; ++i) {
      convert_to_boost_bin_float(&inpvec_b[i], inpvec_a[i]);
    }

    t0 = __rdtsc();
    for (int i = 0; i < VECLEN; ++i)
      outvec_b[i] = inpvec_b[i].backend().exponent() < 128 ? sin(inpvec_b[i]) : 0;
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

  printf("%I64x %I64x %I64x\n", dummy[0], dummy[1], dummy[2]);
}

static const boost_float2048b_t Pi = boost::math::constants::pi<boost_float2048b_t>();
static boost_float2048b_t rem_ulp(const boost_float2048b_t& x, int exp)
{
  return x - ldexp(trunc(ldexp(x, 127-exp)), exp-127);
}

static boost_float2048b_t find_nearest_above_NxPi(int exp, double nPlus)
{
  boost_float2048b_t ulp = ldexp(boost_float2048b_t(1), exp-127);
  boost_float2048b_t n0  = ceil(ldexp(boost_float2048b_t(1), exp) / Pi - nPlus);
  boost_float2048b_t n2  = ceil(ldexp(boost_float2048b_t(2), exp) / Pi - nPlus);
  boost_float2048b_t n   = n0;
  boost_float2048b_t overI  = 1;
  boost_float2048b_t underI = 1;
  for (;;) {
    boost_float2048b_t dx = rem_ulp((n+nPlus)*Pi, exp);
    boost_float2048b_t underV;
    bool done = false;
    do {
      underV = ulp - rem_ulp(underI*Pi, exp);
      if (underV < dx)
        break;
      // update overI and underI
      boost_float2048b_t overV  = rem_ulp(overI*Pi, exp);
      if (underV < overV) {
        overI += underI*trunc(overV/underV);
        overV  = rem_ulp(overI*Pi, exp);
      }
      // here underV > overV
      boost_float2048b_t m = trunc(std::min(underV, underV-dx+overV)/overV);
      underI += overI*m;
      done = (n+underI >= n2);
    } while (!done);

    if (!done) {
      // (underV < dx)
      boost_float2048b_t dn = trunc(dx/underV)*underI;
      if (n+dn >= n2) {
        dn = trunc((n2-1-n)/underI)*underI;
        done = true;
      }
      n += dn;
    }
    if (done)
      break;
  }
  return (n+nPlus)*Pi;
}

static boost_float2048b_t find_nearest_below_NxPi(int exp, double nPlus)
{
  boost_float2048b_t ulp = ldexp(boost_float2048b_t(1), exp-127);
  boost_float2048b_t n0  = ceil(ldexp(boost_float2048b_t(1), exp) / Pi - nPlus);
  boost_float2048b_t n2  = ceil(ldexp(boost_float2048b_t(2), exp) / Pi - nPlus);
  boost_float2048b_t n   = n0;
  boost_float2048b_t overI  = 1;
  boost_float2048b_t underI = 1;
  for (;;) {
    boost_float2048b_t dx = ulp - rem_ulp((n+nPlus)*Pi, exp);
    boost_float2048b_t overV;
    bool done = false;
    do {
      overV = rem_ulp(overI*Pi, exp);
      if (overV < dx)
        break;
      // update overI and underI
      boost_float2048b_t underV = ulp - rem_ulp(underI*Pi, exp);
      if (overV < underV) {
        underI += overI*trunc(underV/overV);
        underV  = ulp - rem_ulp(underI*Pi, exp);
      }
      // here overV > underV
      boost_float2048b_t m = trunc(std::min(overV, overV-dx+underV)/underV);
      overI += underI*m;
      done = (n+overI >= n2);
    } while (!done);

    if (!done) {
      // (underV < dx)
      boost_float2048b_t dn = trunc(dx/overV)*overI;
      if (n+dn >= n2) {
        dn = trunc((n2-1-n)/overI)*overI;
        done = true;
      }
      n += dn;
    }
    if (done)
      break;
  }
  return (n+nPlus)*Pi;
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
  sin(x);

  boost_quadfloat_t b;
  convert_to_boost_bin_float(&b, dbg_xxx[0]); std::cout << "x0 = "; print(b);
  convert_to_boost_bin_float(&b, dbg_xxx[1]); std::cout << "xa = "; print(b);
  convert_to_boost_bin_float(&b, dbg_xxx[2]); std::cout << "xb = "; print(b);

  boost_octafloat_t x0_o;convert_to_boost_bin_float(&x0_o, dbg_xxx[0]);
  boost_octafloat_t xa_o;convert_to_boost_bin_float(&xa_o, dbg_xxx[1]);
  boost_octafloat_t xb_o;convert_to_boost_bin_float(&xb_o, dbg_xxx[2]);
  boost_octafloat_t y_r = sin(x0_o);
  boost_octafloat_t xb_r = x0_o - xa_o;
  std::cout << "xb err = " << xb_o/xb_r - 1 <<  "\n";

  convert_to_boost_bin_float(&b, dbg_xxx[3]); std::cout << "sin(xb) =\n"; print(b);
  boost_octafloat_t ySin_o;convert_to_boost_bin_float(&ySin_o, dbg_xxx[3]);
  boost_octafloat_t ySin_r = sin(xb_r);
  print(ySin_r);
  std::cout << "ySin err = " << (ySin_o-ySin_r)/ulp(ySin_r) << " : " << (ySin_o-ySin_r)/ulp(y_r) << "\n";

  convert_to_boost_bin_float(&b, dbg_xxx[4]); std::cout << "cosn(xb) =\n"; print(b);
  boost_octafloat_t yCosn_o;convert_to_boost_bin_float(&yCosn_o, dbg_xxx[4]);
  boost_decafloat_t yCosn_d = cos(boost_decafloat_t(xb_r)) - 1;
  boost_octafloat_t yCosn_r = boost_octafloat_t(yCosn_d);
  print(yCosn_r);
  std::cout << "yCosn err = " << (yCosn_o-yCosn_r)/ulp(yCosn_r) << " : " << (yCosn_o-yCosn_r)/ulp(y_r) << "\n";

  convert_to_boost_bin_float(&b, dbg_xxx[5]); std::cout << "tabSinA = "; print(b);
  convert_to_boost_bin_float(&b, dbg_xxx[6]); std::cout << "tabCosA = "; print(b);
  boost_octafloat_t tabSinA_o;convert_to_boost_bin_float(&tabSinA_o, dbg_xxx[5]);
  boost_octafloat_t tabCosA_o;convert_to_boost_bin_float(&tabCosA_o, dbg_xxx[6]);

  convert_to_boost_bin_float(&b, dbg_xxx[7]); std::cout << "tabSinB =\n"; print(b);
  boost_octafloat_t tabSinB_o;convert_to_boost_bin_float(&tabSinB_o, dbg_xxx[7]);
  boost_octafloat_t tabSinB_r = sin(xa_o) - tabSinA_o;
  print(tabSinB_r);
  std::cout << "tabSinB err = " << (tabSinB_o-tabSinB_r)/ulp(tabSinB_r) << " : " << (tabSinB_o-tabSinB_r)/ulp(y_r) << "\n";

  convert_to_boost_bin_float(&b, dbg_xxx[8]); std::cout << "tabCosB =\n"; print(b);
  boost_octafloat_t tabCosB_o;convert_to_boost_bin_float(&tabCosB_o, dbg_xxx[8]);
  boost_octafloat_t tabCosB_r = cos(xa_o) - tabCosA_o;
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

