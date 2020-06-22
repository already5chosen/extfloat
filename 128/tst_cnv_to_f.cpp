#include <stdio.h>
#include <stdint.h>
#include <cmath>
#include <random>
#include <functional>           // for std::bind
#include <intrin.h>
#include <iostream>
#include <boost/multiprecision/cpp_bin_float.hpp>

#include "extfloat128.h"
#include "extfloat128_to_bobinf.h"

typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<128, boost::multiprecision::backends::digit_base_2, void, boost::int64_t, INT_MIN, INT_MAX> > boost_quadfloat_t;

static const double MIN_N_ITER  = 1;
static const double MAX_N_ITER  = 1e10;

static void make_random_quadfloat(extfloat128_t* dst, std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr, bool saneExponent = false);
static void speed_test(std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr);


static void make_quadfloat(extfloat128_t* dst, int sign, unsigned exp, uint64_t msw, uint64_t lsw) {
  dst->m_sign           = sign;
  dst->m_exponent      = exp;
  dst->m_significand[0] = lsw;
  dst->m_significand[1] = msw;
}

static bool isEqual(float a, float b) {
  uint32_t ua, ub;
  memcpy(&ua, &a, sizeof(ua));
  memcpy(&ub, &b, sizeof(ub));
  return ua == ub;
}

static void print(const boost_quadfloat_t& a) {
  printf("%d %11I64d", a.backend().sign(), a.backend().exponent());
  // printf("%d %11d", a.backend().sign(), a.backend().exponent());
  const size_t alen = a.backend().bits().size()*sizeof(a.backend().bits().limbs()[0]);
  if (alen == 16) {
    const uint64_t* bits = reinterpret_cast<const uint64_t*>(a.backend().bits().limbs());
    printf(" %016I64x:%016I64x", bits[1], bits[0]);
  }
  else
    printf(" bits.size = %u*%u", a.backend().bits().size(), unsigned(sizeof(a.backend().bits().limbs()[0])));
  std::cout << " " << a << "\n";
}

static void print(float a) {
  uint32_t ua;
  memcpy(&ua, &a, sizeof(ua));
  printf("%08x %.10e\n", ua, double(a));
}

static void print(const extfloat128_t& a) {
  printf("%d %11u %016I64x:%016I64x {%d}\n", a.m_sign, a.m_exponent, a.m_significand[1], a.m_significand[0], a._get_exponent());
}

template <typename T, typename TDst>
void my_eval_convert_to(TDst* res, const T& x)
{
  TDst ret = 0;
  switch (eval_fpclassify(x)) {
    case FP_INFINITE:
      ret = std::numeric_limits<TDst>::infinity();
      break;

    case FP_NAN:
      *res = std::numeric_limits<TDst>::quiet_NaN();
      return;

    case FP_ZERO:
      break;

    default:
    { // == normal, because cpp_bin_float does not support subnormal
      if (x.exponent() >= std::numeric_limits<TDst>::min_exponent - std::numeric_limits<TDst>::digits - 1) {
        if (x.exponent() <= std::numeric_limits<TDst>::max_exponent-1) {
          int e = x.exponent();
          T y = x;
          y.exponent() = 63;
          y.sign()     = 0;

          uint64_t bits;
          eval_convert_to(&bits, y);

          int nhbits = std::numeric_limits<TDst>::digits;
          if (e < std::numeric_limits<TDst>::min_exponent-1)
            nhbits = e - std::numeric_limits<TDst>::min_exponent + 1 + std::numeric_limits<TDst>::digits;
          // round
          uint64_t lbits = bits << nhbits;
          uint64_t hbits = (bits >> 1) >> (63-nhbits); // shift by (64-nhbits) that works for nhbits in range [0..53]
          lbits |= (hbits & 1);                        // assure that tie rounded to even
          const uint64_t BIT63 = uint64_t(1) << 63;
          if (lbits == BIT63) {
            T truncY;
            truncY = bits;
            if (!eval_eq(truncY, y))
              lbits |= 1;
          }
          hbits += (lbits > BIT63);
          // ret = ldexp(double(int64_t(hbits)), e + 1 - nhbits);
          ret = std::ldexp(static_cast<TDst>(static_cast<int64_t>(hbits)), e + 1 - nhbits);
        } else {
          // overflow
          ret = std::numeric_limits<TDst>::infinity();
        }
      }
    }
      break;
  }
  *res = x.sign() ? -ret : ret;
}

template <typename T>
float inline my_convert_to_float(const T& x)
{
  float ret;
  #if 0
  my_eval_convert_to(&ret, x.backend());
  #else
  eval_convert_to(&ret, x.backend());
  #endif
  return ret;
}

namespace boost{ namespace multiprecision{ namespace backends{
double uuu_uuu(double x, int n) {
  // return std::ldexp(x, n);
  uint64_t u;
  memcpy(&u, &x, sizeof(u));
  u += (uint64_t(int64_t(n)) << 52);
  memcpy(&x, &u, sizeof(x));
  return x;
}
}}}


static bool report_mismatch(extfloat128_t a, uint64_t n);
static bool test(extfloat128_t a, int64_t n);

int main(int argz, char** argv)
{
  if (argz < 2)
  {
    fprintf(stderr,
      "tst_cnv_to_f - validate implementation of convert_to_float().\n"
      "Usage:\n"
      "tst_cnv_to_f nIter [s]\n"
      "where\n"
      " s - test with arguments uniformly distributed in subnormal range\n"
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

  bool subnormal = false;
  if (argz > 2)
    subnormal = argv[2][0]=='s';

  if (argz >= 6) {
    extfloat128_t       res;
    int s       = strtol  (argv[2], 0, 0) & 1;
    int e       = strtol  (argv[3], 0, 0);
    uint64_t m1 = strtoull(argv[4], 0, 16);
    uint64_t m0 = strtoull(argv[5], 0, 16);
    m1 |= uint64_t(1) << 63;
    make_quadfloat(&res, s, e, m1, m0);

    // set_div_dbg(1);
    float rres = res.convert_to_float();
    // set_div_dbg(0);

    boost_quadfloat_t ref;
    convert_to_boost_bin_float(&ref, res);
    float rref = ref.convert_to<float>();
    if (!isEqual(rres, rref)) {
      if (report_mismatch(res, 0))
        return 1;
    }
    print(ref);
    print(rref);
  }

  std::mt19937_64 rndGen;
  std::uniform_int_distribution<uint64_t> rndDistr(0, uint64_t(-1));
  //auto rndFunc = std::bind ( rndDistr, rndGen );
  if (!subnormal) {
    for (int64_t n = 0;n < nIter; ++n) {
      extfloat128_t     res;
      make_random_quadfloat(&res, rndGen, rndDistr);

      if (!test(res, n))
        return 1;

      // float rres = res.convert_to_float();

      // boost_quadfloat_t ref;
      // convert_to_boost_bin_float(&ref, res);
      // float rref = ref.convert_to<float>();

      // if (!isEqual(rres, rref)) {
        // if (report_mismatch(res, n))
          // return 1;
      // }
    }
  } else {
    auto rndFunc = std::bind ( rndDistr, rndGen );
    // float mx = 0.0f, mn = FLT_MAX;
    for (int64_t n = 0;n < nIter; ++n) {
      uint64_t msw = rndFunc();
      uint64_t lsw = rndFunc();
      uint64_t exw = rndFunc();
      uint32_t sign = exw & 1;
      uint32_t exp  = (exw & 2) ? (uint32_t(exw) >> 8) % 32 : 0;
      extfloat128_t res =
         extfloat128_t::ldexp(extfloat128_t(msw), -126 -(64*1) - exp)
       + extfloat128_t::ldexp(extfloat128_t(lsw), -126 -(64*2) - exp);
      res.m_sign = sign;
      if (exw & 4) {
        float a0 = res.convert_to_float();
        float a1 = nextafterf(a0, float(FLT_MAX));
        res = (extfloat128_t(a0)+extfloat128_t(a1))*0.5;
        if (exw & 8)
          res.m_significand[0] |= lsw & (uint64_t(-1) >> (39+126+23+res._get_exponent()));
        else if (exw & 16)
          res = nextinout(res, (exw & 32)!=0);
      }

      if (!test(res, n))
        return 1;

      // float rres = res.convert_to_float();

      // extfloat128_t res2 = rres;
      // if (res.m_sign != res2.m_sign) {
        // boost_quadfloat_t b;
        // convert_to_boost_bin_float(&b, res);
        // printf("extfloat128_t fail at iter #%I64d\n", n);
        // print(b);
        // print(rres);
        // print(boost_quadfloat_t(rres));
        // return 1;
      // }
      // if (res != res2) {
        // float nres = nextafterf(rres, float(res2 < res ? FLT_MAX : -FLT_MAX));
        // extfloat128_t dCurr = abs(rres-res);
        // extfloat128_t dNext = abs(nres-res);
        // if (dNext <= dCurr) {
          // bool fail = true;
          // if (dNext == dCurr) {
            // uint32_t u;
            // memcpy(&u, &rres, sizeof(u));
            // fail = (u & 1) != 0;
          // }
          // if (fail) {
            // boost_quadfloat_t b;
            // convert_to_boost_bin_float(&b, res);
            // printf("extfloat128_t fail at iter #%I64d\n", n);
            // print(b);
            // print(rres);
            // print(boost_quadfloat_t(rres));
            // print(nres);
            // print(boost_quadfloat_t(nres));
            // return 1;
          // }
        // }
      // }
      // if (n < 20) {
      // if (uu) {
        // float nres = nextafterf(rres, float(res2 < res ? FLT_MAX : -FLT_MAX));
        // boost_quadfloat_t b;
        // convert_to_boost_bin_float(&b, res);
        // printf("\nextfloat128_t o.k. at iter #%I64d\n", n);
        // print(b);
        // print(rres);
        // print(boost_quadfloat_t(rres));
        // print(nres);
        // print(boost_quadfloat_t(nres));
      // }
      // boost_quadfloat_t ref;
      // convert_to_boost_bin_float(&ref, res);
      // float rref = my_convert_to_float(ref);
      // if (!isEqual(rres, rref)) {
        // printf("boost fail at iter #%I64d\n", n);
        // print(ref);
        // print(rref);
        // print(boost_quadfloat_t(rref));
        // print(rres);
        // print(boost_quadfloat_t(rres));
        // return 1;
      // }

      // if (rres != 0) {
        // float v = fabs(rres);
        // if (mx < v) mx = v;
        // if (mn > v) mn = v;
      // }
    }
    // printf("%e %e %e %e\n", mx, mn, FLT_MIN, nextafterf(float(0), 1.0f));
  }

  printf("o.k.\n");
  speed_test(rndGen, rndDistr);

  return 0;
}

static void report_mismatch1(extfloat128_t a, uint64_t n)
{
  printf("extfloat128_t fail at iter #%I64d\n", n);
  boost_quadfloat_t b;
  convert_to_boost_bin_float(&b, a);
  print(b);
  float da = a.convert_to_float();
  print(da);
  print(boost_quadfloat_t(da));
  // printf("isinf(%f)=%d isnan(%f)=%d\n", da, isinf(da), da, isnan(da));
}

static bool test_x(extfloat128_t a, int64_t n, float ares)
{
  if (isnan(ares)) {
    if (isnan(a))
      return true;
    report_mismatch1(a, n);
    return false;
  }

  extfloat128_t a2 = ares;
  if (a.m_sign != a2.m_sign) {
    report_mismatch1(a, n);
    return false;
  }

  if (a2 == a)
    return true;

  if (isinf(ares)) {
    if (abs(a) > FLT_MAX)
      return true;
    report_mismatch1(a, n);
    return false;
  }

  if (a != a2) {
    float nres = nextafterf(ares, a2 < a ? float(FLT_MAX) : -float(FLT_MAX));
    extfloat128_t dCurr = abs(ares-a);
    extfloat128_t dNext = abs(nres-a);
    if (dNext <= dCurr) {
      bool fail = true;
      if (dNext == dCurr) {
        uint32_t u;
        memcpy(&u, &ares, sizeof(u));
        fail = (u & 1) != 0;
      }
      if (fail) {
        report_mismatch1(a, n);
        print(nres);
        print(boost_quadfloat_t(nres));
        return false;
      }
    }
  }
  return true;
}

static bool test(extfloat128_t a, int64_t n)
{
  float ares = a.convert_to_float();
  if (!test_x(a, n, ares))
    return false;

  boost_quadfloat_t b;
  convert_to_boost_bin_float(&b, a);
  float bres = b.convert_to<float>();
  // float bres = my_convert_to_float(b);

  if (!isEqual(ares, bres)) {
    printf("boost fail at iter #%I64d\n", n);
    print(a);
    print(b);
    print(ares);
    print(bres);
    print(boost_quadfloat_t(bres));
    return false;
  }

  return true;
}


static bool report_mismatch(extfloat128_t a, uint64_t n)
{
  boost_quadfloat_t b;
  convert_to_boost_bin_float(&b, a);
  float ra = a.convert_to_float();
  float rb = b.convert_to<float>();

  if (std::isfinite(ra) && std::isfinite(rb)) {
    if (ra == 0 && rb == 0 && a.m_sign) {
      if (std::signbit(ra))
        return false; // patch. boost returns wrong sign of zero - to be fixed
    } else {
      boost_quadfloat_t da = abs(b-ra);
      boost_quadfloat_t db = abs(b-rb);
      if (da < db)
        return false; // my convertion is better than booost's
    }
  }

  // set_div_dbg(1);
  ra = a.convert_to_float();
  // set_div_dbg(0);

  printf("fail at iteration %I64u.\n", n);
  print(a);
  print(b);
  print(ra);
  print(rb);
  print(boost_quadfloat_t(ra));
  print(boost_quadfloat_t(rb));
  return true;
}

static void make_random_quadfloat(extfloat128_t* dst, std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr, bool saneExponent)
{
  auto rndFunc = std::bind ( rndDistr, std::ref(rndGen) );
  uint64_t msw = rndFunc() | (uint64_t(1) << 63);
  uint64_t lsw = rndFunc();
  uint64_t exw = rndFunc();
  uint32_t expLsw = uint32_t(exw);
  uint32_t expMsw = uint32_t(exw >> 32);
  int sign = expMsw & 1;
  uint32_t biased_exp = expLsw;
  int exp = static_cast<int>(expLsw);
  if (saneExponent) {
    exp = static_cast<int>(expLsw % 256) - 127;
    biased_exp = dst->exponent_bias + exp;
  } else {
    unsigned expSel = (expMsw >> 1) & 63;
    if (expSel != 0) {
      // generate "normal" exponent
      if ((expSel & 2)==2) {
        // choose exponent near float range
        exp = static_cast<int>(expLsw % 512) - 256;
        if (expSel == 2) {
          // choose mantissa in a way that mantissa of result has ~25 significant bits
          lsw &= 3;
          msw &= uint64_t(-1) << 39;
        }
      } else if ((expSel & 1)==1) {
        // choose exponent close to min/max
        if (exp < 0)
          exp = extfloat128_t::max_exponent_val - static_cast<int>(expLsw % 128);
        else
          exp = extfloat128_t::min_exponent_val + static_cast<int>(expLsw % 128);
      }
      // otherwise full range of exponent
      exp = std::min(std::max(exp, extfloat128_t::min_exponent_val), extfloat128_t::max_exponent_val);
      biased_exp = dst->exponent_bias + exp;
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
          msw = dst->qnan_bit;
          break;
      }
    }
  }
  make_quadfloat(dst, sign, biased_exp, msw, lsw);
}

static void speed_test(std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr)
{
  const int VECLEN = 11000;
  const int NITER  = 777;
  extfloat128_t      inpvec_a[VECLEN];
  boost_quadfloat_t  inpvec_b[VECLEN];
  float              outvec_a[VECLEN];
  float              outvec_b[VECLEN];
  uint32_t dummy = 0;
  uint64_t dt_a[NITER];
  uint64_t dt_b[NITER];
  for (int it = 0; it < NITER; ++it) {
    for (int i = 0; i < VECLEN; ++i)
      make_random_quadfloat(&inpvec_a[i], std::ref(rndGen), rndDistr, true);
    uint64_t t0 = __rdtsc();
    for (int i = 0; i < VECLEN; ++i)
      outvec_a[i] = inpvec_a[i].convert_to_float();
    uint64_t t1 = __rdtsc();
    dt_a[it] = t1 - t0;
    for (int i = 0; i < VECLEN; ++i)
      dummy ^= reinterpret_cast<const uint32_t*>(outvec_a)[i];

    for (int i = 0; i < VECLEN; ++i)
      convert_to_boost_bin_float(&inpvec_b[i], inpvec_a[i]);
    t0 = __rdtsc();
    for (int i = 0; i < VECLEN; ++i)
      outvec_b[i] = inpvec_b[i].convert_to<float>();
      // outvec_b[i] = my_convert_to_float(inpvec_b[i]);
    t1 = __rdtsc();
    dt_b[it] = t1 - t0;
    for (int i = 0; i < VECLEN; ++i)
      dummy ^= reinterpret_cast<const uint32_t*>(outvec_b)[i];
  }
  std::nth_element(&dt_a[0], &dt_a[NITER/2], &dt_a[NITER]);
  std::nth_element(&dt_b[0], &dt_b[NITER/2], &dt_b[NITER]);
  printf("my %I64u. boost %I64u\n", dt_a[NITER/2]/VECLEN, dt_b[NITER/2]/VECLEN);

  printf("%x\n", dummy);
}