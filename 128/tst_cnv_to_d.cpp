#include <stdio.h>
#include <stdint.h>
#include <cmath>
#include <random>
#include <functional>           // for std::bind
#include <iostream>
#include <boost/multiprecision/cpp_bin_float.hpp>

#include "extfloat128.h"
#include "extfloat128_to_bobinf.h"

typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<128, boost::multiprecision::backends::digit_base_2> > boost_quadfloat_t;
typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<196, boost::multiprecision::backends::digit_base_2> > boost_quadfloatx_t;
//typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<128, boost::multiprecision::backends::digit_base_2, void, int64_t> > boost_quadfloat_t;

double my_convert_to_double(const boost_quadfloat_t& x)
{
  if (isnormal(x)) {
    if (x.backend().exponent() < -1023+1) {
      // subnormal or zero
      double ret = 0;
      if (x.backend().exponent() >= -1023 - 52) {
        boost_quadfloatx_t tmp = abs(boost_quadfloatx_t(x)) + DBL_MIN;
        ret = tmp.convert_to<double>() - DBL_MIN;
      }
      return x.backend().sign() ? -ret : ret;
    }
  }
  return x.convert_to<double>();
}


static const double MIN_N_ITER  = 1;
static const double MAX_N_ITER  = 1e10;

static void make_random_quadfloat(extfloat128_t* dst, std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr, bool saneExponent = false);
static void speed_test(std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr);

static void make_quadfloat(extfloat128_t* dst, int sign, unsigned exp, uint64_t msw, uint64_t lsw) {
  dst->m_sign           = sign;
  dst->m_exponent       = exp;
  dst->m_significand[0] = lsw;
  dst->m_significand[1] = msw;
}

static bool isEqual(double a, double b) {
  uint64_t ua, ub;
  memcpy(&ua, &a, sizeof(ua));
  memcpy(&ub, &b, sizeof(ub));
  return ua == ub;
}

static void print(const boost_quadfloat_t& a) {
  printf("%d %11d", a.backend().sign(), a.backend().exponent());
  //printf("%d %11I64d", a.backend().sign(), a.backend().exponent());
  const size_t alen = a.backend().bits().size()*sizeof(a.backend().bits().limbs()[0]);
  if (alen == 16) {
    const uint64_t* bits = reinterpret_cast<const uint64_t*>(a.backend().bits().limbs());
    printf(" %016I64x:%016I64x", bits[1], bits[0]);
  }
  else
    printf(" bits.size = %u*%u", a.backend().bits().size(), unsigned(sizeof(a.backend().bits().limbs()[0])));
  std::cout << " " << a << "\n";
}

static void print(double a) {
  uint64_t ua;
  memcpy(&ua, &a, sizeof(ua));
  printf("%016I64x %.20e\n", ua, a);
}

static void print(const extfloat128_t& a) {
  printf("%d %11u %016I64x:%016I64x {%d}\n", a.m_sign, a.m_exponent, a.m_significand[1], a.m_significand[0], a._get_exponent());
}

static bool report_mismatch(extfloat128_t a, uint64_t n);
int main(int argz, char** argv)
{
  if (argz < 2)
  {
    fprintf(stderr,
      "tst_cnv_to_d - validate implementation of convert_to_double().\n"
      "Usage:\n"
      "tst_cnv_to_d nIter [s]\n"
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

    double rres = res.convert_to_double();

    boost_quadfloat_t ref;
    convert_to_boost_bin_float(&ref, res);
    double rref = ref.convert_to<double>();
    print(ref);
    print(rres);
    print(rref);
    if (!isEqual(rres, rref)) {
      if (report_mismatch(res, 0))
        return 1;
    }
    return 0;
  }

  std::mt19937_64 rndGen;
  std::uniform_int_distribution<uint64_t> rndDistr(0, uint64_t(-1));
  //auto rndFunc = std::bind ( rndDistr, rndGen );

  if (!subnormal) {
    for (int64_t n = 0;n < nIter; ++n) {
      extfloat128_t     res;
      make_random_quadfloat(&res, rndGen, rndDistr);
      boost_quadfloat_t ref;
      convert_to_boost_bin_float(&ref, res);

      // double rref = ref.convert_to<double>();
      double rref = my_convert_to_double(ref);
      double rres = res.convert_to_double();
      if (!isEqual(rres, rref)) {
        if (report_mismatch(res, n))
          return 1;
      }
    }
  } else {
    auto rndFunc = std::bind ( rndDistr, rndGen );
    // double mx = 0, mn = DBL_MAX;
    for (int64_t n = 0;n < nIter; ++n) {
      uint64_t msw = rndFunc();
      uint64_t lsw = rndFunc();
      uint64_t exw = rndFunc();
      uint32_t sign = exw & 1;
      uint32_t exp  = (exw & 2) ? (uint32_t(exw) >> 8) % 64 : 0;
      extfloat128_t res =
         extfloat128_t::ldexp(extfloat128_t(msw), -1022 -(64*1) - exp)
       + extfloat128_t::ldexp(extfloat128_t(lsw), -1022 -(64*2) - exp);
      res.m_sign = sign;
      if (exw & 4) {
        double a0 = res.convert_to_double();
        double a1 = nextafter(a0, DBL_MAX);
        res = (extfloat128_t(a0)+extfloat128_t(a1))*0.5;
        if (exw & 8)
          res.m_significand[0] |= lsw & (uint64_t(-1) >> (1022+62+res._get_exponent()));
        else if (exw & 16)
          res = nextinout(res, (exw & 32)!=0);
      }

      double rres = res.convert_to_double();

      extfloat128_t res2 = rres;
      if (res.m_sign != res2.m_sign) {
        boost_quadfloat_t b;
        convert_to_boost_bin_float(&b, res);
        printf("extfloat128_t fail at iter #%I64d\n", n);
        print(b);
        print(rres);
        print(boost_quadfloat_t(rres));
        return 1;
      }
      if (res != res2) {
        double nres = nextafter(rres, res2 < res ? DBL_MAX : -DBL_MAX);
        extfloat128_t dCurr = abs(rres-res);
        extfloat128_t dNext = abs(nres-res);
        if (dNext <= dCurr) {
          bool fail = true;
          if (dNext == dCurr) {
            uint64_t u;
            memcpy(&u, &rres, sizeof(u));
            fail = (u & 1) != 0;
          }
          if (fail) {
            boost_quadfloat_t b;
            convert_to_boost_bin_float(&b, res);
            printf("extfloat128_t fail at iter #%I64d\n", n);
            print(b);
            print(rres);
            print(boost_quadfloat_t(rres));
            print(nres);
            print(boost_quadfloat_t(nres));
            return 1;
          }
        }
      }
      // boost_quadfloat_t ref;
      // convert_to_boost_bin_float(&ref, res);
      // double rref = my_convert_to_double(ref);
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
        // double v = fabs(rres);
        // if (mx < v) mx = v;
        // if (mn > v) mn = v;
      // }
    }
    // printf("%e %e %e %e\n", mx, mn, DBL_MIN, nextafter(double(0), 1.0));
  }

  printf("o.k.\n");
  speed_test(rndGen, rndDistr);

  return 0;
}

static bool report_mismatch(extfloat128_t a, uint64_t n)
{
  boost_quadfloat_t b;
  convert_to_boost_bin_float(&b, a);
  double ra = a.convert_to_double();
  double rb = b.convert_to<double>();

  if (std::isfinite(ra) && std::isfinite(rb)) {
    if (ra == 0 && rb == 0 && a.m_sign) {
      if (std::signbit(ra))
        return false; // patch. boost returns wrong sign of zero - to be fixed
    } else {
      boost_quadfloat_t da = abs(b-ra);
      boost_quadfloat_t db = abs(b-rb);
      if (da <= db) {
        if (da < db)
          return false; // my convertion is better than booost's
        // da==db
        uint64_t ura, urb;
        memcpy(&ura, &ra, sizeof(ura));
        memcpy(&urb, &rb, sizeof(urb));
        if (((ura & 1)==0) && ((urb & 1) != 0))
          return false; // my convertion resolved tie to even - better than booost's which reslved to odd
      }
    }
  }

  ra = a.convert_to_double();

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
    exp = static_cast<int>(expLsw % 4096) - 2048;
    if      (exp < -1024-64) exp += 1024;
    else if (exp > 1024)     exp -= 1024;
    biased_exp = dst->exponent_bias + exp;
  } else {
    unsigned expSel = (expMsw >> 1) & 63;
    if (expSel != 0) {
      // generate "normal" exponent
      if ((expSel & 2)==2) {
        // choose exponent near double range
        exp = static_cast<int>(expLsw % 4096) - 2048;
        if (expSel == 2) {
          // choose mantissa in a way that mantissa of result has ~54 significant bits
          lsw &= 3;
          msw &= uint64_t(-1) << 10;
        }
      } else if ((expSel & 1)==1) {
        // choose exponent close to min/max
        if (exp < 0)
          exp = extfloat128_t::max_exponent_val - static_cast<int>(expLsw % 1024);
        else
          exp = extfloat128_t::min_exponent_val + static_cast<int>(expLsw % 1024);
      }
      // otherwise full range of exponent
      exp = std::min(std::max(exp, extfloat128_t::min_exponent_val), extfloat128_t::max_exponent_val);
      biased_exp = dst->exponent_bias + exp;
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
  extfloat128_t        inpvec_a[VECLEN];
  boost_quadfloat_t  inpvec_b[VECLEN];
  double             outvec_a[VECLEN];
  double             outvec_b[VECLEN];
  uint64_t dummy = 0;
  uint64_t dt_a[NITER];
  uint64_t dt_b[NITER];
  for (int it = 0; it < NITER; ++it) {
    for (int i = 0; i < VECLEN; ++i)
      make_random_quadfloat(&inpvec_a[i], std::ref(rndGen), rndDistr, true);
    uint64_t t0 = __rdtsc();
    for (int i = 0; i < VECLEN; ++i)
      outvec_a[i] = inpvec_a[i].convert_to_double();
    uint64_t t1 = __rdtsc();
    dt_a[it] = t1 - t0;
    for (int i = 0; i < VECLEN; ++i)
      dummy ^= reinterpret_cast<const uint64_t*>(outvec_a)[i];

    for (int i = 0; i < VECLEN; ++i)
      convert_to_boost_bin_float(&inpvec_b[i], inpvec_a[i]);
    t0 = __rdtsc();
    for (int i = 0; i < VECLEN; ++i)
      outvec_b[i] = inpvec_b[i].convert_to<double>();
    t1 = __rdtsc();
    dt_b[it] = t1 - t0;
    for (int i = 0; i < VECLEN; ++i)
      dummy ^= reinterpret_cast<const uint64_t*>(outvec_b)[i];
  }
  std::nth_element(&dt_a[0], &dt_a[NITER/2], &dt_a[NITER]);
  std::nth_element(&dt_b[0], &dt_b[NITER/2], &dt_b[NITER]);
  printf("my %I64u. boost %I64u\n", dt_a[NITER/2]/VECLEN, dt_b[NITER/2]/VECLEN);

  printf("%I64u\n", dummy);
}