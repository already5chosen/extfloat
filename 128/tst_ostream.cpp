#include <stdio.h>
#include <stdint.h>
#include <cmath>
#include <random>
#include <functional>           // for std::bind
#include <iostream>
#include <boost/multiprecision/cpp_bin_float.hpp>

#define USE_EXTFLOAT128_T_OSTREAM
#include "extfloat128.h"
#include "extfloat128_to_bobinf.h"
#include "bobinf_to_extfloat128.h"

typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<128, boost::multiprecision::backends::digit_base_2, void, boost::int64_t, INT_MIN, INT_MAX> > boost_float128_t;
typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<512, boost::multiprecision::backends::digit_base_2, void, boost::int64_t, INT_MIN, INT_MAX> > boost_float512_t;

static const double MIN_N_ITER  = 1;
static const double MAX_N_ITER  = 1e10;

static int make_random_quadfloat(extfloat128_t* dst, std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr, int minExp, int maxExp);
static void speed_test(std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr, int minExp, int maxExp);

static void print(const boost_float128_t& a) {
  printf("%d %11I64d", a.backend().sign(), a.backend().exponent());
  const size_t alen = a.backend().bits().size()*sizeof(a.backend().bits().limbs()[0]);
  if (alen == 8*2) {
    const uint64_t* bits = reinterpret_cast<const uint64_t*>(a.backend().bits().limbs());
    printf(" %016I64x:%016I64x", bits[1], bits[0]);
  }
  else
    printf(" bits.size = %u*%u", a.backend().bits().size(), unsigned(sizeof(a.backend().bits().limbs()[0])));
  std::cout << " " << a << "\n";
}

#if 1
static void print(const boost_float512_t& a) {
  printf("%d %11I64d", a.backend().sign(), a.backend().exponent());
  const size_t alen = a.backend().bits().size()*sizeof(a.backend().bits().limbs()[0]);
  if (alen == 8*8) {
    const uint64_t* bits = reinterpret_cast<const uint64_t*>(a.backend().bits().limbs());
    printf(" %016I64x", bits[7]);
    for (int i = 1; i < 8; ++i)
    printf(":%016I64x", bits[7-i]);
  }
  else
    printf(" bits.size = %u*%u", a.backend().bits().size(), unsigned(sizeof(a.backend().bits().limbs()[0])));
  std::cout << " " << a << "\n";
}
#endif

static void print(const extfloat128_t& a) {
  printf("%d %11u %016I64x:%016I64x {%d}\n", a.m_sign, a.m_exponent, a.m_significand[1], a.m_significand[0], a._get_exponent());
}

static bool report_mismatch(extfloat128_t a, int64_t n, int prec);

int main(int argz, char** argv)
{
  if (argz < 4)
  {
    fprintf(stderr,
      "tst_ostream - validate implementation of ostream& operator<<().\n"
      "Usage:\n"
      "tst_ostream nIter expMin expMax\n"
      "where\n"
      " expMin - minimal value of exponent of test arguments. x = min_exponent_val\n"
      " expMax - maximal value of exponent of test arguments. x = max_exponent_val\n"
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

  int  minExp = extfloat128_t::min_exponent_val;
  int  maxExp = extfloat128_t::max_exponent_val;
  char* arg = argv[2];
  if (arg[0] != 'x') {
    char* endp;
    minExp = strtol(arg, &endp, 0);
    if (endp == arg) {
      fprintf(stderr, "Bad expMin argument '%s'. Not a number\n", arg);
      return 1;
    }
  }
  arg = argv[3];
  if (arg[0] != 'x') {
    char* endp;
    maxExp = strtol(arg, &endp, 0);
    if (endp == arg) {
      fprintf(stderr, "Bad expMax argument '%s'. Not a number\n", arg);
      return 1;
    }
  }

  if (minExp > maxExp) {
    fprintf(stderr, "Bad arguments: expMin > expMax\n");
    return 1;
  }

  printf("exp. range = [%d..%d]\n", minExp, maxExp);

  std::vector<extfloat128_t> specialPoints;
  specialPoints.push_back(extfloat128_t::nan());
  specialPoints.push_back(extfloat128_t(0));
  specialPoints.push_back(-extfloat128_t(0));
  specialPoints.push_back(extfloat128_t::inf());
  specialPoints.push_back(-extfloat128_t::inf());
  for (int exp = -256; exp < 256; ++exp) {
    boost_float512_t b = pow(boost_float512_t(10), exp);
    extfloat128_t x;
    convert_from_boost_bin_float(&x, b);
    specialPoints.push_back(x);
    specialPoints.push_back(nextinout(x, false));
    specialPoints.push_back(nextinout(x, true));
  }

  std::mt19937_64 rndGen;
  std::uniform_int_distribution<uint64_t> rndDistr(0, uint64_t(-1));

  boost_float128_t b;
  extfloat128_t    x;
  int prec;
  std::ostringstream res, ref;
  int64_t n = 0;
  int64_t nSpecialPoints = specialPoints.size();
  nIter += nSpecialPoints;
  res << std::scientific;
  ref << std::scientific;
  while (nIter > 0) {
    if (n < nSpecialPoints) {
      x = specialPoints[n];
      prec = 45;
    } else {
      prec = make_random_quadfloat(&x, rndGen, rndDistr, minExp, maxExp);
    }
    convert_to_boost_bin_float(&b, x);
    res.str(""); res.precision(prec); res << x;
    ref.str(""); ref.precision(prec); ref << b;
    if (res.str() != ref.str()) {
      if (report_mismatch(x, n, prec))
        return 1;
    }
    ++n;
    nIter -= 1;
  }

  printf("o.k.\n");
  speed_test(rndGen, rndDistr, minExp, maxExp);

  return 0;
}

extern int dbg;
static bool report_mismatch(extfloat128_t a, int64_t n, int prec)
{
  boost_float128_t b;
  convert_to_boost_bin_float(&b, a);
  std::ostringstream res, ref;

  res << std::scientific;
  ref << std::scientific;

  res.precision(prec);
  ref.precision(prec);

  res << a;
  ref << b;

  boost_float512_t x = b;
  boost_float512_t xb(ref.str());
  boost_float512_t xa = -1;
  try {
    xa = boost_float512_t(res.str());
    if (res.str().size()==ref.str().size()) {
      if (abs(x-xa) < abs(x-xb)) {
        return false; // mine is better
      }

      extfloat128_t ax;
      convert_from_boost_bin_float(&ax, xa);
      if (ax == a) {
      #if 1
        int eff_prec = std::min(64, prec);
        if (abs(x-xa) < abs(x)*pow(10,-eff_prec)*0.7)
          return false; // mine is good enough
      #endif
      }
    }
  } catch (...) {}
  dbg = 1;
  res.str("");
  res << a;
  dbg = 0;

  if (n >= 0) printf("fail at iteration %I64d.\n", n);
  std::cout << "prec " << prec << "\n";
  print(a);
  print(b);
  print(xa);
  print(xb);
  std::cout << "res " << res.str() << "\n";
  std::cout << "ref " << ref.str() << "\n";

  std::cout << std::scientific;
  std::cout.precision(prec);
  std::cout << "ext " << x << "\n";
  std::cout.precision(prec+4);
  std::cout << "EXT " << x << "\n";
  return true;
}

static int make_random_quadfloat(
  extfloat128_t* dst,
  std::mt19937_64& rndGen,
  std::uniform_int_distribution<uint64_t>& rndDistr,
  int minExp, int maxExp)
{
  auto rndFunc = std::bind ( rndDistr, std::ref(rndGen) );
  uint64_t msw = rndFunc();
  uint64_t lsw = rndFunc();
  uint64_t exw = rndFunc();

  extfloat128_t xx(int64_t(maxExp)+1 - minExp);
  xx *= exw & ((uint64_t(1) << 40)-1);
  int32_t e = minExp + int32_t(trunc(extfloat128_t::ldexp(xx,-40)).convert_to_double());

  int sign = (msw >> 63) & 1;
  msw |= uint64_t(1) << 63;
  *dst = extfloat128_t::construct(sign, e, msw, lsw);
  int precw = exw >> 40;
  return (precw & 15)==0 ? (precw >> 4) & 127 : (precw >> 4) & 63;
}

static void speed_test(std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr, int minExp, int maxExp)
{
  const int VECLEN = 1111;
  const int NITER  = 77;
  int precvec[VECLEN];
  extfloat128_t     inpvec_a[VECLEN];
  boost_float128_t  inpvec_b[VECLEN];
  std::ostringstream out_a;
  std::ostringstream out_b;
  uint64_t dummy = 0;
  uint64_t dt_a[NITER];
  uint64_t dt_b[NITER];
  for (int it = 0; it < NITER; ++it) {
    for (int i = 0; i < VECLEN; ++i)
      precvec[i] = make_random_quadfloat(&inpvec_a[i], std::ref(rndGen), rndDistr, minExp, maxExp);

    out_a.str("");
    out_a << std::scientific;
    uint64_t t0 = __rdtsc();
    for (int i = 0; i < VECLEN; ++i) {
      out_a.precision(precvec[i]);
      out_a << inpvec_a[i];
    }
    uint64_t t1 = __rdtsc();
    dt_a[it] = t1 - t0;

    for (int i = 0; i < VECLEN; ++i) {
      convert_to_boost_bin_float(&inpvec_b[i], inpvec_a[i]);
    }

    out_b.str("");
    out_b << std::scientific;
    t0 = __rdtsc();
    for (int i = 0; i < VECLEN; ++i) {
      out_b.precision(precvec[i]);
      out_b << inpvec_b[i];
    }
    t1 = __rdtsc();
    dt_b[it] = t1 - t0;

    dummy += out_a.str() != out_b.str();
  }
  std::nth_element(&dt_a[0], &dt_a[NITER/2], &dt_a[NITER]);
  std::nth_element(&dt_b[0], &dt_b[NITER/2], &dt_b[NITER]);
  printf("my %I64u. boost %I64u\n", dt_a[NITER/2]/VECLEN, dt_b[NITER/2]/VECLEN);

  printf("%I64x\n", dummy);
}
