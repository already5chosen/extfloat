#include <stdio.h>
#include <stdint.h>
#include <cmath>
#include <random>
#include <functional>           // for std::bind
#include <iostream>
#include <boost/multiprecision/cpp_bin_float.hpp>

#include "extfloat128.h"
#include "extfloat128_to_bobinf.h"

typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<128, boost::multiprecision::backends::digit_base_2, void, boost::int64_t, INT_MIN, INT_MAX> > boost_quadfloat_t;

static const double MIN_N_ITER  = 1;
static const double MAX_N_ITER  = 1e10;

static double make_random_double(std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr, bool saneExponent = false);
static void speed_test(std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr);

static bool isEqual(const extfloat128_t& a, const boost_quadfloat_t& b) {
  if (isnan(a))                                    return isnan(b);
  if ((a.m_sign != 0) != b.backend().sign())       return false;
  if (!isfinite(a))                                return isinf(b);
  if (!isnormal(a)) return b.backend().exponent_zero == b.backend().exponent();
  if (a._get_exponent() != b.backend().exponent()) return false;
  const uint64_t* bits = reinterpret_cast<const uint64_t*>(b.backend().bits().limbs());
  if (a.m_significand[1] != bits[1]) return false;
  if (a.m_significand[0] != bits[0]) return false;
  return true;
}

static void print(double a) {
  uint64_t ua;
  memcpy(&ua, &a, sizeof(ua));
  printf("%016I64x %.20e\n", ua, a);
}

static void print(const boost_quadfloat_t& a) {
  printf("%d %11I64d", a.backend().sign(), a.backend().exponent());
  const size_t alen = a.backend().bits().size()*sizeof(a.backend().bits().limbs()[0]);
  if (alen == 16) {
    const uint64_t* bits = reinterpret_cast<const uint64_t*>(a.backend().bits().limbs());
    printf(" %016I64x:%016I64x", bits[1], bits[0]);
  }
  else {
    printf(" bits.size = %u*%u", a.backend().bits().size(), unsigned(sizeof(a.backend().bits().limbs()[0])));
    const uint8_t* bits = reinterpret_cast<const uint8_t*>(a.backend().bits().limbs());
    for (size_t i = 0; i < alen; ++i)
      printf(" %02x", bits[i]);
  }
  std::cout << " " << a << "\n";
}

static void print(const extfloat128_t& a) {
  printf("%d %11u %016I64x:%016I64x {%d}\n", a.m_sign, a.m_exponent, a.m_significand[1], a.m_significand[0], a._get_exponent());
}

void set_div_dbg(int x);

static bool report_mismatch(double a, uint64_t n);
int main(int argz, char** argv)
{
  if (argz < 2)
  {
    fprintf(stderr,
      "cnv_from_test - validate implementation of from_double().\n"
      "Usage:\n"
      "cnv_from_test nIter\n"
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

  if (argz > 2) {
    double x = strtod(argv[2], 0);
    extfloat128_t     res = x;
    boost_quadfloat_t ref = x;
    uint64_t u;
    memcpy(&u, &x, sizeof(u));
    printf("%.20e %016I64x\n", x, u);
    print(res);
    print(ref);
    printf("%s\n", isEqual(res, ref) ? "o.k." : "fail");
    printf("%d\n", std::isnan(x));
    return 0;
  }

  std::mt19937_64 rndGen;
  std::uniform_int_distribution<uint64_t> rndDistr(0, uint64_t(-1));
  //auto rndFunc = std::bind ( rndDistr, rndGen );
  boost_quadfloat_t ref;
  extfloat128_t       res;
  uint64_t n = 1;
  while (nIter > 0) {
    double x = make_random_double(rndGen, rndDistr);

    ref = x;
    res = x;
    if (!isEqual(res, ref)) {
      if (report_mismatch(x, n))
        return 1;
    }
    ++n;
    nIter -= 1;
  }

  printf("o.k.\n");
  speed_test(rndGen, rndDistr);

  return 0;
}

static bool report_mismatch(double x, uint64_t n)
{
  extfloat128_t     a = x;
  boost_quadfloat_t b = x;

  // set_div_dbg(1);
  a = x;
  // set_div_dbg(0);

  printf("fail at iteration %I64u.\n", n);
  print(x);
  print(a);
  // boost_quadfloat_t ba;
  // convert_to_boost_bin_float(&ba, a);
  // print(ba);
  print(b);
  return true;
}

static double make_random_double(std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr, bool )
{
  auto rndFunc = std::bind ( rndDistr, std::ref(rndGen) );
  uint64_t u = rndFunc();
  double d;
  memcpy(&d, &u, sizeof(d));
  return d;
}

static void speed_test(std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr)
{
  const int VECLEN = 11000;
  const int NITER  = 777;
  double             inpvec[VECLEN];
  extfloat128_t        outvec_a[VECLEN];
  boost_quadfloat_t  outvec_b[VECLEN];
  uint64_t dummy[3] = {0};
  uint64_t dt_a[NITER];
  uint64_t dt_b[NITER];
  for (int it = 0; it < NITER; ++it) {
    for (int i = 0; i < VECLEN; ++i)
      inpvec[i] = make_random_double(std::ref(rndGen), rndDistr, true);
    uint64_t t0 = __rdtsc();
    for (int i = 0; i < VECLEN; ++i)
      outvec_a[i] = inpvec[i];
    uint64_t t1 = __rdtsc();
    dt_a[it] = t1 - t0;
    for (int i = 0; i < VECLEN; ++i) {
      dummy[0] ^= outvec_a[i].m_significand[0];
      dummy[1] ^= outvec_a[i].m_significand[1];
      dummy[2] ^= (uint64_t(outvec_a[i].m_sign) << 32) | outvec_a[i].m_exponent;
    }

    t0 = __rdtsc();
    for (int i = 0; i < VECLEN; ++i)
      outvec_b[i] = inpvec[i];
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