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
typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<256, boost::multiprecision::backends::digit_base_2, void, boost::int64_t, INT_MIN, INT_MAX> > boost_octafloat_t;

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

static void print(const extfloat128_t& a) {
  printf("%d %11u %016I64x:%016I64x {%d}\n", a.m_sign, a.m_exponent, a.m_significand[1], a.m_significand[0], a._get_exponent());
}

static bool report_mismatch(extfloat128_t a, uint64_t n);
int main(int argz, char** argv)
{
  if (argz < 2)
  {
    fprintf(stderr,
      "tst_trunc - validate implementation of trunc().\n"
      "Usage:\n"
      "tst_trunc nIter\n"
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

  std::mt19937_64 rndGen;
  std::uniform_int_distribution<uint64_t> rndDistr(0, uint64_t(-1));
  boost_quadfloat_t ref;
  extfloat128_t       res;
  uint64_t n = 1;
  while (nIter > 0) {
    make_random_quadfloat(&res, rndGen, rndDistr);
    convert_to_boost_bin_float(&ref, res);

    boost_quadfloat_t rref = trunc(ref);
    extfloat128_t     rres = trunc(res);
    if (!isEqual(rres, rref)) {
      if (report_mismatch(res, n))
        return 1;
    }
    ++n;
    nIter -= 1;
  }

  printf("o.k.\n");
  speed_test(rndGen, rndDistr);

  return 0;
}

static bool report_mismatch(extfloat128_t a, uint64_t n)
{
  boost_quadfloat_t b;
  convert_to_boost_bin_float(&b, a);
  extfloat128_t     ra = trunc(a);
  boost_quadfloat_t rb = trunc(b);

  if (ra == 0 && rb == 0 && ra.m_sign != 0 && a.m_sign != 0)
    return false; // patch - boost returns wrong sign of negative zero

  printf("fail at iteration %I64u.\n", n);
  print(a);
  print(b);
  print(ra);
  print(rb);
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
    exp = static_cast<int>(expLsw % 256) - 100;
    biased_exp = dst->exponent_bias + exp;
  } else {
    unsigned expSel = (expMsw >> 1) & 63;
    if (expSel != 0) {
      // generate "normal" exponent
      if ((expSel & 2)==2) {
        // choose exponent near 128-bit integer range
        exp = static_cast<int>(expLsw % 256) - 100;
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
          //biased_exp = dst->inf_nan_biased_exponent; // inf
          biased_exp = dst->zero_biased_exponent; // zero
          break;
        default:
          //biased_exp = dst->inf_nan_biased_exponent; // NaN
          //msw = dst->qnan_bit;
          biased_exp = dst->zero_biased_exponent; // zero
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
  extfloat128_t        outvec_a[VECLEN];
  boost_quadfloat_t  outvec_b[VECLEN];
  uint64_t dummy = 0;
  uint64_t dt_a[NITER];
  uint64_t dt_b[NITER];
  for (int it = 0; it < NITER; ++it) {
    for (int i = 0; i < VECLEN; ++i)
      make_random_quadfloat(&inpvec_a[i], std::ref(rndGen), rndDistr, true);
    uint64_t t0 = __rdtsc();
    for (int i = 0; i < VECLEN; ++i)
      outvec_a[i] = trunc(inpvec_a[i]);
    uint64_t t1 = __rdtsc();
    dt_a[it] = t1 - t0;
    for (int i = 0; i < VECLEN; ++i)
      dummy ^= reinterpret_cast<const uint64_t*>(outvec_a)[i];

    for (int i = 0; i < VECLEN; ++i)
      convert_to_boost_bin_float(&inpvec_b[i], inpvec_a[i]);
    t0 = __rdtsc();
    for (int i = 0; i < VECLEN; ++i)
      outvec_b[i] = trunc(inpvec_b[i]);
    t1 = __rdtsc();
    dt_b[it] = t1 - t0;
    for (int i = 0; i < VECLEN; ++i)
      dummy ^= reinterpret_cast<const uint64_t*>(outvec_b)[i];
  }
  std::nth_element(&dt_a[0], &dt_a[NITER/2], &dt_a[NITER]);
  std::nth_element(&dt_b[0], &dt_b[NITER/2], &dt_b[NITER]);
  printf("my %I64u. boost %I64u\n", dt_a[NITER/2]/VECLEN, dt_b[NITER/2]/VECLEN);

  printf("%I64x\n", dummy);
}
