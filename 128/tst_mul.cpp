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
  //printf("%d %08x", a.backend().sign(), a.backend().exponent());
  const size_t alen = a.backend().bits().size()*sizeof(a.backend().bits().limbs()[0]);
  if (alen == 16) {
    const uint64_t* bits = reinterpret_cast<const uint64_t*>(a.backend().bits().limbs());
    printf(" %016I64x:%016I64x", bits[1], bits[0]);
  }
  else
    printf(" bits.size = %u", a.backend().bits().size());
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

void set_div_dbg(int x);

static void report_mismatch(extfloat128_t a, extfloat128_t b, uint64_t n);
int main(int argz, char** argv)
{
  if (argz < 2)
  {
    fprintf(stderr,
      "mul_test - validate implementation of multiplication.\n"
      "Usage:\n"
      "mul_test nIter\n"
      );
    return 1;
  }

  char* endp;
  double nIter = strtod(argv[1], &endp);
  if (endp == argv[1] || nIter < MIN_N_ITER || nIter >= MAX_N_ITER+1)
  {
    fprintf(stderr, "Illegal nIter argument '%s'. Please specify number in range [%.0f..%.0f]\n",
      argv[1], MIN_N_ITER, MAX_N_ITER);
    return 1;
  }

  if (argz >= 10) {
    boost_quadfloat_t ref[3];
    extfloat128_t       res[3];
    for (int i = 0; i < 2; ++i) {
      int s       = strtol  (argv[2+i*4], 0, 0) & 1;
      int e       = strtol  (argv[3+i*4], 0, 0);
      uint64_t m1 = strtoull(argv[4+i*4], 0, 16);
      uint64_t m0 = strtoull(argv[5+i*4], 0, 16);
      make_quadfloat      (&res[i], s, e, m1, m0);
    }
    // extfloat128_t uu;
    // make_quadfloat(&uu, 0, 0, uint64_t(1) << 63, 1);
    // res[0] *= uu;
    // set_div_dbg(1);
    res[2] = res[0] * res[1];
    // set_div_dbg(0);

    for (int i = 0; i < 2; ++i)
      convert_to_boost_bin_float(&ref[i], res[i]);
    ref[2] = ref[0] * ref[1];
    if (!isEqual(res[2], ref[2])) {
      report_mismatch(res[0], res[1], 0);
      return 1;
    }
  }

  std::mt19937_64 rndGen;
  std::uniform_int_distribution<uint64_t> rndDistr(0, uint64_t(-1));
  boost_quadfloat_t ref1, ref2, ref12;
  extfloat128_t     res1, res2, res12;
  uint64_t n = 1;
  while (nIter > 0) {
    make_random_quadfloat(&res1, rndGen, rndDistr);
    make_random_quadfloat(&res2, rndGen, rndDistr);
    convert_to_boost_bin_float(&ref1, res1);
    convert_to_boost_bin_float(&ref2, res2);

    ref12 = ref1 * ref2;
    res12 = res1 * res2;

    if (!isEqual(res12, ref12)) {
      report_mismatch(res1, res2, n);
      return 1;
    }

    ++n;
    nIter -= 1;
  }

  printf("o.k.\n");
  speed_test(rndGen, rndDistr);

  return 0;
}

static void report_mismatch(extfloat128_t a, extfloat128_t b, uint64_t n)
{
  // set_div_dbg(1);
  extfloat128_t ab = a * b;
  // set_div_dbg(0);
  boost_quadfloat_t ba, bb;
  convert_to_boost_bin_float(&ba, a);
  convert_to_boost_bin_float(&bb, b);
  boost_quadfloat_t bab = ba * bb;

  printf("fail at iteration %I64u\n", n);
  print(a);
  print(b);
  print(ba);
  print(bb);
  print(bab);
  print(ab);
  boost_octafloat_t oa = ba;
  boost_octafloat_t ob = bb;
  boost_octafloat_t oab = oa*ob;
  print(oab);
}

static void make_random_quadfloat(extfloat128_t* dst, std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr, bool saneExponent)
{
  auto rndFunc = std::bind ( rndDistr, std::ref(rndGen) );
  uint64_t msw = rndFunc() | (uint64_t(1) << 63);
  uint64_t lsw = rndFunc();
  uint64_t exw = rndFunc();
  uint32_t expLsw = uint32_t(exw);
  uint32_t expMsw = uint32_t(exw >> 32);
  uint32_t biased_exp = expLsw;
  int sign = expMsw & 1;
  int exp = static_cast<int>(expLsw);
  if (saneExponent) {
    exp >>= 2; // avoid multiplication overflow
    biased_exp = dst->exponent_bias + exp;
  } else {
    int expSel = (expMsw >> 1) & 63;
    if (expSel != 0) {
      // generate "normal" exponent
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
          msw  = dst->qnan_bit;
          sign = 0;
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
  extfloat128_t        inpvec_a[VECLEN+1];
  boost_quadfloat_t  inpvec_b[VECLEN+1];
  extfloat128_t        outvec_a[VECLEN];
  boost_quadfloat_t  outvec_b[VECLEN];
  uint64_t dummy[3] = {0};
  uint64_t dt_a[NITER];
  uint64_t dt_b[NITER];
  for (int it = 0; it < NITER; ++it) {
    for (int i = 0; i < VECLEN+1; ++i)
      make_random_quadfloat(&inpvec_a[i], std::ref(rndGen), rndDistr, true);
    uint64_t t0 = __rdtsc();
    for (int i = 0; i < VECLEN; ++i)
      outvec_a[i] = inpvec_a[i] * inpvec_a[i+1];
    uint64_t t1 = __rdtsc();
    dt_a[it] = t1 - t0;
    for (int i = 0; i < VECLEN; ++i) {
      dummy[0] ^= outvec_a[i].m_significand[0];
      dummy[1] ^= outvec_a[i].m_significand[1];
      dummy[2] ^= (uint64_t(outvec_a[i].m_sign) << 32) | outvec_a[i].m_exponent;
    }

    for (int i = 0; i < VECLEN+1; ++i)
      convert_to_boost_bin_float(&inpvec_b[i], inpvec_a[i]);
    t0 = __rdtsc();
    for (int i = 0; i < VECLEN; ++i)
      outvec_b[i] = inpvec_b[i] * inpvec_b[i+1];
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

  printf("%I64u %I64u %I64x\n", dummy[0], dummy[1], dummy[2]);
}