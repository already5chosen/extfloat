#include <stdio.h>
#include <stdint.h>
#include <cmath>
#include <random>
#include <functional>           // for std::bind
#include <iostream>
#include <intrin.h>
#include <boost/multiprecision/cpp_bin_float.hpp>

#include "extfloat128.h"
#include "extfloat128_to_bobinf.h"

typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<128, boost::multiprecision::backends::digit_base_2, void, boost::int64_t, INT_MIN, INT_MAX> > boost_quadfloat_t;
typedef boost::multiprecision::number<boost::multiprecision::backends::cpp_bin_float<320, boost::multiprecision::backends::digit_base_2, void, boost::int64_t, INT_MIN, INT_MAX>, boost::multiprecision::et_off> boost_decafloat_t;

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

static void print(const boost_decafloat_t& a) {
  printf("%d %11I64d", a.backend().sign(), a.backend().exponent());
  const size_t alen = a.backend().bits().size()*sizeof(a.backend().bits().limbs()[0]);
  if (alen == 40) {
    const uint64_t* bits = reinterpret_cast<const uint64_t*>(a.backend().bits().limbs());
    printf(" %016I64x:%016I64x:%016I64x:%016I64x:%016I64x", bits[4], bits[3], bits[2], bits[1], bits[0]);
  }
  else
    printf(" bits.size = %u", a.backend().bits().size());
  std::cout << " " << a << "\n";
}

static void print(const extfloat128_t& a) {
  printf("%d %11u %016I64x:%016I64x {%d}\n", a.m_sign, a.m_exponent, a.m_significand[1], a.m_significand[0], a._get_exponent());
}

void set_div_dbg(int x);

static bool report_mismatch(extfloat128_t a, uint64_t n);
int main(int argz, char** argv)
{
  if (argz < 2)
  {
    fprintf(stderr,
      "sqrt_test - validate implementation of sqrt().\n"
      "Usage:\n"
      "sqrt_test nIter\n"
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

  if (argz >= 6) {
    extfloat128_t       res;
    int s       = strtol  (argv[2], 0, 0);
    int e       = strtol  (argv[3], 0, 0);
    uint64_t m1 = strtoull(argv[4], 0, 16);
    uint64_t m0 = strtoull(argv[5], 0, 16);
    m1 |= uint64_t(1) << 63;
    make_quadfloat(&res, s & 1, e, m1, m0);

    // set_div_dbg(1);
    extfloat128_t rres = res.sqrt();
    // set_div_dbg(0);

    boost_quadfloat_t ref;
    convert_to_boost_bin_float(&ref, res);
    boost_quadfloat_t rref = sqrt(ref);
    if (!isEqual(rres, rref)) {
      if (report_mismatch(res, 0))
        return 1;
    }
    print(ref);
    print(rref);
    boost_decafloat_t oref(ref);
    boost_decafloat_t orref = sqrt(oref);
    print(orref);
    print(boost_quadfloat_t(orref));
    if (s > 1)
      return 0;
  }

  std::mt19937_64 rndGen;
  std::uniform_int_distribution<uint64_t> rndDistr(0, uint64_t(-1));

  boost_quadfloat_t ref;
  extfloat128_t       res;
  uint64_t n = 1;
  while (nIter > 0) {
    make_random_quadfloat(&res, rndGen, rndDistr);
    convert_to_boost_bin_float(&ref, res);

    boost_quadfloat_t rref = sqrt(ref);
    extfloat128_t       rres = res.sqrt();
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
  extfloat128_t       ra = a.sqrt();
  boost_quadfloat_t rb = sqrt(b);
  boost_decafloat_t o(b);
  boost_decafloat_t ro = sqrt(o);

  boost_quadfloat_t rab;
  convert_to_boost_bin_float(&rab, ra);

  // set_div_dbg(1);
  ra = a.sqrt();
  // set_div_dbg(0);

  printf("fail at iteration %I64u.\n", n);
  print(b);
  print(a);
  print(rb);
  print(ra);
  print(ro);
  // printf("b==0 %d, rab==0 %d, rb.m_sign() %d\n", b==0, rab==0, rb.backend().sign());
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
  uint32_t biased_exp = std::min(std::max(expLsw, extfloat128_t::min_biased_exponent), extfloat128_t::max_biased_exponent);
  if (saneExponent) {
    sign = 0;
  } else {
    unsigned expSel = (expMsw >> 1) & 63;
    if (expSel != 0) {
      if (expSel != 1)
        sign = 0;
      unsigned mskSz = (expMsw >> 7) & 63;
      uint64_t msk = uint64_t(-1) >> mskSz;
      if ((expSel & 3)== 2) {
        msw &= msk;
        lsw &= msk;
      } else if ((expSel & 3) == 3) {
        msw |= ~msk;
        lsw |= ~msk;
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
  }
  make_quadfloat(dst, sign, biased_exp, msw, lsw);
}

static void speed_test(std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr)
{
  const int VECLEN = 111;
  const int NITER  = 777;
  extfloat128_t        inpvec_a[VECLEN];
  boost_quadfloat_t  inpvec_b[VECLEN];
  extfloat128_t        outvec_a[VECLEN];
  boost_quadfloat_t  outvec_b[VECLEN];
  uint64_t dummy[3] = {0};
  uint64_t dt_a[NITER];
  uint64_t dt_b[NITER];
  for (int it = 0; it < NITER; ++it) {
    for (int i = 0; i < VECLEN; ++i)
      make_random_quadfloat(&inpvec_a[i], std::ref(rndGen), rndDistr, true);
    uint64_t t0 = __rdtsc();
    for (int i = 0; i < VECLEN; ++i)
      outvec_a[i] = inpvec_a[i].sqrt();
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
      outvec_b[i] = sqrt(inpvec_b[i]);
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