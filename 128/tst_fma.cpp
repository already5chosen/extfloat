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

typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<128,   boost::multiprecision::backends::digit_base_2, void, boost::int64_t, INT_MIN, INT_MAX> > boost_quadfloat_t;
typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<128*2, boost::multiprecision::backends::digit_base_2, void, boost::int64_t, int64_t(INT_MIN)*3, int64_t(INT_MAX)*3> > boost_8xfloat_t;

static const double MIN_N_ITER  = 1;
static const double MAX_N_ITER  = 1e10;

static void make_random_quadfloat(extfloat128_t* dst, std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr, bool saneExponent = false, const extfloat128_t* pBuddy = 0);
static void speed_test(std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr);

static void print(const boost_8xfloat_t& a);
static boost_quadfloat_t ref_fma(const extfloat128_t& a, const extfloat128_t& b, const extfloat128_t& z)
{
  boost_8xfloat_t tmp[5];
  convert_to_boost_bin_float(&tmp[0], a);
  convert_to_boost_bin_float(&tmp[1], b);
  convert_to_boost_bin_float(&tmp[2], z);
  tmp[3] = tmp[0] * tmp[1];
  tmp[4] = tmp[3] + tmp[2];
  if (fpclassify(tmp[4])==FP_NORMAL) {
    const uint64_t* bits = reinterpret_cast<const uint64_t*>(tmp[4].backend().bits().limbs());
    if (bits[0]==0 && bits[1]==(uint64_t(1)<< 63)) {
      // tie, assure proper rounding
      if (tmp[3].backend().exponent() > tmp[2].backend().exponent()) {
        boost_8xfloat_t d1 = tmp[4] - tmp[3];
        boost_8xfloat_t r  = tmp[2] - d1;
        if (fpclassify(r)==FP_NORMAL) {
          // printf("tie0:");print(tmp[4]);
          // printf("tie1:");print(d1);
          // printf("tie2:");print(r);
          boost_8xfloat_t incr(r.backend().sign() ? -1 : 1);
          // printf("tie3:");print(incr);
          tmp[4] += ldexp(incr, tmp[4].backend().exponent()-255);
          // printf("tie4:");print(tmp[4]);
        }
      }
    }
  }
  return boost_quadfloat_t(tmp[4]);
}

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

static void print(const boost_8xfloat_t& a) {
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

static bool report_mismatch(extfloat128_t a, extfloat128_t b, extfloat128_t z, uint64_t n);
int main(int argz, char** argv)
{
  if (argz < 2)
  {
    fprintf(stderr,
      "add_test - validate implementation of add/sub.\n"
      "Usage:\n"
      "add_test nIter\n"
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

  if (argz >= 14) {
    boost_quadfloat_t ref[4];
    //boost_12xloat_t   xref[4];
    extfloat128_t     res[4];
    bool exitFlag = false;
    for (int i = 0; i < 3; ++i) {
      int s = strtol(argv[2+i*4], 0, 0);
      int e = strtol(argv[3+i*4], 0, 0);
      uint64_t m1 = strtoull(argv[4+i*4], 0, 16);
      uint64_t m0 = strtoull(argv[5+i*4], 0, 16);
      switch (argv[3+i*4][0]) {
        case 'n':
        case 'i':
          make_quadfloat(&res[i], s & 1, res[i].inf_nan_biased_exponent, m1, m0);
          break;
        case 'z':
          make_quadfloat(&res[i], s & 1, res[i].zero_biased_exponent, 0, 0);
          break;
        default:
          make_quadfloat(&res[i], s & 1, res[i].exponent_bias, m1 | (uint64_t(1) << 63), m0);
          res[i]._set_exponent(e);
          break;
      }
      if (s > 1)
        exitFlag = true;
    }
    res[3] = fma(res[0], res[1], res[2]);
    ref[3] = ref_fma(res[0], res[1], res[2]);

    // for (int i = 0; i < 3; ++i) {
      // convert_to_boost_bin_float(&ref[i],  res[i]);
      // convert_to_boost_bin_float(&xref[i], res[i]);
    // }

    if (!isEqual(res[3], ref[3])) {
      report_mismatch(res[0], res[1], res[2], 0);
      return 1;
    }
    if (exitFlag) {
      for (int i = 0; i < 4; ++i)
        print(res[i]);
      print(ref[3]);
      boost_8xfloat_t xref[4];
      for (int i = 0; i < 3; ++i)
        convert_to_boost_bin_float(&xref[i], res[i]);
      xref[3] = xref[0] * xref[1] + xref[2];
      print(xref[3]);
      return 0;
    }
  }

  std::mt19937_64 rndGen;
  std::uniform_int_distribution<uint64_t> rndDistr(0, uint64_t(-1));
  //auto rndFunc = std::bind ( rndDistr, rndGen );
  boost_quadfloat_t ref;
  extfloat128_t     res[4];
  uint64_t n = 1;
  while (nIter > 0) {
    make_random_quadfloat(&res[0], rndGen, rndDistr);
    make_random_quadfloat(&res[1], rndGen, rndDistr);
    extfloat128_t prod = res[0]*res[1];
    make_random_quadfloat(&res[2], rndGen, rndDistr, false, &prod);

    ref    = ref_fma(res[0], res[1], res[2]);
    res[3] = fma    (res[0], res[1], res[2]);
    if (!isEqual(res[3], ref)) {
      if (report_mismatch(res[0], res[1], res[2], n))
        return 1;
    }

    ref    = ref_fma(res[0], res[1], -res[2]);
    res[3] = fma    (res[0], res[1], -res[2]);
    if (!isEqual(res[3], ref)) {
      if (report_mismatch(res[0], res[1], -res[2], n))
        return 1;
    }

    ref    = ref_fma(-res[0], res[1], res[2]);
    res[3] = fma    (-res[0], res[1], res[2]);
    if (!isEqual(res[3], ref)) {
      if (report_mismatch(-res[0], res[1], res[2], n))
        return 1;
    }

    ++n;
    nIter -= 1;
  }

  printf("o.k.\n");
  speed_test(rndGen, rndDistr);

  return 0;
}

static bool report_mismatch(extfloat128_t a, extfloat128_t b, extfloat128_t z, uint64_t n)
{
  extfloat128_t     res = fma(a, b, z);
  boost_quadfloat_t ref = ref_fma(a, b, z);
  boost_8xfloat_t x[4];
  convert_to_boost_bin_float(&x[0], a);
  convert_to_boost_bin_float(&x[1], b);
  convert_to_boost_bin_float(&x[2], z);
  x[3] = x[0] * x[1] + x[2];

  if (ref == 0 && res == res.zero())
    return false; // ignore mismatch in sign of zeros

  printf("fail at iteration %I64u.\n", n);
  print(boost_quadfloat_t(x[0]));
  print(boost_quadfloat_t(x[1]));
  print(boost_quadfloat_t(x[2]));
  print(ref);
  print(res);
  print(x[3]);
  return true;
}

static void make_random_quadfloat(extfloat128_t* dst, std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr, bool saneExponent, const extfloat128_t* pBuddy)
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

  unsigned expSel = (expMsw >> 1) & 63;
  unsigned distE  = saneExponent ? (expMsw >> 7) % 9 : (expMsw >> 7) % 10 ;
  unsigned distM  = 1u << distE;
  int dist = int(expLsw % (distM*2)) - distM;
  if (expSel != 0 || saneExponent) {
    // generate "normal" exponent
    if (saneExponent)
      exp >>= 2; // choose any exponent in "sane" range
    if (pBuddy && isfinite(*pBuddy)) {
      if ((expSel & 2)==2) {
        // choose exponent close to buddy
        int buddyExp = pBuddy->m_exponent == pBuddy->zero_biased_exponent ?
          pBuddy->min_exponent_val : pBuddy->_get_exponent();
        exp = buddyExp + dist;
        if (expSel == 2 && dist > 0 && dist < 128) {
          // choose mantissa in a way that mantissa of result of subtraction is close to 0.5
          uint64_t lsw_b = pBuddy->m_significand[0];
          uint64_t msw_b = pBuddy->m_significand[1];
          msw = (uint64_t(1) << 63);
          lsw = lsw & 1;
          if (dist < 64) {
            msw ^= (msw_b >> dist);
            lsw ^= (lsw_b >> dist) | (msw_b << (64-dist));
          } else {
            lsw ^= (msw_b >> (dist - 64));
          }
        }
      } else if ((expSel & 1)==1) {
        // choose exponent close to min/max
        if (dist <= 0)
          exp = extfloat128_t::max_exponent_val + dist;
        else
          exp = extfloat128_t::min_exponent_val + dist - 1;
      }
      // otherwise full range of exponent
    } else {
      if ((expSel & 7)==3 && !pBuddy && !saneExponent) {
        // reduce number of ls bits in mantissa
        unsigned nBits  = (expMsw >> 7) % 128 + 1;
        if (nBits <= 64) {
          lsw = 0;
          msw &= (uint64_t(-1) << (64-nBits));
        } else {
          lsw &= (uint64_t(-1) << (128-nBits));
        }
      }
    }
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

  make_quadfloat(dst, sign, biased_exp, msw, lsw);
}

static void speed_test(std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr)
{
  const int VECLEN = 5500;
  const int NITER  = 777;
  extfloat128_t      inpvecm_a[VECLEN+1];
  boost_quadfloat_t  inpvecm_b[VECLEN+1];
  extfloat128_t      inpvecz_a[VECLEN];
  boost_quadfloat_t  inpvecz_b[VECLEN];
  extfloat128_t      outvec_a[VECLEN];
  boost_quadfloat_t  outvec_b[VECLEN];
  uint64_t dummy[3] = {0};
  uint64_t dt_a[NITER];
  uint64_t dt_b[NITER];
  for (int it = 0; it < NITER; ++it) {
    for (int i = 0; i < VECLEN+1; ++i)
      make_random_quadfloat(&inpvecm_a[i], std::ref(rndGen), rndDistr, true);
    for (int i = 0; i < VECLEN; ++i) {
      extfloat128_t prod = inpvecm_a[i] * inpvecm_a[i+1];
      make_random_quadfloat(&inpvecz_a[i], std::ref(rndGen), rndDistr, true, &prod);
    }
    uint64_t t0 = __rdtsc();
    for (int i = 0; i < VECLEN; ++i)
      outvec_a[i] = fma(inpvecm_a[i], inpvecm_a[i+1], inpvecz_a[i]);
    uint64_t t1 = __rdtsc();
    dt_a[it] = t1 - t0;
    for (int i = 0; i < VECLEN; ++i) {
      dummy[0] ^= outvec_a[i].m_significand[0];
      dummy[1] ^= outvec_a[i].m_significand[1];
      dummy[2] ^= (uint64_t(outvec_a[i].m_sign) << 32) | outvec_a[i].m_exponent;
    }

    for (int i = 0; i < VECLEN+1; ++i)
      convert_to_boost_bin_float(&inpvecm_b[i], inpvecm_a[i]);
    for (int i = 0; i < VECLEN; ++i)
      convert_to_boost_bin_float(&inpvecz_b[i], inpvecz_a[i]);
    t0 = __rdtsc();
    for (int i = 0; i < VECLEN; ++i) {
      // using boost::multiprecision::default_ops::eval_multiply_add;
      // outvec_b[i] = inpvecz_b[i];
      // outvec_b[i] += inpvecm_b[i] * inpvecm_b[i+1];
      // eval_multiply_add(outvec_b[i].backend(), inpvecm_b[i].backend(), inpvecm_b[i+1].backend());
      outvec_b[i] = inpvecm_b[i] * inpvecm_b[i+1] + inpvecz_b[i];
    }
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