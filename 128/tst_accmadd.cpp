#include <stdio.h>
#include <stdint.h>
#include <cmath>
#include <random>
#include <functional>           // for std::bind
#include <iostream>
#include <boost/multiprecision/cpp_bin_float.hpp>

#include "extfloat128.h"
#include "extfloat128_to_bobinf.h"
#include "bobinf_to_extfloat128.h"

typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<384, boost::multiprecision::backends::digit_base_2, void, boost::int64_t> > boost_float384_t;

static const double MIN_N_ITER  = 1;
static const double MAX_N_ITER  = 1e10;

static void make_random_triade(extfloat128_t::acc_t* dstA, extfloat128_t dstB[2], std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr);

static void print(const boost_float384_t& a, const char* prefix = "") {
  printf("%s%d %11I64d", prefix, a.backend().sign(), a.backend().exponent());
  const size_t alen = a.backend().bits().size()*sizeof(a.backend().bits().limbs()[0]);
  if (alen == 48) {
    const uint64_t* bits = reinterpret_cast<const uint64_t*>(a.backend().bits().limbs());
    printf(" %016I64x:%016I64x:%016I64x:%016I64x:%016I64x:%016I64x", bits[5], bits[4], bits[3], bits[2], bits[1], bits[0]);
  }
  else
    printf(" bits.size = %u", a.backend().bits().size());
  std::cout << " " << a << "\n";
}

static void print(const extfloat128_t::acc_t& a, const char* prefix = "") {
  printf("%s%d %11I64d", prefix, a.m_sign, a.m_exponent-a.exponent_bias);
  const uint64_t* bits = a.m_significand;
  printf(" %016I64x:%016I64x:%016I64x:%016I64x:%016I64x:%016I64x %s\n"
    , bits[5], bits[4], bits[3], bits[2], bits[1], bits[0]
    , a.m_isNormal ? "n" : "");
}

#if 0
static void print(const extfloat128_t& a, const char* prefix = "") {
  printf("%s%d %11u %016I64x:%016I64x {%d}\n", prefix, a.m_sign, a.m_exponent, a.m_significand[1], a.m_significand[0], a._get_exponent());
}
#endif

void convert_to_boost_bin_float(boost_float384_t* pDst, const extfloat128_t::acc_t& src)
{
  if (src.m_exponent == 0) {
    *pDst = 0;
    return;
  }

  *pDst = src.m_significand[0];
  for (int i = 1; i < 5; ++i)
    *pDst =  ldexp(*pDst, -64) + src.m_significand[i];
  *pDst = ldexp(*pDst, -64) + int64_t(src.m_significand[5]);
  int64_t e = src.m_exponent - src.exponent_bias;
  *pDst = ldexp(*pDst, e+1);
  if (src.m_sign)
    *pDst = -*pDst;
}

static bool isLegal(const extfloat128_t::acc_t& a)
{
  if (a.m_exponent == 0)
    return true;
  int64_t y5 = a.m_significand[5];
  const int64_t SMALL_MNT_THR = int64_t(1) << (a.NBITS_MIN-64*4);
  const int64_t LARGE_MNT_THR = int64_t(1) << (a.NBITS_MAX-64*5);
  if (y5 < -LARGE_MNT_THR || y5 >= LARGE_MNT_THR)
    return false;
  int64_t y4 = a.m_significand[4];
  if (y5 == 0) {
    if (y4 >= 0) {
      if (y4 < SMALL_MNT_THR)
        return false;
      if (a.m_isNormal)
        return false;
    }
  } else if (y5 == -1) {
    if (y4 < 0) {
      if (y4 > -SMALL_MNT_THR) {
        printf("(%I64d %I64d)", y4, -SMALL_MNT_THR);
        return false;
      }
      if (a.m_isNormal)
        return false;
    }
  } else if (a.m_isNormal) {
    return false;
  }
  return true;
}

static bool isDiffMajor(const boost_float384_t& ref, const boost_float384_t& res, const boost_float384_t& a)
{
  boost_float384_t diff = res - ref;
  if (diff.backend().exponent() > ref.backend().exponent() - (extfloat128_t::acc_t::NBITS_MIN-0))
    return diff.backend().exponent() > a.backend().exponent() - (extfloat128_t::acc_t::NBITS_MIN-0);
  return false;
}

void set_div_dbg(int x);

static bool report_mismatch(extfloat128_t::acc_t A, const extfloat128_t B[2], uint64_t n, int op);
int main(int argz, char** argv)
{
  if (argz < 2)
  {
    fprintf(stderr,
      "tst_accmadd - validate implementation of acc_t::madd/msub.\n"
      "Usage:\n"
      "tst_accmadd nIter\n"
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

  boost_float384_t refA, refB[2], refR, resR;
  uint64_t n = 1;
  while (nIter > 0) {
    extfloat128_t::acc_t A;
    extfloat128_t        B[2];
    make_random_triade(&A, B, rndGen, rndDistr);

    convert_to_boost_bin_float(&refA, A);
    convert_to_boost_bin_float(&refB[0], B[0]);
    convert_to_boost_bin_float(&refB[1], B[1]);

    refR = refA + refB[0]*refB[1];
    extfloat128_t::acc_t R = A;
    R.madd(B[0], B[1]);
    if (!isLegal(R)) {
      if (report_mismatch(A, B, n, '+'))
        return 1;
    }
    convert_to_boost_bin_float(&resR, R);
    if (resR != refR && isDiffMajor(refR, resR, refA)) {
      if (report_mismatch(A, B, n, '+'))
        return 1;
    }

    refR = refA - refB[0]*refB[1];
    R = A;
    R.msub(B[0], B[1]);
    if (!isLegal(R)) {
      if (report_mismatch(A, B, n, '-'))
        return 1;
    }
    convert_to_boost_bin_float(&resR, R);
    if (resR != refR && isDiffMajor(refR, resR, refA)) {
      if (report_mismatch(A, B, n, '-'))
        return 1;
    }

    #if 0
    convert_from_boost_bin_float(&B, refR);
    if (isfinite(B)) {
      convert_to_boost_bin_float(&refB, B);
      refR = refA + refB;
      A = A0;
      A += B;
      if (!isLegal(A)) {
        if (report_mismatch(A0, B, n, '+'))
          return 1;
      }
      convert_to_boost_bin_float(&resR, A);
      if (resR != refR && isDiffMajor(refR, resR, refA)) {
        if (report_mismatch(A0, B, n, '+'))
          return 1;
      }

      refR = refA - refB;
      A = A0;
      A -= B;
      if (!isLegal(A)) {
        if (report_mismatch(A0, B, n, '+'))
          return 1;
      }
      convert_to_boost_bin_float(&resR, A);
      if (resR != refR && isDiffMajor(refR, resR, refA)) {
        if (report_mismatch(A0, B, n, '-'))
          return 1;
      }
    }
    #endif

    ++n;
    nIter -= 1;
  }

  printf("o.k.\n");
  // speed_test(rndGen, rndDistr, '+');
  // speed_test(rndGen, rndDistr, '-');

  return 0;
}

static bool report_mismatch(extfloat128_t::acc_t A, const extfloat128_t B[2], uint64_t n, int op)
{
  boost_float384_t refA, refB[2], refR, resR;
  convert_to_boost_bin_float(&refA, A);
  for (int k = 0; k < 2; ++k)
    convert_to_boost_bin_float(&refB[k], B[k]);
  refR = op=='-' ? refA - refB[0]*refB[1] : refA + refB[0]*refB[1] ;
  extfloat128_t::acc_t R = A;
  R.maddsub(B[0], B[1], op=='-');
  convert_to_boost_bin_float(&resR, R);

  printf("fail at iteration %I64u. Op = %c\n", n, op);
  print(A,       "a     ");
  print(refA,    "A     ");
  print(refB[0], "B0    ");
  print(refB[1], "B0    ");
  print(refB[0]*refB[1], "B0*B1 ");
  print(refR,  "refR  ");
  print(resR,  "R     ");
  print(R,     "r     ");
  return true;
}

static void make_random_triade(
  extfloat128_t::acc_t* dstA,
  extfloat128_t         dstB[2],
  std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr)
{
  auto rndFunc = std::bind ( rndDistr, std::ref(rndGen) );

  const int asl = sizeof(dstA->m_significand)/sizeof(dstA->m_significand[0]);
  for (int i = 0; i < asl; ++i)
    dstA->m_significand[i] = rndFunc();
  const int bsl = sizeof(dstB->m_significand)/sizeof(dstB->m_significand[0]);

  uint64_t bexp64 = rndFunc();
  uint64_t bExpSum = 0;
  const uint64_t BIT63 = (uint64_t(1) << 63);
  for (int k = 0; k < 2; ++k) {
    for (int i = 0; i < bsl; ++i)
      dstB[k].m_significand[i] = rndFunc();
    dstB[k].m_significand[1] |= BIT63;

    uint32_t bExp = bexp64;
    bExp = std::min(std::max(bExp, extfloat128_t::min_biased_exponent), extfloat128_t::max_biased_exponent);
    dstB[k].m_exponent = bExp;

    bExpSum += bExp;
    bexp64 >>= 32;
  }

  uint64_t aux = rndFunc();

  dstA  ->m_sign = (aux >> 0) & 1;
  dstB[0].m_sign = (aux >> 1) & 1;
  dstB[1].m_sign = (aux >> 2) & 1;

  int deltaE = int((aux >> 3) & 0x7FF) - 0x400; // -1024 to 1023;
  int aNbits = ((aux >> 14) & 127) + dstA->NBITS_MIN+1;
  if (aNbits > dstA->NBITS_MAX)
    aNbits = 320; // normalized case

  // replicate sign bits of significand with accordance to aNbits
  uint64_t hMask = uint64_t(-1);
  uint64_t lMask = 0;
  if (aNbits < 320)
    lMask = uint64_t(-1) << (aNbits-256);
  else
    hMask = uint64_t(-1) << (aNbits-320);
  if ((dstA->m_significand[asl-1] & BIT63) == 0) {
    dstA->m_significand[asl-2] |= (lMask >> 1); // MS bit = '1'
    dstA->m_significand[asl-2] &= ~lMask;
    dstA->m_significand[asl-1] &= ~hMask;
  } else {
    dstA->m_significand[asl-2] &= ~(lMask >> 1); // MS bit = '0'
    dstA->m_significand[asl-2] |= lMask;
    dstA->m_significand[asl-1] |= hMask;
  }

  dstA->m_exponent = (dstA->exponent_bias - dstB->exponent_bias*2 + bExpSum) + deltaE + 320 - aNbits;
  dstA->m_isNormal = false;
}
