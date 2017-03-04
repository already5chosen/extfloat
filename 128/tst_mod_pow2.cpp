#include <stdio.h>
#include <stdint.h>
#include <cmath>
#include <random>
#include <functional>           // for std::bind
#include <iostream>
#include <algorithm>

#include "extfloat128.h"

static const double MIN_N_ITER  = 1;
static const double MAX_N_ITER  = 1e10;


static void make_quadfloat(extfloat128_t* dst, int sign, unsigned exp, uint64_t msw, uint64_t lsw) {
  dst->m_sign           = sign;
  dst->m_exponent      = exp;
  dst->m_significand[0] = lsw;
  dst->m_significand[1] = msw;
}

static bool isEqual(const extfloat128_t& a, const extfloat128_t& b) {
  if (a.m_sign           != b.m_sign)           return false;
  if (a.m_exponent      != b.m_exponent)      return false;
  if (a.m_exponent      != b.m_exponent)      return false;
  if (a.m_significand[1] != b.m_significand[1]) return false;
  if (a.m_significand[0] != b.m_significand[0]) return false;
  return true;
}

static void print(const extfloat128_t& a) {
  printf("%d %11u %016I64x:%016I64x {%d}%s\n", a.m_sign, a.m_exponent, a.m_significand[1], a.m_significand[0], a._get_exponent(), isnan(a) ? " nan" : "");
}

static void make_random_quadfloat(extfloat128_t* dst, std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr);
static uint32_t make_random_exp(std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr, uint32_t xExp);
static bool report_mismatch(extfloat128_t x, int32_t e, uint64_t n);

extfloat128_t ref_mod_pow2(const extfloat128_t& x, int32_t e) {
  extfloat128_t xa = x;
  xa.m_sign = 0;
  if (isinf(xa)) {
    xa.m_significand[1] = extfloat128_t::qnan_bit;
    return xa;
  }
  extfloat128_t xm = xa * extfloat128_t::pow2(-e);
  if (isinf(xm))
    xa = 0;
  else
    xa -= trunc(xm)*extfloat128_t::pow2(e);
  xa.m_sign = x.m_sign;
  return xa;
}

int main(int argz, char** argv)
{
  if (argz < 2)
  {
    fprintf(stderr,
      "tst_mod_pow2 - validate implementation of mod_pow2().\n"
      "Usage:\n"
      "tst_mod_pow2 nIter\n"
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
  uint64_t n = 1;
  while (nIter > 0) {
    extfloat128_t x;
    make_random_quadfloat(&x, rndGen, rndDistr);
    int32_t e = make_random_exp(rndGen, rndDistr, x.m_exponent);

    extfloat128_t ref = ref_mod_pow2(x, e);
    extfloat128_t res = mod_pow2(x, e);
    if (!isEqual(res, ref)) {
      if (report_mismatch(x, e, n))
        return 1;
    }
    ++n;
    nIter -= 1;
  }

  printf("o.k.\n");

  return 0;
}

static bool report_mismatch(extfloat128_t x, int32_t e, uint64_t n)
{
  extfloat128_t ref = ref_mod_pow2(x, e);
  extfloat128_t res = mod_pow2(x, e);

  printf("fail at iteration %I64u.\n", n);
  print(x);
  print(extfloat128_t::pow2(e));
  print(extfloat128_t::pow2(-e)*x);
  print(trunc(extfloat128_t::pow2(-e)*x));
  print(trunc(extfloat128_t::pow2(-e)*x)*extfloat128_t::pow2(e));
  print(ref);
  print(res);
  return true;
}

static void make_random_quadfloat(extfloat128_t* dst, std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr)
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
  if (expSel != 0) {
    // full range of exponent
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
  make_quadfloat(dst, sign, biased_exp, msw, lsw);
}

static uint32_t make_random_exp(std::mt19937_64& rndGen, std::uniform_int_distribution<uint64_t>& rndDistr, uint32_t xExp)
{
  auto rndFunc = std::bind ( rndDistr, std::ref(rndGen) );
  uint64_t exw = rndFunc();
  uint32_t expLsw = uint32_t(exw);
  uint32_t expMsw = uint32_t(exw >> 32);
  unsigned expSel = expMsw & 1;
  if (expSel)
    xExp += (expLsw % 512) - 256; // chose number not far from xExp
  else
    xExp = expLsw; // chose from full range
  return std::min(std::max(xExp, extfloat128_t::min_biased_exponent), extfloat128_t::max_biased_exponent);
}
