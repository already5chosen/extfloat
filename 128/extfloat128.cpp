#include <stdint.h>
#include <string.h>
#include <float.h>
#include <stdio.h>
#include <cmath>
#include <intrin.h>

#include "extfloat128.h"

#if 0
double extfloat128_t::convert_to_double() const
{
  //const uint64_t DBL_MNT_MSK = (uint64_t(1)<<52)-1;
  const int DBL_EXP_MSK = (1<<11)-1;
  const int DBL_EXP_MAX = DBL_EXP_MSK;
  const int DBL_EXP_BIAS = 1023;
  uint32_t exp = m_exponent;
  if (exp >= exponent_bias+1-DBL_EXP_BIAS-53 && exp <= exponent_bias+DBL_EXP_MAX-DBL_EXP_BIAS) {
    // normal
    int64_t lsbits = m_significand[0];
    uint64_t h = m_significand[1];
    lsbits |= h & 1;
    int64_t i = (h >> 1) | (lsbits  != 0);
    if (m_sign)
      i = - i;
    #if (FLT_RADIX==2)
      if (exp >= exponent_bias+1-DBL_EXP_BIAS && exp <= exponent_bias+DBL_EXP_MAX-DBL_EXP_BIAS-1) {
        // normal range
        double d = i;
        uint64_t u;
        memcpy(&u, &d, sizeof(u));
        u += uint64_t(int64_t(exp) - exponent_bias - 62) << 52;
        memcpy(&d, &u, sizeof(u));
        return d;
      } else {
        return scalbn(static_cast<double>(i), int32_t(exp - exponent_bias - 62));
      }
    #else
      return ldexp(static_cast<double>(i), int32_t(exp - exponent_bias - 62));
    #endif
  } else if (exp < exponent_bias) {
    // zero or very small number (underflow);
    return m_sign ? -0.0 : 0.0;
  } else if (exp <= max_biased_exponent) {
    // very big numbers
    return m_sign ? -HUGE_VAL : HUGE_VAL;
  } else {
    // NaN or Inf
    if (m_significand[1] == 0 && m_significand[0] == 0)
      return m_sign ? -HUGE_VAL : HUGE_VAL;
    else
      return ::nan("");
  }
}
#endif

double extfloat128_t::convert_to_double() const
{
  //const uint64_t DBL_MNT_MSK = (uint64_t(1)<<52)-1;
  const int DBL_EXP_MSK = (1<<11)-1;
  const int DBL_EXP_MAX = DBL_EXP_MSK;
  const int DBL_EXP_BIAS = 1023;
  uint32_t exp = m_exponent;
  if (exp >= exponent_bias+1-DBL_EXP_BIAS && exp < exponent_bias+DBL_EXP_MAX-DBL_EXP_BIAS) {
    // normal range
    uint64_t u = m_significand[1];
    uint32_t l = (u & ((1u<<11)-1)) | (m_significand[0] != 0);
    uint64_t h = u >> 11;
    l |= (h & 1);  // tie rounded to even
    h += uint64_t(exp - exponent_bias + DBL_EXP_BIAS - 1) << 52;
    h += (l > (1u << 10));
    h |= uint64_t(m_sign) << 63;
    double d;
    memcpy(&d, &h, sizeof(d));
    return d;
  } else {
    if (exp < exponent_bias+1-DBL_EXP_BIAS-53) {
      // zero or very small number (underflow);
      return m_sign ? -0.0 : 0.0;
    } else if (exp < exponent_bias+DBL_EXP_MAX-DBL_EXP_BIAS) {
      // subnormal
      static extfloat128_t dblMin = DBL_MIN;
      extfloat128_t tmp;
      tmp.m_significand[0] = 0;
      tmp.m_significand[1] = m_significand[1] | (m_significand[0] != 0);
      tmp.m_exponent = exp;
      tmp.m_sign = 0;
      double ret = (tmp + dblMin).convert_to_double() - DBL_MIN;
      return m_sign == 0 ? ret : -ret;
    } else if (exp <= max_biased_exponent) {
      // very big numbers
      return m_sign ? -HUGE_VAL : HUGE_VAL;
    } else {
      // NaN or Inf
      if (m_significand[1] == 0 && m_significand[0] == 0)
        return m_sign ? -HUGE_VAL : HUGE_VAL;
      else
        return ::nan("");
    }
  }
}
#if 0
float extfloat128_t::convert_to_float() const
{
  //const uint64_t DBL_MNT_MSK = (uint64_t(1)<<52)-1;
  const int FLT_EXP_MSK = (1<<8)-1;
  const int FLT_EXP_MAX = FLT_EXP_MSK;
  const int FLT_EXP_BIAS = 127;
  uint32_t exp = m_exponent;
  if (exp >= exponent_bias+1-FLT_EXP_BIAS-24 && exp <= exponent_bias+FLT_EXP_MAX-FLT_EXP_BIAS) {
    // normal
    int64_t lsbits = m_significand[0];
    uint64_t h = m_significand[1];
    lsbits |= h & (uint64_t(-1) >> 31); // combine word 0 with 33 LS bits of word 1
    int32_t i = static_cast<int32_t>(h >> 33) | (lsbits  != 0);
    if (m_sign)
      i = - i;
    #if (FLT_RADIX==2)
      if (exp >= exponent_bias+1-FLT_EXP_BIAS && exp <= exponent_bias+FLT_EXP_MAX-FLT_EXP_BIAS-1) {
        // normal range
        float f = static_cast<float>(i);
        uint32_t u;
        memcpy(&u, &f, sizeof(u));
        u += uint32_t(int32_t(exp) - exponent_bias - 30) << 23;
        memcpy(&f, &u, sizeof(f));
        return f;
      } else {
        //return static_cast<float>(convert_to_double()); // it is really better, but for now I'll remain compatible with boost
        return scalbnf(static_cast<float>(i), exp - exponent_bias - 30);
      }
    #else
      return ldexpf(static_cast<float>(i),  exp - exponent_bias - 30);
    #endif
  } else if (exp < exponent_bias) {
    // zero or very small number (underflow);
    return m_sign ? -0.0f : 0.0f;
  } else if (exp <= max_biased_exponent) {
    // very big numbers
    return m_sign ? -HUGE_VALF : HUGE_VALF;
  } else {
    // NaN or Inf
    if (m_significand[1] == 0 && m_significand[0] == 0)
      return m_sign ? -HUGE_VALF : HUGE_VALF;
    else
      return nanf("");
  }
}
#endif
float extfloat128_t::convert_to_float() const
{
  //const uint64_t DBL_MNT_MSK = (uint64_t(1)<<52)-1;
  const int FLT_EXP_MSK = (1<<8)-1;
  const int FLT_EXP_MAX = FLT_EXP_MSK;
  const int FLT_EXP_BIAS = 127;
  uint32_t exp = m_exponent;
  if (exp >= exponent_bias+1-FLT_EXP_BIAS-24 && exp <= exponent_bias+FLT_EXP_MAX-FLT_EXP_BIAS) {
    // normal or subnormal
    // if (exp >= exponent_bias+1-FLT_EXP_BIAS && exp <= exponent_bias+FLT_EXP_MAX-FLT_EXP_BIAS-2) {
      // // normal, no danger of rounding to inf
      // uint64_t hl = m_significand[1];
      // uint64_t lsbits = m_significand[0] | uint32_t(hl);
      // uint32_t h = uint32_t(hl >> 32) | (lsbits != 0);
      // uint32_t l = h & ((1<<8)-1);
      // h >>= 8;
      // l |= (h & 1); // tie rounded to even
      // h += uint32_t(exp - exponent_bias + FLT_EXP_BIAS - 1) << 23;
      // h += (l > (1<<7));
      // h |= uint32_t(m_sign) << 31;
      // float f;
      // memcpy(&f, &h, sizeof(f));
      // return f;
    // } else {
      // subnormal or too close to FLT_MAX
      const int DBL_EXP_BIAS = 1023;
      uint64_t lsbits = m_significand[0];
      uint64_t h = m_significand[1];
      lsbits |= h & ((1 << 11)-1); // combine word 0 with 11 LS bits of word 1
      h = (h >> 11) | (lsbits != 0);
      h += uint64_t(exp - exponent_bias + DBL_EXP_BIAS - 1) << 52;
      h |= uint64_t(m_sign) << 63;
      double d;
      memcpy(&d, &h, sizeof(h));
      return static_cast<float>(d);
    // }
  } else if (exp < exponent_bias) {
    // zero or very small number (underflow);
    return m_sign ? -0.0f : 0.0f;
  } else if (exp <= max_biased_exponent) {
    // very big numbers
    return m_sign ? -HUGE_VALF : HUGE_VALF;
  } else {
    // NaN or Inf
    if (m_significand[1] == 0 && m_significand[0] == 0)
      return m_sign ? -HUGE_VALF : HUGE_VALF;
    else
      return nanf("");
  }
}

void extfloat128_t::from_double(extfloat128_t* dst, double src) {
  uint64_t usrc;
  memcpy(&usrc, &src, sizeof(usrc));

  const unsigned DBL_EXP_MSK = (1<<11)-1;
  const uint64_t DBL_MNT_MSK = (uint64_t(1)<<52)-1;
  const unsigned DBL_EXP_MAX = DBL_EXP_MSK;
  const unsigned DBL_EXP_BIAS = 1023;
  const uint64_t BIT_63      = uint64_t(1)<<63;
  unsigned exp = (usrc >> 52) & DBL_EXP_MSK;
  uint64_t mnt = usrc & DBL_MNT_MSK;
  if (exp - 1 < DBL_EXP_MAX-1) {
    // normal
    dst->m_significand[0] = 0;
    dst->m_significand[1] = (mnt << 11) + BIT_63;
    dst->_set_exponent(static_cast<int32_t>(exp) - DBL_EXP_BIAS);
    dst->m_sign = static_cast<uint32_t>(usrc >> 63);
  } else if (exp == 0) {
    // subnormal or zero
    if (mnt == 0) {
      // zero
      dst->m_significand[0] = 0;
      dst->m_significand[1] = 0;
      dst->m_exponent = zero_biased_exponent;
      dst->m_sign = static_cast<uint32_t>(usrc >> 63);
    } else {
      // subnormal
      from_double(dst, static_cast<double>(mnt));
      dst->m_exponent -= DBL_EXP_BIAS+51;
      dst->m_sign = static_cast<uint32_t>(usrc >> 63);
    }
  } else {
    dst->m_exponent = inf_nan_biased_exponent;
    if (std::isnan(src)) {
      // NaN
      dst->m_significand[0] = 0;
      dst->m_significand[1] = qnan_bit;
      dst->m_sign = 0;
    } else {
      // inf
      dst->m_significand[0] = 0;
      dst->m_significand[1] = 0;
      dst->m_sign = static_cast<uint32_t>(usrc >> 63);
    }
  }
}

bool extfloat128_t::eval_eq(const extfloat128_t& srcA, const extfloat128_t& srcB)
{
  if (srcA.m_exponent       != srcB.m_exponent)       return false;
  if (srcA.m_significand[0] != srcB.m_significand[0]) return false;
  if (srcA.m_significand[1] != srcB.m_significand[1]) return false;
  if (isnan(srcA))                                    return false;
  if (srcA.m_sign           != srcB.m_sign)
    return srcA.m_exponent == zero_biased_exponent;
  return true;
}

bool extfloat128_t::eval_lg(const extfloat128_t& srcA, const extfloat128_t& srcB)
{
  if (isnan(srcA)) return false;
  if (isnan(srcB)) return false;
  if (srcA.m_exponent       != srcB.m_exponent)       return true;
  if (srcA.m_significand[0] != srcB.m_significand[0]) return true;
  if (srcA.m_significand[1] != srcB.m_significand[1]) return true;
  if (srcA.m_sign           != srcB.m_sign)
    return srcA.m_exponent != zero_biased_exponent;
  return false;
}

bool extfloat128_t::eval_uo(const extfloat128_t& srcA, const extfloat128_t& srcB)
{
  if (isnan(srcA)) return true;
  if (isnan(srcB)) return true;
  return false;
}

bool extfloat128_t::eval_lt(const extfloat128_t& srcA, const extfloat128_t& srcB)
{
  if (isnan(srcA)) return false;
  if (isnan(srcB)) return false;
  // both a and b comparable
  if (srcA.m_sign != srcB.m_sign) {
    if (srcA.m_sign == 0) return false;
    if (srcA.m_exponent != zero_biased_exponent) return true;
    if (srcB.m_exponent != zero_biased_exponent) return true;
    return false; // both a and b are zeros
  }
  // same sign
  if (srcA.m_exponent != srcB.m_exponent)
    return (srcA.m_exponent < srcB.m_exponent) ^ srcA.m_sign;
  // same exponent
  if (srcA.m_significand[1] != srcB.m_significand[1])
    return (srcA.m_significand[1] < srcB.m_significand[1]) ^ srcA.m_sign;
  if (srcA.m_significand[0] != srcB.m_significand[0])
    return (srcA.m_significand[0] < srcB.m_significand[0]) ^ srcA.m_sign;
  return false; // equal
}

bool extfloat128_t::eval_le(const extfloat128_t& srcA, const extfloat128_t& srcB)
{
  if (isnan(srcA)) return false;
  if (isnan(srcB)) return false;
  // both a and b comparable
  if (srcA.m_exponent != srcB.m_exponent) {
    // different exponents
    if (srcA.m_sign != srcB.m_sign) return srcA.m_sign;
    // same sign
    return (srcA.m_exponent < srcB.m_exponent) ^ srcA.m_sign;
  }
  // same exponent
  if (srcA.m_exponent == zero_biased_exponent)
    return true; // both a and b are zeros
  // both non-zeros
  if (srcA.m_sign != srcB.m_sign) return srcA.m_sign;
  // both non-zeros with the same sign and the same exponent
  if (srcA.m_significand[1] != srcB.m_significand[1])
    return (srcA.m_significand[1] < srcB.m_significand[1]) ^ srcA.m_sign;
  if (srcA.m_significand[0] != srcB.m_significand[0])
    return (srcA.m_significand[0] < srcB.m_significand[0]) ^ srcA.m_sign;
  return true; // equal
}

// return
//   -1 for unordered,
//    0 for srcA == srcB,
//    1 for srcA <  srcB,
//    2 for srcA >  srcB
int extfloat128_t::compare_ex(const extfloat128_t& srcA, const extfloat128_t& srcB)
{
  if (isnan(srcA)) return -1;
  if (isnan(srcB)) return -1;
  // both a and b comparable
  if (srcA.m_exponent != srcB.m_exponent) {
    // different exponents
    if (srcA.m_sign != srcB.m_sign) return srcB.m_sign + 1;
    // same sign
    return ((srcA.m_exponent > srcB.m_exponent) ^ srcA.m_sign) + 1;
  }
  // same exponent
  if (srcA.m_exponent == zero_biased_exponent)
    return 0; // both a and b are zeros
  // both non-zeros
  if (srcA.m_sign != srcB.m_sign) return srcB.m_sign + 1;
  // both non-zeros with the same sign and the same exponent
  if (srcA.m_significand[1] != srcB.m_significand[1])
    return ((srcA.m_significand[1] > srcB.m_significand[1]) ^ srcA.m_sign) + 1;
  if (srcA.m_significand[0] != srcB.m_significand[0])
    return ((srcA.m_significand[0] > srcB.m_significand[0]) ^ srcA.m_sign) + 1;
  return 0; // equal
}

namespace {

enum {
  NFO_CLASS_ZERO = 0,
  NFO_CLASS_INF  = 2,
  NFO_CLASS_NAN  = 3,
  NFO_SIGN_POS   = 0 << 2,
  NFO_SIGN_NEG   = 1 << 2,
  NFO_SIGN_A     = 2 << 2,
  NFO_SIGN_B     = 3 << 2,
  NFO_SIGN_MSK   = 3 << 2,
};

}

void extfloat128_t::eval_non_finite(
  extfloat128_t& dst, const extfloat128_t& srcA, const extfloat128_t& srcB, int issub,
  const uint8_t tab[2][2][3])
{
  int classA = 0;
  switch (srcA.m_exponent)
  {
    case zero_biased_exponent:
      classA = 1;
      break;
    case inf_nan_biased_exponent:
      if ((srcA.m_significand[0] | srcA.m_significand[1]) != 0) {
        dst = srcA;
        return;
      }
      classA = 2;
      break;
    default:
      break;
  }
  int classB = 0;
  switch (srcB.m_exponent)
  {
    case zero_biased_exponent:
      classB = 1;
      break;
    case inf_nan_biased_exponent:
      if ((srcB.m_significand[0] | srcB.m_significand[1]) != 0) {
        dst = srcB;
        return;
      }
      classB = 2;
      break;
    default:
      break;
  }
  int signA = srcA.m_sign;
  int signB = srcB.m_sign ^ issub;
  int signX = signA ^ signB;
  int r = (classA == 2) ? tab[0][signX][classB] : tab[1][signX][classA];
  dst.m_significand[0] = 0;
  dst.m_significand[1] = (r & 1) ? qnan_bit : 0;
  dst.m_exponent      = (r & 2) ? inf_nan_biased_exponent : zero_biased_exponent;
  int sign;
  switch (r & NFO_SIGN_MSK)
  {
    case NFO_SIGN_POS: sign = 0;     break;
    case NFO_SIGN_NEG: sign = 1;     break;
    case NFO_SIGN_A:   sign = signA; break;
    case NFO_SIGN_B:
    default:           sign = signB; break;
  }
  dst.m_sign = sign;
}

static const uint8_t multiply_nfo_tab[2][2][3] = {
//  NRM                         ZERO                        INF
 {{ NFO_CLASS_INF+NFO_SIGN_POS, NFO_CLASS_NAN+NFO_SIGN_POS, NFO_CLASS_INF+NFO_SIGN_POS  },  // A=inf, signA==signB
  { NFO_CLASS_INF+NFO_SIGN_NEG, NFO_CLASS_NAN+NFO_SIGN_POS, NFO_CLASS_INF+NFO_SIGN_NEG  }}, // A=inf, signA!=signB
 {{ NFO_CLASS_INF+NFO_SIGN_POS, NFO_CLASS_NAN+NFO_SIGN_POS, NFO_CLASS_INF+NFO_SIGN_POS  },  // B=inf, signA==signB
  { NFO_CLASS_INF+NFO_SIGN_NEG, NFO_CLASS_NAN+NFO_SIGN_POS, NFO_CLASS_INF+NFO_SIGN_NEG  }}, // B=inf, signA!=signB
};
void extfloat128_t::eval_multiply(extfloat128_t& dst, const extfloat128_t& srcA, const extfloat128_t& srcB)
{
  uint32_t expA = srcA.m_exponent;
  uint32_t expB = srcB.m_exponent;
  if (expA != inf_nan_biased_exponent && expB != inf_nan_biased_exponent) {
    if (expA != zero_biased_exponent && expB != zero_biased_exponent) {
      // both numbers are normal
      int64_t expR = static_cast<int64_t>(expA) + expB + 1 - exponent_bias;

      const uint64_t BIT_63      = uint64_t(1) << 63;

      // multiply significands
      uint64_t a0 = srcA.m_significand[0];
      uint64_t a1 = srcA.m_significand[1];
      uint64_t b0 = srcB.m_significand[0];
      uint64_t b1 = srcB.m_significand[1];

      uint64_t x00_h, x00_l = _umul128(a0, b0, &x00_h);
      uint64_t x10_h, x10_l = _umul128(a1, b0, &x10_h);
      uint64_t x0_l = x00_l;
      uint64_t x1_l = x00_h + x10_l;
      uint64_t x2_l = x10_h + (x1_l < x10_l);

      uint64_t x01_h, x01_l = _umul128(a0, b1, &x01_h);
      uint64_t x11_h, x11_l = _umul128(a1, b1, &x11_h);
      uint64_t x0_h = x01_l;
      uint64_t x1_h = x01_h + x11_l;
      uint64_t x2_h = x11_h + (x1_h < x11_l);

      uint64_t x0 = x0_l;
      uint64_t x1 = x1_l + x0_h;
      x2_l  += (x1 < x0_h); // this addition never produces carry, because previous value of x2_l <= 2^64-2
      uint64_t x2 = x2_l + x1_h;
      uint64_t x3 = x2_h + (x2 < x1_h);

      // steaky bits reduction - all x0 bits folded into x1[61..0]
      x1 |= x0 >> 2;
      x1 |= x0 & 3;

      if (x3 < BIT_63) {
        // shift to the left
        x3 = (x3 << 1) | (x2 >> 63);
        x2 = (x2 << 1) | (x1 >> 63);
        x1 = (x1 << 1);
        expR -= 1;
      }

      // round to nearest even
      uint64_t guardsx = x1 | (x2 & 1); // handle tie
      if (guardsx > BIT_63) {
        x2 += 1;
        x3 += (x2 == 0);
        if (x3 == 0) {
          x3 = BIT_63;
          expR += 1;
        }
      }

      if (expR > max_biased_exponent) {
        expR = inf_nan_biased_exponent; // not sure about it
        x3 = x2 = 0;
      }
      if (expR < min_biased_exponent) {
        expR = zero_biased_exponent;
        x3 = x2 = 0;
      }

      dst.m_exponent      = static_cast<uint32_t>(expR);
      dst.m_significand[0] = x2;
      dst.m_significand[1] = x3;
    } else { // a==0 || b==0
      dst.m_exponent      = zero_biased_exponent;
      dst.m_significand[0] = 0;
      dst.m_significand[1] = 0;
    }
    dst.m_sign = srcA.m_sign ^ srcB.m_sign;
  } else {
    eval_non_finite(dst, srcA, srcB, 0, multiply_nfo_tab);
  }
}

void extfloat128_t::eval_multiply_rtz(extfloat128_t& dst, const extfloat128_t& srcA, const extfloat128_t& srcB)
{
  uint32_t expA = srcA.m_exponent;
  uint32_t expB = srcB.m_exponent;
  if (expA != inf_nan_biased_exponent && expB != inf_nan_biased_exponent) {
    if (expA != zero_biased_exponent && expB != zero_biased_exponent) {
      // both numbers are normal
      int64_t expR = static_cast<int64_t>(expA) + expB + 1 - exponent_bias;

      const uint64_t BIT_63      = uint64_t(1) << 63;

      // multiply significands
      uint64_t a0 = srcA.m_significand[0];
      uint64_t a1 = srcA.m_significand[1];
      uint64_t b0 = srcB.m_significand[0];
      uint64_t b1 = srcB.m_significand[1];

      uint64_t x00_h;         _umul128(a0, b0, &x00_h);
      uint64_t x10_h, x10_l = _umul128(a1, b0, &x10_h);
      uint64_t x1_l = x00_h + x10_l;
      uint64_t x2_l = x10_h + (x1_l < x10_l);

      uint64_t x01_h, x01_l = _umul128(a0, b1, &x01_h);
      uint64_t x11_h, x11_l = _umul128(a1, b1, &x11_h);
      uint64_t x0_h = x01_l;
      uint64_t x1_h = x01_h + x11_l;
      uint64_t x2_h = x11_h + (x1_h < x11_l);

      uint64_t x1 = x1_l + x0_h;
      x2_l  += (x1 < x0_h); // this addition never produces carry, because previous value of x2_l <= 2^64-2
      uint64_t x2 = x2_l + x1_h;
      uint64_t x3 = x2_h + (x2 < x1_h);

      if (x3 < BIT_63) {
        // shift to the left
        x3 = (x3 << 1) | (x2 >> 63);
        x2 = (x2 << 1) | (x1 >> 63);
        expR -= 1;
      }

      if (expR > max_biased_exponent) {
        expR = inf_nan_biased_exponent; // not sure about it
        x3 = x2 = 0;
      }
      if (expR < min_biased_exponent) {
        expR = zero_biased_exponent;
        x3 = x2 = 0;
      }

      dst.m_exponent      = static_cast<uint32_t>(expR);
      dst.m_significand[0] = x2;
      dst.m_significand[1] = x3;
    } else { // a==0 || b==0
      dst.m_exponent      = zero_biased_exponent;
      dst.m_significand[0] = 0;
      dst.m_significand[1] = 0;
    }
    dst.m_sign = srcA.m_sign ^ srcB.m_sign;
  } else {
    eval_non_finite(dst, srcA, srcB, 0, multiply_nfo_tab);
  }
}

static const uint8_t reciprocalTab[] = {
 254, 250, 246, 243, 239, 236, 232, 229,
 226, 223, 220, 217, 214, 211, 209, 206,
 204, 201, 199, 196, 194, 192, 189, 187,
 185, 183, 181, 179, 177, 175, 173, 172,
 170, 168, 166, 165, 163, 161, 160, 158,
 157, 155, 154, 152, 151, 150, 148, 147,
 146, 144, 143, 142, 141, 139, 138, 137,
 136, 135, 134, 133, 132, 131, 130, 129,
};

// approximateReciprocal - calculate approximate value of 2^127/(x+1)
// x - argument in range [2^63..2^64-1]
// return - approximate value of 2^127/(x+1).
//          precision: ~53.5 bits.
//          A value guaranteed to be smaller or equal than floor(2^127/(x+1))
static inline uint64_t approximateReciprocal(uint64_t x)
{
  const uint64_t BIT_63 = uint64_t(1)<<63;
  const uint64_t BIT_62 = uint64_t(1)<<62;
  const uint64_t BIT_60 = uint64_t(1)<<60;

  uint64_t y1 = uint64_t(reciprocalTab[(x >> (63-6))-64]);
  uint64_t m2 = (x >> 8)*y1;
  y1 <<= 56;
  uint64_t y2; _umul128(m2, y1, &y2);
  y2 = y1 - y2;
  int64_t  dx = m2 - BIT_63;
  int64_t  dxx; _mul128(dx, dx, &dxx);
  uint64_t uxx = uint64_t(dxx);
  uint64_t y3;  _umul128(BIT_62 + uxx, y2, &y3);
  _umul128(uxx, uxx, &uxx);
  uint64_t y4;  _umul128(BIT_60-32 + uxx, y3-0, &y4);
  return (y4-0) << 7;
}

static const uint8_t divide_nfo_tab[2][2][3] = {
//  NRM                          ZERO                         INF
 {{ NFO_CLASS_INF+NFO_SIGN_POS,  NFO_CLASS_INF+NFO_SIGN_POS,  NFO_CLASS_NAN+NFO_SIGN_POS  },  // A=inf, signA==signB
  { NFO_CLASS_INF+NFO_SIGN_NEG,  NFO_CLASS_INF+NFO_SIGN_NEG,  NFO_CLASS_NAN+NFO_SIGN_POS  }}, // A=inf, signA!=signB
 {{ NFO_CLASS_ZERO+NFO_SIGN_POS, NFO_CLASS_ZERO+NFO_SIGN_POS, NFO_CLASS_NAN+NFO_SIGN_POS  },  // B=inf, signA==signB
  { NFO_CLASS_ZERO+NFO_SIGN_NEG, NFO_CLASS_ZERO+NFO_SIGN_NEG, NFO_CLASS_NAN+NFO_SIGN_POS  }}, // B=inf, signA!=signB
};

#define U64_BITS(first, last) ( (uint64_t(-1) << (first)) & (uint64_t(-1) >> (63-(last))) )
#define U64_BIT(bit_i)        ( uint64_t(1) << (bit_i) )
void extfloat128_t::eval_divide(extfloat128_t& dst, const extfloat128_t& srcA, const extfloat128_t& srcB)
{
  uint32_t expA = srcA.m_exponent;
  uint32_t expB = srcB.m_exponent;
  if (expA != inf_nan_biased_exponent && expB != inf_nan_biased_exponent) {
    if (expA != zero_biased_exponent && expB != zero_biased_exponent) {
      // both numbers are normal
      int64_t expR = static_cast<int64_t>(expA) - expB - 1 + exponent_bias;

      // divide significands
      uint64_t a1 = srcA.m_significand[0];
      uint64_t a2 = srcA.m_significand[1];
      uint64_t bl = srcB.m_significand[0];
      uint64_t bh = srcB.m_significand[1];
      uint64_t a0 = 0;
      if (a2 > bh || (a2 == bh && a1 >= bl)) {
        a0 = (a1 << 63);
        a1 = (a1 >> 1) | (a2 << 63);
        a2 = (a2 >> 1);
        expR += 1;
      }

      // static const double DIV_SCALE = double((int64_t(1) << 53)-2)*double(int64_t(1) << 47)*double(int64_t(1) << 25); // 2^125-eps*2
      // uint64_t invB = int64_t(DIV_SCALE / int64_t(bh>>1)) << 1; // invB = 2^127/bh. Range [2^63..2^64)
      uint64_t invB = approximateReciprocal(bh);

      uint64_t xH, xL, x1, x2;

      uint64_t r0;  _umul128(a2, invB, &r0); r0 *= 2;           // r0 = (a2*2^64)/bh.  Range [2^63..2^64)
      xL  = _umul128(bl, r0, &xH);
      x1  = xH + (a0 < xL);
      a0 -= xL;
      xL  = _umul128(bh, r0, &xH);
      x1 += xL;
      x2  = xH + (x1 < xL);
      x2 += (a1 < x1);
      a1 -= x1;
      a2 -= x2;

      a2 = (a2 << 48) | (a1 >> 16);
      a1 = (a1 << 48) | (a0 >> 16);
      a0 = (a0 << 48);

      uint64_t r1;  _umul128(a2, invB, &r1); r1 *= 2;
      xL  = _umul128(bl, r1, &xH);
      x1  = xH + (a0 < xL);
      a0 -= xL;
      xL  = _umul128(bh, r1, &xH);
      x1 += xL;
      x2  = xH + (x1 < xL);
      x2 += (a1 < x1);
      a1 -= x1;
      a2 -= x2;

      a2 = (a2 << 48) | (a1 >> 16);
      a1 = (a1 << 48) | (a0 >> 16);
      a0 = (a0 << 48);

      uint64_t r2;  _umul128(a2, invB, &r2); r2 *= 2;
      if ((r2 & U64_BITS(20, 31)) == U64_BITS(20, 30)) {
        // Result is close to a middle point between two numbers
        // Let's calculate an exact value of bit r2[31]
        uint64_t r2x = (r2 & U64_BITS(32, 63)) | U64_BIT(31); // bit[31] set, other LS bits cleared

        xL  = _umul128(bl, r2x, &xH);
        x1  = xH + (a0 < xL);
        xL  = _umul128(bh, r2x, &xH);
        x1 += xL;
        x2  = xH + (x1 < xL);
        x2 += (a1 < x1);
        if (a2 >= x2)
          r2 = r2x;
      }

      r2 = (r2 >> 32) + ((r2 >> 31) & 1);
      r1 += (r2 >> 16);
      uint64_t rH = r0 + (r1 >> 48);
      uint64_t rL = (r1 << 16) + (r2 & U64_BITS(0, 15));

      if (expR > max_biased_exponent) {
        expR = inf_nan_biased_exponent; // return inf. not sure about it
        rL = rH = 0;
      }
      if (expR < min_biased_exponent) {
        expR = zero_biased_exponent; // return zero
        rL = rH = 0;
      }

      dst.m_exponent      = expR;
      dst.m_significand[0] = rL;
      dst.m_significand[1] = rH;
    } else if (expA == zero_biased_exponent) {
      dst.m_exponent      = zero_biased_exponent;
      dst.m_significand[0] = 0;
      dst.m_significand[1] = 0;
      if (expB==zero_biased_exponent) { // 0/0
        dst.m_exponent      = inf_nan_biased_exponent;
        dst.m_significand[1] = qnan_bit;
      }
    } else { // (expB == zero_biased_exponent)
      dst.m_exponent      = inf_nan_biased_exponent;
      dst.m_significand[0] = 0;
      dst.m_significand[1] = 0;
    }
    dst.m_sign = srcA.m_sign ^ srcB.m_sign;
  } else {
    eval_non_finite(dst, srcA, srcB, 0, divide_nfo_tab);
  }
}

void extfloat128_t::eval_add(extfloat128_t& dst, const extfloat128_t& srcA, const extfloat128_t& srcB) {
  eval_addsub(dst, srcA, srcB, 0);
}
void extfloat128_t::eval_subtract(extfloat128_t& dst, const extfloat128_t& srcA, const extfloat128_t& srcB) {
  eval_addsub(dst, srcA, srcB, 1);
}

void extfloat128_t::op_int(extfloat128_t* dst, const extfloat128_t* srcA, int srcB, void (*func)(extfloat128_t& dst, const extfloat128_t& srcA, const extfloat128_t& srcB)) {
  op_double(dst, srcA, srcB, func);
}

void extfloat128_t::op_double(extfloat128_t* dst, const extfloat128_t* srcA, double srcB, void (*func)(extfloat128_t& dst, const extfloat128_t& srcA, const extfloat128_t& srcB)) {
  extfloat128_t b(srcB);
  func(*dst, *srcA, b);
}

void extfloat128_t::op_double_rev(extfloat128_t* dst, double srcA, const extfloat128_t* srcB, void (*func)(extfloat128_t& dst, const extfloat128_t& srcA, const extfloat128_t& srcB)) {
  extfloat128_t a(srcA);
  func(*dst, a, *srcB);
}

static const uint8_t add_nfo_tab[2][2][3] = {
//  NRM                          ZERO                         INF
 {{ NFO_CLASS_INF+NFO_SIGN_A, NFO_CLASS_INF+NFO_SIGN_A, NFO_CLASS_INF+NFO_SIGN_A    },  // A=inf, signA==signB
  { NFO_CLASS_INF+NFO_SIGN_A, NFO_CLASS_INF+NFO_SIGN_A, NFO_CLASS_NAN+NFO_SIGN_POS  }}, // A=inf, signA!=signB
 {{ NFO_CLASS_INF+NFO_SIGN_B, NFO_CLASS_INF+NFO_SIGN_B, NFO_CLASS_INF+NFO_SIGN_A    },  // B=inf, signA==signB
  { NFO_CLASS_INF+NFO_SIGN_B, NFO_CLASS_INF+NFO_SIGN_B, NFO_CLASS_NAN+NFO_SIGN_POS  }}, // B=inf, signA!=signB
};

void extfloat128_t::eval_addsub(extfloat128_t& dst, const extfloat128_t& srcA, const extfloat128_t& srcB, int isSub)
{
  const uint64_t BIT_63 = uint64_t(1) << 63;
  const uint64_t BIT_62 = uint64_t(1) << 62;

  uint32_t expA = srcA.m_exponent;
  uint32_t expB = srcB.m_exponent;
  if (expA != inf_nan_biased_exponent && expB != inf_nan_biased_exponent) {
    if (expB != zero_biased_exponent) {
      if (expA != zero_biased_exponent) {
        // both operands normal
        uint64_t uExpDiff = uint64_t(expA) - uint64_t(expB);
        if (uExpDiff + 128 <= 256) {
          // both operands significant
          unsigned signA = srcA.m_sign;
          unsigned sub   = srcB.m_sign ^ isSub ^ signA;
          if (uExpDiff==0) {
            uint64_t a0 = srcA.m_significand[0];
            uint64_t a1 = srcA.m_significand[1];
            uint64_t b0 = srcB.m_significand[0];
            uint64_t b1 = srcB.m_significand[1];
            if (0 == sub) {
              // same exponent : addition
              uint64_t s0 = a0 + b0;
              uint64_t s1 = a1 + b1 + (s0 < b0);
              // shift right
              uint64_t r1 = (s1 >> 1) | BIT_63;
              uint64_t r0a= (s0 >> 1) | (s1 << 63);
              uint64_t r0 = r0a + (r0a & s0 & 1);
              r1 += (r0 < r0a);
              uint32_t expR = expA + 1;
              if (expR > max_biased_exponent) {
                expR = inf_nan_biased_exponent;
                r0 = r1 = 0;
              }
              dst.m_exponent      = expR;
              dst.m_sign           = signA;
              dst.m_significand[0] = r0;
              dst.m_significand[1] = r1;
            } else {
              // same exponent : subtraction
              if (b1 >= a1 && (b1 > a1 || b0 > a0)) {
                uint64_t tmp;
                tmp = a0; a0 = b0; b0 = tmp;
                tmp = a1; a1 = b1; b1 = tmp;
                signA ^= 1;
              }
              uint64_t r0 = a0 - b0;
              uint64_t r1 = a1 - b1 - (a0 < b0);
              int64_t expR = expA;
              // at least 1 MS bits guaranteed to be canceled
              unsigned long shift;
              if (_BitScanReverse64(&shift, r1)) {
                // MS word partially canceled
                r1 = (r1 << (63-shift)) | (r0 >> (shift+1));
                r0 = (r0 << (63-shift));
                expR -= 63-shift;
              }
              else { // MS word is totally canceled
                if (_BitScanReverse64(&shift, r0)) {
                  r1 = r0 << (63-shift);
                  r0 = 0;
                  expR -= 127-shift;
                }
                else {
                  // both words canceled
                  expR  = min_biased_exponent - 1;
                  signA = 0; // According to IEEE rules subtraction of two equal numbers rresults in positive zero
                }
              }
              if (expR < min_biased_exponent) {
                expR = zero_biased_exponent;
                r0 = r1 = 0;
              }
              dst.m_exponent      = expR;
              dst.m_sign           = signA;
              dst.m_significand[0] = r0;
              dst.m_significand[1] = r1;
            }
          } else {
            // expA-expB in range [-128..128]
            int32_t iExpDiff = int32_t(uExpDiff);
            const extfloat128_t* pSrcA = &srcA;
            const extfloat128_t* pSrcB = &srcB;
            int64_t expR = expA;
            if (iExpDiff < 0) {
              expR  -= iExpDiff;
              signA ^= sub;
              iExpDiff = -iExpDiff;
              const extfloat128_t* tmp = pSrcA;
              pSrcA = pSrcB;
              pSrcB = tmp;
            }
            uint64_t a0 = pSrcA->m_significand[0];
            uint64_t a1 = pSrcA->m_significand[1];
            uint64_t b0 = pSrcB->m_significand[0];
            uint64_t b1 = pSrcB->m_significand[1];
            uint64_t r0, r1;
            if (iExpDiff == 1 && sub != 0) {
              // subtract numbers with exponents that differ by 1
              // cancelation of MS bit is 3 times more likely than otherwise, so do it as a default
              a1 = (a1 << 1) | (a0 >> 63);
              a0 = (a0 << 1);
              r0 = a0 - b0;
              r1  = a1 - b1 - (a0 < b0);
              if (a1 < b1 || r1 >= BIT_63) {
                // at least 1 MS bit canceled
                unsigned long shift;
                if (_BitScanReverse64(&shift, r1)) {
                  // MS word partially canceled
                  r1 = (r1 << (63-shift)) | ((r0 >> shift) >> 1);
                  r0 = (r0 << (63-shift));
                  expR -= 64-shift;
                } else { // MS word is totally canceled
                  _BitScanReverse64(&shift, r0); // full cancelation is not possible here
                  r1 = r0 << (63-shift);
                  r0 = 0;
                  expR -= 128-shift;
                }
              } else {
                // no cancelation - round to nearest even
                uint64_t odd = (r0 >> 1) & 1;
                r0 += odd;
                if (r0 < odd)
                  r1 += 1;
                r0 = (r0 >> 1) | (r1 << 63);
                r1 = (r1 >> 1) | BIT_63;
              }
              if (expR < min_biased_exponent) {
                expR = zero_biased_exponent;
                r0 = r1 = 0;
              }
            } else {
              int bshift = iExpDiff & 63;
              uint64_t bL = (b0 << (63-bshift)) << 1;
              b0 = (b0 >> bshift) | ((b1 << (63-bshift)) << 1);
              b1 = (b1 >> bshift);
              if (iExpDiff >= 64) {
                bL = b0 | (bL!=0);
                b0 = b1;
                b1 = 0;
                if (iExpDiff == 128) {
                  bL = b0 | (bL!=0);
                  b0 = 0;
                }
              }
              if (sub == 0) {
                // different exponents : addition
                r0 = a0 + b0;
                b1 += (r0 < b0);
                r1 = a1 + b1;
                if (r1 >= b1) {
                  // no carry out
                  bL |= (r0 & 1);
                  if (bL > BIT_63) {
                    r0 += 1;
                    if (r0 == 0)
                      r1 += 1;
                    if (r1 == 0) {
                      r1 = BIT_63;
                      expR += 1;
                    }
                  }
                } else {
                  // carry out - shift right by 1
                  bL = (bL >> 1) | (bL & 1) | (r0 << 63);
                  r0 = (r0 >> 1) | (r1 << 63);
                  r1 = (r1 >> 1) | BIT_63;
                  bL |= (r0 & 1);
                  if (bL > BIT_63) {
                    r0 += 1;
                    if (r0 == 0)
                      r1 += 1;
                  }
                  expR += 1;
                }
                if (expR > max_biased_exponent) {
                  expR = inf_nan_biased_exponent;
                  r0 = r1 = 0;
                }
              } else {
                // different exponents (diff > 1): subtraction
                r0 = a0 - b0;
                r1 = a1 - b1 - (a0 < b0);
                uint64_t rx = r1 | (r0 != 0);
                if (rx > BIT_63) {
                  // no cancellation
                  // round r0[0]
                  bL |= (r0 & 1);
                  if (bL > BIT_63) {
                    if (r0 == 0)
                      r1 -= 1;
                    r0 -= 1;
                  }
                } else if (rx < BIT_63 || bL > BIT_62) {
                  // MB bit canceled - shift left
                  r1 = (r1 << 1) | (r0 >> 63);
                  r0 = (r0 << 1);
                  // round bL[63]
                  if (bL > BIT_62) {
                    bL >>= 62;
                    bL += 1;
                    bL >>= 1;
                    if (r0 < bL)
                      r1 -= 1;
                    r0 -= bL;
                  }
                  expR -= 1;
                  if (expR < min_biased_exponent) {
                    expR = zero_biased_exponent;
                    r0 = r1 = 0;
                  }
                }
              }
            }
            dst.m_exponent       = expR;
            dst.m_sign           = signA;
            dst.m_significand[0] = r0;
            dst.m_significand[1] = r1;
          }
        } else {
          // one of operands is void
          uint32_t sub = srcA.m_sign ^ srcB.m_sign ^ isSub;
          uint64_t voidSignificand;
          if (expA > expB) {
            voidSignificand = srcB.m_significand[0] | (srcB.m_significand[1] << 1);
            dst = srcA; // b is void
          } else {
            voidSignificand = srcA.m_significand[0] | (srcA.m_significand[1] << 1);
            dst = srcB; // a is void
            dst.m_sign ^= isSub;
          }

          // check for special case 2^N - 2^N-129 - eps
          if (uExpDiff + 129 <= 129*2   &&
              sub != 0                  &&
              voidSignificand != 0      &&
              dst.m_significand[0] == 0 &&
              dst.m_significand[1] == BIT_63) {
            dst.m_exponent       -= 1;
            dst.m_significand[0] = uint64_t(-1);
            dst.m_significand[1] = uint64_t(-1);
          }
        }
      } else {
        // a==0, b != 0
        dst = srcB;
        dst.m_sign ^= isSub;
      }
    } else {
      // b==0
      unsigned signA = srcA.m_sign;
      unsigned signB = srcB.m_sign ^ isSub;
      dst = srcA;
      if (expA == zero_biased_exponent && signA != signB) {
        // two zeros with different signs => +0
        dst.m_sign = 0;
      }
    }
  } else {
    eval_non_finite(dst, srcA, srcB, isSub, add_nfo_tab);
  }
}

void extfloat128_t::eval_fma(extfloat128_t& dst, const extfloat128_t& srcA, const extfloat128_t& srcB, const extfloat128_t& srcZ)
{
  uint32_t expA = srcA.m_exponent;
  uint32_t expB = srcB.m_exponent;
  if (expA != inf_nan_biased_exponent && expB != inf_nan_biased_exponent) {
    // both a and b are finite
    if (expA != zero_biased_exponent && expB != zero_biased_exponent) {
      // both a and b are normal
      int64_t expPr = static_cast<int64_t>(expA) + expB + 1 - exponent_bias;
      const uint64_t BIT_63      = uint64_t(1) << 63;
      const uint64_t BIT_62      = uint64_t(1) << 62;

      // multiply significands
      uint64_t a0 = srcA.m_significand[0];
      uint64_t a1 = srcA.m_significand[1];
      uint64_t b0 = srcB.m_significand[0];
      uint64_t b1 = srcB.m_significand[1];

      uint64_t pr0, pr1, pr2, pr3, la, ha, lb, hb, lc;
      uint8_t c;
      pr0 = _umul128(a0, b0, &pr1);
      la = _umul128(a0, b1, &ha);
      lb = _umul128(a1, b0, &hb);
      pr1 += la; c  = pr1 < la;
      pr1 += lb; c += pr1 < lb;
      pr2 = c;
      pr2 += ha; c  = pr2 < ha;
      pr2 += hb; c += pr2 < hb;
      lc = _umul128(a1, b1, &pr3);
      pr2 += lc; c += pr2 < lc;
      pr3 += c;
      unsigned signPr = srcA.m_sign ^ srcB.m_sign;

      uint32_t expZ32 = srcZ.m_exponent;
      if (expZ32 != inf_nan_biased_exponent) {
        int64_t expZ = (expZ32 == zero_biased_exponent) ? expPr : expZ32;
        if (expZ <= expPr) {
          int64_t dExp64 = expPr - expZ;
          if (dExp64 > 256) dExp64 = 256;
          unsigned dExp = static_cast<unsigned>(dExp64);
          uint64_t z0 = srcZ.m_significand[0];
          uint64_t z1 = srcZ.m_significand[1];
          uint64_t zx = 0;
          unsigned rshift = dExp % 64;
          if (rshift) {
            unsigned lshift = 64-rshift;
            zx = (z0 << lshift);
            z0 = (z0 >> rshift) | (z1 << lshift);
            z1 = (z1 >> rshift);
          }
          unsigned signZ = srcZ.m_sign;
          if (signZ != signPr) {
            // subtraction
            if (dExp < 64*1) {
              uint8_t b, bb;
              b  = pr1 < zx;        pr1 -= zx;
              bb = pr2 < b;         pr2 -= b;
              b  = bb | (pr2 < z0); pr2 -= z0;
                                    pr3 -= b;
              b  = pr3 < z1;        pr3 -= z1;
              if (b) {
                // negate pr3:pr2:pr1:pr0
                pr3 = ~pr3;
                pr2 = ~pr2;
                pr1 = ~pr1;
                pr0 = 0-pr0;
                if (pr0 == 0) {
                  pr1 += 1;
                  if (pr1 == 0) {
                    pr2 += 1;
                    if (pr2 == 0)
                      pr3 += 1;
                  }
                }
                signPr ^= 1;
              }
              if (pr3 < BIT_63) {
                unsigned long shift;
                if (_BitScanReverse64(&shift, pr3)) {
                  // MS word partially canceled
                  expPr -= 64*1-1-shift;
                  unsigned lshift = 63-shift;
                  unsigned rshift = 1+shift;
                  pr3 = (pr3 << lshift) | (pr2 >> rshift);
                  pr2 = (pr2 << lshift) | (pr1 >> rshift);
                  pr1 = (pr1 << lshift) | (pr0 >> rshift);
                  pr0 = (pr0 << lshift);
                } else if (_BitScanReverse64(&shift, pr2)) {
                  // 2nd word not fully canceled
                  expPr -= 64*2-1-shift;
                  unsigned lshift = 63-shift;
                  pr3 = (pr2 << lshift) | ((pr1 >> shift) >> 1);
                  pr2 = (pr1 << lshift) | ((pr0 >> shift) >> 1);
                  pr1 = (pr0 << lshift);
                  pr0 = 0;
                } else if (_BitScanReverse64(&shift, pr1)) {
                  // 3rd word not fully canceled
                  expPr -= 64*3-1-shift;
                  unsigned lshift = 63-shift;
                  pr3 = (pr1 << lshift) | ((pr0 >> shift) >> 1);
                  pr2 = (pr0 << lshift);
                  pr1 = pr0 = 0;
                } else if (_BitScanReverse64(&shift, pr0)) {
                  // 4th word not fully canceled
                  expPr -= 64*4-1-shift;
                  unsigned lshift = 63-shift;
                  pr3 = (pr0 << lshift);
                  pr2 = pr1 = pr0 = 0;
                } else {
                  // full cancelation
                  expPr  = min_biased_exponent - 1;
                  signPr = 0; // According to IEEE rules subtraction of two equal numbers results in positive zero
                  goto assign_pr;
                }
              }
              zx = 0;
            } else {
              if (dExp < 64*2) {
                uint8_t b, bb;
                b  = pr0 < zx;        pr0 -= zx;
                bb = pr1 < b;         pr1 -= b;
                b  = bb | (pr1 < z0); pr1 -= z0;
                bb = pr2 < b;         pr2 -= b;
                b  = bb | (pr2 < z1); pr2 -= z1;
                                      pr3 -= b;
                zx = 0;
              } else if (dExp < 64*3) {
                uint8_t b, bb;
                b  = pr0 < z0;        pr0 -= z0;
                bb = pr1 < b;         pr1 -= b;
                b  = bb | (pr1 < z1); pr1 -= z1;
                bb = pr2 < b;         pr2 -= b; b = bb;
                                      pr3 -= b;
              } else if (dExp < 64*4) {
                uint8_t b, bb;
                b  = pr0 < z1; pr0 -= z1;
                bb = pr1 < b;  pr1 -= b; b = bb;
                bb = pr2 < bb; pr2 -= b; b = bb;
                               pr3 -= b;
                zx = zx | z0;
              } else {
                zx = zx | z0 | z1;
              }
              if (pr3 < BIT_63) {
                if (pr3 < BIT_62) {
                  expPr -= 2;
                  pr3 = (pr3 << 2) | (pr2 >> 62);
                  pr2 = (pr2 << 2) | (pr1 >> 62);
                  pr1 = (pr1 << 2) | (pr0 >> 62);
                  pr0 = (pr0 << 2);
                } else {
                  expPr -= 1;
                  pr3 = (pr3 << 1) | (pr2 >> 63);
                  pr2 = (pr2 << 1) | (pr1 >> 63);
                  pr1 = (pr1 << 1) | (pr0 >> 63);
                  pr0 = (pr0 << 1);
                }
              }
            }
            // prepare for rounding
            pr1 |= (pr0 != 0);
            zx = (zx != 0) & (pr1==BIT_63);
            pr1 -= zx;
          } else {
            // addition
            uint8_t c;
            if (dExp < 64*1) {
              pr1 += zx; c  = pr1 < zx;
              pr2 += c;  c  = pr2 < c;
              pr2 += z0; c |= pr2 < z0;
              pr3 += c;  c  = pr3 < c;
              pr3 += z1; c |= pr3 < z1;
              zx = 0;
            } else if (dExp < 64*2) {
              pr0 += zx; c  = pr0 < zx;
              pr1 += c;  c  = pr1 < c;
              pr1 += z0; c |= pr1 < z0;
              pr2 += c;  c  = pr2 < c;
              pr2 += z1; c |= pr2 < z1;
              pr3 += c;  c  = pr3 < c;
              zx = 0;
            } else if (dExp < 64*3) {
              pr0 += z0; c  = pr0 < z0;
              pr1 += c;  c  = pr1 < c;
              pr1 += z1; c |= pr1 < z1;
              pr2 += c;  c  = pr2 < c;
              pr3 += c;  c  = pr3 < c;
            } else if (dExp < 64*4) {
              pr0 += z1; c  = pr0 < z1;
              pr1 += c;  c  = pr1 < c;
              pr2 += c;  c  = pr2 < c;
              pr3 += c;  c  = pr3 < c;
              pr0 |= ((zx|z0) != 0);
              zx = zx | z0;
            } else {
              zx = zx | z0 | z1;
              c = 0;
            }
            if (c) {
              expPr += 1;
              zx |= (pr0 & 1);
              pr0 = (pr0 >> 1) | (pr1 << 63);
              pr1 = (pr1 >> 1) | (pr2 << 63);
              pr2 = (pr2 >> 1) | (pr3 << 63);
              pr3 = (pr3 >> 1) | BIT_63;
            } else if (pr3 < BIT_63) {
              expPr -= 1;
              pr3 = (pr3 << 1) | (pr2 >> 63);
              pr2 = (pr2 << 1) | (pr1 >> 63);
              pr1 = (pr1 << 1) | (pr0 >> 63);
              pr0 = (pr0 << 1);
            }
            // prepare for rounding
            pr1 |= ((pr0|zx)!=0);
          }
          // round
          {
          pr1 |= (pr2 & 1);
          uint64_t incr = (pr1 > BIT_63);
          pr2 += incr;
          pr3 += (pr2 < incr);
          }
          if (pr3 == 0) {
            pr3 = BIT_63;
            expPr += 1;
          }
          // check range of exponent
          // printf("%I64d *\n", expPr);
          if (expPr <= max_biased_exponent) {
            if (expPr >= min_biased_exponent) {
            } else {
              expPr = zero_biased_exponent;
              pr2 = pr3 = 0;
            }
          } else {
            expPr = inf_nan_biased_exponent;
            pr2 = pr3 = 0;
          }
          assign_pr:
          dst.m_significand[0] = pr2;
          dst.m_significand[1] = pr3;
          dst.m_exponent      = expPr;
          dst.m_sign           = signPr;
        } else { // expZ > expPr
          uint64_t z1;
          uint64_t z2 = srcZ.m_significand[0];
          uint64_t z3 = srcZ.m_significand[1];
          unsigned signZ = srcZ.m_sign;
          int64_t dExp64 = expZ - expPr;
          if (dExp64 <= 129) {
            unsigned dExp = static_cast<unsigned>(dExp64);
            if (dExp == 1 && signZ != signPr) {
              // subtract numbers with exponents that differ by 1
              // cancellation of MS bit is 3 times more likely than otherwise, so do it as a default
              z3 = (z3 << 1) | (z2 >> 63);
              z2 = (z2 << 1);
              uint8_t b0 = 0  < pr0; uint64_t z0 = 0 - pr0;
              uint8_t b1 = 0  < pr1; z1 = 0 - pr1;
              uint8_t b2 = z2 < pr2; z2 -= pr2;
              uint8_t b3 = z3 < pr3; z3 -= pr3;
              b1 |= z1 < b0; z1 -= b0;
              b2 |= z2 < b1; z2 -= b1;
              b3 |= z3 < b2; z3 -= b2;
              if (b3) {
                // at least 1 MS bit canceled
                unsigned long shift;
                if (_BitScanReverse64(&shift, z3)) {
                  // MS word partially canceled
                  unsigned lshift = 63 - shift;
                  if (lshift != 0) {
                    unsigned rshift = shift + 1;
                    z3 = (z3 << lshift) | (z2 >> rshift);
                    z2 = (z2 << lshift) | (z1 >> rshift);
                    z1 = (z1 << lshift) | (z0 >> rshift);
                    z0 = (z0 << lshift);
                  }
                  expZ -= 64-shift;
                } else { // MS word is totally canceled
                  z3 = z2; z2 = z1 ; z1 = z0; z0 = 0;
                  _BitScanReverse64(&shift, z3); // full cancellation is not possible here
                  unsigned lshift = 63 - shift;
                  if (lshift != 0) {
                    unsigned rshift = shift + 1;
                    z3 = (z3 << lshift) | (z2 >> rshift);
                    z2 = (z2 << lshift) | (z1 >> rshift);
                    z1 = (z1 << lshift);
                  }
                  expZ -= 128-shift;
                }
                z1 |= (z0 != 0);
              } else {
                // no cancelation - round to nearest even
                uint64_t odd = (z2 >> 1) & 1;
                z2 += odd | ((z1|z0)!=0);
                if (z2 < odd)
                  z3 += 1;
                z2 = (z2 >> 1) | (z3 << 63);
                z3 = (z3 >> 1) | BIT_63;
                z1 = z0 = 0;
              }
              //printf("%d:%I64x:%I64x:%I64x *%d\n", signZ, z3, z2, z1, b3);
            } else {
              pr1 |= (pr0 != 0);
              unsigned rshift = dExp % 64;
              if (rshift) {
                unsigned lshift = 64-rshift;
                uint64_t prx = (pr1 << lshift);
                pr1 = (pr1 >> rshift) | (pr2 << lshift);
                pr2 = (pr2 >> rshift) | (pr3 << lshift);
                pr3 = (pr3 >> rshift);
                pr1 |= (prx != 0);
              }
              if (signZ == signPr) {
                // addition
                uint8_t c;
                if        (dExp < 64*1) {
                  z2 += pr2; c  = z2 < pr2;
                  z3 += c;   c  = z3 < c;
                  z3 += pr3; c |= z3 < pr3;
                } else if (dExp < 64*2) {
                  z2 += pr3; c  = z2 < pr3;
                  z3 += c;   c  = z3 < c;
                  pr1 = pr2 | (pr1 != 0);
                } else {
                  c = 0;
                  pr1 = pr3 | ((pr2 | pr1) != 0);
                }
                if (c) {
                  expZ += 1;
                  uint64_t prx = (pr1 & 1);
                  pr1 = (pr1 >> 1) | (z2 << 63);
                  z2  = (z2  >> 1) | (z3 << 63);
                  z3  = (z3  >> 1) | BIT_63;
                  pr1 |= prx;
                }
                z1 = pr1;
              } else {
                // different exponents (diff > 1): subtraction
                // at most 1 bit can be canceled
                if        (dExp < 64*1) {
                  z3 -= (z2 < pr2);
                  z2 -= pr2;
                  z3 -= pr3;
                } else if (dExp < 64*2) {
                  z3 -= (z2 < pr3);
                  z2 -= pr3;
                  pr1 = pr2 | (pr1 != 0);
                } else {
                  pr1 = pr3 | ((pr2 | pr1) != 0);
                }
                uint64_t decr = (pr1 != 0);
                z3 -= (z2 < decr);
                z2 -= decr;
                z1 = 0 - pr1;
                if (z3 < BIT_63) {
                  expZ -= 1;
                  z3 = (z3 << 1) | (z2 >> 63);
                  z2 = (z2 << 1) | (z1 >> 63);
                  z1 = (z1 << 1);
                }
                // printf("%d:%I64x:%I64x:%I64x\n", signZ, z3, z2, z1);
              }
            }
            // rounding
            z1 |= (z2 & 1); // tie rounded to nearest even
            uint64_t incr = (z1 > BIT_63);
            z2 += incr;
            z3 += (z2 < incr);
            if (z3 == 0) {
              z3 = BIT_63;
              expZ += 1;
            }
            // check range of exponent
            if (expZ <= max_biased_exponent) {
              if (expPr >= min_biased_exponent) {
              } else {
                expZ = zero_biased_exponent;
                z2 = z3 = 0;
              }
            } else {
              expZ = inf_nan_biased_exponent;
              z2 = z3 = 0;
            }
          }
          // printf("%u **\n", expZ);
          dst.m_significand[0] = z2;
          dst.m_significand[1] = z3;
          dst.m_exponent       = expZ;
          dst.m_sign           = signZ;
        }
      } else { // z non-finite
        extfloat128_t prod;
        prod.m_exponent = exponent_bias;
        prod.m_sign      = signPr;
        eval_non_finite(dst, prod, srcZ, 0, add_nfo_tab);
      }
    } else { // a==0 || b==0
      uint32_t expZ = srcZ.m_exponent;
      if (expZ != zero_biased_exponent) {
        dst = srcZ;
      } else {
        // z == 0
        dst.m_exponent      = zero_biased_exponent;
        dst.m_significand[0] = 0;
        dst.m_significand[1] = 0;
        unsigned signPr = srcA.m_sign ^ srcB.m_sign;
        unsigned signZ = srcZ.m_sign;
        dst.m_sign = signPr == signZ ? signZ : 0;
      }
    }
  } else {
    extfloat128_t prod;
    eval_non_finite(prod, srcA, srcB, 0, multiply_nfo_tab);
    eval_non_finite(dst, prod, srcZ, 0, add_nfo_tab);
  }
}


#define U32_BITS(first, last) ( (uint32_t(-1) << (first)) & (uint32_t(-1) >> (31-(last))) )
#define U32_BIT(bit_i)        ( uint32_t(1) << (bit_i) )

// calculate approximation of 2^63*sqrt(2^63)/sqrt(m*2^e) with precision of 49+ bits.
// result guaranteed to be < 2^63
static int64_t rsqrt_u64(uint64_t m, int e)
{
  enum {
    RSQRT_TAB_ABITS = 6,
    RSQRT_TAB_SZ    = 1 << RSQRT_TAB_ABITS,
  };
  static uint16_t rsqrtTab[RSQRT_TAB_SZ+1] = {
    65535,  64523,  63534,  62567,  61622,  60697,  59792,  58907,
    58039,  57190,  56358,  55543,  54743,  53960,  53191,  52438,
    51698,  50972,  50259,  49560,  48873,  48198,  47534,  46883,
    46242,  45612,  44993,  44384,  43785,  43196,  42616,  42045,
    41483,  40930,  40386,  39849,  39321,  38801,  38288,  37783,
    37285,  36794,  36310,  35833,  35363,  34899,  34441,  33990,
    33545,  33105,  32672,  32244,  31821,  31404,  30993,  30586,
    30185,  29789,  29397,  29010,  28628,  28251,  27878,  27510,
    27145,
  };
  // do NR iteration with second-order polinomial P= ((A-x)^2+B)*C
  // where
  //   A = 1.6667058478613405725
  //   B = 2.2223267059585696259
  //   C = 0.37497796157994628397
  const uint64_t COEFF_A = uint64_t((uint64_t(1) << 61)*1.6667058478613405725+0.5);
  //const uint64_t COEFF_B = uint64_t((uint64_t(1) << 62)*2.2223267059585696259+0.5);
  const uint64_t COEFF_B = uint64_t((uint64_t(1) << 62)*2.2223267059585696259+0.5)+2455;
  const uint64_t COEFF_C0 = uint64_t((uint64_t(1) << 63)*2.0*0.37497796157994628397+0.5);
  const uint64_t COEFF_C1 = uint64_t((uint64_t(1) << 63)*2.0*0.37497796157994628397*0.70710678118654752440084436210485+0.5);
  static const uint64_t c_tab[] = {COEFF_C0, COEFF_C1};
  int tabI = (m>>(63-RSQRT_TAB_ABITS)) & (RSQRT_TAB_SZ-1); // MSB is not a part of an index
  uint32_t t0 = rsqrtTab[tabI+0];
  uint32_t t1 = rsqrtTab[tabI+1];
  uint32_t interp = (m>>(63-RSQRT_TAB_ABITS-22)) & ((1u<<22)-1);
  uint64_t y0 = ((uint64_t(t0) << 22) - (t0 - t1)*interp + (uint64_t(0x10001) << 22))<<24;
  // return y0 << 1;
  uint64_t yy; _umul128(y0, y0, &yy);
  uint64_t y1; _umul128(y0, c_tab[e], &y1);
  uint64_t mx; _umul128(m, yy, &mx);
  uint64_t amx = (COEFF_A - mx)<<2;
  uint64_t amxx; _umul128(amx, amx, &amxx);
  _umul128(y1<<2, amxx+COEFF_B, &y1);
  return y1;
}

void extfloat128_t::eval_sqrt(extfloat128_t& dst, const extfloat128_t& src)
{
  uint32_t expA = src.m_exponent - 1;
  if (expA < inf_nan_biased_exponent-1) {
    if (src.m_sign == 0) {
      int e = expA & 1;
      uint64_t x3 = src.m_significand[1];
      uint64_t x2 = src.m_significand[0];
      int64_t  yw = rsqrt_u64(x3, e);
      uint64_t r1; _umul128(x3, uint64_t(yw) << 1, &r1);

      r1 <<= e;
      uint64_t rr3, rr2, rr1, rr0;
      rr2 = _umul128(r1, r1, &rr3);
      uint64_t rr = (rr3 << 44) | (rr2 >> 20);
      uint64_t x  = (x3  << 43) | (x2  >> 21);
      int64_t  d = rr - (x << e);
      // if (div_dbg) printf("%016I64x d\n", d);
      _mul128(d, yw, &d);
      // if (div_dbg) printf("%016I64x d\n", d);
      // if (div_dbg) printf("%016I64x r1\n", r1);
      uint64_t dL = uint64_t(d) << 21;
      uint64_t dH = d >> 43;
      uint64_t r0 = uint64_t(-1) - dL;
      r1 -=  dH + 1;
      // if (div_dbg) printf("%016I64x:%016I64x r\n", r1, r0);

      // multiply 128b x 128b => 256b. Calculate only bits[191..64]
      uint64_t rr00_h;          _umul128(r0, r0, &rr00_h);
      uint64_t rr01_h, rr01_l = _umul128(r0, r1, &rr01_h);
      uint64_t rrL_1 = rr00_h + rr01_l;
      uint64_t rrL_2 = rr01_h + (rrL_1 < rr01_l);
      uint64_t rr11_l = r1 * r1;
      uint64_t rrH_0 = rr01_l;
      uint64_t rrH_1 = rr01_h + rr11_l;
      rr1 = rrL_1 + rrH_0;
      rr2 = rrL_2 + (rr1 < rrH_0);
      rr2 = rr2 + rrH_1;

      // if (div_dbg) printf("%016I64x:%016I64x rr\n", rr2, rr1);
      // if (div_dbg) printf("%016I64x:%016I64x x\n", x2, x1);
      rr = (rr2 << 30) | (rr1 >> 34);
      x  = (x2  << 29) << e;
      d = rr - x;
      // if (div_dbg) printf("%016I64x rr\n", rr);
      // if (div_dbg) printf("%016I64x x\n",  x);
      // if (div_dbg) printf("%016I64x d\n",  d);
      _mul128(d, yw, &d);
      dH = d >> 29;
      r1 -= r0 < dH;
      r0 -= dH;
      r1 -= (d >> 63);
      // if (div_dbg) printf("%016I64x:%016I64x r\n", r1, r0);
      // if (div_dbg) printf("%016I64x:%016I64x - %08x\n", r1, r0, uint32_t(d) & U32_BITS(0, 30));
      dL = (d >> 28) & 1;
      uint32_t dRem = uint32_t(d) & U32_BITS(12, 28);
      if (dRem >= U32_BITS(12, 27) && dRem <= U32_BIT(28)) {
        uint64_t rr00_h, rr00_l = _umul128(r0, r0, &rr00_h);
        uint64_t rr01_l = r0 * r1;
        rr0 = rr00_l;
        rr1 = rr01_l+rr01_l+rr00_h;
        rr1 -= rr0 < r0;
        // rr0 -= r0;
        rr1 -= r1;
        // if (div_dbg) printf("%016I64x:%016I64x rr\n", rr1, rr0);
        // if (div_dbg) printf("%016I64x:%016I64x x1:x0\n", x1, uint64_t(0));
        uint64_t x1 = (x2 & ~e) << 63;
        dL = (int64_t(rr1 - x1) >= 0) ? 1 : 0;
      }
      r1 -= (r0 < dL);
      r0 -= dL;
      // if (div_dbg) printf("%016I64x:%016I64x r\n", r1, r0);
      dst.m_significand[1] = r1;
      dst.m_significand[0] = r0;
      dst.m_sign           = 0;
      dst.m_exponent      = (expA >> 1) + ((exponent_bias+1) >> 1);
    } else {
      dst.m_exponent      = inf_nan_biased_exponent;
      dst.m_sign           = 0;
      dst.m_significand[0] = 0;
      dst.m_significand[1] = qnan_bit;
    }
  } else {
    if (expA != inf_nan_biased_exponent-1) { // src.exp==zero_biased_exponent
      dst = src; // zero. Yes as strange as it sounds, according to IEEE sqrt(-0)=-0
    } else { // exp==inf_nan_biased_exponent
      if ((src.m_significand[0] | src.m_significand[1]) != 0) {
        dst = src; // NaN
      } else { // exponent_infinity
        dst.m_exponent      = inf_nan_biased_exponent;
        dst.m_sign           = 0;
        dst.m_significand[0] = 0;
        dst.m_significand[1] = src.m_sign ? qnan_bit : 0;
      }
    }
  }
}

static void eval_rsqrt_tie_break(extfloat128_t& dst, uint64_t borrow, uint64_t x0, uint64_t x1)
{
  uint64_t r0 = dst.m_significand[0];
  uint64_t r1 = dst.m_significand[1];
  if (!borrow) {
    r1 -= r0 == 0;
    r0 -= 1;
  }
  // calculate (r1:r0+0.5)^2, ignore fractional part of result
  uint64_t rr0,rr1,rr2,rr3;
  uint64_t h, l;
  uint8_t c;
  l = _umul128(r0, r0, &h);
  rr0  = l;
  rr0 += r0; c = rr0 < r0;
  rr1 = c;
  rr1 += h;  c  = rr1 < h;
  rr1 += r1; c += rr1 < r1;
  l = _umul128(r0, r1, &h);
  rr1 += l;  c += rr1 < l;
  rr1 += l;  c += rr1 < l;
  rr2 = c;
  rr2 += h;  c  = rr2 < h;
  rr2 += h;  c += rr2 < h;
  l = _umul128(r1, r1, &h);
  rr2 += l;  c += rr2 < l;
  rr3 = h + c;

  // m = (r1:r0+0.5)^2 * (x1:x0)
  // we are interested in bit[255]. In orede to find we calculate words 1, 2 and 3.
  uint64_t m1,m2,m3;
  uint64_t ha, la, hb, lb;
  _umul128(rr0, x0, &m1);
  la = _umul128(rr0, x1, &ha);
  lb = _umul128(rr1, x0, &hb);
  m1 += la; c  = m1 < la;
  m1 += lb; c += m1 < lb;
  m2 = c;
  m2 += ha; c  = m2 < ha;
  m2 += hb; c += m2 < hb;
  la = _umul128(rr1, x1, &ha);
  lb = _umul128(rr2, x0, &hb);
  m2 += la; c += m2 < la;
  m2 += lb; c += m2 < lb;
  m3 = c;
  m3 += ha;
  m3 += hb;
  m3 += rr2 * x1;
  m3 += rr3 * x0;

  uint64_t inc = m3 >> 63;
  r0 += inc;
  r1 += r0 < inc;
  dst.m_significand[0] = r0;
  dst.m_significand[1] = r1;
}

void extfloat128_t::eval_rsqrt(extfloat128_t& dst, const extfloat128_t& src)
{
  uint32_t expA = src.m_exponent - 1;
  if (expA < inf_nan_biased_exponent-1) {
    if (src.m_sign == 0) {
      const uint64_t BIT_63      = uint64_t(1) << 63;
      int e = expA & 1;
      uint64_t x2 = src.m_significand[0];
      uint64_t x3 = src.m_significand[1];
      if (x2 == 0 && x3 == BIT_63 && e == 0) {
        dst.m_significand[0] = x2;
        dst.m_significand[1] = x3;
        dst.m_sign           = 0;
        dst.m_exponent      = (exponent_bias + ((exponent_bias-1) >> 1)) - (expA >> 1);
      } else {
        int64_t  yw = rsqrt_u64(x3, e);
        uint64_t r1 = uint64_t(yw) << 1;
        uint64_t rr3, rr2, rr1;
        rr2 = _umul128(r1, r1, &rr3);
        uint64_t m1, m2, m3;
        m2 = _umul128(rr3, x3, &m3);
        uint64_t m23_h; _umul128(rr2, x3, &m23_h);
        m2 += m23_h;
        m3 += m2 < m23_h;
        uint64_t m32_h; _umul128(rr3, x2, &m32_h);
        m2 += m32_h;
        m3 += m2 < m32_h;
        int64_t m = int64_t((m3 << 44)|(m2>>20));
        int64_t d; _mul128(m<<e, yw, &d);
        uint64_t d3 = (d >> 43);
        uint64_t d2 = uint64_t(d) << (64-43);
        uint64_t r0 = 0 - d2;
        r1 -= (d2 > 0);
        r1 -= d3;

        uint64_t rr00_h, rr01_h, rr01_l, rr11_h, rr11_l, carry;
        rr01_l = _umul128(r0, r1, &rr01_h);
        rr3 = rr01_h >> 63;
        rr2 = (rr01_h << 1) | (rr01_l >> 63);
        rr1 = (rr01_l << 1);
        _umul128(r0, r0, &rr00_h);
        rr1 += rr00_h; carry = rr1 < rr00_h;
        rr2 += carry;  carry = rr2 < carry;
        rr3 += carry;
        rr11_l = _umul128(r1, r1, &rr11_h);
        rr2 += rr11_l; carry = rr2 < rr11_l;
        rr3 += rr11_h + carry;

        // m3   m2   m1
        // 33_h 33_l
        //      32_h 32_l
        //           31_h 31_l
        //      23_h 23_l
        //           22_h 22_l
        //                21_h 21_l
        uint64_t m22_h, m23_l, m31_h, m32_l, m33_l;
                _umul128(x2, rr2, &m22_h);
        m23_l = _umul128(x2, rr3, &m23_h);
        uint64_t m1a = m22_h + m23_l; carry = m1a < m23_l;
        uint64_t m2a = m23_h + carry;

                _umul128(x3, rr1, &m31_h);
        m32_l = _umul128(x3, rr2, &m32_h);
        uint64_t m1b = m31_h + m32_l; carry = m1b < m32_l;
        uint64_t m2b = m32_h + carry;
        m33_l = x3 * rr3;
        m2b += m33_l;
        m1 = m1a + m1b; carry = m1 < m1b;
        m2 = m2a + m2b + carry;
        m = int64_t((m2 << 28)|(m1 >> 36));
        _mul128(m<<e, r1>>1, &d);
        d3 = (d >> 63); // sign
        d2 = (d >> 27);
        r1 -= r0 < d2;
        r0 -= d2;
        r1 -= d3;
        uint64_t borrow = (d >> 26) & 1;
        uint32_t dRem = uint32_t(d) & U32_BITS(11, 26);
        r1 -= r0 < borrow;
        r0 -= borrow;

        dst.m_significand[1] = r1;
        dst.m_significand[0] = r0;
        dst.m_sign           = 0;
        dst.m_exponent      = (exponent_bias + ((exponent_bias-1) >> 1) - 1) - (expA >> 1);
        if (dRem >= U32_BITS(11, 25) && dRem <= U32_BIT(26))
          eval_rsqrt_tie_break(dst, borrow, x2, x3);
      }
    } else {
      dst.m_exponent      = inf_nan_biased_exponent;
      dst.m_sign           = 0;
      dst.m_significand[0] = 0;
      dst.m_significand[1] = qnan_bit;
    }
  } else {
    if (expA != inf_nan_biased_exponent-1) { // src.exp==zero_biased_exponent
      dst.m_exponent      = inf_nan_biased_exponent;
      dst.m_sign           = src.m_sign; // i decided that rsqrt(-0)=-0. It does not make much sense, but consistent with IEEE definition of sqrt(-0)
      dst.m_significand[0] = 0;
      dst.m_significand[1] = 0;
    } else { // exp==inf_nan_biased_exponent
      if ((src.m_significand[0] | src.m_significand[1]) != 0) {
        dst = src; // NaN
      } else { // exponent_infinity
        if (src.m_sign) {
          // rsqrt(-inf)=>Nan
          dst.m_exponent      = inf_nan_biased_exponent;
          dst.m_sign           = 0;
          dst.m_significand[0] = 0;
          dst.m_significand[1] = qnan_bit;
        } else {
          // rsqrt(+inf)=>zero
          dst.m_exponent      = zero_biased_exponent;
          dst.m_sign           = 0;
          dst.m_significand[0] = 0;
          dst.m_significand[1] = 0;
        }
      }
    }
  }
}

// return sign(a)*mod(abs(a), 2**exp)
extfloat128_t mod_pow2(const extfloat128_t& a, int32_t exp)
{
  uint32_t expA = a.m_exponent;
  if (expA != extfloat128_t::inf_nan_biased_exponent) {
    if (expA != extfloat128_t::zero_biased_exponent) {
      int32_t ea = expA - extfloat128_t::exponent_bias;
      if (ea < exp) {
        return a;
      } else {
        // ea >= exp
        extfloat128_t ret;
        uint32_t dExp = uint32_t(ea) - uint32_t(exp);
        if (dExp < 127) {
          uint64_t m0, m1;
          if (dExp < 63) {
            m0 = a.m_significand[0];
            m1 = a.m_significand[1] & (uint64_t(-1) >> (dExp + 1));
          } else {
            m0 = a.m_significand[0] & (uint64_t(-1) >> (dExp - 63));
            m1 = 0;
          }
          unsigned long shift;
          int64_t expA64 = expA;
          if (_BitScanReverse64(&shift, m1)) {
            m1 = (m1 << (63 - shift)) | (m0 >> (1+shift));
            m0 = (m0 << (63 - shift));
            expA64 -= 63-shift;
          } else if (_BitScanReverse64(&shift, m0)) {
            m1 = (m0 << (63 - shift));
            m0 = 0;
            expA64 -= 127-shift;
          } else {
            expA64 = extfloat128_t::zero_biased_exponent;
          }
          if (expA64 <= extfloat128_t::zero_biased_exponent) {
            expA64 = extfloat128_t::zero_biased_exponent;
            m0 = m1 = 0;
          }
          ret.m_significand[0] = m0;
          ret.m_significand[1] = m1;
          ret.m_exponent = uint32_t(expA64);
        } else {
          // zero
          ret.m_significand[0] =
          ret.m_significand[1] = 0;
          ret.m_exponent = extfloat128_t::zero_biased_exponent;
        }
        ret.m_sign = a.m_sign;
        return ret;
      }
    } else { // zero
      return a;
    }
  } else { // NaN or Inf
    extfloat128_t ret = a;
    if ((a.m_significand[0] | a.m_significand[1]) == 0) {
      ret.m_significand[1] = extfloat128_t::qnan_bit; // mod_pow2(inf, exp) = NaN
      ret.m_sign = 0;
    }
    return ret;
  }
}

// trunc
// Rounds x toward zero, returning the nearest integral value that is not larger in magnitude than x
extfloat128_t trunc(const extfloat128_t& x)
{
  extfloat128_t ret = x;
  uint32_t exp = x.m_exponent;
  if (exp >= extfloat128_t::exponent_bias) {
    uint32_t dExp = exp - extfloat128_t::exponent_bias;
    if (dExp < 127) {
      if (dExp < 64) {
        ret.m_significand[1] &= uint64_t(-1) << (63 - dExp);
        ret.m_significand[0] = 0;
      } else {
        ret.m_significand[0] &= uint64_t(-1) << (127 - dExp);
      }
    }
  } else {
    ret.m_significand[0] = ret.m_significand[1] = 0;
    ret.m_exponent = extfloat128_t::zero_biased_exponent;
  }
  return ret;
}


// round
// Returns the integral value that is nearest to x, with halfway cases rounded away from zero
extfloat128_t round(const extfloat128_t& x)
{
  extfloat128_t ret = x;
  uint32_t exp = x.m_exponent;
  if (exp >= extfloat128_t::exponent_bias - 1) {
    // abs(x) >= 0.5
    if (exp < extfloat128_t::exponent_bias + 127) {
      int32_t dExp = exp - extfloat128_t::exponent_bias;
      if (dExp >= 0) {
        ret += extfloat128_t::pow2(-1, x.m_sign);
        dExp = ret.m_exponent - extfloat128_t::exponent_bias;
        if (dExp < 64) {
          ret.m_significand[1] &= uint64_t(-1) << (63 - dExp);
          ret.m_significand[0] = 0;
        } else {
          ret.m_significand[0] &= uint64_t(-1) << (127 - dExp);
        }
      } else { // 0.5 <= abs(x) < 1
        ret = extfloat128_t::one(x.m_sign);
      }
    }
  } else {
    ret.m_significand[0] = ret.m_significand[1] = 0;
    ret.m_exponent = extfloat128_t::zero_biased_exponent;
  }
  return ret;
}


void extfloat128_t::to_17bytes(uint8_t* dst) const
{
  memcpy(dst, m_significand, sizeof(m_significand));
  dst[sizeof(m_significand)] = uint8_t(((uint32_t(_get_exponent())+64) << 1) | m_sign);
}

void extfloat128_t::from_17bytes(const uint8_t* src)
{
  memcpy(m_significand, src, sizeof(m_significand));
  uint8_t sexp = src[sizeof(m_significand)];
  m_sign = sexp & 1;
  _set_exponent(int32_t(sexp>>1)-64);
}

 extfloat128_t extfloat128_t::ulp() const
 {
   if (isfinite(*this)) {
    uint64_t exp = static_cast<uint64_t>(m_exponent) - 127;
    if (exp >= min_biased_exponent)
      return pow2(static_cast<int32_t>(exp-exponent_bias));
    else
      return pow2(min_exponent_val);
   } else {
     return *this;
   }
 }

void extfloat128_t::from_17bytes_fix(extfloat128_t dst[2], const uint8_t* src)
{
  dst[0].m_significand[0] = dst[0].m_significand[1] = 0;
  memcpy(&dst[0].m_significand[0], src, 15);
  dst[1].m_significand[0] = dst[1].m_significand[1] = 0;
  memcpy(&dst[1].m_significand[1], src+15, 2);
  for (int i = 0; i < 2; ++i) {
    uint64_t m0 = dst[i].m_significand[0];
    uint64_t m1 = dst[i].m_significand[1];
    unsigned long shift;
    if (_BitScanReverse64(&shift, m1)) {
      m1 = (m1 << (63 - shift)) | ((m0 >> 1) >> shift);
      m0 = (m0 << (63 - shift));
      dst[i].m_exponent = exponent_bias - 9 + i*56 - 63 + shift;
    } else if (_BitScanReverse64(&shift, m0)) {
      m1 = (m0 << (63 - shift));
      m0 = 0;
      dst[i].m_exponent = exponent_bias - 9 + i*56 - 127 + shift;
    } else {
      dst[i].m_exponent = zero_biased_exponent;
    }
    dst[i].m_significand[0] = m0;
    dst[i].m_significand[1] = m1;
    dst[i].m_sign = 0;
  }
}

void extfloat128_t::to_18bytes_fp(uint8_t* dst, const extfloat128_t src[2])
{
  src[0].to_17bytes(&dst[1]);
  uint64_t LSB = 0;
  uint32_t dExp = src[0].m_exponent - src[1].m_exponent;
  if (dExp < 138) {
    LSB = src[1].m_significand[1] >> 55;
    LSB = LSB >> (dExp - 128);
    LSB = (LSB+1) >> 1;
    if (LSB > 127)
      LSB = 127;
  }
  dst[0] = static_cast<uint8_t>((LSB << 1) | (src[1].m_sign & 1));
}

void extfloat128_t::from_18bytes_fp(extfloat128_t dst[2], const uint8_t* src)
{
  dst[0].from_17bytes(&src[1]);
  uint8_t b1 = src[0];
  int v1 = b1 >> 1;
  dst[1] = v1;
  dst[1].m_sign = b1 & 1;
  if (v1)
    dst[1].m_exponent += dst[0]._get_exponent() - 7 - 128;
}

extfloat128_t extfloat128_t::ldexp(const extfloat128_t& x, int exp)
{
  uint32_t srcExp = x.m_exponent;
  if (srcExp != zero_biased_exponent && srcExp != inf_nan_biased_exponent) {
    int64_t dstExp = int64_t(srcExp) + exp;
    extfloat128_t ret;
    if (dstExp <= max_biased_exponent) {
      if (dstExp >= min_biased_exponent) {
        ret.m_significand[0] = x.m_significand[0];
        ret.m_significand[1] = x.m_significand[1];
        ret.m_exponent = dstExp;
      } else {
        ret = zero();
      }
    } else {
      ret = inf();
    }
    ret.m_sign = x.m_sign;
    return ret;
  }
  return x;
}

void extfloat128_t::from_uint64(extfloat128_t* dst, uint64_t src)
{
  *dst = zero();
  unsigned long shift;
  if (_BitScanReverse64(&shift, src)) {
    dst->m_significand[1] = src << (63 - shift);
    dst->m_exponent = exponent_bias + shift;
  }
}

// *this should be non-nan
// when *this == 0 and out = false return x
// when (*this == +inf or *this == -inf) and out = true return x
// otherwise
// when out = false return next representable value in the direction of 0
// when out = true  return next representable value in the  direction away of zero
extfloat128_t extfloat128_t::nextinout_core(uint32_t out) const
{
  extfloat128_t ret = *this;
  uint32_t exp = m_exponent;
  if (exp != inf_nan_biased_exponent && exp != zero_biased_exponent) {
    if (out == 0) {
      // toward zero
      ret.m_significand[0] = m_significand[0] - 1;
      if (m_significand[0] == 0) {
        ret.m_significand[1] = m_significand[1] - 1;
        if (m_significand[1] == (uint64_t(1)<<63)) {
          ret.m_significand[1] = uint64_t(-1);
          ret.m_exponent = exp - 1;
          if (exp == 1) {
            ret.m_significand[0] = ret.m_significand[1] = 0; // zero
          }
        }
      }
    } else {
      // away from zero
      ret.m_significand[0] = m_significand[0] + 1;
      if (m_significand[0] == uint64_t(-1)) {
        ret.m_significand[1] = m_significand[1] + 1;
        if (m_significand[1] == uint64_t(-1)) {
          ret.m_significand[1] = (uint64_t(1)<<63);
          ret.m_exponent = exp + 1;
          if (exp == max_biased_exponent) {
            ret.m_significand[0] = ret.m_significand[1] = 0; // inf
          }
        }
      }
    }
  } else if (exp == zero_biased_exponent) {
    // +0 or -0
    if (out != 0) {
      // away from zero
      ret.m_exponent = min_biased_exponent;
      ret.m_significand[1] = (uint64_t(1)<<63);
    }
  } else {
    // +inf or -inf
    if (out == 0) {
      // toward zero
      ret.m_exponent = max_biased_exponent;
      ret.m_significand[0] = ret.m_significand[1] = uint64_t(-1);
    }
  }
  return ret;
}

// *this should be non-nan
// when *this == +inf and down = false return x
// when *this == -inf and down = true  return x
// otherwise
// when down = false return next representable value in the direction of +inf
// when down = true  return next representable value in the direction of -inf
extfloat128_t extfloat128_t::nextupdown_core(uint32_t down) const
{
  uint32_t exp = m_exponent;
  if (exp != zero_biased_exponent)
    return nextinout_core(m_sign ^ down ^ 1);
  // *this == 0
  extfloat128_t ret;
  ret.m_significand[0] = 0;
  ret.m_significand[1] = (uint64_t(1)<<63);
  ret.m_exponent = min_biased_exponent;
  ret.m_sign = down;
  return ret;
}

extfloat128_t nextafter(const extfloat128_t& x, const extfloat128_t& y)
{
  switch (extfloat128_t::compare_ex(x, y)) {
    case 0: return x;                    // y==x
    case 1: return x.nextupdown_core(0); // y > x, next up
    case 2: return x.nextupdown_core(1); // y < x, next down
    default:
      break;
  }
  return isnan(x) ? x : y;
}

// when x==NAN return x
// when x == +inf and down = false return x
// when x == -inf and down = true  return x
// otherwise
// when down = false return next representable value in the direction of +inf
// when down = true  return next representable value in the direction of -inf
extfloat128_t nextupdown(const extfloat128_t& x, bool down)
{
  if (isnan(x))
    return x;
  return x.nextupdown_core(down);
}

// when x==NAN return x
// when x == 0 and out = false return x
// when (x == +inf or x == -inf) and out = true return x
// otherwise
// when out = false return next representable value in the direction of 0
// when out = true  return next representable value in the  direction away of zero
extfloat128_t nextinout(const extfloat128_t& x, bool out)
{
  if (isnan(x))
    return x;
  return x.nextinout_core(out);
}
