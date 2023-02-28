#include <stdint.h>
#include <string.h>

#include "mulq_test_generator.h"

#ifndef __FLOAT_WORD_ORDER__
 #error "__FLOAT_WORD_ORDER__ undefined"
#endif

#if (__FLOAT_WORD_ORDER__ != __BYTE_ORDER__)
 #error "__FLOAT_WORD_ORDER__ has to be equal to __BYTE_ORDER__"
#endif

#if (__FLOAT_WORD_ORDER__ != __ORDER_LITTLE_ENDIAN__) && (__FLOAT_WORD_ORDER__ != __ORDER_BIG_ENDIAN__)
 #error "__FLOAT_WORD_ORDER__ has to be either Little or Big Endian"
#endif

static
void set_float128(__float128* dst, uint64_t lo, uint64_t hi)
{
  uint64_t w[2];
#if (__FLOAT_WORD_ORDER__ == __ORDER_LITTLE_ENDIAN__)
  w[0] = lo;
  w[1] = hi;
#else
  w[1] = lo;
  w[0] = hi;
#endif
  memcpy(dst, w, sizeof(*dst));
}

static
void make_float128(__float128* dst, const uint64_t ra[2], int fp_class)
{
  const uint64_t BIT_1  = (uint64_t)1;
  const uint64_t BIT_63 = BIT_1 << 63;
  const uint64_t BIT_48 = BIT_1 << 48;
  const uint64_t MSK_48 = BIT_48 - 1;
  const uint64_t MSK_15 = (uint64_t)(-1) >> (64-15);
  uint64_t hi = ra[1];
  uint64_t lo = ra[0];
  if (fp_class >= 4) {    // 252 out of 256 are fully random, so, mostly normal
    if (fp_class < 224) { // 220 out of 252 has reduced range to prevent over/under flow
      uint64_t mant_and_sign = hi & (MSK_48 | BIT_63);
      uint64_t exp = (hi >> 48) & MSK_15;
      exp = exp/2 + MSK_15/4;
      hi = mant_and_sign | (exp << 48);
    }
  } else {
    uint64_t sign = hi & BIT_63;
    if (fp_class == 0)        { // zero
      hi = lo  = 0;
    } else if (fp_class == 1) { // subnormal
      uint64_t exp = (ra[1] >> 48) & MSK_15;
      int nbits = ((exp*112)>>15)+1; // [1:112]
      if (nbits <= 64) {
        lo = (lo >> (64-nbits)) | (BIT_1 << (nbits-1));
        hi = 0;
      } else { // [65:112]
        hi = (hi & (MSK_48 >> (113-nbits))) | (BIT_1 << (nbits-65));
      }
    } else if (fp_class == 3) { // infinity
      lo  = 0;
      hi  = (MSK_15 << 48);
    } else                    { // NaN
      hi |= (MSK_15 << 48);
    }
    hi |= sign;
  }
  set_float128(dst, lo, hi);
}

static void clear_lsbits(uint64_t* pHi, uint64_t* pLo, int nbits)
{
  if (nbits <= 48) {
    *pHi &= (uint64_t)(-1) << (48-nbits);
    *pLo = 0;
  } else { // [49:112]
    *pHi &= (uint64_t)(-1) << (112-nbits);
  }
}

void make_test_values_for_mulq(__float128 dst[2], const uint64_t random_words[5])
{
  const uint64_t BIT_1  = (uint64_t)1;
  const uint64_t BIT_63 = BIT_1 << 63;
  const uint64_t BIT_48 = BIT_1 << 48;
  const uint64_t MSK_48 = BIT_48 - 1;
  const uint64_t MSK_15 = (uint64_t)(-1) >> (64-15);
  const uint64_t MSK_14 = (uint64_t)(-1) >> (64-14);
  const uint64_t EXP_BIAS  = MSK_14;

  uint64_t class_w = random_words[0];
  uint8_t  r_class = class_w >> (8*0);
  uint8_t  x_class = class_w >> (8*1);
  uint8_t  y_class = class_w >> (8*2);
  if (r_class >= 128) { // half of the time x and y chosen independently
    make_float128(&dst[0], &random_words[1], x_class);
    make_float128(&dst[1], &random_words[3], y_class);
  } else { // x and y generated together for desired properties of result
    uint32_t prm1 = (uint16_t)class_w >> (8*1);
    uint32_t prm2 = (uint16_t)class_w >> (8*3);
    uint64_t x_lo = random_words[1];
    uint64_t x_hi = random_words[2];
    uint64_t y_lo = random_words[3];
    uint64_t y_hi = random_words[4];
    if (r_class > 0) { // reduce number of significant bits in product of mantissa
      uint32_t r_nbits = prm1*225 >> 16;
      uint32_t x_nbits = r_nbits*prm2 >> 16;
      if (x_nbits > 112) x_nbits = 112;
      uint32_t y_nbits = r_nbits - x_nbits;
      if (y_nbits > 112) y_nbits = 112;
      clear_lsbits(&x_hi, &x_lo, x_nbits);
      clear_lsbits(&y_hi, &y_lo, y_nbits);
    }
    if (r_class < 2) { // exponent of result near subnormal range
      uint32_t x_exp = (x_hi >> 48) & MSK_15;
      uint32_t y_exp = (y_hi >> 48) & MSK_15;
      uint32_t r_exp = (((x_exp + y_exp) * 128) >> 16) + EXP_BIAS - 120;
      x_exp = (r_exp * (uint64_t)x_exp) >> 15;
      y_exp = r_exp - x_exp;
      x_hi = (x_hi & (MSK_48 | BIT_63)) | ((uint64_t)x_exp << 48);
      y_hi = (y_hi & (MSK_48 | BIT_63)) | ((uint64_t)y_exp << 48);
    }
    set_float128(&dst[0], x_lo, x_hi);
    set_float128(&dst[1], y_lo, y_hi);
  }
}
