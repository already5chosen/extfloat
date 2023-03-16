#include <stdint.h>
#include <string.h>
#include <fenv.h>
#include "private/rounding_modes.h"

#ifndef __FLOAT_WORD_ORDER__
 #error "__FLOAT_WORD_ORDER__ undefined"
#endif

#if (__FLOAT_WORD_ORDER__ != __BYTE_ORDER__)
 #error "__FLOAT_WORD_ORDER__ has to be equal to __BYTE_ORDER__"
#endif

#if (__FLOAT_WORD_ORDER__ != __ORDER_LITTLE_ENDIAN__) && (__FLOAT_WORD_ORDER__ != __ORDER_BIG_ENDIAN__)
 #error "__FLOAT_WORD_ORDER__ has to be either Little or Big Endian"
#endif

enum {
#if (__FLOAT_WORD_ORDER__ == __ORDER_LITTLE_ENDIAN__)
  DWORD0_OFFSET =  0*sizeof(uint32_t),
  DWORD1_OFFSET =  1*sizeof(uint32_t),
  DWORD2_OFFSET =  2*sizeof(uint32_t),
  DWORD3_OFFSET =  3*sizeof(uint32_t),
#else
  DWORD0_OFFSET =  3*sizeof(uint32_t),
  DWORD1_OFFSET =  2*sizeof(uint32_t),
  DWORD2_OFFSET =  1*sizeof(uint32_t),
  DWORD3_OFFSET =  0*sizeof(uint32_t),
#endif
};

static __inline void set4_f128(__float128* dst, uint32_t w3, uint32_t w2, uint32_t w1, uint32_t w0)
{
  memcpy((char*)dst + DWORD0_OFFSET, &w0, sizeof(w0));
  memcpy((char*)dst + DWORD1_OFFSET, &w1, sizeof(w1));
  memcpy((char*)dst + DWORD2_OFFSET, &w2, sizeof(w2));
  memcpy((char*)dst + DWORD3_OFFSET, &w3, sizeof(w3));
}

static __inline uint32_t get_f128_w0(const __float128* src) {
  uint32_t u;
  memcpy(&u, (const char*)src + DWORD0_OFFSET, sizeof(u));
  return u;
}
static __inline uint32_t get_f128_w1(const __float128* src) {
  uint32_t u;
  memcpy(&u, (const char*)src + DWORD1_OFFSET, sizeof(u));
  return u;
}
static __inline uint32_t get_f128_w2(const __float128* src) {
  uint32_t u;
  memcpy(&u, (const char*)src + DWORD2_OFFSET, sizeof(u));
  return u;
}
static __inline uint32_t get_f128_w3(const __float128* src) {
  uint32_t u;
  memcpy(&u, (const char*)src + DWORD3_OFFSET, sizeof(u));
  return u;
}

static __inline
__float128 mk4_f128(uint32_t w3, uint32_t w2, uint32_t w1, uint32_t w0)
{
  __float128 f;
  set4_f128(&f, w3, w2, w1, w0);
  return f;
}

__float128 __multf3(__float128 srcx, __float128 srcy)
{
  // return srcx * srcy;
  const int      EXP_BIAS     = 0x3FFF;
  const int      EXP_NAN_INF  = 0x7FFF;
  const uint32_t BIT_31       = (uint32_t)1 << 31;
  const uint32_t BIT_30       = (uint32_t)1 << 30;
  const uint32_t BIT_16       = (uint32_t)1 << 16;
  const uint32_t BIT_15       = (uint32_t)1 << 15;
  const uint32_t MSK_16       = BIT_16 - 1;
  // MS Word definitions
  const uint32_t SIGN_BIT     = BIT_31;
  const uint32_t INF_MSW      = (uint32_t)EXP_NAN_INF << 16;
  const uint32_t QNAN_BIT     = BIT_15;
  const uint32_t QNAN_MSW     = INF_MSW | QNAN_BIT;
  const uint32_t SNAN_MSW     = INF_MSW;
  const uint32_t MIN_NORMAL_MSW = BIT_16;

  uint32_t xHi = get_f128_w3(&srcx);
  uint32_t yHi = get_f128_w3(&srcy);
  int xBiasedExp  = (xHi*2) >> 17;
  int yBiasedExp  = (yHi*2) >> 17;
  const __float128 *pSrcx = &srcx;
  const __float128 *pSrcy = &srcy;
  if (__builtin_expect(xBiasedExp==EXP_NAN_INF,0)) {
    // exchange x and y
    pSrcx = &srcy;
    pSrcy = &srcx;
    xBiasedExp = yBiasedExp;
    yBiasedExp = EXP_NAN_INF;
  }
  if (__builtin_expect(yBiasedExp==EXP_NAN_INF,0)) {
    // y is NaN or Inf
    uint32_t x0  = get_f128_w0(pSrcx);
    uint32_t x1  = get_f128_w1(pSrcx);
    uint32_t x2  = get_f128_w2(pSrcx);
    uint32_t xHi = get_f128_w3(pSrcx);

    uint32_t y0  = get_f128_w0(pSrcy);
    uint32_t y1  = get_f128_w1(pSrcy);
    uint32_t y2  = get_f128_w2(pSrcy);
    uint32_t yHi = get_f128_w3(pSrcy);

    if (((yHi<<16)|y0|y1|y2) != 0)  { // y is NaN
      if ((yHi & QNAN_BIT)==0)      { // y is SNaN
        feraiseexcept(FE_INVALID);    // raise invalid operand exception
        yHi |= QNAN_BIT;              // turn y into QNaN
      } else {                        // y is QNaN
        if ((xHi & QNAN_MSW)==SNAN_MSW) { // x is SNaN or Inf
          if (((xHi<<16)|x0|x1|x2) != 0)  // x is SNaN
            feraiseexcept(FE_INVALID);    // raise invalid operand exception
        }
      }
      return mk4_f128(yHi, y2, y1, y0); // return y
    }

    // y is inf
    if (xBiasedExp==EXP_NAN_INF) {
      if (((xHi<<16)|x0|x1|x2) != 0)  { // x is NaN
        if ((xHi & QNAN_BIT)==0)      { // x is SNaN
          feraiseexcept(FE_INVALID);    // raise invalid operand exception
          xHi |= QNAN_BIT;              // turn x into QNaN
        }
        return mk4_f128(xHi, x2, x1, x0); // return x
      }
    }

    if (((xHi<<1)|x0|x1|x2) == 0) {    // x is zero
      return mk4_f128(QNAN_MSW,0,0,0); // return inf*zero => QNaN
    }

    // x is neither zero nor NaN
    uint32_t msw = ((yHi^xHi) & SIGN_BIT) | INF_MSW;
    return mk4_f128(msw, 0, 0, 0);     // return INF with combined sign
  }

  __float128 normalizedSrc;
  if (__builtin_expect(xBiasedExp==0,0)) {
    // exchange x and y
    pSrcx = &srcy;
    pSrcy = &srcx;
    xBiasedExp = yBiasedExp;
    yBiasedExp = 0;
  }
  if (__builtin_expect(yBiasedExp==0,0)) {
    // y is subnormal or zero
    uint32_t y0  = get_f128_w0(pSrcy);
    uint32_t y1  = get_f128_w1(pSrcy);
    uint32_t y2  = get_f128_w2(pSrcy);
    uint32_t yHi = get_f128_w3(pSrcy);
    uint32_t ySign = yHi & SIGN_BIT;
    yHi &= ~SIGN_BIT;
    if (yHi == 0) {
      if (y2 == 0) {
        if ((y0 | y1)==0) { // zero
          uint32_t msw = (get_f128_w3(pSrcx) & SIGN_BIT) ^ ySign;
          return mk4_f128(msw, 0, 0, 0); // return zero with combined sign
        }
        // non-zero
        do {
          y2  = y1;
          y1  = y0;
          y0  = 0;
          yBiasedExp -= 32;
        } while (y2 == 0);
      }
      // y2 != 0
      unsigned lz = __builtin_clz(y2);
      yBiasedExp -= lz + 16;
      if (lz > 0) {
        y2  = (y2  << lz) | (y1  >> (32-lz));
        y1  = (y1  << lz) | (y0  >> (32-lz));
        y0  = (y0  << lz);
      }
      yHi =               (y2  >> (32-17));
      y2  = (y2  << 17) | (y1  >> (32-17));
      y1  = (y1  << 17) | (y0  >> (32-17));
      y0  = (y0  << 17);
    } else { // (yHi != 0)
      unsigned lz = __builtin_clz(yHi) - 15; // lz > 0
      yHi = (yHi << lz) | (y2  >> (32-lz));
      y2  = (y2  << lz) | (y1  >> (32-lz));
      y1  = (y1  << lz) | (y0  >> (32-lz));
      y0  = (y0  << lz);
      yBiasedExp -= lz-1;
    }
    yHi |= ySign; // restore sign
    set4_f128(&normalizedSrc, yHi,y2,y1,y0);
    pSrcy = &normalizedSrc;

    if (xBiasedExp + yBiasedExp < EXP_BIAS - 115) {
      // Exponent underflow
      xHi = get_f128_w3(pSrcx);
      if ((xHi<<1) == 0) {
        uint32_t x0  = get_f128_w0(pSrcx);
        uint32_t x1  = get_f128_w1(pSrcx);
        uint32_t x2  = get_f128_w2(pSrcx);
        if ((x0|x1|x2) == 0) {    // x is zero
          uint32_t msw = (xHi & SIGN_BIT) ^ ySign;
          return mk4_f128(msw, 0, 0, 0); // return zero with combined sign
        }
      }
      xBiasedExp = 127; // cause detection of Exponent underflow down the road
    }
  }

  int resBiasedExp = xBiasedExp + yBiasedExp - EXP_BIAS - 1;
  xHi = get_f128_w3(pSrcx);
  yHi = get_f128_w3(pSrcy);
  uint32_t xySign = (xHi ^ yHi) & SIGN_BIT;
  uint32_t x3 = (xHi & MSK_16) | BIT_16; // isolate mantissa and set hidden bit
  uint32_t y3 = (yHi & MSK_16) | BIT_16;

  uint32_t x0  = get_f128_w0(pSrcx);
  uint32_t x1  = get_f128_w1(pSrcx);
  uint32_t x2  = get_f128_w2(pSrcx);

  uint32_t y0  = get_f128_w0(pSrcy);
  uint32_t y1  = get_f128_w1(pSrcy);
  uint32_t y2  = get_f128_w2(pSrcy);

  // multiply x by y
  uint64_t x0y0 = (uint64_t)x0 * y0;
  uint32_t xy0  = (uint32_t)x0y0;

  uint64_t x1y0 = (uint64_t)x1 * y0;
  uint64_t x0y1 = (uint64_t)x0 * y1;
  x1y0 += (uint32_t)(x0y0 >> 32);
  x0y1 += (uint32_t)(x1y0);
  uint32_t xy1  = (uint32_t)x0y1;
  uint64_t acc2 = (uint64_t)(uint32_t)(x1y0 >> 32) + (uint32_t)(x0y1 >> 32);

  uint64_t x2y0 = (uint64_t)x2 * y0;
  uint64_t x1y1 = (uint64_t)x1 * y1;
  uint64_t x0y2 = (uint64_t)x0 * y2;
  acc2 += x2y0;
  x1y1 += (uint32_t)acc2;
  uint64_t acc3 = (uint64_t)(uint32_t)(x1y1 >> 32) + (uint32_t)(acc2 >> 32);
  x0y2 += (uint32_t)x1y1;
  uint32_t xy2  = (uint32_t)x0y2;
  acc3 += (uint32_t)(x0y2 >> 32);

  uint64_t x3y0 = (uint64_t)x3 * y0; // < 2**49
  uint64_t x2y1 = (uint64_t)x2 * y1;
  uint64_t x1y2 = (uint64_t)x1 * y2;
  uint64_t x0y3 = (uint64_t)x0 * y3; // < 2**49
  acc3 += x3y0;
  acc3 += x0y3;
  x2y1 += (uint32_t)acc3;
  uint64_t acc4 = (uint64_t)(uint32_t)(x2y1 >> 32) + (uint32_t)(acc3 >> 32);
  x1y2 += (uint32_t)x2y1;
  uint32_t xy3  = (uint32_t)x1y2;
  acc4 += (uint32_t)(x1y2 >> 32);

  uint64_t x3y1 = (uint64_t)x3 * y1; // < 2**49
  uint64_t x2y2 = (uint64_t)x2 * y2;
  uint64_t x1y3 = (uint64_t)x1 * y3; // < 2**49
  acc4 += x3y1;
  acc4 += x1y3;
  x2y2 += (uint32_t)acc4;
  uint32_t xy4  = (uint32_t)x2y2;
  uint64_t acc5 = (uint64_t)(uint32_t)(x2y2 >> 32) + (uint32_t)(acc4 >> 32);

  uint64_t x3y2 = (uint64_t)x3 * y2; // < 2**49
  uint64_t x2y3 = (uint64_t)x2 * y3; // < 2**49
  acc5 += x3y2;
  acc5 += x2y3;
  uint32_t xy5  = (uint32_t)acc5;

  uint64_t x3y3 = (uint64_t)x3 * y3; // < 2**34
  x3y3 += (uint32_t)(acc5 >> 32);
  uint32_t xy6  = (uint32_t)x3y3;
  uint32_t xy7  = (uint32_t)(x3y3 >> 32);
  // printf("%08x\n", xy7);
  // printf("%08x\n", xy6);
  // printf("%08x\n", xy5);
  // printf("%08x\n", xy4);
  // printf("%08x\n", xy3);
  // printf("%08x\n", xy2);
  // printf("%08x\n", xy1);
  // printf("%08x\n", xy0);
  xy1 |= xy0; // sticky bits
  xy2 |= xy1; // sticky bits

  // normalize and round to nearest
  unsigned msbit = (xy7 >> 1);
  resBiasedExp += (int)msbit;
  int rm = fast_fegetround();
  if (rm != fast_FE_TONEAREST) {
    if (xySign) {
      if (rm == fast_FE_DOWNWARD)
        rm = fast_FE_UPWARD;
      else if (rm == fast_FE_UPWARD)
        rm = fast_FE_DOWNWARD;
    }
    // At this point rm==fast_FE_UPWARD means 'away from zero', rm==fast_FE_DOWNWARD means 'toward zero'
  }
  uint32_t rnd_incr = (uint32_t)(msbit + 1) << 15;
  uint32_t rnd_msk = rnd_incr*2 - 1;

  // conditionally raise Inexact exception
  #if 0
  if ((xy2 | (xy3 & rnd_msk)) != 0)
    feraiseexcept(FE_INEXACT); // raise Inexact exception
  #else
  static const double inexact_exception_tab[] = { 0, 0x1.0p-64 };
  volatile double dummy = inexact_exception_tab[(xy2 | (xy3 & rnd_msk)) != 0] + 1.0;
  (void)dummy;
  #endif

  if (rm != fast_FE_TONEAREST) {
    xy3 |= (xy2 != 0);
    xy2 = (xy3 & rnd_msk);
    rnd_incr = (rm == fast_FE_UPWARD) ? rnd_msk : 0;
  }
  xy3 += rnd_incr;
  if (__builtin_expect(xy3 < rnd_incr, 0)) { // carry out of xy3
    xy4 += 1;
    if (xy4 < 1) {
      xy5 += 1;
      if (xy5 < 1) {
        xy6 += 1;
        if (xy6 < 1) {
          xy7 += 1;
        }
      }
    }
  }

  // Compose and return result
  int lshift = 16 - msbit;
  uint32_t res3 = (xy7 << lshift) | (xy6 >> (32-lshift));
  uint32_t res2 = (xy6 << lshift) | (xy5 >> (32-lshift));
  uint32_t res1 = (xy5 << lshift) | (xy4 >> (32-lshift));
  uint32_t res0 = (xy4 << lshift) | (xy3 >> (32-lshift));
  res3 += ((uint32_t)(unsigned)(resBiasedExp) << 16);

  if (__builtin_expect_with_probability(res3-MIN_NORMAL_MSW >= INF_MSW-MIN_NORMAL_MSW, 0, 1.0)) {
    // exponent overflow or underflow
    res3 -= MIN_NORMAL_MSW;
    if (res3 >= BIT_31+BIT_30) { // result is subnormal or zero
      unsigned rshift = (1<<16) - (unsigned)(res3 >> 16);
      res3 = (res3 & MSK_16) | BIT_16;
      // Undo "normal" rounding
      if (rnd_incr) {
        uint32_t rnd_msk = (rnd_incr & 1) ? rnd_incr : rnd_incr * 2 - 1;
        if ((xy3 & rnd_msk) < rnd_incr) {
          if (res0 == 0) {
            if (res1 == 0) {
              if (res2 == 0) {
                res3 -= 1;
              }
              res2 -= 1;
            }
            res1 -= 1;
          }
          res0 -= 1;
        }
        xy2 |= (xy3 - rnd_incr) & rnd_msk; // sticky bits
      }

      // de-normalize and round again
      if (rshift >= 32) {
        if (rshift >= 115) // underflow
          rshift = 115;
        while (rshift > 32) {
          xy2  |= res0;        // sticky bits
          res0  = res1;
          res1  = res2;
          res2  = res3;
          res3  = 0;
          rshift -= 32;
        }
        if (rshift == 32) {
          xy2  |= (res0 << 1); // sticky bits
          res0  = (res0 >> 31) | (res1 << 1);
          res1  = (res1 >> 31) | (res2 << 1);
          res2  = (res2 >> 31) | (res3 << 1);
          res3  = (res3 >> 31);
          rshift = 1;
        }
      }
      // rshift in [1:31]
      rnd_incr = (uint32_t)1 << (rshift-1);
      uint32_t rnd_msk = rnd_incr * 2 - 1;
      if (rm != fast_FE_TONEAREST) {
        res0 |= (xy2 != 0);
        rnd_incr = (rm==fast_FE_UPWARD) ? rnd_msk: 0;
      }
      if (__builtin_add_overflow(res0, rnd_incr, &res0)) {
        if (  __builtin_add_overflow(res1, 1, &res1))
          if (__builtin_add_overflow(res2, 1, &res2))
            res3 += 1;
      }
      xy3 = res0;
      res0 = (res0 >> rshift) | (res1 << (32-rshift));
      res1 = (res1 >> rshift) | (res2 << (32-rshift));
      res2 = (res2 >> rshift) | (res3 << (32-rshift));
      res3 = (res3 >> rshift);
      if (xy2 | ((xy3-rnd_incr) & rnd_msk))
        feraiseexcept(FE_UNDERFLOW | FE_INEXACT); // raise Underflow+Inexact exception
    } else { // Overflow
      feraiseexcept(FE_OVERFLOW | FE_INEXACT); // raise Overflow+Inexact exception
      res3 = INF_MSW;
      res0 = res1 = res2 = 0;
      if (rm != fast_FE_TONEAREST && rm != fast_FE_UPWARD) {
        // in case of rounding toward zero return maximal normal value
        res0 = res1 = res2 = (uint32_t)-1;
        res3 += res0;
      }
    }
  }
  res3 |= xySign;
  if (rm == fast_FE_TONEAREST) {
    uint32_t rnd_msk = rnd_incr * 2 - 1;
    if (__builtin_expect_with_probability((xy3 & rnd_msk)==0, 0, 1.0)) { // possibly a tie
      if (xy2 == 0)           // indeed, a tie
        res0 &= ~(uint32_t)1; // break tie to even
    }
  }
  return mk4_f128(res3, res2, res1, res0);
}
