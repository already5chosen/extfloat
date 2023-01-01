#include <stdint.h>
#include <string.h>

#ifndef __FLOAT_WORD_ORDER__
 #error "__FLOAT_WORD_ORDER__ undefined"
#endif

#if (__FLOAT_WORD_ORDER__ != __BYTE_ORDER__)
 #error "__FLOAT_WORD_ORDER__ has to be equal to __BYTE_ORDER__"
#endif

#if (__FLOAT_WORD_ORDER__ != __ORDER_LITTLE_ENDIAN__) && (__FLOAT_WORD_ORDER__ != __ORDER_BIG_ENDIAN__)
 #error "__FLOAT_WORD_ORDER__ has to be either Little or Big Endian"
#endif

#ifdef __amd64
#ifdef __WIN64

// Under x86-64 Windows ABI _Float128 is just another 16-byte structure
typedef struct {
  char b[sizeof(_Float128)];
} soft__float128;

#define __float128 soft__float128

#else /* SYSV x86-64 ABI */

#ifdef __SSE2__
#include <x86intrin.h>
#define  CLEVER_SSE_F128_TO_INTEGER 1
#endif

#endif /* __WIN64 */
#endif /* __amd64 */

#include "mulsq.h"

static __inline
unsigned __int128 f128_to_u128(__float128 f)
{
  #ifdef CLEVER_SSE_F128_TO_INTEGER
  __m128i v;
  memcpy(&v, &f, sizeof(v));
  #ifdef  __SSE4_1__
  uint64_t lo = _mm_extract_epi64(v, 0);
  uint64_t hi = _mm_extract_epi64(v, 1);
  #else
  uint64_t lo, hi;
  memcpy(&lo, &v, sizeof(lo));
  v = _mm_unpackhi_epi64(v, v);
  memcpy(&hi, &v, sizeof(hi));
  #endif
  return ((unsigned __int128)hi << 64) | lo;
  #else
  unsigned __int128 u;
  memcpy(&u, &f, sizeof(u));
  return u;
  #endif
}

enum {
  LO_QWORD_OFFSET =
#if (__FLOAT_WORD_ORDER__ == __ORDER_LITTLE_ENDIAN__)
  0,
#else
  sizeof(uint64_t),
#endif
  HI_QWORD_OFFSET = sizeof(uint64_t) - LO_QWORD_OFFSET,
};

static __inline void set_f128(__float128* dst, uint64_t hi, uint64_t lo)
{
  memcpy((char*)dst + LO_QWORD_OFFSET, &lo, sizeof(uint64_t));
  memcpy((char*)dst + HI_QWORD_OFFSET, &hi, sizeof(uint64_t));
}

static __inline uint64_t get_f128_lo(const __float128* src)
{
  uint64_t u;
  memcpy(&u, (const char*)src + LO_QWORD_OFFSET, sizeof(u));
  return u;
}

static __inline uint64_t get_f128_hi(const __float128* src)
{
  uint64_t u;
  memcpy(&u, (const char*)src + HI_QWORD_OFFSET, sizeof(u));
  return u;
}

void mulSq(__float128* dst, const __float128* src, __float128 factor, int srclen)
{
  typedef unsigned __int128 __uint128;
  if (srclen <= 0)
    return;

  const uint64_t BIT_47      = (uint64_t)1 << 47;
  const uint64_t BIT_48      = (uint64_t)1 << 48;
  const uint64_t BIT_63      = (uint64_t)1 << 63;
  const uint64_t SIGN_BIT    = BIT_63;
  const int      EXP_BIAS    = 0x3FFF;
  const int      EXP_NAN_INF = 0x7FFF;
  const uint64_t MSK_48      = BIT_48 - 1;
  const uint64_t INF_MSW     = (uint64_t)EXP_NAN_INF << 48;
  const uint64_t QNAN_MSW    = INF_MSW | BIT_47;

  // preprocess factor
  __uint128 u_x = f128_to_u128(factor);
  uint64_t xHi = (uint64_t)(u_x >> 64);
  uint64_t xLo = (uint64_t)u_x;
  int xBiasedExp = (xHi*2) >> 49;
  if (__builtin_expect(xBiasedExp==EXP_NAN_INF,0)) {
    // x is NaN or Inf
    if (((xHi<<16)|xLo) != 0) { // x is NaN
      // dst[] = factor
      do {
        set_f128(dst, xHi, xLo);
        ++dst;
      } while (--srclen > 0);
    } else { // x is inf
      do {
        uint64_t yHi = get_f128_hi(src);
        uint64_t resHi = xHi ^ (yHi & SIGN_BIT); // common case - result=Inf with proper sign
        uint64_t resLo = 0;
        int yBiasedExp = (yHi*2) >> 49;
        if (yBiasedExp == 0) {
          uint64_t yLo = get_f128_lo(src);
          if (((yHi<<16)|yLo) == 0) { // y is zero
            resHi = QNAN_MSW;         // inf*zero => QNaN
          }
        } else if (yBiasedExp == EXP_NAN_INF) {
          uint64_t yLo = get_f128_lo(src);
          if (((yHi<<16)|yLo) != 0) { // y is NaN
            resHi = yHi; resLo = yLo; // result = y
          }
        }
        set_f128(dst, resHi, resLo);
        ++src;
        ++dst;
      } while (--srclen > 0);
    }
    return;
  }

  uint64_t xSign = xHi & SIGN_BIT;
  xHi  &= MSK_48; // isolate mantissa
  if (__builtin_expect(xBiasedExp==0,0)) { // x is subnormal or zero
    if (xHi == 0) {
      if (xLo == 0) { // x is zero
        do {
          uint64_t yHi = get_f128_hi(src);
          uint64_t resHi = (xSign ^ yHi) & SIGN_BIT; // common case - result=zero with proper sign
          uint64_t resLo = 0;
          int yBiasedExp = (yHi*2) >> 49;
          if (yBiasedExp == EXP_NAN_INF) {
            resHi = yHi;
            resLo = get_f128_lo(src);
            if (((resHi<<16)|resLo) == 0) { // y is Inf
              resHi = QNAN_MSW;             // zero*inf => QNaN
            }
          }
          set_f128(dst, resHi, resLo);
          ++src;
          ++dst;
        } while (--srclen > 0);
        return;
      }
      do {
        xHi = xLo >> 16;
        xLo = xLo << 48;
        xBiasedExp -= 48;
      } while (xHi == 0);
    }
    // x is not zero - normalize
    unsigned lz = __builtin_clzll(xHi) - 15; // lz > 0
    xHi = (xHi << lz) | (xLo >> (64-lz));
    xLo = (xLo << lz);
    xBiasedExp -= lz-1;
  }
  xHi |= BIT_48; // set hidden bit

  // Common case - factor is a finite number != 0
  do {
    uint64_t yLo = get_f128_lo(src);
    uint64_t yHi = get_f128_hi(src);
    int yBiasedExp = (yHi*2) >> 49;
    uint64_t xySign = (xSign ^ yHi) & SIGN_BIT;
    yHi &= MSK_48; // isolate mantissa
    if (__builtin_expect(yBiasedExp==EXP_NAN_INF,0)) {
      uint64_t resHi = get_f128_hi(src);
      if ((yHi | yLo) == 0) // y == Inf
        resHi ^= xSign;     // set proper sign
      set_f128(dst, resHi, yLo);
      goto increment_pointers;
    }

    // y is finite
    if (__builtin_expect(yBiasedExp==0,0)) { // y is subnormal or zero
      if (yHi == 0) {
        if (yLo == 0) { // y is zero
          set_f128(dst, xySign, 0);
          goto increment_pointers;
        }
        do {
          yHi = yLo >> 16;
          yLo = yLo << 48;
          yBiasedExp -= 48;
        } while (yHi == 0);
      }
      // x is not zero - normalize
      unsigned lz = __builtin_clzll(yHi) - 15; // lz > 0
      yHi = (yHi << lz) | (yLo >> (64-lz));
      yLo = (yLo << lz);
      yBiasedExp -= lz-1;
    }
    yHi |= BIT_48; // set hidden bit
    int resBiasedExp = xBiasedExp + yBiasedExp - EXP_BIAS - 1;

    // multiply x by y
    __uint128 mxy0 = (__uint128)xLo * yLo;
    uint64_t xy0 = (uint64_t)mxy0;
    uint64_t xy1 = (uint64_t)(mxy0>>64);
    __uint128 mxy1 = (__uint128)xHi * yLo + xy1;
     mxy1         += (__uint128)xLo * yHi; // overflow here is impossible because xH and yHi are < 2*49
     xy1          =  (uint64_t)mxy1;
     uint64_t xy2 = (uint64_t)(mxy1>>64);
    __uint128 mxy2 = (__uint128)xHi * yHi + xy2;
     xy2          =  (uint64_t)mxy2;
     uint64_t xy3 = (uint64_t)(mxy2>>64);

    // normalize and round to nearest
    unsigned msbit = (xy3 >> 33);
    resBiasedExp += (int)msbit;
    if (__builtin_expect(resBiasedExp >= EXP_NAN_INF-1, 0)) {   // overflow
      set_f128(dst, xySign | ((uint64_t)EXP_NAN_INF << 48), 0); // return Inf
      goto increment_pointers;
    }

    if (__builtin_expect(resBiasedExp < 0, 0)) {
      // result is subnormal or underflow (zero)
      unsigned rshift =  msbit - resBiasedExp;
      if (rshift >= 64) {
        if (rshift > 114) { // underflow
          set_f128(dst, xySign, 0); // return zero
          goto increment_pointers;
        }
        xy0 |= xy1;
        xy1  = xy2;
        xy2  = xy3;
        xy3  = 0;
        rshift -= 64;
      }
      if (rshift > 0) {
        xy0 =  xy0            | (xy1 << (64-rshift));
        xy1 = (xy1 >> rshift) | (xy2 << (64-rshift));
        xy2 = (xy2 >> rshift) | (xy3 << (64-rshift));
        xy3 = (xy3 >> rshift);
      }
      resBiasedExp = 0;
      msbit = 0;
    }

    uint64_t rnd_incr = (uint64_t)(msbit + 1) << 47;
    xy1 += rnd_incr;
    if (__builtin_expect(xy1 < rnd_incr, 0)) { // carry out of xy1
      xy2 += 1;
      if (xy2 < 1)
        xy3 += 1;
    }

    // Compose and return result
    int lshift = 16 - msbit;
    uint64_t resHi = (xy3 << lshift) | (xy2 >> (64-lshift));
    uint64_t resLo = (xy2 << lshift) | (xy1 >> (64-lshift));
    resHi += ((uint64_t)(unsigned)(resBiasedExp) << 48);
    resHi |= xySign;
    uint64_t rnd_msk = rnd_incr * 2 - 1;
    if (__builtin_expect_with_probability((xy1 & rnd_msk)==0, 0, 1.0)) { // possibly a tie
      if (xy0 == 0)            // indeed, a tie
        resLo &= ~(uint64_t)1; // break tie to even
    }
    set_f128(dst, resHi, resLo);
    increment_pointers:
    ++src;
    ++dst;
  } while (--srclen > 0);
}
