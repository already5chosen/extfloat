#include <stdint.h>
#include <string.h>
#include <fenv.h>
#include "private/rounding_modes.h"

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

static __inline
__float128 mk_f128(uint64_t hi, uint64_t lo)
{
  #ifdef CLEVER_SSE_F128_TO_INTEGER
  __m128i u = _mm_set_epi64x(hi, lo);
  #else
  unsigned __int128 u = ((unsigned __int128)hi << 64) | lo;
  #endif
  __float128 f;
  memcpy(&f, &u, sizeof(f));
  return f;
}

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

__float128 __multf3(__float128 srcx, __float128 srcy)
{
  // return srcx * srcy;
  typedef unsigned __int128 __uint128;
  __uint128 u_x = f128_to_u128(srcx);
  __uint128 u_y = f128_to_u128(srcy);
  uint64_t xHi = (uint64_t)(u_x >> 64);
  uint64_t xLo = (uint64_t)u_x;
  uint64_t yHi = (uint64_t)(u_y >> 64);
  uint64_t yLo = (uint64_t)u_y;

  const uint64_t BIT_47       = (uint64_t)1 << 47;
  const uint64_t BIT_48       = (uint64_t)1 << 48;
  const uint64_t BIT_62       = (uint64_t)1 << 62;
  const uint64_t BIT_63       = (uint64_t)1 << 63;
  const uint64_t SIGN_BIT     = BIT_63;
  const int      EXP_BIAS     = 0x3FFF;
  const int      EXP_NAN_INF  = 0x7FFF;
  const uint64_t MSK_48       = BIT_48 - 1;
  const uint64_t INF_MSW      = (uint64_t)EXP_NAN_INF << 48;
  const uint64_t QNAN_BIT     = BIT_47;
  const uint64_t QNAN_MSW     = INF_MSW | QNAN_BIT;
  const uint64_t SNAN_MSW     = INF_MSW;
  const uint64_t MIN_NORMAL_MSW = BIT_48;

  int xBiasedExp  = (xHi*2) >> 49;
  int yBiasedExp  = (yHi*2) >> 49;
  if (__builtin_expect(xBiasedExp==EXP_NAN_INF,0)) {
    // exchange x and y
    { uint64_t tmp = xHi; xHi = yHi; yHi = tmp; }
    { uint64_t tmp = xLo; xLo = yLo; yLo = tmp; }
    xBiasedExp = yBiasedExp; yBiasedExp = EXP_NAN_INF;
  }
  if (__builtin_expect(yBiasedExp==EXP_NAN_INF,0)) {
    // y is NaN or Inf
    if (((yHi<<16)|yLo) != 0)  { // y is NaN
      if ((yHi & QNAN_BIT)==0) { // y is SNaN
        feraiseexcept(FE_INVALID); // raise invalid operand exception
        yHi |= QNAN_BIT;         // turn y into QNaN
      } else {                   // y is QNaN
        if ((xHi & QNAN_MSW)==SNAN_MSW) { // x is SNaN or Inf
          if (((xHi<<16)|xLo) != 0)       // x is SNaN
            feraiseexcept(FE_INVALID);    // raise invalid operand exception
        }
      }
      return mk_f128(yHi, yLo); // return y
    }
    // y is inf
    if (xBiasedExp==EXP_NAN_INF) {
      if (((xHi<<16)|xLo) != 0)  { // x is NaN
        if ((xHi & QNAN_BIT)==0) { // x is SNaN
          feraiseexcept(FE_INVALID); // raise invalid operand exception
          xHi |= QNAN_BIT;         // turn x into QNaN
        }
        return mk_f128(xHi, xLo); // return x
      }
    }
    if (((xHi<<1)|xLo) == 0) {    // x is zero
      return mk_f128(QNAN_MSW, yLo); // return inf*zero => QNaN
    }

    // x is normal, subnormal or Inf
    // return Inf with correct sign
    return mk_f128(yHi ^ (xHi & SIGN_BIT), yLo);
  }

  if (__builtin_expect(xBiasedExp==0,0)) {
    // exchange x and y
    { uint64_t tmp = xHi; xHi = yHi; yHi = tmp; }
    { uint64_t tmp = xLo; xLo = yLo; yLo = tmp; }
    xBiasedExp = yBiasedExp; yBiasedExp = 0;
  }
  uint64_t xySign = (xHi ^ yHi) & SIGN_BIT;
  if (__builtin_expect(yBiasedExp==0,0)) {
    // y is subnormal or zero
    uint64_t ySign = yHi >> 63;
    yHi &= ~SIGN_BIT;
    if (yHi == 0) {
      if (yLo == 0) { // y is zero
        return mk_f128(xySign, yLo);  // return zero with correct sign
      }
      do {
        yHi = yLo >> 16;
        yLo = yLo << 48;
        xBiasedExp -= 48;
      } while (yHi == 0);
    }
    unsigned lz = __builtin_clzll(yHi) - 15; // lz > 0
    yHi = (yHi << lz) | (yLo >> (64-lz));
    yLo = (yLo << lz);
    yBiasedExp -= lz-1;
    yHi |= ySign << 63;

    if (xBiasedExp + yBiasedExp < EXP_BIAS - 115) {
      // Exponent underflow
      if (((xHi*2) | xLo) == 0)      // x is zero
        return mk_f128(xySign, xLo); // return zero with correct sign

      xBiasedExp = 127; // cause detection of Exponent underflow down the road
    }
  }

  int resBiasedExp = xBiasedExp + yBiasedExp - EXP_BIAS - 1;
  xHi = (xHi & MSK_48) | BIT_48; // isolate mantissa and set hidden bit
  yHi = (yHi & MSK_48) | BIT_48;

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
  uint64_t rnd_incr = (uint64_t)(msbit + 1) << 47;
  if (rm != fast_FE_TONEAREST) {
    uint64_t rnd_msk = rnd_incr*2 - 1;
    xy1 |= (xy0 != 0);
    xy0 = (xy1 & rnd_msk);
    rnd_incr = (rm == fast_FE_UPWARD) ? rnd_msk : 0;
  }
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
  if (__builtin_expect_with_probability(resHi-MIN_NORMAL_MSW >= INF_MSW-MIN_NORMAL_MSW, 0, 1.0)) {
    // exponent overflow or underflow
    resHi -= MIN_NORMAL_MSW;
    if (resHi >= BIT_63+BIT_62) { // result is subnormal or zero
      unsigned rshift = (1<<16) - (unsigned)(resHi >> 48);
      resHi = (resHi & MSK_48) | BIT_48;
      // Undo "normal" rounding
      if (rnd_incr) {
        uint64_t rnd_msk = (rnd_incr & 1) ? rnd_incr : rnd_incr * 2 - 1;
        if ((xy1 & rnd_msk) < rnd_incr) {
          if (resLo == 0)
            resHi -= 1;
          resLo -= 1;
        }
        xy0 |= (xy1 - rnd_incr) & rnd_msk; // sticky bits
      }

      // de-normalize and round again
      if (rshift >= 64) {
        if (rshift >= 115)
          rshift = 115;
        xy0  |= (resLo << 1);        // sticky bits
        resLo = (resLo >> 63) | (resHi << 1);
        resHi = (resHi >> 63);
        rshift -= 63; // [64:114] => [1:51]
      }
      // rshift in [1:63]
      rnd_incr = (uint64_t)1 << (rshift-1);
      uint64_t rnd_msk = rnd_incr * 2 - 1;
      if (rm != fast_FE_TONEAREST) {
        if (rm==fast_FE_UPWARD) {
          resLo |= (xy0 != 0);
          rnd_incr = rnd_msk;
        } else {
          rnd_incr = 0;
        }
      }
      if (__builtin_add_overflow(resLo, rnd_incr, &resLo))
        resHi += 1;
      xy1 = resLo;
      resLo = (resLo >> rshift) | (resHi << (64-rshift));
      resHi = (resHi >> rshift);
      if (xy0 | ((xy1-rnd_incr) & rnd_msk))
        feraiseexcept(FE_UNDERFLOW | FE_INEXACT); // raise Underflow+Inexact exception
    } else { // Overflow
      feraiseexcept(FE_OVERFLOW | FE_INEXACT); // raise Overflow+Inexact exception
      resHi = INF_MSW;
      resLo = 0;
      if (rm != fast_FE_TONEAREST && rm != fast_FE_UPWARD) {
        // in case of rounding toward zero return maximal normal value
        resLo = (uint64_t)-1;
        resHi += resLo;
      }
    }
  }
  if (rm == fast_FE_TONEAREST) {
    // test for tie
    uint64_t rnd_msk = rnd_incr * 2 - 1;
    if (__builtin_expect_with_probability((xy1 & rnd_msk)==0, 0, 1.0)) { // possibly a tie
      if (xy0 == 0)            // indeed, a tie
        resLo &= ~(uint64_t)1; // break tie to even
    }
  }
  resHi |= xySign;
  return mk_f128(resHi, resLo);
}
