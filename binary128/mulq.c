#include <stdint.h>
#include <string.h>

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
  uint64_t lo = _mm_extract_epi64(v, 0);
  uint64_t hi = _mm_extract_epi64(v, 1);
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

  const uint64_t BIT_48       = (uint64_t)1 << 48;
  const uint64_t BIT_63       = (uint64_t)1 << 63;
  const uint64_t SIGN_BIT     = BIT_63;
  const int      EXP_BIAS     = 0x3FFF;
  const int      EXP_NAN_INF  = 0x7FFF;
  const uint64_t MSK_48       = BIT_48 - 1;

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
    if (((yHi<<16)|yLo) != 0) { // y is NaN
      return mk_f128(yHi, yLo); // return y
    }
    // y is inf
    if (xBiasedExp==EXP_NAN_INF) {
      if (((xHi<<16)|xLo) != 0)   // x is NaN
        return mk_f128(xHi, xLo); // return x
    }
    if (((xHi<<1)|xLo) == 0) {    // x is zero
      return mk_f128(yHi|1, yLo); // return inf*zero => NaN
    }
    xBiasedExp = yBiasedExp;      // cause overflow detection down the road
  }

  if (__builtin_expect(xBiasedExp==0,0)) {
    // exchange x and y
    { uint64_t tmp = xHi; xHi = yHi; yHi = tmp; }
    { uint64_t tmp = xLo; xLo = yLo; yLo = tmp; }
    xBiasedExp = yBiasedExp; yBiasedExp = 0;
  }
  if (__builtin_expect(yBiasedExp==0,0)) {
    // y is subnormal or zero
    uint64_t ySign = yHi >> 63;
    yHi &= ~SIGN_BIT;
    if (yHi == 0) {
      if (yLo == 0) {
        yLo = (uint64_t)-1;
        yBiasedExp = -(1 << 17);  // cause underflow detection down the road
      }
      do {
        yHi = yLo >> 16;
        yLo = yLo << 48;
        xBiasedExp -= 48;
      } while (yHi == 0);
    }
    if (yHi != 0) { // y is not zero - normalize
      unsigned lz = __builtin_clzll(yHi) - 15; // lz > 0
      yHi = (yHi << lz) | (yLo >> (64-lz));
      yLo = (yLo << lz);
      yBiasedExp -= lz-1;
    }
    yHi |= ySign << 63;
    // we do not care about the case when both x and y
    // are subnormals or zero
    // This case will be handled by exponent underflow
  }

  int resBiasedExp = xBiasedExp + yBiasedExp - EXP_BIAS - 1;
  uint64_t xySign = (xHi ^ yHi) & SIGN_BIT;
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
  if (__builtin_expect(resBiasedExp >= EXP_NAN_INF-1, 0)) { // overflow
    xy3 = xy2 = xy1 = 0;
    resBiasedExp = EXP_NAN_INF; // Inf
  }
  if (__builtin_expect(resBiasedExp < 0, 0)) {
    // result is subnormal or underflow (zero)
    unsigned rshift =  msbit - resBiasedExp;
    if (rshift >= 64) {
      if (rshift > 114) { // underflow
        return mk_f128(xySign, 0); // return zero
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
  return mk_f128(resHi, resLo);
}
