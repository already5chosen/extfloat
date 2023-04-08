#include <stdint.h>
#include <string.h>
#include <fenv.h>
#ifdef __SSE2__
#include <x86intrin.h>
#endif

#if defined(__amd64) && defined(__WIN64)
 #define X64_WINABI_QUIRCK 1
#else
 #define X64_WINABI_QUIRCK 0
#endif


static __inline
unsigned __int128 mk_u128(uint64_t hi, uint64_t lo)
{
  return ((unsigned __int128)hi << 64) | lo;
}

static __inline
__float128 mk_f128(uint64_t hi, uint64_t lo)
{
  #ifdef __SSE2__
  __m128i u = _mm_set_epi64x(hi, lo);
  #else
  unsigned __int128 u = ((unsigned __int128)hi << 64) | lo;
  #endif
  __float128 f;
  memcpy(&f, &u, sizeof(f));
  return f;
}

#if X64_WINABI_QUIRCK

static __inline _Float128* set_f128(_Float128* pRet, uint64_t hi, uint64_t lo)
{
  memcpy(pRet, &lo, sizeof(lo));
  memcpy((char*)pRet + sizeof(lo), &hi, sizeof(hi));
  return pRet;
}

static __inline uint64_t get_f128_lo(const _Float128* src)
{
  uint64_t u;
  memcpy(&u, src, sizeof(u));
  return u;
}

static __inline uint64_t get_f128_hi(const _Float128* src)
{
  uint64_t u;
  memcpy(&u, (const char*)src + sizeof(u), sizeof(u));
  return u;
}

#define ADDQ_CORE_RETURN(pRet, hi, lo) return set_f128(pRet, hi, lo)

#else
#define ADDQ_CORE_RETURN(pRet, hi, lo) return mk_f128(hi, lo)
#endif


static __inline
unsigned __int128 f128_to_u128(__float128 f)
{
  #ifdef __SSE2__
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

static __inline
uint64_t shrdq(uint64_t hi, uint64_t lo, unsigned rshift)
{
  return (uint64_t)(mk_u128(hi, lo) >> (rshift % 64));
}

static __inline
uint64_t shldq(uint64_t hi, uint64_t lo, unsigned lshift) // lshift in range [0:63]
{
  return (hi << lshift) | ((lo >> (lshift ^ 63)) >> 1);
  // return (uint64_t)((mk_u128(hi, lo) << (lshift % 64)) >> 64);
}

static __inline
uint64_t unsafe_shldq(uint64_t hi, uint64_t lo, unsigned lshift)  // lshift in range [1:63]
{
  return (hi << lshift) | (lo >> (64-lshift));
  // return (uint64_t)((mk_u128(hi, lo) << (lshift % 64)) >> 64);
}

static __inline uint64_t d2u(double d) {
  uint64_t u;
  memcpy(&u, &d, sizeof(u));
  return u;
}

static __inline double u2d(uint64_t u) {
  double d;
  memcpy(&d, &u, sizeof(d));
  return d;
}

static
#if X64_WINABI_QUIRCK
_Float128*
addq_core(_Float128* ret, const _Float128* pX, uint64_t  yLo, uint64_t  yHi)
#else
_Float128
addq_core(unsigned __int128 u_x, uint64_t yLo, uint64_t yHi)
#endif
{
  #if X64_WINABI_QUIRCK
  uint64_t xHi = get_f128_hi(pX);
  uint64_t xLo = get_f128_lo(pX);
  #else
  uint64_t xHi = (uint64_t)(u_x >> 64);
  uint64_t xLo = (uint64_t)(u_x);
  #endif
  if (mk_u128(xHi*2, xLo) < mk_u128(yHi*2, yLo)) {
    { uint64_t tmp = xLo; xLo = yLo; yLo = tmp; }
    { uint64_t tmp = xHi; xHi = yHi; yHi = tmp; }
  }

  const uint64_t MANT_H_MSK = (uint64_t)-1 >> 16;
  const uint64_t BIT_15     = (uint64_t)1  << 15;
  const uint64_t BIT_16     = (uint64_t)1  << 16;
  const uint64_t BIT_47     = (uint64_t)1  << 47;
  const uint64_t BIT_48     = (uint64_t)1  << 48;
  const uint64_t BIT_63     = (uint64_t)1  << 63;
  const uint64_t MSK_15     = BIT_15 - 1;
  const uint64_t MSK_16     = BIT_16 - 1;
  const uint64_t INF_MSW    = MSK_15 << 48;
  const uint64_t QNAN_BIT   = BIT_47;
  const uint64_t QNAN_MSW   = INF_MSW | QNAN_BIT;
  const uint64_t DBL_ONE_MSW = (uint64_t)0x3ff << 52;
  const uint64_t DBL_INF_MSW = (uint64_t)0x7ff << 52;

  unsigned exp_x = (xHi >> 48) & 0x7FFF;
  unsigned exp_y = (yHi >> 48) & 0x7FFF;
  uint64_t sub   = (xHi ^ yHi) & BIT_63;
  uint64_t sign_x = xHi & BIT_63;
  uint64_t xHiWord = xHi;
  xHi &= MANT_H_MSK;
  yHi &= MANT_H_MSK;
  if (__builtin_expect(exp_x == 0x7FFF, 0)) { // x is Inf or NaN
    if ((xLo | xHi)==0)      { // x is Inf
      if (exp_y == 0x7FFF)   { // y is Inf
        if (sub) {
          xHiWord = QNAN_MSW;  // Inf-Inf => QNaN
          feraiseexcept(FE_INVALID); // raise invalid operand exception
        }
      }
    } else {  // x is NaN
      if ((xHiWord & QNAN_BIT)==0) { // x is SNaN
        feraiseexcept(FE_INVALID);   // raise invalid operand exception
        xHiWord |= QNAN_BIT;         // turn x into QNaN
      } else {                       // x is QNaN
        if (exp_y == 0x7FFF) {       // y is Inf or NaN
          if (yLo | yHi)     {       // y is NaN
            if ((yHi & QNAN_BIT)==0) // y is SNaN
              feraiseexcept(FE_INVALID);   // raise invalid operand exception
          }
        }
      }
    }
    ADDQ_CORE_RETURN(ret, xHiWord, xLo);
  }

  xHi |= BIT_48;          // add hidden mantissa bit
  yHi |= BIT_48;
  if (__builtin_expect(exp_y == 0, 0)) { // y is subnormal or zero
    yHi &= MANT_H_MSK;    // remove hidden mantissa bit
    exp_y = 1;            // adjust exponent
    if (exp_x == 0)     { // x is subnormal or zero
      xHi &= MANT_H_MSK;
      exp_x = 1;
    }
    if ((yHi | yLo)==0)   { // y is Zero
      if ((xHi | xLo)==0) { // x is Zero
        // Sum of zeros is zero, but the sign of zero depends
        // both on signs of operand and, for opposite signs, on rounding mode
        // Doing the same operation in binary64 (a.k.a double) will produce a correct sign
        xHiWord = d2u(u2d(sign_x)+u2d(sign_x^sub));
      }
      ADDQ_CORE_RETURN(ret, xHiWord, xLo);
    }
  }

  unsigned delta_exp = exp_x - exp_y;
  if (delta_exp > 14) {
    // align mantissa of y with x * 2**-14
    unsigned rshift_x = delta_exp - 14;
    unsigned rshift = rshift_x % 64;
    uint64_t yG = (yLo << (rshift ^ 63)) << 1; // fraction of LS bit
    yLo = (yLo >> rshift) | ((yHi << (rshift ^ 63)) << 1);
    yHi = (yHi >> rshift);
    if (__builtin_expect((rshift_x >= 64),0)) {
      if (rshift_x > 100)
        yHi = 1;
      yG  |= yLo;
      yLo  = yHi;
      yHi  = 0;
    }
    yLo |= (yG != 0); // sticky bit
    delta_exp = 14;
    exp_y     = exp_x - 14;
  }

  // if (sub) y = - y; Implemented in almost branchless code
  sub = (int64_t)sub >> 63; // convert to all '1' or all '0'. This code is not portable by c rules, but guaranteed to work in gcc
  yLo ^= sub;
  yHi ^= sub;
  yLo -= sub;
  if (__builtin_expect_with_probability(yLo == 0,0, 1.0))
    yHi -= sub;

  uint64_t resHi, resLo;
  unsigned exp_res;

  // align mantissa of x with y
  xHi =  shldq(xHi, xLo, delta_exp);
  xLo =  xLo << delta_exp;
  // add mantissa
  xLo += yLo;
  xHi += yHi + (xLo < yLo);

  // Normalize
  if (__builtin_expect(xHi == 0, 0)) { // MS word is fully canceled
    if (xLo == 0) { // full cancellation
      // Result of full cancellation is zero, but the sign of zero
      // depends on rounding mode
      // Doing the same operation in binary64 (a.k.a double) will produce a correct sign
      resHi = d2u(u2d(DBL_ONE_MSW)+u2d(DBL_ONE_MSW^(sub&BIT_63)));
      ADDQ_CORE_RETURN(ret, resHi, xLo);
    }
    unsigned lz = __builtin_clzll(xLo);
    unsigned lshift = lz + 49;
    if (__builtin_expect((exp_y <= lshift), 0)) { // subnormal result
      lshift = exp_y - 1;
      if (lshift >= 64) {
        resHi = xLo << (lshift % 64);
        resLo = 0;
      } else {
        resHi = shldq(xHi, xLo, lshift);
        resLo = xLo << lshift;
      }
      exp_res = 0;
    } else {
      xLo <<= lz;
      resHi = xLo >> 15;
      resLo = xLo << 49;
      exp_res = exp_y - lshift - 1;
    }
  } else {
    unsigned lz = __builtin_clzll(xHi);
    if (lz <= 15) {
      // No shift or right shift. Rounding required. Subnormal result impossible
      // This case expected to be the most common in real world
      uint64_t rnd_bits = (xLo << lz) & MSK_16; // isolate bits related to rounding. LS bit of result in bit 15
      unsigned rshift = 15 - lz; // rshift in range [0:15]
      resLo = (xLo >> rshift) | ((xHi << (rshift^63)) << 1);
      resHi = xHi >> rshift;

      // Rounding
      // Round by mocking binary128 rounding with binary64 (a.k.a. double) rounding
      // Doing it this way we avoid costly reading of current rounding mode
      // It also has a desirable side effect of correctly setting Inexact flag
      const uint64_t DBL_BIAS = 2;
      uint64_t rnd_u1 = rnd_bits | (DBL_BIAS << 52) | sign_x;
      uint64_t rnd_u2 = ((DBL_BIAS+15) << 52) | sign_x;
      uint64_t rnd_sum_u = d2u(u2d(rnd_u1) + u2d(rnd_u2));
      // After binary64 rounding, LS bit of result resides in bit 0 of the rnd_sum_u
      uint64_t rnd_incr = (resLo ^ rnd_sum_u) & 1;
      if (__builtin_add_overflow(resLo, rnd_incr, &resLo))
        resHi += 1;
      // Finish rounding
      exp_res = exp_y + rshift - 1;
    } else {
      // Left shift. No need for rounding. Subnormal result possible
      unsigned lshift = lz - 15; // lz = 1 to 48
      exp_res = exp_y - lshift - 1;
      if (__builtin_expect((exp_y <= lshift), 0)) { // subnormal result
        lshift = exp_y - 1;
        exp_res = 0;
        resHi = shldq(xHi, xLo, lshift);
        resLo = xLo << lshift;
      } else {
        resHi = unsafe_shldq(xHi, xLo, lshift);
        resLo = xLo << lshift;
      }
    }
  }

  // combine sign+exponent+mantissa
  resHi += (uint64_t)exp_res << 48;
  if (__builtin_expect(resHi >= INF_MSW, 0)) { // overflow
    // Mock binary64 overflow
    uint64_t mocked_u = d2u(u2d((DBL_INF_MSW-1) | sign_x) * 2.0);
    uint64_t round_to_max_normal = mocked_u & 1;
    resHi = INF_MSW - round_to_max_normal;
    resLo = 0 -  round_to_max_normal;
  }
  ADDQ_CORE_RETURN(ret, resHi | sign_x, resLo);
}

#if X64_WINABI_QUIRCK

_Float128* __addtf3(_Float128* pRet, const _Float128* pX, const _Float128* pY)
{
  uint64_t yLo = get_f128_lo(pY);
  uint64_t yHi = get_f128_hi(pY);
  return addq_core(pRet, pX, yLo, yHi);
}

_Float128* __subtf3(_Float128* pRet, const _Float128* pX, const _Float128* pY)
{
  uint64_t yHi = get_f128_hi(pY) ^ ((uint64_t)1 << 63);
  uint64_t yLo = get_f128_lo(pY);
  return addq_core(pRet, pX, yLo, yHi);
}

#else

__float128 __addtf3(__float128 x, __float128 y)
{
  unsigned __int128 ux = f128_to_u128(x);
  unsigned __int128 uy = f128_to_u128(y);
  uint64_t yHi = (uint64_t)(uy >> 64);
  uint64_t yLo = (uint64_t)(uy);
  return addq_core(ux, yLo, yHi);
}

__float128 __subtf3(__float128 x, __float128 y)
{
  unsigned __int128 ux = f128_to_u128(x);
  unsigned __int128 uy = f128_to_u128(y);
  uint64_t yHi = (uint64_t)(uy >> 64) ^ ((uint64_t)1 << 63);
  uint64_t yLo = (uint64_t)(uy);
  return addq_core(ux, yLo, yHi);
}

#endif
