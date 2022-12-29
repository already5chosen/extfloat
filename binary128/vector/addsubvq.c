#include <stdint.h>
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

#include "addsubvq.h"

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

static __inline
unsigned __int128 mk_u128(uint64_t hi, uint64_t lo)
{
  return ((unsigned __int128)hi << 64) | lo;
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
}

static __inline
uint64_t unsafe_shldq(uint64_t hi, uint64_t lo, unsigned lshift)  // lshift in range [1:63]
{
  return (hi << lshift) | (lo >> (64-lshift));
}

static void addsubvq(__float128* dst, const __float128* src1, const __float128* src2, int srclen, int opsub)
{
  if (srclen <= 0)
    return;

  uint64_t subq = (uint64_t)opsub << 63;

  const uint64_t MANT_H_MSK = (uint64_t)-1 >> 16;
  const uint64_t BIT_14     = (uint64_t)1  << 14;
  const uint64_t BIT_15     = (uint64_t)1  << 15;
  const uint64_t BIT_16     = (uint64_t)1  << 16;
  const uint64_t BIT_48     = (uint64_t)1  << 48;
  const uint64_t BIT_63     = (uint64_t)1  << 63;
  const uint64_t MSK_15     = BIT_15 - 1;
  const uint64_t INF_MSW    = MSK_15 << 48;

  do {
    uint64_t xHi = get_f128_hi(src1);
    uint64_t xLo = get_f128_lo(src1);
    uint64_t yHi = get_f128_hi(src2);
    uint64_t yLo = get_f128_lo(src2);
    uint64_t xHi2 = xHi*2;
    uint64_t yHi2 = yHi*2;
    yHi ^= subq;
    if (mk_u128(xHi2, xLo) < mk_u128(yHi2, yLo)) {
      { uint64_t tmp = xLo; xLo = yLo; yLo = tmp; }
      { uint64_t tmp = xHi; xHi = yHi; yHi = tmp; }
    }
    // abs(x) >= abs(y)

    unsigned exp_x = (xHi >> 48) & 0x7FFF;
    unsigned exp_y = (yHi >> 48) & 0x7FFF;
    uint64_t sub   = (xHi ^ yHi) & BIT_63;
    uint64_t sign_x = xHi & BIT_63;
    uint64_t xHiWord = xHi;
    xHi &= MANT_H_MSK;
    yHi &= MANT_H_MSK;
    if (__builtin_expect(exp_x == 0x7FFF, 0)) { // x is Inf or NaN
      if ((xLo | xHi)==0)      { // x is Inf
        if (exp_y == 0x7FFF)   { // y is Inf. It can't be NaN, because abs(NaN) > abs(Inf)
          if (sub)
            xLo = 1;             // Inf-Inf => NaN
        }
      } else {  // x is NaN
      }
      set_f128(dst, xHiWord, xLo);
      goto increment_pointers;
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
          if (sub)
            xHiWord = 0;      // sum of zeros with different signed => +zero
        }
        set_f128(dst, xHiWord, xLo);
        goto increment_pointers;
      }
    }

    // if (sub) y = - y; Implemented in almost branchless code
    sub = (int64_t)sub >> 63; // convert to all '1' or all '0'. This code is not portable by c rules, but guaranteed to work in gcc
    yLo ^= sub;
    yHi ^= sub;
    yLo -= sub;
    if (__builtin_expect_with_probability(yLo == 0,0, 1.0))
      yHi -= sub;

    unsigned delta_exp = exp_x - exp_y;
    uint64_t resHi, resLo;
    unsigned exp_res;
    if (delta_exp <= 14) {
      // align mantissa of x with y
      xHi =  shldq(xHi, xLo, delta_exp);
      xLo =  xLo << delta_exp;
      // add mantissa
      xLo += yLo;
      xHi += yHi + (xLo < yLo);

      // Normalize
      if (__builtin_expect(xHi == 0, 0)) { // MS word is fully canceled
        if (xLo == 0) { // full cancellation, result == +0
          set_f128(dst, xHi, xLo);
          goto increment_pointers;
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
          unsigned rshift = 15 - lz; // rshift in range [0:15]
          uint64_t rnd_incr = ((uint64_t)(1) << rshift) >> 1;
          xLo += rnd_incr;
          // if (__builtin_expect(xLo < rnd_incr,0))
            // xHi += 1;
          xHi += xLo < rnd_incr;
          // resLo = shrdq(xHi, xLo, rshift);
          resLo = (xLo >> rshift) | ((xHi << (rshift^63)) << 1);
          resHi = xHi >> rshift;
          uint64_t rnd_msk = rnd_incr * 2 - 1;
          if (__builtin_expect_with_probability((xLo & rnd_msk)==0, 0, 1.0)) { // a tie
            resLo &= ~(uint64_t)1; // break tie to even
          }
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
    } else {
      // align mantissa of y with x
      unsigned rshift = delta_exp % 64;
      uint64_t yG = (yLo << (rshift ^ 63)) << 1; // fraction of LS bit
      yLo = (yLo >> rshift) | ((yHi << (rshift ^ 63)) << 1);
      yHi = (int64_t)yHi >> rshift;
      if (__builtin_expect((delta_exp >= 64),0)) {
        if (__builtin_expect((delta_exp > 114),1)) {
          // abs(y) < 1/4 LSbit(abs(x)) => result = x
          set_f128(dst, xHiWord, xLo);
          goto increment_pointers;
        }
        // 64 <= delta_exp <= 114
        // Let's do a trick:
        //	Since delta_exp in [64:114] delat_exp % 64 is in [0:50]. It means that 14 LS bits of yG==0
        // On the other hand, since delta_exp > 1 we are guranteed to have no moore than 1 MS bit of mnt_hi(x) canceled during summation
        // which means that we need at most 2 exact MS bits in yG and the rest are sticky bits.
        // Then we can calculate yG = (yG >> 8) | yLo
        yG  = (yG >> 8) | yLo;
        yLo = yHi;
        yHi = (int64_t)yHi >> 63; // expend sign bit
      }
      // add mantissa
      xLo += yLo;
      xHi += yHi + (xLo < yLo);

      // Normalize
      unsigned lz = __builtin_clzll(xHi); // lz = 14 to 16
      resHi = unsafe_shldq(xHi, xLo, lz);
      resLo = unsafe_shldq(xLo, yG,  lz);
      uint64_t resG  = yG << lz;
      exp_res = exp_x + 15 - lz;

      resHi &= ~BIT_63; // clear hidden bit
      // round to nearest
      resLo += BIT_14;
      if (__builtin_expect(resLo < BIT_14,0))
        resHi += 1;
      // resHi += resLo < BIT_14;
      unsigned rnd_rem = resLo & MSK_15;
      if (__builtin_expect(rnd_rem == 0, 0)) { // a possible tie
        if (resG == 0)      // a tie, indeed
          resLo &= -BIT_16; // break tie to even
      }
      resLo = shrdq(resHi, resLo, 15);
      resHi = resHi >> 15;
    }
    // combine sign+exponent+mantissa
    resHi += (uint64_t)exp_res << 48;
    if (__builtin_expect(resHi >= INF_MSW, 0)) { // overflow
      resHi = INF_MSW; // Inf
      resLo = 0;       // Inf
    }
    set_f128(dst, resHi | sign_x, resLo);

    increment_pointers:
    ++src1;
    ++src2;
    ++dst;
  } while (--srclen);
}

void addvq(__float128* dst, const __float128* src1, const __float128* src2, int srclen)
{
  addsubvq(dst, src1, src2, srclen, 0);
}

void subvq(__float128* dst, const __float128* src1, const __float128* src2, int srclen)
{
  addsubvq(dst, src1, src2, srclen, 1);
}
