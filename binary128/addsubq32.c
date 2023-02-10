#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <fenv.h>
#ifdef __i386__
#include <x86intrin.h>
#endif


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

  QWORD0_OFFSET =  0*sizeof(uint64_t),
  QWORD1_OFFSET =  1*sizeof(uint64_t),
#else
  DWORD0_OFFSET =  3*sizeof(uint32_t),
  DWORD1_OFFSET =  2*sizeof(uint32_t),
  DWORD2_OFFSET =  1*sizeof(uint32_t),
  DWORD3_OFFSET =  0*sizeof(uint32_t),

  QWORD0_OFFSET =  1*sizeof(uint64_t),
  QWORD1_OFFSET =  0*sizeof(uint64_t),
#endif
};

static __inline
__float128*
set4_f128(__float128* dst, uint32_t w3, uint32_t w2, uint32_t w1, uint32_t w0)
{
  memcpy((char*)dst + DWORD0_OFFSET, &w0, sizeof(w0));
  memcpy((char*)dst + DWORD1_OFFSET, &w1, sizeof(w1));
  memcpy((char*)dst + DWORD2_OFFSET, &w2, sizeof(w2));
  memcpy((char*)dst + DWORD3_OFFSET, &w3, sizeof(w3));
  return dst;
}

static __inline
__float128 mk4_f128(uint32_t w3, uint32_t w2, uint32_t w1, uint32_t w0)
{
  __float128 f;
  set4_f128(&f, w3, w2, w1, w0);
  return f;
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
uint32_t shrdd(uint32_t hi, uint32_t lo, unsigned rshift)  // rshift in range [0:31]
{
  return (lo >> rshift) | ((hi << (rshift ^ 31)) << 1);
}

static __inline
uint32_t shldd(uint32_t hi, uint32_t lo, unsigned lshift) // lshift in range [0:31]
{
  return (hi << lshift) | ((lo >> (lshift ^ 31)) >> 1);
}

static __inline
uint32_t unsafe_shrdd(uint32_t hi, uint32_t lo, unsigned rshift)  // rshift in range [1:31]
{
  return (lo >> rshift) | (hi << (32-rshift));
}

static __inline
uint32_t unsafe_shldd(uint32_t hi, uint32_t lo, unsigned lshift)  // lshift in range [1:31]
{
  return (hi << lshift) | (lo >> (32-lshift));
}

static __inline
void add128(
  uint32_t* r3, uint32_t *r2, uint32_t* r1, uint32_t* r0,
  uint32_t a3, uint32_t a2, uint32_t a1, uint32_t a0,
  uint32_t b3, uint32_t b2, uint32_t b1, uint32_t b0)
{
#ifdef __i386__
  unsigned char c;
  c = _addcarry_u32(0, a0, b0, r0);
  c = _addcarry_u32(c, a1, b1, r1);
  c = _addcarry_u32(c, a2, b2, r2);
      _addcarry_u32(c, a3, b3, r3);
#else
  uint64_t aHi = ((uint64_t)a3 << 32) | a2;
  uint64_t aLo = ((uint64_t)a1 << 32) | a0;
  uint64_t bHi = ((uint64_t)b3 << 32) | b2;
  uint64_t bLo = ((uint64_t)b1 << 32) | b0;
  aLo += bLo;
  aHi += bHi + (aLo < bLo);
  *r0 = (uint32_t)aLo;
  *r1 = (uint32_t)(aLo >> 32);
  *r2 = (uint32_t)aHi;
  *r3 = (uint32_t)(aHi >> 32);
#endif
}

static __inline
void incr128by32(uint32_t* a3, uint32_t *a2, uint32_t* a1, uint32_t* a0, uint32_t b0)
{
  if (__builtin_expect_with_probability(__builtin_add_overflow(*a0, b0, a0),0, 1.0)) {
#ifdef __i386__
    unsigned char c;
    c = _addcarry_u32(0, *a1,  1, a1);
    c = _addcarry_u32(c, *a2,  0, a2);
        _addcarry_u32(c, *a3,  0, a3);
#else
    if (  __builtin_add_overflow(*a1, 1,  a1)) {
      if (__builtin_add_overflow(*a2, 1,  a2)) {
        *a3 += 1;
      }
    }
#endif
  }
}

__float128 __addtf3(__float128 x, __float128 y)
{
  const uint32_t MANT_H_MSK = (uint32_t)-1 >> 16;
  const uint32_t BIT_14     = (uint32_t)1  << 14;
  const uint32_t BIT_15     = (uint32_t)1  << 15;
  const uint32_t BIT_16     = (uint32_t)1  << 16;
  const uint32_t BIT_31     = (uint32_t)1  << 31;
  const uint32_t MSK_15     = BIT_15 - 1;
  const uint32_t EXP_MSK    = MSK_15 << 16;
  const uint32_t INF_MSW    = EXP_MSK;
  const uint32_t QNAN_BIT   = BIT_15;
  const uint32_t QNAN_MSW   = INF_MSW | QNAN_BIT;

  uint32_t yHiw = get_f128_w3(&y);
  uint32_t xHiw = get_f128_w3(&x);
  const _Float128* pX = &x;
  const _Float128* pY = &y;
  if (xHiw*2 < yHiw*2) {
    { const _Float128* tmp = pX;   pX   = pY;   pY   = tmp; }
    { uint32_t   tmp = xHiw; xHiw = yHiw; yHiw = tmp; }
  }

  if (__builtin_expect((xHiw & EXP_MSK) == INF_MSW, 0)) { // x is Inf or NaN
    uint32_t x0  = get_f128_w0(pX);
    uint32_t x1  = get_f128_w1(pX);
    uint32_t x2  = get_f128_w2(pX);
    uint32_t x3  = xHiw & MANT_H_MSK;
    if ((x0|x1|x2|x3)==0)    { // x is Inf
      if ((yHiw & EXP_MSK) == INF_MSW) { // y is Inf or NaN
        uint32_t y0  = get_f128_w0(pY);
        uint32_t y1  = get_f128_w1(pY);
        uint32_t y2  = get_f128_w2(pY);
        if (y0|y1|y2)        { // y is NaN
          if ((yHiw & QNAN_BIT)==0)    // y is SNaN
            feraiseexcept(FE_INVALID); // raise invalid operand exception
          x0 = y0;
          x1 = y1;
          x2 = y2;
          xHiw = yHiw;
        } else {  // y is Inf
          uint32_t sub = (xHiw ^ yHiw) & BIT_31;
          if (sub)
            xHiw = QNAN_MSW;  // Inf-Inf => QNaN
        }
      }
    } else {  // x is NaN
      if ((xHiw & QNAN_BIT)==0) {   // x is SNaN
        feraiseexcept(FE_INVALID);  // raise invalid operand exception
        xHiw |= QNAN_BIT;           // turn x into QNaN
      } else {                      // x is QNaN
        if ((yHiw & EXP_MSK) == INF_MSW) { // y is Inf or NaN
          uint32_t y0  = get_f128_w0(pY);
          uint32_t y1  = get_f128_w1(pY);
          uint32_t y2  = get_f128_w2(pY);
          uint32_t y3  = yHiw & MANT_H_MSK;
          if (y0|y1|y2|y3)       {       // y is NaN
            if ((yHiw & QNAN_BIT)==0)    // y is SNaN
              feraiseexcept(FE_INVALID); // raise invalid operand exception
          }
        }
      }
    }
    return mk4_f128(xHiw, x2, x1, x0);
  }

  if (__builtin_expect(xHiw*2 == yHiw*2, 0)) {
    bool swap = false;
    uint32_t x2  = get_f128_w2(pX);
    uint32_t y2  = get_f128_w2(pY);
    if (x2 < y2) {
      swap = true;
    } else if (x2 == y2) {
      uint32_t x1  = get_f128_w1(pX);
      uint32_t y1  = get_f128_w1(pY);
      if (x1 < y1) {
        swap = true;
      } else if (x1 == y1) {
        uint32_t x0  = get_f128_w0(pX);
        uint32_t y0  = get_f128_w0(pY);
        if (x0 == y0) {       // abs(x) == abs(y)
          if (xHiw != yHiw) { // x == -y => subtraction
            return mk4_f128(0, 0, 0, 0); // +0
          }
        } else if (x0 < y0) {
          swap = true;
        }
      }
    }
    if (swap) {
      { const _Float128* tmp = pX;   pX   = pY;   pY   = tmp; }
      { uint32_t   tmp = xHiw; xHiw = yHiw; yHiw = tmp; }
    }
  }

  // at this point either abs(x) > abs(y) or x==y
  unsigned exp_x = (xHiw >> 16) & MSK_15;
  unsigned exp_y = (yHiw >> 16) & MSK_15;
  uint32_t sub = (xHiw ^ yHiw) & BIT_31;
  uint32_t sign_x = xHiw & BIT_31;
  uint32_t x3 = xHiw & MANT_H_MSK;
  uint32_t y3 = yHiw & MANT_H_MSK;
  x3 |= BIT_16;        // add hidden mantissa bit
  y3 |= BIT_16;
  if (__builtin_expect(exp_y == 0, 0)) { // y is subnormal or zero
    y3 &= MANT_H_MSK;  // remove hidden mantissa bit
    exp_y = 1;         // adjust exponent
    if (exp_x == 0)  { // x is subnormal or zero
      x3 &= MANT_H_MSK;
      exp_x = 1;
    }
    uint32_t y0  = get_f128_w0(pY);
    uint32_t y1  = get_f128_w1(pY);
    uint32_t y2  = get_f128_w2(pY);
    if ((y0|y1|y2|y3)==0) { // y is Zero
      uint32_t x0 = get_f128_w0(pX);
      uint32_t x1 = get_f128_w1(pX);
      uint32_t x2 = get_f128_w2(pX);
      return mk4_f128(xHiw, x2, x1, x0);
    }
  }

  uint32_t y0  = get_f128_w0(pY);
  uint32_t y1  = get_f128_w1(pY);
  uint32_t y2  = get_f128_w2(pY);
  // if (sub) y = - y; Implemented in almost branchless code
  sub = (int32_t)sub >> 31; // convert to all '1' or all '0'. This code is not portable by C rules, but guaranteed to work in gcc
  y0 ^= sub;
  y1 ^= sub;
  y2 ^= sub;
  y3 ^= sub;
  incr128by32(&y3,&y2,&y1,&y0, sub & 1);

  uint32_t x0  = get_f128_w0(pX);
  uint32_t x1  = get_f128_w1(pX);
  uint32_t x2  = get_f128_w2(pX);

  unsigned delta_exp = exp_x - exp_y;
  unsigned exp_res;
  uint32_t res3, res2, res1, res0;
  if (delta_exp <= 14) {
    // align mantissa of x with y
    x3 =  shldd(x3, x2, delta_exp);
    x2 =  shldd(x2, x1, delta_exp);
    x1 =  shldd(x1, x0, delta_exp);
    x0 =  x0 << delta_exp;

    // add mantissa
    add128(&x3,&x2,&x1,&x0, x3,x2,x1,x0, y3,y2,y1,y0);

    // Normalize
    if (__builtin_expect(x3 == 0, 0)) { // MS word is fully canceled
      // Left shift. Result is exact (no need for rounding). Subnormal result possible
      unsigned lzw = 0;
      do {
        x3 = x2;
        x2 = x1;
        x1 = x0;
        x0 = 0;
        lzw += 32;
      } while (x3 == 0);
      unsigned lzb = __builtin_clz(x3);
      // shift MS bit into x3[31]
      x3 = shldd(x3, x2, lzb);
      x2 = shldd(x2, x1, lzb);
      x1 = x1 << lzb;
      // shift MS bit into x3[16]
      res0 = x1 << 17;
      res1 = unsafe_shrdd(x2, x1, 15);
      res2 = unsafe_shrdd(x3, x2, 15);
      res3 = x3 >> 15;
      unsigned lshift = lzw + lzb - 15;
      exp_res = exp_y - lshift - 1;
      if (__builtin_expect((exp_y <= lshift), 0)) { // subnormal result
        // de-normalize
        unsigned rshift = lshift + 1 - exp_y;
        while (rshift >= 32) {
          res0 = res1;
          res1 = res2;
          res2 = res3;
          res3 = 0;
          rshift -= 32;
        }
        res0 = shrdd(res1, res0, rshift);
        res1 = shrdd(res2, res1, rshift);
        res2 = shrdd(res3, res2, rshift);
        res3 = res3 >> rshift;
        exp_res = 0;
      }
    } else {
      unsigned lz = __builtin_clz(x3);
      if (lz <= 15) {
        // No shift or right shift. Rounding required. Subnormal result impossible
        // This case expected to be the most common in real world
        unsigned rshift = 15 - lz; // rshift in range [0:15]
        uint32_t rnd_incr = ((uint32_t)(1) << rshift) >> 1;
        incr128by32(&x3,&x2,&x1,&x0, rnd_incr);
        // Rounding done, except for ties
        res0 = shrdd(x1, x0, rshift);
        res1 = shrdd(x2, x1, rshift);
        res2 = shrdd(x3, x2, rshift);
        res3 = x3 >> rshift;
        uint32_t rnd_msk = rnd_incr * 2 - 1;
        if (__builtin_expect_with_probability((x0 & rnd_msk)==0, 0, 1.0)) { // a tie
          res0 &= ~(uint32_t)1; // break tie to even
        }
        exp_res = exp_y + rshift - 1;
      } else { // 16 <= lz <= 31
        // Left shift. Result is exact (no need for rounding). Subnormal result possible
        unsigned lshift = lz - 15; // lz = 1 to 16
        exp_res = exp_y - lshift - 1;
        if (__builtin_expect((exp_y <= lshift), 0)) { // subnormal result
          lshift = exp_y - 1;
          exp_res = 0;
          res3 = shldd(x3, x2, lshift);
          res2 = shldd(x2, x1, lshift);
          res1 = shldd(x1, x0, lshift);
          res0 = x0 << lshift;
        } else {
          res3 = unsafe_shldd(x3, x2, lshift);
          res2 = unsafe_shldd(x2, x1, lshift);
          res1 = unsafe_shldd(x1, x0, lshift);
          res0 = x0 << lshift;
        }
      }
    }
  } else { // delta_exp > 14
    if (delta_exp > 114) {
      return mk4_f128(xHiw, x2, x1, x0);
    }

    // align mantissa of y with x
    unsigned rshift = delta_exp % 32;
    uint32_t yG = (y0 << (rshift ^ 31)) << 1; // fraction of LS bit
    y0 =  shrdd(y1, y0, rshift);
    y1 =  shrdd(y2, y1, rshift);
    y2 =  shrdd(y3, y2, rshift);
    y3 =  (int32_t)y3 >> rshift;
    if (__builtin_expect((delta_exp >= 32),0)) {
      uint32_t signExt = (int32_t)y3 >> 31; // extend sign bit
      uint32_t sticky = yG;
      if (delta_exp >= 96) {
        sticky |= y1 | y0;
        yG = y2;
        y0 = y3;
        y1 = y2 = signExt;
      } else if (delta_exp >= 64) {
        sticky |= y0;
        yG = y1;
        y0 = y2;
        y1 = y3;
        y2 = signExt;
      } else { // 32 <= delta_exp < 64
        yG = y0;
        y0 = y1;
        y1 = y2;
        y2 = y3;
      }
      y3  = signExt;
      yG |= (sticky != 0);
    }

    // add mantissa
    add128(&x3,&x2,&x1,&x0, x3,x2,x1,x0, y3,y2,y1,y0);

    // Normalize
    unsigned lz = __builtin_clz(x3); // lz = 14 to 16
    res3 = unsafe_shldd(x3, x2, lz);
    res2 = unsafe_shldd(x2, x1, lz);
    res1 = unsafe_shldd(x1, x0, lz);
    res0 = unsafe_shldd(x0, yG, lz);
    uint32_t resG  = yG << lz;
    exp_res = exp_x + 15 - lz;

    res3 &= ~BIT_31; // clear hidden bit
    // round to nearest
    incr128by32(&res3,&res2,&res1,&res0, BIT_14);
    unsigned rnd_rem = res0 & MSK_15;
    if (__builtin_expect(rnd_rem == 0, 0)) { // a possible tie
      if (resG == 0)     // a tie, indeed
        res0 &= -BIT_16; // break tie to even
    }
    res0 = unsafe_shrdd(res1, res0, 15);
    res1 = unsafe_shrdd(res2, res1, 15);
    res2 = unsafe_shrdd(res3, res2, 15);
    res3 = res3 >> 15;
  }
  // combine sign+exponent+mantissa
  res3 += (uint32_t)exp_res << 16;
  if (__builtin_expect(res3 >= INF_MSW, 0)) { // overflow
    res3 = INF_MSW;         // Inf
    res2 = res1 = res0 = 0; // Inf
  }
  return mk4_f128(res3 | sign_x, res2, res1, res0);
}

__float128 __subtf3(__float128 x, __float128 y)
{
#if defined(__i386__)
  asm("xorb $128, 15+%0;\n"
  #ifdef _WIN32
      "jmp	___addtf3;"
  #else
      "jmp	__addtf3;"
  #endif
      : "=m" (y));
  __builtin_unreachable();
  return 0;
#else
  uint32_t yw3 = get_f128_w3(&y) ^ ((uint32_t)1 << 31);
  memcpy((char*)&y + DWORD3_OFFSET, &yw3, sizeof(yw3));
  return __addtf3(x, y);
#endif
}
