#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <intrin.h>

#include "extfloat128.h"


// extended precision type - mostly for trigs and other special functions
// The type intended for intermediate values, mostly accumulators
// It has very wide exponent range and 2-3 times wider mantissa than extfloat128_t
// It does not support NaN.
// It does not support signed zeros except for export of extfloat128_t
// It does not support overflow except for export of extfloat128_t
// Overflow on arithmetic operations is mostly not checked

// extfloat128_t::acc_t initialization
void extfloat128_t::acc_t::doNormalizeSign()
{
  const uint64_t BIT_63 = uint64_t(1)<<63;
  uint64_t y5 = m_significand[5];
  if (y5 & BIT_63) {
    // significand is negative - negate
    int i;
    for (i = 0; m_significand[i] == 0; ++i) ;
    m_significand[i] = 0-m_significand[i];
    ++i;
    for (; i < 6; ++i)
      m_significand[i] = ~m_significand[i];
    m_sign ^= 1;
  }
}
void extfloat128_t::acc_t::doNormalize()
{
  m_isNormal = true;
  if (m_exponent != 0) {
    doNormalizeSign();
    const uint64_t BIT_63 = uint64_t(1)<<63;
    uint64_t y5 = m_significand[5];
    if (y5 != 0) {
      // shift to the right
      uint64_t y4 = m_significand[4];
      uint64_t y3 = m_significand[3];
      uint64_t y2 = m_significand[2];
      uint64_t y1 = m_significand[1];
      uint64_t y0 = m_significand[0];
      unsigned long shift;
      _BitScanReverse64(&shift, y5);
      m_significand[4] = (y5 << (63-shift)) | (y4 >> (shift+1));
      m_significand[3] = (y4 << (63-shift)) | (y3 >> (shift+1));
      m_significand[2] = (y3 << (63-shift)) | (y2 >> (shift+1));
      m_significand[1] = (y2 << (63-shift)) | (y1 >> (shift+1));
      m_significand[0] = (y1 << (63-shift)) | (y0 >> (shift+1));
      m_significand[5] = 0;
      m_exponent += shift + 1;
    } else {
      uint64_t y4 = m_significand[4];
      if ((y4 & BIT_63)==0) {
        // shift to the left
        uint64_t y3 = m_significand[3];
        uint64_t y2 = m_significand[2];
        uint64_t y1 = m_significand[1];
        uint64_t y0 = m_significand[0];
        uint64_t e  = m_exponent;
        // word-level shifts
        if (y4 == 0) {
          if ((y3 | y2 | y1 | y0)==0) {
            m_exponent = 0;
            return;
          }
          do {
            y4 = y3;
            y3 = y2;
            y2 = y1;
            y1 = y0;
            y0 = 0;
            e -= 64;
          } while (y4==0);
        }
        if ((y4 & BIT_63)==0) {
          unsigned long shift;
          _BitScanReverse64(&shift, y4);
          y4 = (y4 << (63-shift)) | (y3 >> (shift+1));
          y3 = (y3 << (63-shift)) | (y2 >> (shift+1));
          y2 = (y2 << (63-shift)) | (y1 >> (shift+1));
          y1 = (y1 << (63-shift)) | (y0 >> (shift+1));
          y0 = (y0 << (63-shift));
          e -= 63-shift;
        }
        m_exponent = e;
        m_significand[4] = y4;
        m_significand[3] = y3;
        m_significand[2] = y2;
        m_significand[1] = y1;
        m_significand[0] = y0;
      }
    }
  }
}

void extfloat128_t::acc_t::from_17bytes_fix(const uint8_t* src)
{
  m_significand[0] = m_significand[1] = m_significand[2] = 0;
  memcpy(&m_significand[3], src, 16);
  m_significand[5] = src[16];
  m_sign = 0;
  m_exponent = exponent_bias - 9;
  m_isNormal = false;
}

// extfloat128_t::acc_t floating point arithmetic

void extfloat128_t::acc_t::mulx(const extfloat128_t& a, const extfloat128_t& b)
{
  uint32_t s  = a.m_sign ^ b.m_sign;
  uint32_t ea = a.m_exponent;
  uint32_t eb = b.m_exponent;
  if (ea != extfloat128_t::zero_biased_exponent && eb != extfloat128_t::zero_biased_exponent) {
    uint64_t e = exponent_bias - a.exponent_bias - b.exponent_bias + 1 + ea + eb;
    m_exponent = e;
    m_sign     = s;
    m_isNormal = false;

    // multiply significands
    uint64_t a0 = a.m_significand[0];
    uint64_t a1 = a.m_significand[1];
    uint64_t b0 = b.m_significand[0];
    uint64_t b1 = b.m_significand[1];

    uint64_t pr0, pr1, pr2, pr3, la, ha, lb, hb, lc;
    uint8_t c;
    pr0 = _umul128(a0, b0, &pr1);
    la  = _umul128(a0, b1, &ha);
    lb  = _umul128(a1, b0, &hb);
    pr1 += la; c  = pr1 < la;
    pr1 += lb; c += pr1 < lb;
    pr2 = c;
    pr2 += ha; c  = pr2 < ha;
    pr2 += hb; c += pr2 < hb;
    lc = _umul128(a1, b1, &pr3);
    pr2 += lc; c += pr2 < lc;
    pr3 += c;

    m_significand[1] = pr0;
    m_significand[2] = pr1;
    m_significand[3] = pr2;
    m_significand[4] = pr3;
    m_significand[0] = m_significand[5] = 0;
  } else {
    clear();
  }
}

void extfloat128_t::acc_t::PartialNormalize()
{
  m_isNormal = false;
  const uint64_t SMALL_MNT_THR = uint64_t(1) << (NBITS_MIN-64*4);
  const uint64_t LARGE_MNT_THR = uint64_t(1) << (NBITS_MAX-64*5);
  uint64_t y5 = m_significand[5];
  if (y5 + 1 < 2) {
    // y5 contains only sign bit
    // check if significand is not too small
    uint64_t y4 = m_significand[4];
    if ((y4 ^ y5) < SMALL_MNT_THR) {
      // partial of full normalization is necessary
      uint64_t y3 = m_significand[3];
      uint64_t y2 = m_significand[2];
      uint64_t y1 = m_significand[1];
      uint64_t y0 = m_significand[0];
      if ((y4 ^ y5) != 0 || (y3 ^ y5) >= SMALL_MNT_THR) {
        m_significand[0] = 0;
        m_significand[1] = y0;
        m_significand[2] = y1;
        m_significand[3] = y2;
        m_significand[4] = y3;
        m_significand[5] = y4;
        m_exponent -= 64;
      } else {
        // y5:y4 contains only sign bits and y3 is also small
        m_significand[0] = 0;
        m_significand[1] = 0;
        m_significand[2] = y0;
        m_significand[3] = y1;
        m_significand[4] = y2;
        m_significand[5] = y3;
        m_exponent -= 128;
        if ((y3 ^ y5) == 0 && (y2 ^ y5) < SMALL_MNT_THR) {
          // significand is still too small
          // hopefully, this case is rare, let's normalize
          doNormalize();
        }
      }
    }
  } else if (y5+LARGE_MNT_THR > LARGE_MNT_THR*2) {
    m_significand[0] = m_significand[1];
    m_significand[1] = m_significand[2];
    m_significand[2] = m_significand[3];
    m_significand[3] = m_significand[4];
    m_significand[4] = y5;
    m_significand[5] = int64_t(y5) < 0 ? uint64_t(-1) : 0;
    m_exponent += 64;
  }
}

static void addsub_significand_5to2(uint64_t* y, uint64_t x2, uint64_t x3, uint64_t x4, int de, bool sub) {
  uint64_t y2 = y[2];
  uint64_t y3 = y[3];
  uint64_t y4 = y[4];
  if (sub) {
    // subtraction. Y -= X
    uint64_t c;
    if (de < 0) {
      // Y -= (X << 64)
      c = (y3 < x2); y3 -= x2;
      x3 += c; c = (x3 < c) | (y4 < x3); y4 -= x3;
      y[5] -= x4;
    } else if (de < 64*1) {
      // Y -= X
      c = (y2 < x2); y2 -= x2;
      x3 += c; c = (x3 < c) | (y3 < x3); y3 -= x3;
      x4 += c; c = (x4 < c) | (y4 < x4); y4 -= x4;
    } else {
      uint64_t cc;
      uint64_t y1 = y[1];
      if (de < 64*2) {
        // Y -= (X >> 64)
        c = (y1 < x2); y1 -= x2;
        x3 += c; c = (x3 < c) | (y2 < x3); y2 -= x3;
        x4 += c; c = (x4 < c) | (y3 < x4); y3 -= x4;
      } else {
        uint64_t y0 = y[0];
        if (de < 64*3) {
          // Y -= (X >> 64*2)
          c = (y0 < x2); y0 -= x2;
          x3 += c; c = (x3 < c) | (y1 < x3); y1 -= x3;
          x4 += c; c = (x4 < c) | (y2 < x4); y2 -= x4;
        } else {
          if (de < 64*4) {
            // Y -= (X >> 64*3)
            c = (y0 < x3); y0 -= x3;
            x4 += c; c = (x4 < c) | (y1 < x4); y1 -= x4;
          } else { // (64*4 <= de < 64*5)
            // Y -= (X >> 64*4)
            c = (y0 < x4); y0 -= x4;
            cc = c; c = (y1 < c); y1 -= cc;
          }
          cc = c; c = (y2 < c); y2 -= cc;
        }
        cc = c; c = (y3 < c); y3 -= cc;
        y[0] = y0;
      }
      cc = c; c = (y4 < c); y4 -= cc;
      y[1] = y1;
    }
    y[5] -= c;
  } else {
    // addition
    uint64_t c;
    if (de < 0) {
      // Y += (X << 64)
      y3 += x2; c = (y3 < x2);
      x3 += c; y4 += x3; c = (x3 < c) | (y4 < x3);
      y[5] += x4;
    } else if (de < 64*1) {
      // Y += X
      y2 += x2; c = (y2 < x2);
      x3 += c; y3 += x3; c = (x3 < c) | (y3 < x3);
      x4 += c; y4 += x4; c = (x4 < c) | (y4 < x4);
    } else {
      uint64_t y1 = y[1];
      if (de < 64*2) {
        // Y += (X >> 64)
        y1 += x2; c = (y1 < x2);
        x3 += c; y2 += x3; c = (x3 < c) | (y2 < x3);
        x4 += c; y3 += x4; c = (x4 < c) | (y3 < x4);
      } else {
        uint64_t y0 = y[0];
        if (de < 64*3) {
          // Y += (X >> 64*2)
          y0 += x2; c = (y0 < x2);
          x3 += c; y1 += x3; c = (x3 < c) | (y1 < x3);
          x4 += c; y2 += x4; c = (x4 < c) | (y2 < x4);
        } else {
          if (de < 64*4) {
            // Y += (X >> 64*3)
            y0 += x3; c = (y0 < x3);
            x4 += c; y1 += x4; c = (x4 < c) | (y1 < x4);
          } else { // (64*4 <= de < 64*5)
            // Y += (X >> 64*4)
            y0 += x4; c = (y0 < x4);
            y1 += c; c = (y1 < c);
          }
          y2 += c; c = (y2 < c);
        }
        y[0] = y0;
        y3 += c; c = (y3 < c);
      }
      y[1] = y1;
      y4 += c; c = (y4 < c);
    }
    y[5] += c;
  }
  y[2] = y2;
  y[3] = y3;
  y[4] = y4;
}

void extfloat128_t::acc_t::addsub(const extfloat128_t& x, bool sub)
{
  uint32_t xe32 = x.m_exponent;
  if (xe32 != extfloat128_t::zero_biased_exponent) {
    uint64_t xe = (exponent_bias - x.exponent_bias) + xe32;
    uint64_t ye = m_exponent;
    if (xe + 64*5 > ye) {
      // x is not void
      uint32_t xs = x.m_sign ^ uint32_t(sub);
      uint64_t ys = m_sign;
      if (xe == ye) {
        // that is a very common case, so it handled separately
        uint64_t y3 = m_significand[3];
        uint64_t y4 = m_significand[4];
        uint64_t y5 = m_significand[5];
        uint64_t x3 = x.m_significand[0];
        uint64_t x4 = x.m_significand[1];
        if (xs != ys) {
          // subtraction. Y -= X
          uint64_t c = (y3 < x3); y3 -= x3;
          y5 -= (y4 < x4); y4 -= x4;
          y5 -= (y4 < c);  y4 -= c;
          if ((((y4^y5)|(y3^y5))==0) && (y5+1 < 2)) {
            // 3 MS words contain nothing, but sign bit - very common case of subtraction of result of round/trunc
            y4 = m_significand[2];
            y3 = m_significand[1];
            m_significand[2] = m_significand[0];
            m_significand[1] = m_significand[0] = 0;
            m_exponent = ye - 128;
          }
        } else {
          // addition. Y += X
          y3 += x3;
          uint64_t c = y3 < x3;
          y4 += x4; y5 += (y4 < x4);
          y4 += c;  y5 += (y4 < c);
        }
        m_significand[5] = y5;
        m_significand[4] = y4;
        m_significand[3] = y3;
      } else {
        if (ye < xe) {
          uint64_t wshift = (xe - ye + (64*6-NBITS_MAX)+1) / 64;
          if (wshift < 6) {
            // *this is not void
            unsigned i0 = wshift;
            if (i0 > 0) {
              m_exponent = ye = (ye + i0*64);
              // shift significand i0 words to the right
              unsigned i;
              for (i = 0; i0 < 6; ++i, ++i0)
                m_significand[i] = m_significand[i0];
              uint64_t msw = int64_t(m_significand[5]) < 0 ? uint64_t(-1) : 0;
              for (; i < 6; ++i)
                m_significand[i] = msw; // replicate sign word
            }
          } else {
            // *this is void
            *this = x;
            m_sign ^= uint8_t(sub);
            return;
          }
        }
        // shift X to the right
        int de = int(ye - xe); // de in range [-59..64*5-1]
        uint64_t x2 = 0;
        uint64_t x3 = x.m_significand[0];
        uint64_t x4 = x.m_significand[1];
        unsigned rshift = unsigned(de) % 64;
        if (rshift != 0) {
          unsigned lshift = 64 - rshift;
          x2 =                  (x3 << lshift);
          x3 = (x3 >> rshift) | (x4 << lshift);
          x4 = (x4 >> rshift);
        }
        addsub_significand_5to2(m_significand, x2, x3, x4, de, xs != ys);
      }
      PartialNormalize();
    }
  }
}
#if 0
static void print(const extfloat128_t::acc_t& a, const char* prefix = "") {
  printf("%s%d %11I64d", prefix, a.m_sign, a.m_exponent-a.exponent_bias);
  const uint64_t* bits = a.m_significand;
  printf(" %016I64x:%016I64x:%016I64x:%016I64x:%016I64x:%016I64x %s\n"
    , bits[5], bits[4], bits[3], bits[2], bits[1], bits[0]
    , a.m_isNormal ? "n" : "");
}
#endif

void extfloat128_t::acc_t::maddsub(const extfloat128_t& a, const extfloat128_t& b, bool sub)
{
  uint32_t ae = a.m_exponent;
  uint32_t be = b.m_exponent;
  if (ae != extfloat128_t::zero_biased_exponent && be != extfloat128_t::zero_biased_exponent) {
    uint64_t xe = exponent_bias - a.exponent_bias - b.exponent_bias + 1 + ae + be;
    uint64_t ye = m_exponent;

    if (xe + 64*5 > ye) {
      // product X=a*b is not void
      m_isNormal = false;

      // multiply significands
      uint64_t a0 = a.m_significand[0];
      uint64_t a1 = a.m_significand[1];
      uint64_t b0 = b.m_significand[0];
      uint64_t b1 = b.m_significand[1];

      uint64_t x1, x2, x3, x4, la, ha, lb, hb, lc;
      uint8_t c;
      x1 = _umul128(a0, b0, &x2);
      la = _umul128(a0, b1, &ha);
      lb = _umul128(a1, b0, &hb);
      x2 += la; c  = x2 < la;
      x2 += lb; c += x2 < lb;
      x3 = c;
      x3 += ha; c  = x3 < ha;
      x3 += hb; c += x3 < hb;
      lc = _umul128(a1, b1, &x4);
      x3 += lc; c += x3 < lc;
      x4 += c;
      uint64_t x0 = 0;

      uint32_t xs = a.m_sign ^ b.m_sign ^ uint32_t(sub);
      if (ye < xe) {
        uint64_t wshift = (xe - ye + (64*6-NBITS_MAX)+1) / 64;
        if (wshift < 6) {
          // accumulator Y is not void
          unsigned i0 = wshift;
          if (i0 > 0) {
            m_exponent = ye = (ye + i0*64);
            // shift significand i0 words to the right
            unsigned i;
            for (i = 0; i0 < 6; ++i, ++i0)
              m_significand[i] = m_significand[i0];
            uint64_t msw = int64_t(m_significand[5]) < 0 ? uint64_t(-1) : 0;
            for (; i < 6; ++i)
              m_significand[i] = msw; // replicate sign word
          }
        } else {
          // accumulator Y is void
          m_significand[0] = 0;
          m_significand[1] = x1;
          m_significand[2] = x2;
          m_significand[3] = x3;
          m_significand[4] = x4;
          m_significand[5] = 0;
          m_exponent = xe;
          m_sign     = xs;
          return;
        }
      }

      // print(*this, "a ");
      // shift product X to the right
      unsigned de = ye + 64 - xe; // de in range [5..64*6-1]
      unsigned rshift = de % 64;
      if (rshift != 0) {
        unsigned lshift = 64 - rshift;
        x0 =                  (x1 << lshift);
        x1 = (x1 >> rshift) | (x2 << lshift);
        x2 = (x2 >> rshift) | (x3 << lshift);
        x3 = (x3 >> rshift) | (x4 << lshift);
        x4 = (x4 >> rshift);
      }
      unsigned wshift = de / 64; // wshift ==0 means x4 aligned with m_significand[5]
      while (wshift > 1) {
        x0 = x1;
        x1 = x2;
        x2 = x3;
        x3 = x4;
        x4 = 0;
        wshift -= 1;
      }
      uint64_t* y =  &m_significand[1 - wshift];

      // printf("xe=%I64d, de=%u\n", xe-exponent_bias, de, wshift, i0);
      uint32_t ys = m_sign;
      if (ys == xs) {
        // addition
        uint64_t c, xx, yy;
        yy = y[0];     xx = x0 + yy; c = (xx < yy);            y[0] = xx;
        yy = y[1] + c; xx = x1 + yy; c = (xx < yy) | (yy < c); y[1] = xx;
        yy = y[2] + c; xx = x2 + yy; c = (xx < yy) | (yy < c); y[2] = xx;
        yy = y[3] + c; xx = x3 + yy; c = (xx < yy) | (yy < c); y[3] = xx;
        yy = y[4] + c; xx = x4 + yy; c = (xx < yy) | (yy < c); y[4] = xx;
        if (wshift != 0)
          m_significand[5] += c;
      } else {
        // subtraction
        uint64_t c, xx, yy;
        yy = y[0];              c = (yy < x0);            yy -= x0; y[0] = yy;
        yy = y[1]; xx = x1 + c; c = (xx < c) | (yy < xx); yy -= xx; y[1] = yy;
        yy = y[2]; xx = x2 + c; c = (xx < c) | (yy < xx); yy -= xx; y[2] = yy;
        yy = y[3]; xx = x3 + c; c = (xx < c) | (yy < xx); yy -= xx; y[3] = yy;
        yy = y[4]; xx = x4 + c; c = (xx < c) | (yy < xx); yy -= xx; y[4] = yy;
        if (wshift != 0)
          m_significand[5] -= c;
      }
      // print(*this, "b ");
      PartialNormalize();
    }
  }
}

 // return value: -1 if rounded toward 0, 1 if rounded away from zero, 0 if unchanged
int extfloat128_t::acc_t::round_to_nearest_tie_to_even()
{
  int ret = 0;
  uint64_t e = m_exponent;
  if (e < exponent_bias+319) {
    if (e > exponent_bias-(NBITS_MAX-320)) {
      if (!m_isNormal)
        doNormalizeSign();
      unsigned halfBitPos = exponent_bias + 318 - e;
      unsigned halfWordPos   = halfBitPos / 64;
      unsigned halfInWordPos = halfBitPos % 64;
      unsigned oneBitPos  = halfBitPos + 1;
      uint64_t lsbits = 0;
      for (unsigned i = 0; i < halfWordPos; ++i) {
        lsbits |= m_significand[i];
        m_significand[i] = 0;
      }
      uint64_t hw = m_significand[halfWordPos];
      uint64_t hg = (hw << (63-halfInWordPos)) | (lsbits != 0);
      hw &= (uint64_t(-2) << halfInWordPos);
      m_significand[halfWordPos] = hw;
      ret = (hg == 0) ? 0 : -1; // rounded toward zero
      uint64_t odd = (m_significand[oneBitPos/64] >> (oneBitPos%64)) & 1;
      hg |= odd; // in order to break tie toward even
      if (hg > (uint64_t(1) << 63)) {
        // add 1
        uint64_t inc = (uint64_t(1) << (oneBitPos%64));
        unsigned i = oneBitPos/64;
        m_significand[i] += inc;
        // propagate carry
        while (m_significand[i]==0) {
          ++i;
          m_significand[i] += 1;
        }
        m_isNormal = false;
        if (i == 5)
          doNormalize();
        ret = 1; // rounded away from zero
      } else if ((m_significand[5] | m_significand[4])==0) {
        clear();
      }
    } else {
      ret = (e == 0) ? 0 : -1;
      clear();
    }
  }
  return ret;
}

static inline uint64_t load_u64(const uint32_t* src) {
  uint64_t ret;
  memcpy(&ret, src, sizeof(ret));
  return ret;
}

// [0] : hi(a[1]*b[0]) + lo(a[1]*b[1]) +  hi(a[0]*b[1]) + lo(a[0]*b[2])
// [1] : hi(a[1]*b[1]) + lo(a[1]*b[2]) +  hi(a[0]*b[2]) + lo(a[0]*b[3])
// [2] : hi(a[1]*b[2]) + lo(a[1]*b[3]) +  hi(a[0]*b[3]) + lo(a[0]*b[4])
// [3] : hi(a[1]*b[3]) + lo(a[1]*b[4]) +  hi(a[0]*b[4]) + lo(a[0]*b[5])
// [4] : hi(a[1]*b[4]) + lo(a[1]*b[5]) +  hi(a[0]*b[5]) + lo(a[0]*b[6])
// extfloat128_t::acc_t fix-point arithmetic
void extfloat128_t::acc_t::fix_mulx5(const extfloat128_t& a, const uint32_t b[7*2], int b_exponent)
{
  uint32_t ea = a.m_exponent;
  if (ea != extfloat128_t::zero_biased_exponent) {
    uint64_t e = exponent_bias - a.exponent_bias + 1 + ea + b_exponent;
    m_exponent = e;
    m_sign     = a.m_sign;
    m_isNormal = false;

    uint64_t a0 = a.m_significand[0];
    uint64_t a1 = a.m_significand[1];
    uint64_t bx, r00, r10;
    bx = load_u64(b); b += 2;
    _umul128(a1, bx, &r10);
    bx = load_u64(b); b += 2;
    _umul128(a0, bx, &r00);
    uint64_t c = 0;
    for (int i = 0; i < 4; ++i) {
      uint64_t r01, r11, ml;
      ml = _umul128(a1, bx, &r11);
      r10 += ml;
      r11 += (r10 < ml);

      bx = load_u64(b); b += 2;
      ml = _umul128(a0, bx, &r01);
      r00 += ml;
      r01 += (r00 < ml);

      uint64_t s = r00 + r10;
      uint64_t c0 = s < r10;
      s += c;
      c = c0 + (s < c);
      m_significand[i] = s;
      r00 = r01;
      r10 = r11;
    }
    r10 += a1 * bx;
    bx = load_u64(b);
    r00 += a0 * bx;
    m_significand[4] = r10 + r00 + c;
    m_significand[5] = 0;
  } else {
    clear();
  }
}

// [0] : hi(a[1]*b[0]) + lo(a[1]*b[1]) +  hi(a[0]*b[1]) + lo(a[0]*b[2])
// [1] : hi(a[1]*b[1]) + lo(a[1]*b[2]) +  hi(a[0]*b[2]) + lo(a[0]*b[3])
// [2] : hi(a[1]*b[2]) + lo(a[1]*b[3]) +  hi(a[0]*b[3]) + lo(a[0]*b[4])
void extfloat128_t::acc_t::fix_mulx3(const extfloat128_t& a, const uint32_t b[5*2], int b_exponent)
{
  uint32_t ea = a.m_exponent;
  if (ea != extfloat128_t::zero_biased_exponent) {
    uint64_t e = exponent_bias - a.exponent_bias + 1 + ea + b_exponent;
    m_exponent = e;
    m_sign     = a.m_sign;
    m_isNormal = false;
    m_significand[0] = 0;
    m_significand[1] = 0;

    uint64_t a0 = a.m_significand[0];
    uint64_t a1 = a.m_significand[1];
    uint64_t bx, r00, r10;
    bx = load_u64(b); b += 2;
    _umul128(a1, bx, &r10);
    bx = load_u64(b); b += 2;
    _umul128(a0, bx, &r00);
    uint64_t c = 0;
    for (int i = 0; i < 2; ++i) {
      uint64_t r01, r11, ml;
      ml = _umul128(a1, bx, &r11);
      r10 += ml;
      r11 += (r10 < ml);

      bx = load_u64(b); b += 2;
      ml = _umul128(a0, bx, &r01);
      r00 += ml;
      r01 += (r00 < ml);

      uint64_t s = r00 + r10;
      uint64_t c0 = s < r10;
      s += c;
      c = c0 + (s < c);
      m_significand[i+2] = s;
      r00 = r01;
      r10 = r11;
    }
    r10 += a1 * bx;
    bx = load_u64(b);
    r00 += a0 * bx;
    m_significand[4] = r10 + r00 + c;
    m_significand[5] = 0;
  } else {
    clear();
  }
}


// extfloat128_t::ex_t export

// if (bRound) round to nearest otherwise truncate
extfloat128_t extfloat128_t::acc_t::to_extfloat128_t(bool bRound)
{
  normalize();
  uint64_t e  = m_exponent;
  uint64_t m0 = m_significand[3];
  uint64_t m1 = m_significand[4];
  if (bRound) {
    uint64_t inc = m_significand[2] >> 63;
    m0 += inc; inc = m0 < inc;
    m1 += inc; inc = m1 < inc;
    e  += inc;
    m1 |= (uint64_t(1) << 63);
  }
  uint32_t rete = uint32_t(e - (exponent_bias - extfloat128_t::exponent_bias));
  if (e < exponent_bias + extfloat128_t::min_exponent_val ||
      e > exponent_bias + extfloat128_t::max_exponent_val) {
    rete = (e < exponent_bias) ? extfloat128_t::zero_biased_exponent : extfloat128_t::inf_nan_biased_exponent;
    m1 = m0 = 0;
  }
  extfloat128_t ret;
  ret.m_significand[0] = m0;
  ret.m_significand[1] = m1;
  ret.m_exponent       = rete;
  ret.m_sign           = m_sign;
  return ret;
}

bool extfloat128_t::acc_t::is_odd() const
{
  uint64_t e = m_exponent;
  if (e <= exponent_bias+319 && e > exponent_bias-(NBITS_MAX-320)) {
    unsigned oneBitPos = exponent_bias + 319 - e;
    return (m_significand[oneBitPos/64] >> (oneBitPos%64)) & 1;
  }
  return 0;
}
