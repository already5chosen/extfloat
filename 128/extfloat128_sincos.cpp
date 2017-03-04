#include <stdio.h>
#include "extfloat128.h"


// sinPi coefficients
static double sinPi_coeff_l[3] = {
  7.95205400147551260000e-007,
 -2.19153534478302170000e-005,
  4.66302805767612550000e-004,
};
static uint8_t sinPi_coeff_h[7][17] = {
{  4,103,176,102, 20,118, 64,222, 17,251, 56, 68,239,167,131,241,113,}, // -7.37043094571435040000e-003
{ 52,165,199, 86,239, 34, 99, 61,200, 13, 60,247, 67, 26, 60,168,120,}, // 8.21458866111282330000e-002
{ 43,177,224, 50,178,191,112, 47,243, 45,236, 21,115,102,105,153,127,}, // -5.99264529320792110000e-001
{ 22,134,209,  3, 79, 34, 52, 63,146, 14, 87,173, 59,227, 53,163,130,}, // 2.55016403987734550000e+000
{ 38, 72,242, 42,113, 47,199, 93,245,149,242, 45, 49,231, 93,165,133,}, // -5.16771278004997030000e+000
{100,251,173, 93,198, 71,254,232, 58,115,118,150,123, 46,189,187, 81,}, // -8.74227800037248510000e-008
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,219, 15,201,130,}, // 3.14159274101257320000e+000
};
// end of sinPi coefficients

static const uint8_t sin_coeff[9][17] = {
{  3,243, 97, 45, 27,  4, 27,209,  3,  7,119,194, 29,145,147,202, 30,}, //  2.81131275270721464813e-15
{ 48,195,158,150,205, 70,119,150,  3,116, 68,158, 52,159, 63,215, 47,}, // -7.64716372123615355900e-13
{ 54,184,110,115, 61, 84,196, 85, 46, 62, 99, 67,157, 48,146,176, 62,}, //  1.60590438368211971656e-10
{245,127,137,155, 58, 51,127, 69,121, 28, 39,170, 63, 43, 50,215, 77,}, // -2.50521083854417202239e-08
{ 69,253, 28,247,214, 36, 10, 86,125,156, 57,182, 42, 29,239,184, 90,}, //  2.75573192239858925110e-06
{174,107,  3,155,253, 12,208,  0, 13,208,  0, 13,208,  0, 13,208,103,}, // -1.98412698412698412526e-04
{165, 72,135,136,136,136,136,136,136,136,136,136,136,136,136,136,114,}, //  8.33333333333333321769e-03
{170,170,170,170,170,170,170,170,170,170,170,170,170,170,170,170,123,}, // -1.66666666666666657415e-01
{192,157, 21,138, 43,  4,169,110,221,169,138,255,238,132,167,151, 93,}, // -1.32820649256081252783e-44 * 2**128
};

static const uint8_t cos_coeff[9][17] = {
{116,120,181,135,195,232, 29, 50, 75, 85, 44,176,227,116, 60,215, 38,}, //  4.77920277998764881261e-14
{217,146,179,242,238,130, 19, 85,133,101,130, 20, 64,165,203,201, 55,}, // -1.14707455776208345692e-11
{101,191, 10,186,219, 75, 45,166, 98,189,190,198,127,199,118,143, 70,}, //  2.08767569878673060932e-09
{136,149, 85, 64,104,108,105,101,144,227,250,196,187,125,242,147, 85,}, // -2.75573192239858882758e-07
{ 61,247, 51, 61,242, 65,198,  0, 13,208,  0, 13,208,  0, 13,208, 96,}, //  2.48015873015873015658e-05
{172,227, 22, 24, 89, 11,182, 96, 11,182, 96, 11,182, 96, 11,182,109,}, // -1.38888888888888894189e-03
{ 52,187,164,170,170,170,170,170,170,170,170,170,170,170,170,170,118,}, //  4.16666666666666643537e-02
{249,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,125,}, // -5.00000000000000000000e-01
{160, 79,159,251,255, 69, 28,162, 23,104,181,174, 21,180, 22,180,101,}, // -2.52357998257430706828e-43 * 2**128
};

// sin table
const static uint8_t sin_tab[32][17] = { // sin(x*pi*0.5/32) where x[0]=1.25, x[1]=1.75, x[2..29] = i+0.5, x[30]=3.25, x[31]=30.75
{251, 52,199,  0,160, 54, 14, 10,143,246,111, 16,252, 60,183,178, 15,},
{163,133, 41,203, 38,129, 95,186,127,193, 24,164,154, 10,208,246, 21,},
{226, 70,241,218,202, 53,229,  8, 78, 52, 14,115,169, 86, 78, 86, 31,},
{ 35, 67,124, 44,214,238,188,251,143,169,140,127, 22,137, 40,196, 43,},
{ 44,243,160,178,121,107, 24,175, 67,249,197,158,252,212,  4, 23, 56,},
{247,130, 33,154,169,122,171, 32, 76,130,221,217,199,138, 73, 71, 68,},
{123,  5,230,129, 81,162,175,196,151, 12,  5,152, 93, 80,114, 77, 80,},
{239,219,241, 54, 58,155,247, 97, 78,187,122, 22,233,195, 20, 34, 92,},
{ 23,114,177, 86,221, 94, 72,  9,212,182, 40,182,163, 14,229,189,103,},
{ 90,230,104, 61,254, 13,103, 57, 20, 90,120, 17,199,100,186, 25,115,},
{181,104, 92, 72,137, 72, 62,224, 19,237,231,106,226,111,147, 46,126,},
{ 30, 76,200,216, 97,131,202,146,184, 33, 20, 89,218,160,154,245,136,},
{ 69,110,193,113, 30,145, 33,120,177, 68,245,150,232,102, 42,104,147,},
{135,103,174, 56,150, 84, 57, 94,226,227,201,133,  2, 73,209,127,157,},
{155,162,104, 52,153,117,110,159, 20,158, 72, 47, 31,223, 85, 54,167,},
{166,209,188,148, 92, 13,208,202,228,218,246,102,233,168,186,133,176,},
{ 68,219,244,139, 43,178,227,225,157, 26,178,252,127,191, 65,104,185,},
{216,214, 11,178,184, 92,113,240,189,144,110,187,252, 95,112,216,193,},
{ 76, 37, 50,105,225, 27, 36, 53,131,122,218, 31,147, 76, 18,209,201,},
{228,145,239,171,232, 32,234, 79,116,237, 14, 60, 49,  2, 61, 77,209,},
{221,163,116, 53,102,  0,254,185, 36,219,252, 15,168,192, 82, 72,216,},
{213,108,174,195,254,113,170, 25, 75,251, 76,169,124, 99,  5,190,222,},
{135, 92, 81,245,103,174,133, 39, 18,180,167,143,160,  9, 89,170,228,},
{ 61, 35,232,181,161,180, 69,173, 21, 98,205, 73,110,138,166,  9,234,},
{252, 26,176,100,182, 15,194,248,134,  7,227, 17,102,182,157,216,238,},
{109,162, 86, 93, 69,220,107, 72,165,116,143,  8, 71, 98, 71, 20,243,},
{176, 34,252, 31,204, 89, 31,121, 44,232, 25, 75, 66, 59,  7,186,246,},
{122,128, 41,155,209,  8,174,  4, 69, 40, 70, 44, 39, 99,157,199,249,},
{ 36,105, 23, 92,203,255,120,103, 37,171, 73, 93,138,211, 39, 59,252,},
{205,238,116,182,157,149, 29,205,144, 61,154,254, 12,135, 35, 19,254,},
{ 95,248,195,113, 66, 84,  1, 25, 31,209, 72,216,234,229, 87, 14,255,},
{ 61,107,174,204,115, 73, 49,156, 81,  3,106,141,115, 44,171,132,255,},
};
// end of sin table

static const uint32_t invPi_tab[64+4] = {
 0xcaf27f1d, 0x9f3a1f35,
 0x6b1e5ef8, 0xc33d26ef,
 0x98327dbb, 0x32c2de4f,
 0x3f7e33e8, 0xa5ff0705,
 0x5719053e, 0xddaf44d1,
 0x8b961ca6, 0x8359c476,
 0xdce8092a, 0x19c367cd,
 0x8c6b47c4, 0x60e27bc0,
 0xca73a8c9, 0x06061556,
 0x4d732731, 0x8dffd880,
 0x14a06840, 0x6599855f,
 0x5ee61b08, 0xa9e39161,
 0x9af4361d, 0xf0cfbc20,
 0xfc7b6bab, 0x56033046,
 0x1f8d5d08, 0x6bfb5fb1,
 0x8a5292ea, 0x3d0739f7,
 0xebe5f17b, 0x7527bac7,
 0x9e5fea2d, 0x4f463f66,
 0x27cb09b7, 0x6d367ecf,
 0x5a0a6d1f, 0xef2f118b,
 0xde05980f, 0x1ff897ff,
 0xbdf9283b, 0x9c845f8b,
 0x835339f4, 0x3991d639,
 0xb45f7e41, 0xe99c7026,
 0x2ebb4484, 0xe88235f5,
 0xb129a73e, 0xfe1deb1c,
 0x09d1921c, 0x06492eea,
 0x424dd2e0, 0xb7246e3a,
 0xdebbc561, 0xfe5163ab,
 0x3c439041, 0xdb629599,
 0xf534ddc0, 0xfc2757d1,
 0x4e441529, 0xa2f9836e, // (2/pi)*2**64
 0, 0, 0, 0,
};


struct cossinPi_core_ctrl_t {
  int tabi;
  int y_sign;
  int x_sign;
};
static void cossinPi_core_bodyA(extfloat128_t& x, int isSin, cossinPi_core_ctrl_t* pCtrl)
{ // x in range [0..2)
  int ix  = static_cast<int>(x.scale_and_trunc(7)); // trunc(abs(x)*128); ix in range [0..255]
  int ixs = ix + 64 - (isSin << 6);     // shift cos by +pi/2
  pCtrl->y_sign = ((ixs & 128) != 0);   // quadrants 2 & 3 are same as 0 & 1, but with opposite sign
  int x_sign = (ixs & 64) ? -1 : 0;
  ixs = (ixs ^ x_sign) & 63;            // quadrant 1 => 0
  pCtrl->x_sign = x_sign & 1;

  int tabi = ixs >> 1;
  int dx = 2;
  if (((ix + 4) & 63) < 8) {
    static const uint8_t dx_tab[8] = {  0,  0, 1, 3,  1,  3,  4, 4 };
    static const int8_t  ti_tab[8] = { -1, -1, 0, 1, 30, 31, 32, 32};
    dx   = dx_tab[ix & 7];
    tabi = ti_tab[ixs & 7];
  }
  pCtrl->tabi = tabi;
  dx += ((ix << 1) & -4);
  if (dx)
    x -= dx * (1.0/256);
}

static extfloat128_t cossinPi_core_bodyB(const extfloat128_t& x, cossinPi_core_ctrl_t* pCtrl, extfloat128_t::acc_t* xEx)
{
  // calculate ySin = sin(x)
  extfloat128_t xx = x * x;
  double xx_d = xx.convert_to_double();
  double ySin_p_d = sinPi_coeff_l[0];
  for (int i = 1; i < 3; ++i)
    ySin_p_d = ySin_p_d * xx_d + sinPi_coeff_l[i];
  extfloat128_t ySin_p(ySin_p_d), coef;
  for (int i = 0; i < 6; ++i) {
    coef.from_17bytes(sinPi_coeff_h[i]);
    ySin_p = fma(ySin_p, xx, coef); // use fma() instead of * and + not because it's more precise, but because it's a little faster
  }
  coef.from_17bytes(sinPi_coeff_h[6]);
  extfloat128_t::acc_t ySinEx;
  ySinEx.mulx(x, coef);
  ySinEx.madd(x, ySin_p);
  int tabi = pCtrl->tabi;
  if (xEx && tabi < 11) {
    *xEx -= x;
    ySinEx.madd(xEx->trunc(), coef);
  }
  extfloat128_t ySin = ySinEx.round();

  if (tabi >= 0) {
    // calculate yCosn = cos(x) - 1
    extfloat128_t::acc_t yCosnEx;
    yCosnEx.clear();
    if (!is_zero(ySin)) {
      yCosnEx.mulx(ySin, ySin);
      yCosnEx.m_sign = true; // -ySin^2
      extfloat128_t::acc_t yCosSqr = yCosnEx;
      yCosSqr += extfloat128_t::one();
      extfloat128_t yCosn = yCosSqr.round().sqrt() - extfloat128_t::one();
      // improve precision by iterative formula cos(x)-1 = (-sin(x)^2 - (cos(x)-1)^2)/2
      yCosnEx.msub(yCosn, yCosn);
      yCosnEx.m_exponent -= 1; // /= 2
    }

    if (tabi < 32) {
      // res = tabSin*yCos + tabCos*ySin
      //     = tabSin*(yCosn+1) + tabCos*ySin = tabSin + tabSin*yCosn + tabCos*ySin
      //     = tabSinA + tabSinB + (tabSinA + tabSinB)*yCosn + tabCosA*ySin + tabCosB*ySin
      extfloat128_t yCosn = yCosnEx.round();
      extfloat128_t tabCos[2]; // [0] - fine part, [1] - coarse part
      extfloat128_t::from_17bytes_fix(tabCos, sin_tab[31-tabi]);
      extfloat128_t::acc_t acc;
      acc.from_17bytes_fix(sin_tab[tabi]);
      acc.madd(acc.round(), yCosn);
      if (tabi < 11) {
        ySinEx -= ySin;
        acc.madd(tabCos[1], ySinEx.trunc());
      }
      acc.madd(tabCos[0], ySin);
      acc.madd(tabCos[1], ySin);
      ySin = acc.round();
    } else {
      yCosnEx += extfloat128_t::one();
      ySin = yCosnEx.round();
    }
  }
  ySin.m_sign = pCtrl->y_sign;
  return ySin;
}

// cossinPi_core
// return isSin ? sin(a*pi) : cos(a*pi)
extfloat128_t extfloat128_t::cossinPi_core(const extfloat128_t& a, int isSin) {
  extfloat128_t x = mod_pow2(a, 1); // reduce to range (-2..+2)
  if (isfinite(x)) {
    uint32_t y_sign = x.m_sign & isSin;
    x.m_sign = 0;                                     // x in range [0..2)
    cossinPi_core_ctrl_t ctrl;
    cossinPi_core_bodyA(x, isSin, &ctrl);
    x.m_sign    ^= ctrl.x_sign;
    ctrl.y_sign ^= y_sign;
    return cossinPi_core_bodyB(x, &ctrl, 0);
  }
  return x;
}

static extfloat128_t sin45deg = extfloat128_t::pow2(-1).sqrt();

// static int64_t convert_to_int64(const extfloat128_t& x) {
  // unsigned rshift = (extfloat128_t::exponent_bias + 63 - x.m_exponent) % 64;
  // int64_t ret = x.m_significand[1] >> rshift;
  // return x.m_sign ? -ret : ret;
// }

#if 0
void pp(const extfloat128_t::acc_t& x, const char* prefix)
{
  printf("%s %d %016I64x:%016I64x:%016I64x:%016I64x:%016I64x:%016I64x %I64d\n"
  ,prefix
  ,x.m_sign
  ,x.m_significand[5]
  ,x.m_significand[4]
  ,x.m_significand[3]
  ,x.m_significand[2]
  ,x.m_significand[1]
  ,x.m_significand[0]
  ,int64_t(x.m_exponent-x.exponent_bias)
  );
}
void pp(const extfloat128_t& x, const char* prefix)
{
  printf("%s %d %016I64x:%016I64x %d\n"
  ,prefix
  ,x.m_sign
  ,x.m_significand[1]
  ,x.m_significand[0]
  ,int32_t(x.m_exponent-x.exponent_bias)
  );
}
#endif

// cossin_core
// return isSin ? sin(a) : cos(a)
extfloat128_t extfloat128_t::cossin_core(const extfloat128_t& a, int isSin) {
  uint32_t e = a.m_exponent;
  if (e < extfloat128_t::exponent_bias + 1024) {
    uint32_t change_y_sign = 0;
    extfloat128_t af = a;
    af.m_sign = 0;
    int a_sign = a.m_sign & isSin;
    if (e >= extfloat128_t::exponent_bias - 4) {
      extfloat128_t::acc_t xAcc;
      uint32_t bit_offset = 0;
      if (e > af.exponent_bias + 1)
        bit_offset = e - (af.exponent_bias + 1);
      uint32_t dword_offset = bit_offset/32;
      const uint32_t* pInvPi = &invPi_tab[64-10-dword_offset];
      xAcc.fix_mulx3(af, pInvPi+4, -2-dword_offset*32);
      if (((xAcc.m_significand[4] + 1) & 0xFFFF) < 2)
        xAcc.fix_mulx5(af, pInvPi, -2-dword_offset*32);
      xAcc.m_significand[4] &= uint64_t(-1) >> (bit_offset % 32); // reduce to range [0..+2)
      extfloat128_t x = xAcc.trunc();
      extfloat128_t xMod2 = x;
      cossinPi_core_ctrl_t ctrl;
      cossinPi_core_bodyA(xMod2, isSin, &ctrl);
      xAcc -= x;
      xAcc += xMod2;
      xAcc.m_sign ^= ctrl.x_sign;
      ctrl.y_sign ^= a_sign;
      return cossinPi_core_bodyB(xAcc.round(), &ctrl, &xAcc);
    }

    // abs(af) < 2**(-4)
    extfloat128_t y;
    if (e >= extfloat128_t::exponent_bias - 128) {
      // 2**(-128) <= abs(af) < 2**(-4)
      // calculate sin/cos with polynomial over af**2
      extfloat128_t aa = af * af;
      extfloat128_t coef;
      if (isSin) {
        y.from_17bytes(sin_coeff[0]);
        for (int i = 1; i < 8; ++i) {
          coef.from_17bytes(sin_coeff[i]);
          y = fma(y, aa, coef); // use fma() instead of * and + not because it's more precise, but because it's a little faster
        }
        // coef.from_17bytes(sin_coeff[8]);
        // coef.m_exponent -= 128; // last coefficient requires additional scaling
        // y = fma(y, aa, coef);
        y *= aa; // last coefficient is so small that we can ignore it
        y = fma(af, y, af); // Last step. Here fma() is used for better precision
      } else {
        y.from_17bytes(cos_coeff[0]);
        for (int i = 1; i < 8; ++i) {
          coef.from_17bytes(cos_coeff[i]);
          y = fma(y, aa, coef); // use fma() instead of * and + not because it's more precise, but because it's a little faster
        }
        // coef.from_17bytes(cos_coeff[8]);
        // coef.m_exponent -= 128; // last coefficient requires additional scaling
        // y = fma(y, aa, coef);
        // last coefficient is so small that we can ignore it
        y = fma(aa, y, extfloat128_t::one()); // Last step. Here fma() is used for better precision
      }
    } else {
      // very small absolute value
      y = isSin ? af : one();
    }
    y.m_sign ^= change_y_sign ^ a_sign;
    return y;

  } else {
    if (isfinite(a)) {
      return sin45deg;
    } else if (isnan(a)) {
      return a;
    } else {
      return nan();
    }
  }
}

