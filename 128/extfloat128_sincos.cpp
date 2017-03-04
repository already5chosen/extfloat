//#include <cmath>
//#include <stdio.h>
#include "extfloat128.h"


// sin coefficients
static double sin_coeff_l[3] = {
  7.95205400147551260000e-007,
 -2.19153534478302170000e-005,
  4.66302805767612550000e-004,
};
static uint8_t sin_coeff_h[7][17] = {
{  4,103,176,102, 20,118, 64,222, 17,251, 56, 68,239,167,131,241,113,}, // -7.37043094571435040000e-003
{ 52,165,199, 86,239, 34, 99, 61,200, 13, 60,247, 67, 26, 60,168,120,}, // 8.21458866111282330000e-002
{ 43,177,224, 50,178,191,112, 47,243, 45,236, 21,115,102,105,153,127,}, // -5.99264529320792110000e-001
{ 22,134,209,  3, 79, 34, 52, 63,146, 14, 87,173, 59,227, 53,163,130,}, // 2.55016403987734550000e+000
{ 38, 72,242, 42,113, 47,199, 93,245,149,242, 45, 49,231, 93,165,133,}, // -5.16771278004997030000e+000
{100,251,173, 93,198, 71,254,232, 58,115,118,150,123, 46,189,187, 81,}, // -8.74227800037248510000e-008
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,219, 15,201,130,}, // 3.14159274101257320000e+000
};
// end of sin coefficients


// sin table
static uint8_t sin_tab[31][17] = { // sin(x*pi)-4*x*(1-x) where x = (i+1)*0.5/32
{ 18, 23,197,183,178, 84,121,161,118,243,  9,236,134,248,178,143, 12,},
{ 33,113,238,119, 41,145, 45,231,216,225, 43,180, 41,188,166, 23, 25,},
{160,127,240,129, 71,104,162,  2,192, 68,116,194, 28,221, 32,144, 37,},
{243, 19, 54,243,147,  3, 48, 50,151,108, 21, 76,211,120,112,241, 49,},
{ 85, 99,150,220,233, 18,119,136,144, 94, 53,190, 66,246,242, 51, 62,},
{199,169,251, 53,237,185,211, 37,215,162, 22,124, 86,187, 24, 80, 74,},
{ 26, 31,252,119,194,156,202,143, 64,248,115,127,172,214,105, 62, 86,},
{190, 27,238,126,207, 82, 97,145,152, 70,139,165,186,154,138,247, 97,},
{  4, 24, 24, 23, 12,186, 92, 22,155,173,  0,115,133, 39, 64,116,109,},
{214,142,135,  5,161,241,116, 36, 54,120,236,216, 27,224,116,173,120,},
{215,190, 10,136,242, 23,151,215,191,180,108,255, 23,201, 60,156,131,},
{125,103, 72,245,191,236,207,164,187,100, 67, 70,115,205,217, 57,142,},
{107,198,165,135,172, 69,200,206, 25,  8,167,129, 11,231,191,127,152,},
{196, 77, 35,154, 53,238,103, 81, 59,192,176,238, 72, 40,153,103,162,},
{ 48, 76,169,193,137,128,218,236, 27, 21,253,100,103,164, 73,235,171,},
{ 29,159,190, 74,117,179,137,125, 89,132,100,222,249, 51,243,  4,181,},
{156, 92,109,141, 82, 15, 50,133,172,240,118,125, 85, 19,249,174,189,},
{212,106, 50,103,146,217, 37,218, 67,167,  5,186,168, 88,  3,228,197,},
{157,169,164,121,113,219, 49,175, 35,158,  5, 58,156, 63,  2,159,205,},
{180,104, 62,200,218,182,232, 48,246, 25, 24, 13,117, 72, 49,219,212,},
{182,251, 67,218, 83, 50,182, 25, 44,135,236,113,203, 40, 26,148,219,},
{254, 15, 88,140, 47, 55,168,232,244,145,134,237,  5,140,151,197,225,},
{ 47, 47,169, 72,157, 82, 41, 81, 18,134,151, 59,230,161,215,107,231,},
{241,128, 97, 29,172, 49,  2, 97,126, 69, 49,106,148,121, 94,131,236,},
{138, 21,115, 40, 91,179,125, 18,103,253, 37, 55,180, 39,  8,  9,241,},
{ 95,197, 24,183,  3,127, 92, 60, 22,236,210,110, 49,182, 10,250,244,},
{150,186,145,136,152,180,198,173,199, 82,185,134,145,220,247, 83,248,},
{198, 46,114, 42,253, 97,163,114, 33, 86,129,229,186,127,190, 20,251,},
{162,189,171, 29, 28,149,189,230,234, 11,181, 40, 69,248,171, 58,253,},
{ 90,239,233, 38,199,253, 14, 57, 65,240, 44, 41,137, 30,109,196,254,},
{207, 61, 69, 89,175,218,142, 30, 66, 29,239,107,203, 27, 15,177,255,},
};
// end of sin table

// cossinPi_core
// return isSin ? sin(a*pi) : sin(a*pi)
extfloat128_t extfloat128_t::cossinPi_core(const extfloat128_t& a, int isSin) {
  extfloat128_t x = mod_pow2(a, 1); // reduce to range (-2..+2)
  if (isfinite(x)) {
    uint32_t y_sign = x.m_sign & isSin;
    x.m_sign = 0;
    int ix  = static_cast<int>(x.scale_and_trunc(7)); // trunc(abs(x)*128)
    int ixs = ix + 64 - (isSin << 6); // shift cos by +pi/2
    y_sign ^= ((ixs & 128) != 0); // quadrants 2 & 3 are same as 0 & 1, but with opposit sign
    int x_sign = (ixs & 64) ? 63 : 0;
    ixs = (ixs ^ x_sign) & 63;    // qaudrant 1 => 0
    int longPoly = (ixs == 1);    // special handling for range [0.5/64..0.5/32) - don't use tab. Use longer poly instead
    int ia = ((ixs ^ longPoly)+1) >> 1;

    int dx = ((ix ^ longPoly)+1) >> 1;
    if (dx)
      x -=  dx * (1.0/64);
    x.m_sign ^= (x_sign & 1);

    // calculate ySin = sin(x)
    extfloat128_t xx = x * x;
    double xx_d = xx.convert_to_double();
    double ySin_p_d = sin_coeff_l[0];
    for (int i = 1; i < 3; ++i)
      ySin_p_d = ySin_p_d * xx_d + sin_coeff_l[i];
    extfloat128_t ySin_p(ySin_p_d), coef;
    for (int i = 0; i < 6; ++i) {
      coef.from_17bytes(sin_coeff_h[i]);
      ySin_p = fma(ySin_p, xx, coef); // use fma() instead of * and + not because it's more precise, but because it's a little faster
    }
    ySin_p *= x;
    coef.from_17bytes(sin_coeff_h[6]);
    extfloat128_t ySin = fma(x, coef, ySin_p); // Last step. Here fma() is used for better precision

    if (ia != 0) {
      // calculate yCosn = cos(x) - 1
      extfloat128_t yCosn = (-fma(ySin, ySin, -extfloat128_t::one())).sqrt() - extfloat128_t::one();
      // improve precision
      yCosn = -fma(ySin, ySin, yCosn*yCosn);
      // /= 2
      if (yCosn.m_exponent > yCosn.min_biased_exponent)
        yCosn.m_exponent -= 1;
      else
        yCosn = extfloat128_t::zero();

      if (ia < 32) {
        // res = tabSin*yCos + tabCos*ySin
        //     = tabSin*(yCosn+1) + tabCos*ySin = tabSin + tabSin*yCosn + tabCos*ySin
        //     = tabSinA + tabSinB + (tabSinA + tabSinB)*yCosn + tabCosA*ySin + tabCosB*ySin
        // where
        // t = ia * (0.5/32) = ia/64;
        extfloat128_t tabSin[2]; // [0] - fine part, [1] - coarse part
        extfloat128_t tabCos[2];
        from_17bytes_fix(tabSin, sin_tab[ia-1]);
        from_17bytes_fix(tabCos, sin_tab[31-ia]);
        extfloat128_t acc = fma(tabSin[1] + tabSin[0], yCosn, tabSin[0]);
        extfloat128_t ySin_l = fma(x, coef, -ySin) + ySin_p;
        acc = fma(tabCos[1], ySin_l, acc);
        acc = fma(tabCos[0], ySin,   acc);
        acc = fma(tabCos[1], ySin,   acc);
        ySin = acc + tabSin[1];
      } else {
        ySin = yCosn + extfloat128_t::one();
      }
    }
    ySin.m_sign = y_sign;
    return ySin;
  }
  return x;
}
