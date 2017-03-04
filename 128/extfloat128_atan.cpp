#include "extfloat128.h"
#include <stdio.h>

static const double atan_estimate_coeff[2] = {
   -4.55189784406698834118e+01,
    1.62644927809243569072e+02,
};
static const uint8_t tan_tab[40][18] = {
 {143, 74,187, 96,216, 23,164,108,161,147,253,211, 95,142, 64,226,116,198,}, // 2.76186733959772959257e-02 = tan( 0.5625*pi/64) atan() err = -1.18e-41, rel  -58.58*2**(-135)
 { 36, 32,185, 75,152,227,100,178,208,110,  3,159,253, 86, 72,138,118,171,}, // 3.37603948664203373786e-02 = tan( 0.6875*pi/64) atan() err = -2.12e-41, rel  -85.80*2**(-135)
 {215,211,105,  2, 59,222,141, 65,150,231,142,130,252, 17,115,163,118,  1,}, // 3.99046614258258061647e-02 = tan( 0.8125*pi/64) atan() err = -1.13e-43, rel   -0.39*2**(-135)
 {  4, 12, 80,108,247, 70,115,213,  3,214, 38,175,224,244,160,188,118, 29,}, // 4.60519376310430592070e-02 = tan( 0.9375*pi/64) atan() err = 9.83e-42, rel   29.23*2**(-135)
 {198,246,237, 94,244, 99,251,244,174, 82,233,230,106,197,108,226,118, 30,}, // 5.52795135139894783238e-02 = tan( 1.1250*pi/64) atan() err = 1.23e-41, rel   30.43*2**(-135)
 { 44, 90,212, 28, 40, 42,166,137,239,164,255, 43,224,187,112,138,120, 45,}, // 6.75978353297067152683e-02 = tan( 1.3750*pi/64) atan() err = 2.23e-41, rel   45.22*2**(-135)
 {  8, 16,245, 99,119,118, 98,206,156,239, 86,168,130,206,181,163,120,213,}, // 7.99366124747779849269e-02 = tan( 1.6250*pi/64) atan() err = -2.57e-41, rel  -44.01*2**(-135)
 { 97,170,199,205,244, 89, 59,229,112,188, 15, 21,162,149,  7,189,120, 67,}, // 9.22996225441535811385e-02 = tan( 1.8750*pi/64) atan() err = 4.51e-41, rel   67.07*2**(-135)
 {105,216, 64,205, 10,150,164, 10,120, 13, 64,254,187,113, 30,227,120, 33,}, // 1.10897911595913029048e-01 = tan( 2.2500*pi/64) atan() err = 2.63e-41, rel   32.58*2**(-135)
 { 56,184,162, 66,172,114,  3,132,255,126,243, 72, 45,108, 19,139,122,201,}, // 1.35816278709387727730e-01 = tan( 2.7500*pi/64) atan() err = -5.56e-41, rel  -56.31*2**(-135)
 {107,250, 20, 46,239,218,185,101, 53,187,114,122,166, 83,195,164,122, 42,}, // 1.60901362453489155113e-01 = tan( 3.2500*pi/64) atan() err = 4.94e-41, rel   42.34*2**(-135)
 {127, 26,162,179,  8,121, 36,201,253, 64,203,218,167, 98,167,190,122,244,}, // 1.86185399527583728130e-01 = tan( 3.7500*pi/64) atan() err = -1.76e-41, rel  -13.09*2**(-135)
 { 78,234,112,158,196, 88,198,183,148,146,103, 90, 90,103,242,229,122, 21,}, // 2.24557509317129311288e-01 = tan( 4.5000*pi/64) atan() err = 3.39e-41, rel   21.02*2**(-135)
 { 49,  2,  7,138,101, 99, 79,218,203, 73,191,172,233,129,176,141,124,182,}, // 2.76737270140414326480e-01 = tan( 5.5000*pi/64) atan() err = -1.48e-40, rel  -74.90*2**(-135)
 {129, 51, 38, 56, 10, 15,199, 46,125,187,216,128, 11, 87, 36,169,124, 52,}, // 3.30355377344333900336e-01 = tan( 6.5000*pi/64) atan() err = 1.20e-40, rel   51.63*2**(-135)
 { 32,127, 46,214,237,177, 89,216,140, 14,129,158,181, 12,128,197,124,225,}, // 3.85742566271121245514e-01 = tan( 7.5000*pi/64) atan() err = -8.73e-41, rel  -32.46*2**(-135)
 {230,203,188,176, 59,133,254,213,127,176, 61,243,194, 56,244,226,124,  2,}, // 4.43269513890864330641e-01 = tan( 8.5000*pi/64) atan() err = 6.82e-42, rel    2.24*2**(-135)
 {196, 97,234,242,208,119,221,249, 23,186,231,211,218, 12,220,128,126,218,}, // 5.03357699799294233678e-01 = tan( 9.5000*pi/64) atan() err = -1.35e-40, rel  -39.49*2**(-135)
 { 45,231,200,204,195,228,190,181,107, 40,193, 35,120,175,  5,145,126,185,}, // 5.66493002730343975237e-01 = tan(10.5000*pi/64) atan() err = -2.71e-40, rel  -71.95*2**(-135)
 {128, 41, 67,241, 53, 98, 80, 87,176,212, 46,231,220, 54, 28,162,126,234,}, // 6.33243016177569173486e-01 = tan(11.5000*pi/64) atan() err = -9.51e-41, rel  -23.05*2**(-135)
 {203, 41,193,191, 62,212,194,155, 67,177,237,168,163,168, 75,180,126, 40,}, // 7.04279460865044226736e-01 = tan(12.5000*pi/64) atan() err = 1.77e-40, rel   39.53*2**(-135)
 {  7,239,131,140,173, 86,142, 64,135, 86,247,194,223,203,200,199,126, 47,}, // 7.80407659653943652778e-01 = tan(13.5000*pi/64) atan() err = 2.29e-40, rel   47.19*2**(-135)
 { 44,180,171,122,213,120,229,199,111,145,215, 96, 12,190,211,220,126,253,}, // 8.62605932256739871278e-01 = tan(14.5000*pi/64) atan() err = -1.93e-41, rel   -3.71*2**(-135)
 { 61,238, 83, 25,228,111, 87,189, 34, 51,184, 72,126,117,187,243,126, 40,}, // 9.52079146700925305069e-01 = tan(15.5000*pi/64) atan() err = 2.25e-40, rel   40.49*2**(-135)
 {207,192, 99,199,224,244, 76, 12,147,202,228, 65,132, 78,113,134,128,  8,}, // 1.05033284623985978534e+00 = tan(16.5000*pi/64) atan() err = 4.51e-41, rel    7.62*2**(-135)
 {199,171,  7,150,134, 89,201, 34,  3, 73,120,124,237, 55, 99,148,128,202,}, // 1.15927790733343472063e+00 = tan(17.5000*pi/64) atan() err = -3.43e-40, rel  -54.70*2**(-135)
 {115, 33,123,123,156, 73,253,214,240, 99, 22,250,197, 79,  4,164,128, 44,}, // 1.28138158003655444617e+00 = tan(18.5000*pi/64) atan() err = 2.90e-40, rel   43.64*2**(-135)
 { 49, 17,195,233, 59, 88,132, 66, 94,253,143, 50, 49,252,190,181,128, 21,}, // 1.41989090349409252667e+00 = tan(19.5000*pi/64) atan() err = 1.44e-40, rel   20.60*2**(-135)
 { 21,  3,245, 52,  0, 71,215,125,153, 83,145, 16,163, 83, 34,202,128,255,}, // 1.57917256796020888387e+00 = tan(20.5000*pi/64) atan() err = -1.46e-41, rel   -1.99*2**(-135)
 {125,110,106, 71,141, 53, 91,114, 20,122,224, 53,  4,156,243,225,128,  8,}, // 1.76524687009419145589e+00 = tan(21.5000*pi/64) atan() err = 5.85e-41, rel    7.59*2**(-135)
 {147,189,185, 70, 75,243,124, 18,183, 81,122,182,214,213, 74,254,128,235,}, // 1.98665879234336495429e+00 = tan(22.5000*pi/64) atan() err = -1.76e-40, rel  -21.84*2**(-135)
 {108, 24, 11,244,215,255, 43,224,214,110,156, 63, 53,182, 97,144,130,232,}, // 2.25596385192915871443e+00 = tan(23.5000*pi/64) atan() err = -2.13e-40, rel  -25.22*2**(-135)
 {130,144,161,174, 16,163,140,116,250,124, 50,240, 63,236,233,165,130,  8,}, // 2.59240251773807273139e+00 = tan(24.5000*pi/64) atan() err = 6.60e-41, rel    7.51*2**(-135)
 { 67, 13,106,246,164, 29, 55,174,217,192,237,135,107, 19,187,193,130,250,}, // 3.02704320431777418321e+00 = tan(25.5000*pi/64) atan() err = -6.08e-41, rel   -6.64*2**(-135)
 {152,  9, 41,246, 53, 42,142,253, 93, 56, 63,137, 41, 43, 68,231,130,  5,}, // 3.61353568130742841547e+00 = tan(26.5000*pi/64) atan() err = 4.50e-41, rel    4.73*2**(-135)
 {  4, 24, 27,188,  2,121,137,140, 17,246,178,138,243,161,128,142,132, 11,}, // 4.45320222441441071481e+00 = tan(27.5000*pi/64) atan() err = 1.05e-40, rel   10.61*2**(-135)
 {135, 43,121, 49, 32,149,244, 58,222,133, 13, 70,200,168,107,184,132, 10,}, // 5.76314200511880958544e+00 = tan(28.5000*pi/64) atan() err = 1.01e-40, rel    9.84*2**(-135)
 { 61,194, 17,125,171,151,137, 92, 43, 14, 63, 92,155,125,185,129,134,  9,}, // 8.10778580367690793196e+00 = tan(29.5000*pi/64) atan() err = 9.79e-41, rel    9.25*2**(-135)
 { 50, 19,123, 63, 94,242, 24,229,103,249,123,233,  1, 30,232,216,134,  3,}, // 1.35566692423524255418e+01 = tan(30.5000*pi/64) atan() err = 3.62e-41, rel    3.31*2**(-135)
 {116, 80, 53, 47,202,130,123, 68,226, 57, 85, 37,175, 34,241,162,138,  1,}, // 4.07354838720833001275e+01 = tan(31.5000*pi/64) atan() err = 1.23e-41, rel    1.09*2**(-135)
};
static const double atan_coeff_l[3] = {
    1.51086894537297775076e-02,
   -1.67530771282955841250e-02,
    1.87241108873836895055e-02,
};
static const uint8_t atan_coeff_h[9][17] = {
 {182,223, 23, 29, 31,219,141,  2,227, 22,187,201,151,242,214,173,117,}, // -2.12206590788847176e-02
 {183, 81, 99,210,  0,189,104,192,246, 96,103,155, 17,142,149,200,116,}, //  2.44853758602915778e-02
 { 92,172,177, 19,197,199,216,106,222,123, 52,230,113,214, 13,237,117,}, // -2.89372623803446048e-02
 {175,183,172, 43, 88,177,141,252,207, 18, 32, 41, 41,202,221,144,118,}, //  3.53677651315322944e-02
 {217,176,215,  3, 93,173,154,251, 47, 24,224, 52,199,186, 65,186,119,}, // -4.54728408833986672e-02
 { 95,140,145, 93, 14, 19,185,201, 84,119, 54,216,241, 53, 97,130,120,}, //  6.36619772367581355e-02
 {204,124, 70,156,194, 31,223,250, 55, 28,176,189, 61,175, 76,217,121,}, // -1.06103295394596897e-01
 {138,192,182,129,187,105,234,163,175, 78,248, 83, 42,136,156,220, 74,}, //  1.28412766334518303e-08
 {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,131,249,162,124,}, //  3.18309873342514038e-01
};


enum {
  xFlag_none        = 0,
  xFlag_x_present   = 1,
  xFlag_x_negative  = 2,
  xFlag_dxx_present = 4,
  xFlag_dxx_acos    = 8,
  xFlag_positive = xFlag_x_present,
  xFlag_negative = xFlag_x_present + xFlag_x_negative,
  xFlag_asin     = xFlag_x_present + xFlag_dxx_present,
  xFlag_acos_p   = xFlag_x_present + xFlag_dxx_present + xFlag_dxx_acos,
  xFlag_acos_n   = xFlag_x_present + xFlag_dxx_present + xFlag_dxx_acos + xFlag_x_negative,
};
// atan(x)/pi
// 2**(-127) <= y/x < 2**132
extfloat128_t atanPi_core(extfloat128_t y, const extfloat128_t& x, int xFlag, double dxx)
{
  const double dInvPi = 1.0/3.1415926535897932384626433832795;
  // estimate atan(y)*256/pi with sufficient precision for lookup in 32-entry table
  double da = y.convert_to_double();
  double dx = 1.0;
  double dxy = dxx;
  if (xFlag != xFlag_none) {
    dx = x.convert_to_double();
    da /= dx;
    if (xFlag & xFlag_dxx_present) {
      dxy *= (dInvPi * 0.5);
      if (xFlag == xFlag_asin)
        dxy *= da;
      else
        dxy /= da;
    }
  }
  static const double tan30Deg = 0.57735026918962576450914878050196;
  static const double tan60Deg = 1.7320508075688772935274463415059;
  double adj = 0;
  if (da > tan30Deg) {
    da = (da - tan60Deg)/(1 + da*tan60Deg);
    adj = 512.0/3;
  }
  double estAtan = (da*da*atan_estimate_coeff[0]+atan_estimate_coeff[1])*da + adj;

  // use estimate to peek entry in tan_tab, with value which is close to y
  int ix = int(estAtan);
  double dAdj  = 0;
  double atan0 = 0;
  if (ix > 3) {
    unsigned iTab;
    if (ix < 8) {
      iTab = ix - 4;
      ix   = iTab*2 + 9;
    } else if (ix < 16) {
      iTab = unsigned(ix)/2;
      ix   = iTab*4 + 2;
    } else if (ix < 32) {
      iTab = unsigned(ix)/4 + 4;
      ix   = iTab*8 - 28;
    } else {
      ix = (ix <= 255) ? ix : 255;
      iTab = unsigned(ix)/8 + 8;
      ix   = iTab*16 - 120;
    }
    extfloat128_t tabVal;
    tabVal.from_17bytes(tan_tab[iTab]);
    int iAdj = static_cast<int8_t>(tan_tab[iTab][17]);
    atan0 = (1.0/1024)*ix;
    const double TWO_POW_MIN_135 = 2.2958874039497802890014385492622e-41; // 1/(128*2^128)
    dAdj  = iAdj*TWO_POW_MIN_135*atan0;
    #if 0
    y = (y - tabVal)/fma(y, tabVal, extfloat128_t::one()); // y = tan(atan(y)-atan(tabVal)),
                                                           // where atan(tabVal)=(ix/512)*(pi/2)
    #else
    extfloat128_t num = (xFlag != xFlag_none) ? fma(tabVal, -x, y) : y - tabVal;
    extfloat128_t den = fma(y, tabVal, x);
    extfloat128_t den_lo;
    if (ix < 144)
      den_lo = fma(y, tabVal, x-den);
    y = num/den;                                         // y = tan(atan(y)-atan(tabVal)),
                                                         // where atan(tabVal)=(ix/512)*(pi/2)
    if (ix < 256) {
      extfloat128_t xxErr = fma(y, den, -num);
      if (ix < 144)
        xxErr = fma(y, den_lo, xxErr);
      double xErr = (xxErr/den).convert_to_double();
      dAdj -= xErr*dInvPi;
    }
    #endif
  } else if (xFlag != xFlag_none) {
    extfloat128_t num = y;
    y /= x;
    dAdj = fma(y, x, -num).convert_to_double()*(-dInvPi)/dx;
  }
  if (xFlag & xFlag_dxx_present)
    dAdj += dxy;

  // now y is in range (-0.02542851375449105972..0.02542851375449105972)
  // i.e. range slightly wider than [-tan(pi/128)..tan(pi/128)]
  //
  // Calculate atan(y)/pi by polynomial
  extfloat128_t yy = y * y;
  double yy_d = yy.convert_to_double();
  double yAtan_p_d = atan_coeff_l[0];
  for (int i = 1; i < 3; ++i)
    yAtan_p_d = yAtan_p_d * yy_d + atan_coeff_l[i];
  extfloat128_t yAtan_p(yAtan_p_d), coef;
  for (int i = 0; i < 8; ++i) {
    coef.from_17bytes(atan_coeff_h[i]);
    yAtan_p = fma(yAtan_p, yy, coef); // use fma() instead of * and + not because it's more precise, but because it's a little faster
  }
  coef.from_17bytes(atan_coeff_h[8]);
  // Last step. Here fma() is used for better precision
  extfloat128_t ret = (ix > 3 || xFlag != xFlag_none) ?
    fma(y, yAtan_p, dAdj) :
    y*yAtan_p;
  ret = fma(y, coef, ret);
  if (ix > 3) {
    return (xFlag & xFlag_x_negative) ? (1-atan0) - ret : atan0 + ret;
  } else {
    return (xFlag & xFlag_x_negative) ? extfloat128_t::one() - ret : ret;
  }
}

extfloat128_t atanPi(const extfloat128_t& a) // atan(a)/pi
{
  uint32_t aExp = a.m_exponent;
  if (aExp > extfloat128_t::exponent_bias - 128 && aExp < extfloat128_t::exponent_bias + 132) {
    // non-trivial case
    extfloat128_t ret = atanPi_core(abs(a), extfloat128_t::one(), xFlag_none, 0);
    ret.m_sign = a.m_sign;
    return ret;
  } else {
    if (aExp < extfloat128_t::exponent_bias) {
      return fma(a, extfloat128_t::invPi(), a*extfloat128_t::invPi_lo());
    } else if (!isnan(a)) {
      extfloat128_t ret = 0.5;
      ret.m_sign = a.m_sign;
      return ret;
    } else {
      return a;
    }
  }
}

extfloat128_t asinPi(const extfloat128_t& a) // asin(a)/pi
{
  uint32_t aExp = a.m_exponent;
  if (aExp > extfloat128_t::exponent_bias - 66 && aExp < extfloat128_t::exponent_bias) {
    // non-trivial case
    extfloat128_t xx = abs(fma(a, a, extfloat128_t::one(1)));
    extfloat128_t x  = xx.sqrt();
    double dxx = 0;
    int xFlag = xFlag_positive;
    if (aExp < extfloat128_t::exponent_bias - 1 || a.m_significand[1] < (uint64_t(255) << 56)) {
      xx.m_significand[0] = 0;
      dxx = (fma(x, x, -xx) + fma(a, a, xx - extfloat128_t::one())).convert_to_double();
      // xx -= extfloat128_t::one();
      // dxx = (fma(x, x, extfloat128_t::one(1) - xx) + fma(a, a, xx)).convert_to_double();
      xFlag = xFlag_asin;
    }
    extfloat128_t ret = atanPi_core(abs(a), x, xFlag, dxx);
    ret.m_sign = a.m_sign;
    return ret;
  } else {
    if (aExp < extfloat128_t::exponent_bias) {
      return fma(a, extfloat128_t::invPi(), a*extfloat128_t::invPi_lo());
    } else {
      if (aExp != extfloat128_t::exponent_bias ||
          a.m_significand[0] != 0 ||
          a.m_significand[1] != (uint64_t(1) << 63) ) {
        if (isnan(a))
          return a;
        else
          return extfloat128_t::nan(0);
      }
      // a = +1 or -1
      return extfloat128_t::pow2(-1, a.m_sign); // 0.5*sign(a)
    }
  }
}

extfloat128_t acosPi(const extfloat128_t& a) // acos(a)/pi
{
  uint32_t aExp = a.m_exponent;
  if (aExp > extfloat128_t::exponent_bias - 45 && aExp < extfloat128_t::exponent_bias) {
    //return extfloat128_t::pow2(-1) - asinPi(a);
    // non-trivial case
    extfloat128_t yy = abs(fma(a, a, extfloat128_t::one(1)));
    extfloat128_t y  = yy.sqrt();
    yy -= extfloat128_t::one();
    double dxx = -((fma(y, y, extfloat128_t::one(1) - yy) + fma(a, a, yy)).convert_to_double());
    int xFlag = a.m_sign==0 ? xFlag_acos_p : xFlag_acos_n;
    extfloat128_t ret = atanPi_core(y, abs(a), xFlag, dxx);
    return ret;
  } else {
    if (aExp < extfloat128_t::exponent_bias) {
      return fma(a, -extfloat128_t::invPi(), extfloat128_t::pow2(-1)); // 0.5 - a/pi
    } else {
      if (aExp != extfloat128_t::exponent_bias ||
          a.m_significand[0] != 0 ||
          a.m_significand[1] != (uint64_t(1) << 63) ) {
        if (isnan(a))
          return a;
        else
          return extfloat128_t::nan(0);
      }
      // a = +1 or -1
      return a.m_sign==0 ? extfloat128_t::zero() : extfloat128_t::one();
    }
  }
}

extfloat128_t atan2Pi(const extfloat128_t& y, const extfloat128_t& x) // atan(y,x)/pi
{
  extfloat128_t ret;
  if (isfinite(y)) {
    if (isfinite(x)) {
      uint32_t yExp = y.m_exponent;
      uint32_t xExp = x.m_exponent;
      if (yExp != extfloat128_t::zero_biased_exponent && xExp != extfloat128_t::zero_biased_exponent) {
        int64_t dExp = int64_t(yExp) - xExp;
        if (dExp < 132) {
          extfloat128_t num = y; num.m_sign = 0;
          extfloat128_t den = x; den.m_sign = 0;
          if (dExp > -66) {
            // non-trivial case
            num.m_exponent = extfloat128_t::exponent_bias + uint32_t(dExp);
            den.m_exponent = extfloat128_t::exponent_bias;
            ret = atanPi_core(num, den, x.m_sign == 0 ? xFlag_positive : xFlag_negative, 0);
          } else {
            // very close to X-axis
            extfloat128_t rat = num/den;
            extfloat128_t ratErr = fma(rat, -den, num)/den;
            ret = fma(rat, extfloat128_t::invPi(), fma(rat,extfloat128_t::invPi_lo(),ratErr*extfloat128_t::invPi()));
            if (x.m_sign != 0)
              ret = extfloat128_t::one() - ret;
          }
        } else {
          ret = 0.5; // on Y-axis
        }
      } else { // y==0 || x==0
        if (yExp != extfloat128_t::zero_biased_exponent)
          ret = 0.5; // y != 0, x==0
        else if (xExp != extfloat128_t::zero_biased_exponent)
          ret = (x.m_sign==0) ? extfloat128_t::zero() : extfloat128_t::one(); // y==0, x != 0
        else  // both y and x == 0
          return extfloat128_t::nan(0);
      }
    } else {
      if (isnan(x))
        return x;
      // x = inf
      ret = (x.m_sign==0) ? extfloat128_t::zero() : extfloat128_t::one();
    }
  } else {
    if (isnan(y))
      return y;
    // y = inf
    if (isfinite(x)) {
      ret = 0.5;
    } else {
      if (isnan(x))
        return x;
      // both y and x = inf
      return extfloat128_t::nan(0);
    }
  }
  ret.m_sign = y.m_sign;
  return ret;
}
