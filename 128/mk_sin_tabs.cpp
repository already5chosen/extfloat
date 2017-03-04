#include <stdio.h>
#include <iostream>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/math/constants/constants.hpp>

#include "extfloat128.h"
#include "extfloat128_to_bobinf.h"

template<typename T>
void convert_from_boost_bin_float(extfloat128_t* pDst, const T& src)
{
  typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<128, boost::multiprecision::backends::digit_base_2, void, boost::int64_t, INT_MIN, INT_MAX> > boost_quadfloat_t;
  boost_quadfloat_t qsrc(src);
  pDst->_set_exponent(qsrc.backend().exponent());
  pDst->m_significand[0] = reinterpret_cast<uint64_t*>(qsrc.backend().bits().limbs())[0];
  pDst->m_significand[1] = reinterpret_cast<uint64_t*>(qsrc.backend().bits().limbs())[1];
  pDst->m_sign = qsrc.backend().sign();
}

typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<320, boost::multiprecision::backends::digit_base_2, void, boost::int64_t, INT_MIN, INT_MAX> > boost_decafloat_t;
typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<136, boost::multiprecision::backends::digit_base_2, void, boost::int64_t, INT_MIN, INT_MAX> > boost_float136_t;


static void print_tab(const boost_decafloat_t* src, int n) {
  for (int i = 0; i < n; ++i) {
    extfloat128_t x;
    convert_from_boost_bin_float(&x, src[i]);
    uint8_t b[17];
    x.to_17bytes(b);
    printf("{");
    for (int k = 0; k < 17; ++k)
      printf("%3d,", b[k]);
    printf("}, // %.20e\n", x.convert_to_double());
  }
}

static void print_tab_fx(const boost_decafloat_t* src, int n) {
  for (int i = 0; i < n; ++i) {
    boost_float136_t x = boost_float136_t(src[i]);
    uint8_t b[17];
    for (int k = 0; k < 17; ++k) {
      x *= 256;
      boost_float136_t dig = trunc(x);
      x -= dig;
      b[16-k] = uint8_t(dig);
    }
    printf("{");
    for (int k = 0; k < 17; ++k)
      printf("%3d,", b[k]);
    printf("},\n");
  }
}

// // sin coefficients
// uint8_t tt[][17] = {
// { 61,185,157,118,229, 80, 80,229, 79,201,154,201, 87, 25,118,213, 86,},
// {144, 85,243,101,150, 94,  7,141,138, 28,186,170,248,220,214,183, 97,},
// {150,171, 50, 38,113, 37, 43,104,148, 25,107, 12,104, 26,122,244,104,},
// {  4,103,176,102, 20,118, 64,222, 17,251, 56, 68,239,167,131,241,113,},
// { 52,165,199, 86,239, 34, 99, 61,200, 13, 60,247, 67, 26, 60,168,120,},
// { 43,177,224, 50,178,191,112, 47,243, 45,236, 21,115,102,105,153,127,},
// { 22,134,209,  3, 79, 34, 52, 63,146, 14, 87,173, 59,227, 53,163,130,},
// { 38, 72,242, 42,113, 47,199, 93,245,149,242, 45, 49,231, 93,165,133,},
// {209, 28,220,128,139, 98,198,196, 52,194,104, 33,162,218, 15,201,130,},
// };
// // end of sin coefficients

int main()
{
  boost_decafloat_t x    = boost::math::constants::pi<boost_decafloat_t>();
  boost_decafloat_t xx   = -x * x;
  boost_decafloat_t fact = 1;
  const int N_SIN_COEFF = 9;
  boost_decafloat_t sin_coef[N_SIN_COEFF+1];
  for (int i = 0; i < N_SIN_COEFF; ++i) {
    sin_coef[N_SIN_COEFF-1-i] = x / fact;
    x *= xx;
    fact *= (i*2+2)*(i*2+3);
  }
  sin_coef[N_SIN_COEFF] = sin_coef[N_SIN_COEFF-1].convert_to<float>();
  sin_coef[N_SIN_COEFF-1] -= sin_coef[N_SIN_COEFF];
  printf("// sin coefficients\n");
  print_tab(sin_coef, N_SIN_COEFF+1);
  printf("// end of sin coefficients\n");

  // for (int i = 0; i < N_SIN_COEFF; ++i) {
    // extfloat128_t x;
    // x.from_17bytes(tt[i]);
    // boost_decafloat_t dx;
    // convert_to_boost_bin_float(&dx,x);
    // std::cout
      // << std::setprecision(36)
      // << sin_coef[i]
      // << " : "
      // << dx
      // << "\n";
  // }

  #if 0

  const int SIN_TAB_SZ = 32;
  double x_tab[SIN_TAB_SZ];
  for (int i = 0; i < 2; ++i)
    x_tab[i] = (1.25+0.5*i)*0.5/SIN_TAB_SZ;
  for (int i = 2; i < SIN_TAB_SZ-2; ++i)
    x_tab[i] = (i+0.5)*0.5/SIN_TAB_SZ;
  for (int i = SIN_TAB_SZ-2; i < SIN_TAB_SZ; ++i)
    x_tab[i] = 0.5 - x_tab[SIN_TAB_SZ-1-i];

  boost_decafloat_t sin_tab[SIN_TAB_SZ];
  for (int i = 0; i < SIN_TAB_SZ; ++i)
    sin_tab[i] = sin(boost::math::constants::pi<boost_decafloat_t>()*x_tab[i]);

  printf("// sin table\n");
  print_tab_fx(sin_tab, SIN_TAB_SZ);
  printf("// end of sin table\n");

  #else

  const int SIN_TAB_SZ = 32;
  double x_tab[SIN_TAB_SZ];
  for (int i = 0; i < SIN_TAB_SZ; ++i)
    x_tab[i] = (i+0.5)*0.5/SIN_TAB_SZ;

  boost_decafloat_t sin_tab[SIN_TAB_SZ];
  for (int i = 0; i < SIN_TAB_SZ; ++i)
    sin_tab[i] = sin(boost::math::constants::pi<boost_decafloat_t>()*x_tab[i]);

  printf("// sin table\n");
  print_tab_fx(sin_tab, SIN_TAB_SZ);
  printf("// end of sin table\n");

  #endif

  return 0;
}