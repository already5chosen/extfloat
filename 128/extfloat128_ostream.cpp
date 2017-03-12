#define USE_EXTFLOAT128_T_OSTREAM
#include "extfloat128.h"
#include <cmath>
#include <cfloat>
#include <cstdio>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include "extfloat128_to_bobinf.h"
#include "bobinf_to_extfloat128.h"

typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<192+48*1, boost::multiprecision::backends::digit_base_2> > boost_float192_t;

class extfloat128_ostream_init_t {
public:
  extfloat128_t oneE18[2];
  extfloat128_t oneE36[2];
  extfloat128_t oneExx[6][2]; // 1e1, 1e2, 1e4, 1e8, 1e16, 1e32
  int64_t       oneEyy[18];   // 1e1, 1e2, 1e3, ... 1e18
  extfloat128_ostream_init_t() {
    extfloat128_t x = 10;
    for (int i = 0; i < 6; ++i) {
      oneExx[i][0] = x;
      oneExx[i][1] = inv(x);
      x *= x;
    }
    oneE18[0] = oneExx[1][0]*oneExx[4][0];
    oneE36[0] = oneExx[2][0]*oneExx[5][0];
    oneE18[1] = inv(oneE18[0]);
    oneE36[1] = inv(oneE36[0]);
    int64_t y = 10;
    for (int i = 0; i < 18; ++i) {
      oneEyy[i] = y;
      y *= 10;
    }
  }
  extfloat128_t inv(const extfloat128_t& x) {
    extfloat128_t y = extfloat128_t::one() / x;
    if (fma(x, y, extfloat128_t::one(1)) > extfloat128_t::zero())
      y = nextafter(y, 0);
    return y;
  }
};
static extfloat128_ostream_init_t tab;

static void pow5(extfloat128_t dst[2], int e)
{
  boost_float192_t b = pow(boost_float192_t(5), e);
  convert_from_boost_bin_float(&dst[0], b);
  boost_float192_t b1;
  convert_to_boost_bin_float(&b1, dst[0]);
  b -= b1;
  convert_from_boost_bin_float(&dst[1], b);
}

// input:  x < 1e72
// result: remdiv[0] = x % 1e36, remdiv[1] = floor(x/1e36)
static void remdivE36(extfloat128_t remdiv[2], extfloat128_t::acc_t& x)
{
  remdiv[1] = extfloat128_t::zero();
  for (;;) {
    remdiv[0] = x.trunc();
    if (remdiv[0] < tab.oneE36[0])
      break;
    extfloat128_t t = trunc(multiply_rtz(tab.oneE36[1], remdiv[0]));
    if (!is_zero(t)) {
      remdiv[1] += t;
      x.msub(tab.oneE36[0], t);
    } else {
      remdiv[1] += extfloat128_t::one();
      x -= tab.oneE36[0];
    }
  }
}

// input:  x < 1.1e36
// result: remdiv[0] = x % 1e18, remdiv[1] = floor(x/1e18)
static void remdivE18(int64_t remdiv[2], extfloat128_t& x)
{
  remdiv[1] = 0;
  while (x >= tab.oneE18[0]) {
    extfloat128_t t = trunc(multiply_rtz(tab.oneE18[1], x));
    if (!is_zero(t)) {
      remdiv[1] += t.convert_to_int64();
      x = fma(tab.oneE18[0], -t, x);
    } else {
      remdiv[1] += 1;
      x -= tab.oneE18[0];
    }
  }
  remdiv[0] = x.convert_to_int64();
}

#if 0
static int64_t toDecimalStep(extfloat128_t::acc_t& acc, const extfloat128_t& scale, const extfloat128_t& invScale)
{
  extfloat128_t::acc_t tmp;
  tmp.clear();
  extfloat128_t y = acc.trunc();
  while (y >= scale) {
    extfloat128_t::eval_multiply_rtz(y, y, invScale);
    y = trunc(y);
    if (is_zero(y))
      y = extfloat128_t::one();
    tmp += y;
    acc.msub(y, scale);
    y = acc.trunc();
  }
  acc = tmp;
  return y.convert_to_int64();
}
#endif

// return floor(log10(abs(x))) or floor(log10(abs(x)))-1
static int32_t log10_estimate(const extfloat128_t& x)
{
  extfloat128_t mx = x;
  mx.m_sign = 0;
  mx.m_exponent = mx.exponent_bias;
  double md = mx.convert_to_double();
  static const double LOG10_2 = 0.30102999566398119521373889472449;
  return static_cast<int32_t>(floor(((md-1 - 1e-5) + x._get_exponent()) * LOG10_2));
}

int64_t extfloat128_t::convert_to_int64() const
{
  boost_float192_t b;
  convert_to_boost_bin_float(&b, *this);
  return b.convert_to<int64_t>();
}

#if 0
static void pp(extfloat128_t::acc_t x)
{
  printf("%d %lld %016llx %016llx %016llx %016llx %016llx %016llx = %g\n"
  , x.m_sign, x.m_exponent-x.exponent_bias
  ,x.m_significand[5]
  ,x.m_significand[4]
  ,x.m_significand[3]
  ,x.m_significand[2]
  ,x.m_significand[1]
  ,x.m_significand[0]
  ,x.trunc().convert_to_double()
  );
}
#endif

static int64_t divBy10andRound(int64_t x, int roundDir)
{
  int64_t xi = x / 10;
  int xr = x - xi*10;
  int xr4 = xr*4 - roundDir*2 + (xi & 1);
  return xi + (xr4 > 20);
}

int dbg = 0;

static int scalebyPow10(extfloat128_t::acc_t* dst, const extfloat128_t& x, int32_t scale10)
{
  extfloat128_t mult[2];
  pow5(mult, scale10);
  dst->mulx(mult[0], x);
  dst->madd(mult[1], x);
  dst->normalize();
  dst->m_exponent += scale10;
  dst->m_sign = 0;
  if (dbg) {
    extfloat128_t xx = dst->trunc();
    extfloat128_t xi = trunc(xx);
    extfloat128_t xr = xx - xi;
    printf("%18e %.15f scale10=%d\n", xi.convert_to_double(), xr.convert_to_double(), scale10);
  }
  return dst->round_to_nearest_tie_to_even();
}

// scale10 in range [-38..-1]
// return value: -1 if rounded toward 0, 1 if rounded away from zero, 0 if unchanged
static int scalebyPow10_smallNeg(extfloat128_t::acc_t* dst, extfloat128_t x, int scale10)
{
  x.m_sign = 0;
  *dst = x;
  int roundDir = dst->round_to_nearest_tie_to_even();
  x = dst->trunc();
  for (int powlog = 5; scale10 < 0; --powlog) {
    int pow = 1 << powlog;
    if (pow + scale10 <= 0) {
// if (dbg) printf("scale10=%d, x=%e * %e / %e\n", scale10, x.convert_to_double(), tab.oneExx[powlog][0].convert_to_double(), tab.oneExx[powlog][1].convert_to_double());
      scale10 += pow;
      // divide
      extfloat128_t div = extfloat128_t::zero();
      extfloat128_t rem = x;
      while (rem >= tab.oneExx[powlog][0]) {
        extfloat128_t inc = trunc(multiply_rtz(rem, tab.oneExx[powlog][1]));
        if (is_zero(inc))
          inc = extfloat128_t::one();
        div += inc;
        rem = fma(-div, tab.oneExx[powlog][0], x);
      }
      rem = round(rem);
if (dbg) printf(" x %e; div %e; rem %e %d %d scale %e: %d\n", x.convert_to_double(), div.convert_to_double(), rem.convert_to_double(), is_zero(rem), rem._get_exponent(), tab.oneExx[powlog][0].convert_to_double(), roundDir);
      if (!is_zero(rem)) {
        extfloat128_t rem2 = extfloat128_t::ldexp(rem, 1);
        if (rem2 >= tab.oneExx[powlog][0]) {
          if (rem2 > tab.oneExx[powlog][0]) {
            roundDir = 1;
          } else if (roundDir == 0) {
            // a tie
            roundDir = is_zero(mod_pow2(div, 1)) ? -1 : 1; // round to even
          } else {
            roundDir = -roundDir;
          }
          if (roundDir > 0)
            div += extfloat128_t::one();
        } else {
          roundDir = -1;
        }
      }
if (dbg) printf(" x %e; div %e; rem %e : %d\n", x.convert_to_double(), div.convert_to_double(), rem.convert_to_double(), roundDir);
      x = div;
    }
  }
if (dbg) printf("scale10=%d, x=%e\n", scale10, x.convert_to_double());
  *dst = x;
  return roundDir;
}

std::ostream& operator<<(std::ostream& os, const extfloat128_t& x)
{
  const int def_prec   = 41;
  const int max_prec   = 256;
  const int max_digits = 64;
  std::streamsize prec = os.precision();
  if (prec <= 0)
    prec = def_prec;
  if (prec > max_prec)
    prec = max_prec;

  int digits = (prec <= max_digits) ? prec : max_digits;
  char buf[max_prec + 32];
  if (isnormal(x)) {
    int32_t exp10 = log10_estimate(x);
    int32_t scale10 = digits - exp10; // sometimes we'll have one more digit than requested
    extfloat128_t::acc_t acc;
//    int roundDir = (scale10 < 0 && scale10 > -39 && exp10 < 39) ?
    int roundDir = (scale10 < 0 && scale10 > -39 && exp10 < 40 && digits < 38) ?
      scalebyPow10_smallNeg(&acc, x, scale10) :
      scalebyPow10         (&acc, x, scale10);
    if (dbg) {
      extfloat128_t::acc_t tmp = acc;
      printf("roundDir %d. acc=%e\n", roundDir, tmp.trunc().convert_to_double());
    }
    extfloat128_t digE36[2];
    int nDigE18 = digits/18 + 1;
    int nDigE36 = (nDigE18+1)/2;
    if (nDigE36 > 1)
      remdivE36(digE36, acc);
    else
      digE36[0] = acc.trunc();

    int64_t digE18[4];
    nDigE36 = nDigE18/2;
    for (int i = 0; i < nDigE36; ++i)
      remdivE18(&digE18[i*2], digE36[i]);
    if (nDigE36*2 < nDigE18)
      digE18[nDigE36*2] = digE36[nDigE36].convert_to_int64();
    if (dbg) printf("digits %d. nDigE18 %d: %lld %lld %lld %lld scale10=%d exp10=%d\n", digits, nDigE18, digE18[0], digE18[1], digE18[2], digE18[3], scale10, exp10);

    int msDigits = digits + 1 - (nDigE18-1)*18;
    bool extraDigit = (digE18[nDigE18-1] >= tab.oneEyy[msDigits-1]); // 1 MS word contains one more decimal digit than requested
    if (extraDigit) {
      exp10 += 1;
      digE18[0] = divBy10andRound(digE18[0], roundDir);
      if (nDigE18 > 1) {
        const int64_t oneE17 = tab.oneEyy[17-1];
        if (digE18[0] >= oneE17) {
          digE18[0] -= oneE17;
          // propagate carry
          digE18[1] += 1;
          const int64_t oneE18 = tab.oneEyy[18-1];
          for (int i = 1; i < nDigE18-1 && digE18[i] >= oneE18; ++i) {
            digE18[i]   -= oneE18;
            digE18[i+1] += 1;
          }
        }
      }
    }
    if (dbg) printf("Digits %d. nDigE18 %d: %lld %lld %lld %lld scale10=%d exp10=%d\n", digits, nDigE18, digE18[0], digE18[1], digE18[2], digE18[3], scale10, exp10);

    char* p = &buf[2];
    p += sprintf(p, "%lld", digE18[nDigE18-1]);
    for (int i = nDigE18-2; i >= 0; --i)
      p += sprintf(p, ((i == 0 && extraDigit) ? "%017lld" :  "%018lld"), digE18[i]);
    if (digits < prec) {
      memset(p, '0', prec-digits);
      p += prec-digits;
    }
    *p++ = 'e';
    *p++ = (exp10 < 0) ? '-' : '+';
    sprintf(p, "%02d", abs(exp10));
    buf[1] = buf[2];
    buf[2] = '.';
    buf[0] = '-';
    os << (x.m_sign ? &buf[0] : &buf[1]);
  } else {
    if (isnan(x)) {
      os << "nan";
    } else if (is_zero(x)) {
      char* p = buf;
      if (x.m_sign)
        *p++ = '-';
      *p++ = '0';
      if (prec > 0) {
        *p++ = '.';
        memset(p, '0', prec);
        p += prec;
      }
      strcpy(p, "e+00");
      os << buf;
    } else {
      os << (x.m_sign ? "-inf" : "inf");
    }
  }
  return os;
}
