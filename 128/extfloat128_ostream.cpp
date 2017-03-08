#define USE_EXTFLOAT128_T_OSTREAM
#include "extfloat128.h"
#include <cmath>
#include <cfloat>
#include <cstdio>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include "extfloat128_to_bobinf.h"
#include "bobinf_to_extfloat128.h"

typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<192+64, boost::multiprecision::backends::digit_base_2> > boost_float192_t;

static void pow5(extfloat128_t dst[2], int e)
{
  boost_float192_t b = pow(boost_float192_t(5), e);
  convert_from_boost_bin_float(&dst[0], b);
  boost_float192_t b1;
  convert_to_boost_bin_float(&b1, dst[0]);
  b -= b1;
  convert_from_boost_bin_float(&dst[1], b);
}

static int64_t toDecimalStep(extfloat128_t::acc_t& acc, const extfloat128_t& scale, const extfloat128_t& invScale)
{
  extfloat128_t::acc_t tmp;
  tmp.clear();
  extfloat128_t y = acc.trunc();
  while (y >= scale) {
    y = trunc(y * invScale);
    if (is_zero(y))
      y = extfloat128_t::one();
    tmp += y;
    acc.msub(y, scale);
    y = acc.trunc();
  }
  acc = tmp;
  return y.convert_to_int64();
}

#if 0
static int32_t toDecimalStep(extfloat128_t::acc_t& acc, const extfloat128_t scale)
{
  extfloat128_t y = acc.trunc();
  acc -= y;
  extfloat128_t::acc_t tmp;
  tmp.mulx(y, scale);
  extfloat128_t yt = trunc(tmp.trunc());
  tmp -= yt;
  double dt = yt.convert_to_double();
  y = acc.round();
  tmp.madd(y, scale);
  yt = trunc(tmp.trunc());
  tmp -= yt;
  dt += yt.convert_to_double();
  acc = tmp;
  return uint32_t((int64_t(dt)));
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

static const extfloat128_t oneE18    = extfloat128_t(int64_t(1000000)*int64_t(1000000)*int64_t(1000000));
static const extfloat128_t invOneE18 = extfloat128_t::construct(0, -1, uint64_t(-1), uint64_t(-1))/oneE18;
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
// printf("invOneE18=%g\n", invOneE18.convert_to_double());
    int32_t exp10 = log10_estimate(x);
    int32_t scale10 = digits - exp10; // sometimes we'll have one more digit than requested
    extfloat128_t mult[2];
    pow5(mult, scale10);
    extfloat128_t::acc_t acc;
    acc.mulx(mult[0], x);
    acc.madd(mult[1], x);
    acc.normalize();
    acc.m_exponent += scale10;
    acc.m_sign = 0;
    acc.round_to_nearest_tie_to_even();
    int64_t digE18[max_digits/18+2];
    int nIter = (digits+1)/18;
    for (int i = 0; i < nIter; ++i) {
      digE18[i] = toDecimalStep(acc, oneE18, invOneE18);
    }
    if (nIter*18 < digits+1) {
      digE18[nIter] = acc.trunc().convert_to_int64();
      ++nIter;
    }
    // printf("[%d %d %d] %lld %lld %lld\n", exp10, scale10, nIter, digE18[0], digE18[1], digE18[2]);

    char* p = &buf[2];
    p += sprintf(p, "%lld", digE18[nIter-1]);
    if (nIter > 1) {
      for (int i = 2; i < nIter; ++i)
        p += sprintf(p, "%018lld", digE18[nIter-i]);
      int bufDigits = p - &buf[2];
      int remDigits = digits + 1 - bufDigits;
      if (remDigits == 18) {
        p += sprintf(p, "%018lld", digE18[0]);
      } else {
        p += sprintf(p, "%017lld", (digE18[0]+5)/10);
        exp10 += 1;
      }
    } else {
      int bufDigits = p - &buf[2];
      if (bufDigits > digits + 1) {
        p += sprintf(p, "%lld", (digE18[0]+5)/10);
        exp10 += 1;
      }
    }
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
