#define USE_EXTFLOAT128_T_OSTREAM
#include "extfloat128.h"
#include <cmath>
#include <cfloat>
#include <cstdio>

class extfloat128_ostream_init_t {
public:
  extfloat128_t oneE18[2];
  extfloat128_t oneE36[2];
  extfloat128_t oneExx[7*2][2]; // 1e1, 1e2, 1e3,..., 1e7,  1e8, 1e16, 1e24, ..., 1e56
  int64_t       oneEyy[18];     // 1e1, 1e2, 1e3, ... 1e18
  extfloat128_t ln10[2];        // log(10) with > 256 bit precision
  extfloat128_t ln2[2];         // log(2)  with > 256 bit precision
  extfloat128_t expTab[8][7][2];
  extfloat128_t expPoly[5];     // 1/8!, 1/7!, 1/6!, 1/5!, 1/4!
  extfloat128_t::acc_t expPoly3; // 1/3!

  extfloat128_ostream_init_t() {
    extfloat128_t x = 10;
    for (int k = 0; k < 2; ++k) {
      extfloat128_t m = x;
      for (int i = 0; i < 7; ++i) {
        oneExx[k*7+i][0] = m;
        oneExx[k*7+i][1] = inv(m);
        m *= x;
      }
      x = m;
    }
    oneE18[0] = oneExx[7*0+1][0]*oneExx[7*1+1][0];
    oneE36[0] = oneE18[0]*oneE18[0];
    oneE18[1] = inv(oneE18[0]);
    oneE36[1] = inv(oneE36[0]);
    int64_t y = 10;
    for (int i = 0; i < 18; ++i) {
      oneEyy[i] = y;
      y *= 10;
    }

    // log_(ln10, 10);
    // log_(ln2,    2.0);
    // printf("%d %d %016llx %016llx\n", ln10[0].m_sign, ln10[0]._get_exponent(), ln10[0].m_significand[1], ln10[0].m_significand[0]);
    // printf("%d %d %016llx %016llx\n", ln10[1].m_sign, ln10[1]._get_exponent(), ln10[1].m_significand[1], ln10[1].m_significand[0]);
    // printf("%d %d %016llx %016llx\n", ln2[0].m_sign, ln2[0]._get_exponent(), ln2[0].m_significand[1], ln2[0].m_significand[0]);
    // printf("%d %d %016llx %016llx\n", ln2[1].m_sign, ln2[1]._get_exponent(), ln2[1].m_significand[1], ln2[1].m_significand[0]);
    ln10[0] = ln10[0].construct(0, 1,    0x935d8dddaaa8ac16, 0xea56d62b82d30a29);
    ln10[1] = ln10[1].construct(1, -130, 0xeb8098312d10378b, 0xe1cf0beff07e9490);
    ln2[0]  = ln2 [0].construct(0, -1,   0xb17217f7d1cf79ab, 0xc9e3b39803f2f6af);
    ln2[1]  = ln2 [1].construct(0, -130, 0x81e6864ce5316c5b, 0x141a2eb71755f458);

    extfloat128_t::acc_t ex0 = extfloat128_ostream_init_t::exp(extfloat128_t::pow2(-25));
    extfloat128_t::acc_t ex[7];
    for (int k = 7; k >= 0; --k) {
      ex[1-1] = mult(ex0, ex0);
      ex[2-1] = mult(ex[1-1], ex[1-1]);
      ex[3-1] = mult(ex[2-1], ex[1-1]);
      ex[4-1] = mult(ex[2-1], ex[2-1]);
      ex[5-1] = mult(ex[4-1], ex[1-1]);
      ex[6-1] = mult(ex[4-1], ex[2-1]);
      ex[7-1] = mult(ex[4-1], ex[3-1]);
      ex0 = ex[4-1];
      for (int i = 0; i < 7; ++i) {
        expTab[k][i][0] = ex[i].round();
        ex[i] -= expTab[k][i][0];
        ex[i] /= expTab[k][i][0];  // ex[i]/expTab[k][i][0] - 1, calculated in more precise way
        expTab[k][i][1] = ex[i].round();
      }
    }

    int fact = 24;
    for (int k = 0; k < 5; ++k) {
      expPoly[4-k] = extfloat128_t::one() / fact;
      fact *= k + 5;
    }
    expPoly3 = extfloat128_t::one();
    expPoly3 /= 6;
  }
  extfloat128_t inv(const extfloat128_t& x) {
    extfloat128_t y = extfloat128_t::one() / x;
    if (fma(x, y, extfloat128_t::one(1)) > extfloat128_t::zero())
      y = nextafter(y, 0);
    return y;
  }

  static void split(extfloat128_t dst[3], extfloat128_t::acc_t x)
  {
    dst[2] = x.trunc();
    x -= dst[2];
    dst[1] = x.trunc();
    x -= dst[1];
    dst[0] = x.trunc();
  }

  static extfloat128_t::acc_t mult(const extfloat128_t::acc_t& a, const extfloat128_t::acc_t& b)
  {
    extfloat128_t x[2][3];
    split(x[0], a);
    split(x[1], b);
    extfloat128_t::acc_t prod;
    prod.clear();
    for (int i = 0; i < 3; ++i) {
      for (int k = 0; k < 3; ++k) {
        prod.madd(x[0][i], x[1][k]);
      }
    }
    return prod;
  }

  // x must be small
  static extfloat128_t::acc_t exp(const extfloat128_t& x)
  {
    extfloat128_t::acc_t res;
    res = 1;
    for (int k = 12; k > 0; --k) {
      res *= x;
      res /= k;
      res += 1;
    }
    return res;
  }

  // void log_(extfloat128_t dst[2], double x)
  // {
    // boost_float256_t b = log(boost_float256_t(x));
    // convert_from_boost_bin_float(&dst[0], b);
    // boost_float256_t b1;
    // convert_to_boost_bin_float(&b1, dst[0]);
    // b -= b1;
    // convert_from_boost_bin_float(&dst[1], b);
  // }
};
static extfloat128_ostream_init_t tab;

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
      remdiv[1] += t.convert_to_uint64();
      x = fma(tab.oneE18[0], -t, x);
    } else {
      remdiv[1] += 1;
      x -= tab.oneE18[0];
    }
  }
  remdiv[0] = x.convert_to_uint64();
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
  return y.convert_to_uint64();
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

// src in range [0..log(2)+eps]
// return exp(src) with precision of 216 bits or better (65 decimal digits)
static extfloat128_t::acc_t exp216(extfloat128_t::acc_t& src)
{
  extfloat128_t x = src.trunc();
  src -= x;
  extfloat128_t adj = src.trunc();

  extfloat128_t xh = trunc(extfloat128_t::ldexp(x, 24));
  unsigned msbits = xh.convert_to_uint64();
  if (msbits > 0)
    x -= extfloat128_t::ldexp(xh, -24);

  extfloat128_t xy = tab.expPoly[0];
  for (int i = 1; i < 5; ++i)
    xy = fma(xy, x, tab.expPoly[i]);

  extfloat128_t::acc_t res = tab.expPoly3;
  res.madd(xy, x);
  res *= x;
  res += extfloat128_t::pow2(-1); // += 0.5
  res *= x;
  res += extfloat128_t::one();
  res *= x;
  res += extfloat128_t::one();

  if (msbits > 0) {
    for (int k = 0; k < 8; ++k) {
      unsigned idx = msbits & 7;
      if (idx != 0) {
        res *= tab.expTab[7-k][idx-1][0];
        adj += tab.expTab[7-k][idx-1][1];
      }
      msbits >>= 3;
    }
  }
  // res *= (1+adj) where adj is very small
  res.madd(adj, res.round());
  return res;
}

static int scalebyPow10(extfloat128_t::acc_t* dst, const extfloat128_t& x, int32_t scale10)
{
  static const double LOG2_10 = 3.3219280948873623478703194294894;
  int64_t scale2 = int64_t(floor(scale10*LOG2_10-1e-5));
  extfloat128_t::acc_t ln;
  ln.mulx(tab.ln10[0], scale10);
  ln.madd(tab.ln10[1], scale10);
  ln.msub(tab.ln2[0],  scale2);
  ln.msub(tab.ln2[1],  scale2);
  *dst = exp216(ln);
  *dst *= x;
  dst->m_exponent += scale2;
  dst->normalize();
  dst->m_sign = 0;
  // if (dbg) {
    // extfloat128_t xx = dst->trunc();
    // extfloat128_t xi = trunc(xx);
    // extfloat128_t xr = xx - xi;
    // printf("%18e %.15f scale10=%d\n", xi.convert_to_double(), xr.convert_to_double(), scale10);
  // }
  return dst->round_to_nearest_tie_to_even();
}

// scale10 in range [-38..-1]
// return value: -1 if rounded toward 0, 1 if rounded away from zero, 0 if unchanged
static int scalebyPow10_small(extfloat128_t::acc_t* dst, extfloat128_t x, int scale10)
{
  x.m_sign = 0;
  // if (dbg) {
  //  std::cout << "scalebyPow10_small("; pp(x); std::cout << "   , " << std::dec << scale10 << ");\n";
  // }
  if (scale10 <= 0) {
    // divide
    *dst = x;
    int roundDir = dst->round_to_nearest_tie_to_even();
    unsigned scale10n = -scale10;
    int scaleIx = scale10n / 8;
    for (int k = 7; k >= 0; k -= 7) {
      if (scaleIx != 0) {
        int ix = k + scaleIx - 1; // index in oneExx
        // divide
        extfloat128_t::acc_t div; div.clear();
        extfloat128_t rem = dst->trunc();
        // int uu = 0;
        while (rem >= tab.oneExx[ix][0]) {
          extfloat128_t inc = trunc(multiply_rtz(rem, tab.oneExx[ix][1]));
          if (!is_zero(inc)) {
            dst->msub(inc, tab.oneExx[ix][0]);
            div += inc;
          } else {
            *dst -= tab.oneExx[ix][0];
            div  += extfloat128_t::one();
          }
          rem = dst->trunc();
          // ++uu;
          // if (dbg) {
            // std::cout << "uu=" << uu << "\n";
            // std::cout << "div "; pp(div.trunc()); std::cout << "\n";
            // std::cout << "rem "; pp(rem); std::cout << "\n";
          // }
          // if (uu == 10)
            // break;
        }
        dst->normalize();
        if (!dst->is_zero()) {
          *dst -= extfloat128_t::ldexp(tab.oneExx[ix][0], -1);
          dst->normalize();
          if (!dst->is_zero()) {
            roundDir = dst->m_sign ? -1 : 1;
          } else if (roundDir == 0) {
            // a tie
            roundDir = div.is_odd() ? 1 : -1; // round to even
          } else {
            roundDir = -roundDir;
          }
          if (roundDir > 0)
            div += extfloat128_t::one();
        }
        *dst = div;
      }
      scaleIx = scale10n % 8;
    }
    return roundDir;
  } else {
    // multiply
    unsigned scale10p = scale10;
    int scaleIL = scale10p % 8;
    int scaleIH = scale10p / 8;
    if (scaleIL != 0) {
      int ix = scaleIL - 1;   // index in oneExx
      dst->mulx(x, tab.oneExx[ix][0]);
      if (scaleIH != 0) {
        ix = 7 + scaleIH - 1; // index in oneExx
        extfloat128_t a = dst->trunc();
        *dst -= a;
        extfloat128_t b = dst->trunc();
        dst->mulx(a, tab.oneExx[ix][0]);
        dst->madd(b, tab.oneExx[ix][0]);
      }
    } else {
      int ix = 7 + scaleIH - 1; // index in oneExx
      dst->mulx(x, tab.oneExx[ix][0]);
    }
    return dst->round_to_nearest_tie_to_even();
  }
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
    int roundDir = (scale10 > -56 && scale10 < 56) ?
      scalebyPow10_small(&acc, x, scale10) :
      scalebyPow10      (&acc, x, scale10);
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
      digE18[nDigE36*2] = digE36[nDigE36].convert_to_uint64();
    // if (dbg) printf("digits %d. nDigE18 %d: %lld %lld %lld %lld scale10=%d exp10=%d\n", digits, nDigE18, digE18[0], digE18[1], digE18[2], digE18[3], scale10, exp10);

    int msDigits = digits + 1 - (nDigE18-1)*18;
    bool extraDigit = (digE18[nDigE18-1] >= tab.oneEyy[msDigits-1]); // 1 MS word contains one more decimal digit than requested
    if (extraDigit) {
      exp10 += 1;
      digE18[0] = divBy10andRound(digE18[0], roundDir);
      if (nDigE18 > 1) {
        if (msDigits == 18) {
          // keep digE18[] values in range [0..1E18-1]
          digE18[nDigE18-1] -= tab.oneEyy[18-1];
          digE18[nDigE18] = 1;
          nDigE18 += 1;
        }
        int64_t maxVal = tab.oneEyy[17-1]; // 1e17
        for (int i = 0; digE18[i] >= maxVal; ++i) {
          digE18[i]   -= maxVal;
          digE18[i+1] += 1;          // propagate carry
          maxVal = tab.oneEyy[18-1]; // for all, but the least significant digit, the limit = 1e18
        }
      }
    }
    // if (dbg) printf("Digits %d. nDigE18 %d: %lld %lld %lld %lld scale10=%d exp10=%d\n", digits, nDigE18, digE18[0], digE18[1], digE18[2], digE18[3], scale10, exp10);

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
