#pragma once

#include <stdint.h>
#include <limits.h>
#include <string.h>
#ifdef USE_EXTFLOAT128_T_OSTREAM
#include <ostream>
#endif

class extfloat128_t {
public:
  extfloat128_t (){}
  extfloat128_t (int a) {
    from_double(this, a);
  }
  extfloat128_t (double a) {
    from_double(this, a);
  }
  extfloat128_t (uint64_t a) {
    from_uint64(this, a);
  }
  extfloat128_t (int64_t a) {
    from_int64(this, a);
  }

  #ifdef USE_EXTFLOAT128_T_OSTREAM
  friend std::ostream& operator<<(std::ostream& os, const extfloat128_t& x);
  #endif

  //extfloat128_t operator=(double a);
  //extfloat128_t operator=(int a);

  extfloat128_t operator+=(const extfloat128_t& a) {
    extfloat128_t::eval_addsub(*this, *this, a, 0);
    return *this;
  }
  extfloat128_t operator-=(const extfloat128_t& a) {
    extfloat128_t::eval_addsub(*this, *this, a, 1);
    return *this;
  }
  extfloat128_t operator*=(const extfloat128_t& a) {
    extfloat128_t::eval_multiply(*this, *this, a);
    return *this;
  }
  extfloat128_t operator/=(const extfloat128_t& a) {
    extfloat128_t::eval_divide(*this, *this, a);
    return *this;
  }
  extfloat128_t operator-() const {
    extfloat128_t ret = *this;
    ret.m_sign = m_sign ^ 1;
    return ret;
  }

  double convert_to_double() const;
  float  convert_to_float() const;
  int64_t convert_to_int64() const;
  extfloat128_t sqrt() const {
    extfloat128_t res;
    eval_sqrt(res, *this);
    return res;
  }
  extfloat128_t rsqrt() const {
    extfloat128_t res;
    eval_rsqrt(res, *this);
    return res;
  }
  extfloat128_t ulp() const;

  static bool eval_eq      (const extfloat128_t& srcA, const extfloat128_t& srcB);
  static bool eval_lt      (const extfloat128_t& srcA, const extfloat128_t& srcB);
  static bool eval_gt      (const extfloat128_t& srcA, const extfloat128_t& srcB) { return eval_lt(srcB, srcA); }
  static bool eval_le      (const extfloat128_t& srcA, const extfloat128_t& srcB);
  static bool eval_ge      (const extfloat128_t& srcA, const extfloat128_t& srcB) { return eval_le(srcB, srcA); }
  static bool eval_lg      (const extfloat128_t& srcA, const extfloat128_t& srcB);
  static bool eval_uo      (const extfloat128_t& srcA, const extfloat128_t& srcB);
  static void eval_addsub  (extfloat128_t& dst, const extfloat128_t& srcA, const extfloat128_t& srcB, int issub);
  static void eval_multiply(extfloat128_t& dst, const extfloat128_t& srcA, const extfloat128_t& srcB);
  static void eval_divide  (extfloat128_t& dst, const extfloat128_t& srcA, const extfloat128_t& srcB);
  static void eval_add     (extfloat128_t& dst, const extfloat128_t& srcA, const extfloat128_t& srcB);
  static void eval_subtract(extfloat128_t& dst, const extfloat128_t& srcA, const extfloat128_t& srcB);
  static void eval_fma     (extfloat128_t& dst, const extfloat128_t& srcA, const extfloat128_t& srcB, const extfloat128_t& srcZ);

  static void from_double(extfloat128_t* dst, double src);
  static void from_uint64(extfloat128_t* dst, uint64_t src);
  static void from_int64(extfloat128_t* dst, int64_t src) {
    if (src < 0) {
      from_uint64(dst, -src);
      dst->m_sign = 1;
    } else {
      from_uint64(dst, src);
    }
  }
  static extfloat128_t ldexp(const extfloat128_t& x, int exp);

  static extfloat128_t zero(uint32_t sign = 0) {
    extfloat128_t ret;
    ret.m_significand[0] = ret.m_significand[1] = 0;
    ret.m_exponent = zero_biased_exponent;
    ret.m_sign = sign;
    return ret;
  }
  static extfloat128_t one(uint32_t sign = 0) { return pow2(0, sign); }
  static extfloat128_t inf() {
    extfloat128_t ret;
    ret.m_significand[0] = 0;
    ret.m_significand[1] = 0;
    ret.m_exponent       = inf_nan_biased_exponent;
    ret.m_sign           = 0;
    return ret;
  }
  static extfloat128_t nan(const char* tagp=0) {
    extfloat128_t ret;
    ret.m_significand[0] = 0;
    ret.m_significand[1] = qnan_bit;
    ret.m_exponent      = inf_nan_biased_exponent;
    ret.m_sign           = 0;
    if (tagp)
      strncpy(reinterpret_cast<char*>(ret.m_significand), tagp, sizeof(ret.m_significand[0]));
    return ret;
  }

  friend extfloat128_t sqrt(const extfloat128_t& a);
  friend extfloat128_t sinPi(const extfloat128_t& a) { return cossinPi_core(a, 1); }; // sin(a*pi)
  friend extfloat128_t cosPi(const extfloat128_t& a) { return cossinPi_core(a, 0); }; // cos(a*pi)
  friend extfloat128_t atanPi(const extfloat128_t& a);                          // atan(a)/pi
  friend extfloat128_t atan2Pi(const extfloat128_t& a, const extfloat128_t& b); // atan(a,b)/pi
  friend extfloat128_t asinPi(const extfloat128_t& a);                          // asin(a)/pi
  friend extfloat128_t acosPi(const extfloat128_t& a);                          // acos(a)/pi

  friend extfloat128_t sin(const extfloat128_t& a) { return cossin_core(a, 1); }; // sin(a) - not optimized as well as sinPi()
  friend extfloat128_t cos(const extfloat128_t& a) { return cossin_core(a, 0); }; // cos(a) - not optimized as well as cosPi()

  static extfloat128_t construct(uint32_t sign, int32_t exp, uint64_t significand_hi,  uint64_t significand_lo) {
    extfloat128_t ret;
    ret.m_significand[0] = significand_lo;
    ret.m_significand[1] = significand_hi;
    ret.m_exponent       = exp + exponent_bias;
    ret.m_sign           = sign;
    return ret;
  }
  static extfloat128_t pow2(int32_t exp, uint32_t sign = 0) { return construct(sign, exp, uint64_t(1) << 63, 0); }
  static extfloat128_t pi()       { return construct(0,    1, 0xc90fdaa22168c234, 0xc4c6628b80dc1cd1); }
  static extfloat128_t pi_lo()    { return construct(0, -129, 0xa4093822299f31d0, 0x082efa98ec4e6c89); }
  static extfloat128_t pi_lo2()   { return construct(0, -258, 0x8a5043cc71a026ef, 0x7ca8cd9e69d218da); }
  static extfloat128_t pi_lo3()   { return construct(0, -387, 0xfd4f5920da0ebc8b, 0x01eca9292ae3dba2); }
  static extfloat128_t invPi()    { return construct(0,   -2, 0xa2f9836e4e441529, 0xfc2757d1f534ddc1); }
  static extfloat128_t invPi_lo() { return construct(1, -132, 0x9275a99b0ef1bef8, 0x06ba71508510ea79); }
  static extfloat128_t invPi_lo2(){ return construct(1, -262, 0x8db91c5bdb22d1ff,0x9b6d115f62e6de30); }
  static extfloat128_t invPi_lo3(){ return construct(1, -393, 0xf10a71a76b2c608b,0xbee50568a25dbd8b); }
  friend extfloat128_t mod_pow2(const extfloat128_t& a, int32_t exp); // return sign(a)*mod(abs(a), 2**exp)

  friend extfloat128_t ceil (const extfloat128_t& x); // Rounds x upward, returning the smallest integral value that is not less than x
  friend extfloat128_t floor(const extfloat128_t& x); // Rounds x downward, returning the largest integral value that is not greater than x
  friend extfloat128_t round(const extfloat128_t& x); // Returns the integral value that is nearest to x, with halfway cases rounded away from zero
  friend extfloat128_t trunc(const extfloat128_t& x); // Rounds x toward zero, returning the nearest integral value that is not larger in magnitude than x
  static int           compare_ex(const extfloat128_t& x, const extfloat128_t& y);
                                                      // return -1 for unordered, 0 for x == y, 1 for x < y, 2 for x > y

  friend extfloat128_t nextafter(const extfloat128_t& x, const extfloat128_t& y);
  friend extfloat128_t nextupdown(const extfloat128_t& x, bool down); // when x==NAN return x
                                                                      // when x == +inf and down = false return x
                                                                      // when x == -inf and down = true  return x
                                                                      // otherwise
                                                                      // when down = false return next representable value in the direction of +inf
                                                                      // when down = true  return next representable value in the direction of -inf
  friend extfloat128_t nextinout(const extfloat128_t& x, bool out);   // when x==NAN return x
                                                                      // when x == 0 and out = false return x
                                                                      // when (x == +inf or x == -inf) and out = true return x
                                                                      // otherwise
                                                                      // when out = false return next representable value in the direction of 0
                                                                      // when out = true  return next representable value in the  direction away of zero


  static void op_int       (extfloat128_t* dst, const extfloat128_t* srcA, int srcB,    void (*func)(extfloat128_t& dst, const extfloat128_t& srcA, const extfloat128_t& srcB));
  static void op_double    (extfloat128_t* dst, const extfloat128_t* srcA, double srcB, void (*func)(extfloat128_t& dst, const extfloat128_t& srcA, const extfloat128_t& srcB));
  static void op_double_rev(extfloat128_t* dst, double srcA, const extfloat128_t* srcB, void (*func)(extfloat128_t& dst, const extfloat128_t& srcA, const extfloat128_t& srcB));

  static const uint32_t zero_biased_exponent    = 0;
  static const uint32_t min_biased_exponent     = 1;
  static const uint32_t max_biased_exponent     = uint32_t(-1) - 1;
  static const uint32_t inf_nan_biased_exponent = uint32_t(-1);
  static const uint32_t exponent_bias           = (uint32_t(1)<<31)-1;
  static const uint64_t qnan_bit                = (uint64_t(1)<<63);

  static const int32_t min_exponent_val = min_biased_exponent - exponent_bias;
  static const int32_t max_exponent_val = max_biased_exponent - exponent_bias;

  int32_t _get_exponent() const     { return m_exponent - exponent_bias; }
  void    _set_exponent(int32_t e)  { m_exponent = e + exponent_bias;    }

  void to_17bytes(uint8_t* dst) const;
  void from_17bytes(const uint8_t* src);
  static void from_17bytes_fix(extfloat128_t dst[2], const uint8_t* src);
  static void to_18bytes_fp  (uint8_t* dst, const extfloat128_t src[2]);
  static void from_18bytes_fp(extfloat128_t dst[2], const uint8_t* src);

//private:
  uint64_t m_significand[2];
  uint32_t m_exponent;
  uint32_t m_sign;
  // static void   boost_addsub(extfloat128_t* dst, const extfloat128_t* srcA, const extfloat128_t* srcB, int issub);
  // static void   boost_mul(extfloat128_t* dst,    const extfloat128_t* srcA, const extfloat128_t* srcB);
  // static void   to_boost_quadfloat(void* dst,  const extfloat128_t* src);
private:
  // static void   boost_from_double(extfloat128_t* dst, double src);
  // static double boost_convert_to_double(const extfloat128_t* src);
  static void   addsub_core(extfloat128_t* dst, const extfloat128_t* srcA, const extfloat128_t* srcB, int invA, int invB);
  static int    compare_abs_core(const extfloat128_t* srcA, const extfloat128_t* srcB);
  static void   eval_sqrt(extfloat128_t& dst, const extfloat128_t& src);
  static void   eval_rsqrt(extfloat128_t& dst, const extfloat128_t& src);
  static void   eval_non_finite(extfloat128_t& dst, const extfloat128_t& srcA, const extfloat128_t& srcB, int issub, const uint8_t tab[2][2][3]);
  static extfloat128_t cossinPi_core(const extfloat128_t& a, int isSin); // isSin ? sin(a*pi) : cos(a*pi)
  static extfloat128_t cossin_core  (const extfloat128_t& a, int isSin); // isSin ? sin(a)    : cos(a)

  static extfloat128_t cossinPi_core_body(extfloat128_t& x, int isSin, extfloat128_t*); // isSin ? sin(a*pi) : cos(a*pi)

public:
  uint64_t scale_and_trunc(int e) const // return trunc(abs(*this)*2^e)
  { // don't check for Inf/Nan, on overflow return 0
    uint64_t rshift = int64_t(exponent_bias + 63) - e - m_exponent;
    return rshift < 64 ? m_significand[1] >> rshift : 0;
  }

private:
  extfloat128_t nextinout_core(uint32_t out) const;
  // *this should be non-nan
  // when *this == 0 and out = false return x
  // when (*this == +inf or *this == -inf) and out = true return x
  // otherwise
  // when out = false return next representable value in the direction of 0
  // when out = true  return next representable value in the  direction away of zero

  extfloat128_t nextupdown_core(uint32_t down) const;
  // *this should be non-nan
  // when *this == +inf and down = false return x
  // when *this == -inf and down = true  return x
  // otherwise
  // when down = false return next representable value in the direction of +inf
  // when down = true  return next representable value in the direction of -inf

public:
  // extended accumulator type - mostly for trigs and other special functions
  // The type intended for intermediate values, mostly accumulators
  // It has very wide exponent range and 2-3 times wider mantissa (272 to 380 bits) than extfloat128_t
  // It does not support NaN.
  // It does not support overflow except for export of extfloat128_t
  // It does not support signed zeros except for export of extfloat128_t
  // Overflow on arithmetic operations is mostly not checked
  struct acc_t {
    // initialize and import
    void clear() {
      m_exponent = 0;
    }
    void operator=(const extfloat128_t& a) {
      m_significand[0] = m_significand[1] = m_significand[2] = 0;
      m_significand[3] = a.m_significand[0];
      m_significand[4] = a.m_significand[1];
      m_significand[5] = 0;
      m_sign     = a.m_sign;
      m_exponent = (a.m_exponent == a.zero_biased_exponent) ?
        0 : (exponent_bias-a.exponent_bias)+a.m_exponent;
      m_isNormal = true;
    }
    void normalize() {
      if (!m_isNormal)
        doNormalize();
    }
    void from_17bytes_fix(const uint8_t* src);

    // floating point arithmetic
    void mulx  (const extfloat128_t& a, const extfloat128_t& b);
    void addsub(const extfloat128_t& a, bool sub);
    void operator+=(const extfloat128_t& a) {
      addsub(a, false);
    }
    void operator-=(const extfloat128_t& a) {
      addsub(a, true);
    }
    void maddsub(const extfloat128_t& a, const extfloat128_t& b, bool sub);
    void madd(const extfloat128_t& a, const extfloat128_t& b) {
      maddsub(a, b, false);
    }
    void msub(const extfloat128_t& a, const extfloat128_t& b) {
      maddsub(a, b, true);
    }
    void negate() {
      m_sign ^= 1;
    }
    void round_to_nearest_tie_to_even();
    // fix-point arithmetic
    void fix_mulx5(const extfloat128_t& a, const uint32_t b[7*2], int b_exponent);
    void fix_mulx3(const extfloat128_t& a, const uint32_t b[5*2], int b_exponent);

    // export
    extfloat128_t to_extfloat128_t(bool bRound) ; // if (bRound) round to nearest otherwise truncate
    extfloat128_t round() { return to_extfloat128_t(true);  }
    extfloat128_t trunc() { return to_extfloat128_t(false); }

    // representation
    static const uint64_t exponent_bias = (uint64_t(1)<<63)-1;
    static const int NBITS_MIN = 256 + 24;
    static const int NBITS_MAX = 320 + 60;
    uint64_t m_significand[6]; // LS word first. In normal form m_significand[5]==0, bit 63 of m_significand[4]==1
                               // m_significand[5] holds signed (int64_t) value
    uint64_t m_exponent;
    uint8_t  m_sign;
    bool     m_isNormal; // for zero this member undefined

    private:
    void doNormalizeSign();
    void doNormalize();
    void PartialNormalize();
  };
};

inline bool operator==(const extfloat128_t& l, const extfloat128_t& r) {
  return extfloat128_t::eval_eq(l, r);
}
inline bool operator!=(const extfloat128_t& l, const extfloat128_t& r) {
  return !extfloat128_t::eval_eq(l, r);
}
inline extfloat128_t operator+(const extfloat128_t& l, const extfloat128_t& r) {
  extfloat128_t ret;
  extfloat128_t::eval_addsub(ret, l, r, 0);
  return ret;
}
inline extfloat128_t operator-(const extfloat128_t& l, const extfloat128_t& r) {
  extfloat128_t ret;
  extfloat128_t::eval_addsub(ret, l, r, 1);
  return ret;
}
inline extfloat128_t operator*(const extfloat128_t& l, const extfloat128_t& r) {
  extfloat128_t ret;
  extfloat128_t::eval_multiply(ret, l, r);
  return ret;
}
inline extfloat128_t operator/(const extfloat128_t& l, const extfloat128_t& r) {
  extfloat128_t ret;
  extfloat128_t::eval_divide(ret, l, r);
  return ret;
}

inline extfloat128_t operator+(const extfloat128_t& l, int r) {
  extfloat128_t ret;
  extfloat128_t::op_int(&ret, &l, r, extfloat128_t::eval_add);
  return ret;
}
inline extfloat128_t operator-(const extfloat128_t& l, int r) {
  extfloat128_t ret;
  extfloat128_t::op_int(&ret, &l, r, extfloat128_t::eval_subtract);
  return ret;
}
inline extfloat128_t operator*(const extfloat128_t& l, int r) {
  extfloat128_t ret;
  extfloat128_t::op_int(&ret, &l, r, extfloat128_t::eval_multiply);
  return ret;
}
inline extfloat128_t operator/(const extfloat128_t& l, int r) {
  extfloat128_t ret;
  extfloat128_t::op_int(&ret, &l, r, extfloat128_t::eval_divide);
  return ret;
}

inline extfloat128_t operator+(const extfloat128_t& l, double r) {
  extfloat128_t ret;
  extfloat128_t::op_double(&ret, &l, r, extfloat128_t::eval_add);
  return ret;
}
inline extfloat128_t operator-(const extfloat128_t& l, double r) {
  extfloat128_t ret;
  extfloat128_t::op_double(&ret, &l, r, extfloat128_t::eval_subtract);
  return ret;
}
inline extfloat128_t operator*(const extfloat128_t& l, double r) {
  extfloat128_t ret;
  extfloat128_t::op_double(&ret, &l, r, extfloat128_t::eval_multiply);
  return ret;
}
inline extfloat128_t operator/(const extfloat128_t& l, double r) {
  extfloat128_t ret;
  extfloat128_t::op_double(&ret, &l, r, extfloat128_t::eval_divide);
  return ret;
}

inline extfloat128_t operator+(double l, const extfloat128_t& r) {
  extfloat128_t ret;
  extfloat128_t::op_double_rev(&ret, l, &r, extfloat128_t::eval_add);
  return ret;
}
inline extfloat128_t operator-(double l, const extfloat128_t& r) {
  extfloat128_t ret;
  extfloat128_t::op_double_rev(&ret, l, &r, extfloat128_t::eval_subtract);
  return ret;
}
inline extfloat128_t operator*(double l, const extfloat128_t& r) {
  extfloat128_t ret;
  extfloat128_t::op_double_rev(&ret, l, &r, extfloat128_t::eval_multiply);
  return ret;
}
inline extfloat128_t operator/(double l, const extfloat128_t& r) {
  extfloat128_t ret;
  extfloat128_t::op_double_rev(&ret, l, &r, extfloat128_t::eval_divide);
  return ret;
}

inline extfloat128_t abs(const extfloat128_t& a)
{
  extfloat128_t ret = a;
  ret.m_sign = 0;
  return ret;
}

inline bool isfinite(const extfloat128_t& a)
{
  return a.m_exponent != a.inf_nan_biased_exponent;
}

inline bool isnormal(const extfloat128_t& a)
{
  return uint32_t(a.m_exponent-a.min_biased_exponent) <= (a.max_biased_exponent-a.min_biased_exponent);
}

inline bool isnan(const extfloat128_t& a)
{
  return
    (a.m_exponent == a.inf_nan_biased_exponent) &&
    ((a.m_significand[0] | a.m_significand[1]) != 0);
}

inline bool isinf(const extfloat128_t& a)
{
  return
    (a.m_exponent == a.inf_nan_biased_exponent) &&
    ((a.m_significand[0] | a.m_significand[1]) == 0);
}

inline bool is_zero(const extfloat128_t& a)
{
  return a.m_exponent == a.zero_biased_exponent;
}

inline extfloat128_t fma(const extfloat128_t& x, const extfloat128_t& y, const extfloat128_t& z)
{
  extfloat128_t res;
  extfloat128_t::eval_fma(res, x, y, z);
  return res;
}

inline bool isgreater      (const extfloat128_t&  x, const extfloat128_t&  y) { return extfloat128_t::eval_gt(x, y); }
inline bool isgreaterequal (const extfloat128_t&  x, const extfloat128_t&  y) { return extfloat128_t::eval_ge(x, y); }
inline bool isless         (const extfloat128_t&  x, const extfloat128_t&  y) { return extfloat128_t::eval_lt(x, y); }
inline bool islessequal    (const extfloat128_t&  x, const extfloat128_t&  y) { return extfloat128_t::eval_le(x, y); }
inline bool islessgreater  (const extfloat128_t&  x, const extfloat128_t&  y) { return extfloat128_t::eval_lg(x, y); }
inline bool isunordered    (const extfloat128_t&  x, const extfloat128_t&  y) { return extfloat128_t::eval_uo(x, y); }

inline bool operator>      (const extfloat128_t&  x, const extfloat128_t&  y) { return extfloat128_t::eval_gt(x, y); }
inline bool operator>=     (const extfloat128_t&  x, const extfloat128_t&  y) { return extfloat128_t::eval_ge(x, y); }
inline bool operator<      (const extfloat128_t&  x, const extfloat128_t&  y) { return extfloat128_t::eval_lt(x, y); }
inline bool operator<=     (const extfloat128_t&  x, const extfloat128_t&  y) { return extfloat128_t::eval_le(x, y); }
