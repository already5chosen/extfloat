#include <stdio.h>
#include <stdint.h>
#include <cmath>
#include <iostream>
#include <boost/multiprecision/cpp_bin_float.hpp>

#include "extfloat128.h"
#include "extfloat128_to_bobinf.h"

typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<128, boost::multiprecision::backends::digit_base_2, void, boost::int64_t, INT_MIN, INT_MAX> > boost_quadfloat_t;

static void print(const boost_quadfloat_t& a) {
  printf("%d %11I64d", a.backend().sign(), a.backend().exponent());
  const size_t alen = a.backend().bits().size()*sizeof(a.backend().bits().limbs()[0]);
  if (alen == 16) {
    const uint64_t* bits = reinterpret_cast<const uint64_t*>(a.backend().bits().limbs());
    printf(" %016I64x:%016I64x", bits[1], bits[0]);
  }
  else
    printf(" bits.size = %u*%u", a.backend().bits().size(), unsigned(sizeof(a.backend().bits().limbs()[0])));
  std::cout << " " << a << "\n";
}

static void print(const extfloat128_t& a) {
  printf("%d %11u %016I64x:%016I64x {%d}\n", a.m_sign, a.m_exponent, a.m_significand[1], a.m_significand[0], a._get_exponent());
}

static const unsigned N_VALS = (2*2*2+2)*2 + 1;

static extfloat128_t getTestValue(unsigned i)
{
  if (i < N_VALS-1) {
    // non NaN
    extfloat128_t ret;
    uint32_t sign = i & 1;
    i >>= 1;
    switch (i) {
      case 2*2*2+0: ret = extfloat128_t::zero();
      case 2*2*2+1: ret = extfloat128_t::inf();
      default:
        ret.m_significand[0] = i & 1; i >>= 1;
        ret.m_significand[1] = (i & 1) + (uint64_t(1) << 63); i >>= 1;
        ret.m_exponent       = i + extfloat128_t::exponent_bias;
        break;
    }
    ret.m_sign = sign;
    return ret;
  }
  return extfloat128_t::nan();
}

int main(int argz, char** argv)
{
  printf("operator==");
  for (unsigned xi = 0; xi < N_VALS; ++xi) {
    for (unsigned yi = 0; yi < N_VALS; ++yi) {
      extfloat128_t x = getTestValue(xi);
      extfloat128_t y = getTestValue(yi);
      bool res = (x == y);
      boost_quadfloat_t bx; convert_to_boost_bin_float(&bx, x);
      boost_quadfloat_t by; convert_to_boost_bin_float(&by, y);
      bool ref = (bx == by);
      if (res != ref) {
        printf(" fail. %d != %d\n", res, ref);
        print(x);
        print(bx);
        print(y);
        print(by);
        return 1;
      }
    }
  }
  printf(" o.k.\n");

  printf("operator!=");
  for (unsigned xi = 0; xi < N_VALS; ++xi) {
    for (unsigned yi = 0; yi < N_VALS; ++yi) {
      extfloat128_t x = getTestValue(xi);
      extfloat128_t y = getTestValue(yi);
      bool res = (x != y);
      boost_quadfloat_t bx; convert_to_boost_bin_float(&bx, x);
      boost_quadfloat_t by; convert_to_boost_bin_float(&by, y);
      bool ref = (bx != by);
      if (res != ref) {
        printf(" fail. %d != %d\n", res, ref);
        print(x);
        print(bx);
        print(y);
        print(by);
        return 1;
      }
    }
  }
  printf(" o.k.\n");

  printf("isless()");
  for (unsigned xi = 0; xi < N_VALS; ++xi) {
    for (unsigned yi = 0; yi < N_VALS; ++yi) {
      extfloat128_t x = getTestValue(xi);
      extfloat128_t y = getTestValue(yi);
      bool res = isless(x, y);
      boost_quadfloat_t bx; convert_to_boost_bin_float(&bx, x);
      boost_quadfloat_t by; convert_to_boost_bin_float(&by, y);
      bool ref = (bx < by);
      if (res != ref) {
        printf(" fail. %d != %d\n", res, ref);
        print(x);
        print(bx);
        print(y);
        print(by);
        return 1;
      }
    }
  }
  printf(" o.k.\n");

  printf("islessequal()");
  for (unsigned xi = 0; xi < N_VALS; ++xi) {
    for (unsigned yi = 0; yi < N_VALS; ++yi) {
      extfloat128_t x = getTestValue(xi);
      extfloat128_t y = getTestValue(yi);
      bool res = islessequal(x, y);
      boost_quadfloat_t bx; convert_to_boost_bin_float(&bx, x);
      boost_quadfloat_t by; convert_to_boost_bin_float(&by, y);
      bool ref = (bx <= by);
      if (res != ref) {
        printf(" fail. %d != %d\n", res, ref);
        print(x);
        print(bx);
        print(y);
        print(by);
        return 1;
      }
    }
  }
  printf(" o.k.\n");

  printf("isgreater()");
  for (unsigned xi = 0; xi < N_VALS; ++xi) {
    for (unsigned yi = 0; yi < N_VALS; ++yi) {
      extfloat128_t x = getTestValue(xi);
      extfloat128_t y = getTestValue(yi);
      bool res = isgreater(x, y);
      boost_quadfloat_t bx; convert_to_boost_bin_float(&bx, x);
      boost_quadfloat_t by; convert_to_boost_bin_float(&by, y);
      bool ref = (bx > by);
      if (res != ref) {
        printf(" fail. %d != %d\n", res, ref);
        print(x);
        print(bx);
        print(y);
        print(by);
        return 1;
      }
    }
  }
  printf(" o.k.\n");

  printf("isgreaterequal()");
  for (unsigned xi = 0; xi < N_VALS; ++xi) {
    for (unsigned yi = 0; yi < N_VALS; ++yi) {
      extfloat128_t x = getTestValue(xi);
      extfloat128_t y = getTestValue(yi);
      bool res = isgreaterequal(x, y);
      boost_quadfloat_t bx; convert_to_boost_bin_float(&bx, x);
      boost_quadfloat_t by; convert_to_boost_bin_float(&by, y);
      bool ref = (bx >= by);
      if (res != ref) {
        printf(" fail. %d != %d\n", res, ref);
        print(x);
        print(bx);
        print(y);
        print(by);
        return 1;
      }
    }
  }
  printf(" o.k.\n");

  printf("islessgreater()");
  for (unsigned xi = 0; xi < N_VALS; ++xi) {
    for (unsigned yi = 0; yi < N_VALS; ++yi) {
      extfloat128_t x = getTestValue(xi);
      extfloat128_t y = getTestValue(yi);
      bool res = islessgreater(x, y);
      boost_quadfloat_t bx; convert_to_boost_bin_float(&bx, x);
      boost_quadfloat_t by; convert_to_boost_bin_float(&by, y);
      bool ref = (bx < by) || (bx > by);
      if (res != ref) {
        printf(" fail. %d != %d\n", res, ref);
        print(x);
        print(bx);
        print(y);
        print(by);
        return 1;
      }
    }
  }
  printf(" o.k.\n");

  printf("isunordered()");
  for (unsigned xi = 0; xi < N_VALS; ++xi) {
    for (unsigned yi = 0; yi < N_VALS; ++yi) {
      extfloat128_t x = getTestValue(xi);
      extfloat128_t y = getTestValue(yi);
      bool res = isunordered(x, y);
      boost_quadfloat_t bx; convert_to_boost_bin_float(&bx, x);
      boost_quadfloat_t by; convert_to_boost_bin_float(&by, y);
      bool ref = !(bx <= by) && !(bx > by);
      if (res != ref) {
        printf(" fail. %d != %d\n", res, ref);
        print(x);
        print(bx);
        print(y);
        print(by);
        return 1;
      }
    }
  }
  printf(" o.k.\n");

  return 0;
}
