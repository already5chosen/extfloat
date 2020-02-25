#include <iostream>
#include <random>
#include <functional>
#include "sqrt448.h"


static void find_tip_points(const cpp_bin_float_132& src, cpp_bin_float_132 dst[2]);

int main()
{
  for (int i=0; i < 10000;i++) {
    for (int wi = 0; wi < 7; ++wi) {
      cpp_bin_float_132 x(1);
      reinterpret_cast<uint64_t*>(&x.backend().bits())[wi] += i;
      for (int e = 0; e < 2; ++e) {
        cpp_bin_float_132 res = my_sqrt(x);
        cpp_bin_float_132 ref = sqrt(x);
        if (ref != res) {
          std::cout
            << "i=" << i << ":" << wi << ":" << e << std::endl
            << std::setprecision(150)
            << "x=" << x << std::endl
            << "res=" << res << std::endl
            << "ref=" << ref << std::endl;
          return 1;
        }
        x += x;
      }
    }
  }
  std::cout << "[0] o.k\n";

  for (int i=0; i < 100000;i++) {
    cpp_bin_float_132 x = i;
    cpp_bin_float_132 res = my_sqrt(x);
    cpp_bin_float_132 ref = sqrt(x);
    if (ref != res) {
      std::cout
        << "i=" << i << std::endl
        << std::setprecision(150)
        << "x=" << x << std::endl
        << "res=" << res << std::endl
        << "ref=" << ref << std::endl;
      return 1;
    }
  }
  std::cout << "[1] o.k\n";

  cpp_bin_float_132 y = 1;
  for (int i=0; i < 10000;i++) {
    for (int m = 1; m < 3; ++m) {
      cpp_bin_float_132 x = y*m;
      cpp_bin_float_132 res = my_sqrt(x);
      cpp_bin_float_132 ref = sqrt(x);
      if (ref != res) {
        std::cout
          << "i=" << i << std::endl
          << std::setprecision(150)
          << "x=" << x << std::endl
          << "res=" << res << std::endl
          << "ref=" << ref << std::endl;
        return 1;
      }
    }
    y = nextafter(y, cpp_bin_float_132(2));
  }
  std::cout << "[2] o.k\n";

  y = 1;
  for (int i=0; i < 10000;i++) {
    for (int m = 1; m < 3; ++m) {
      cpp_bin_float_132 x = y*m;
      cpp_bin_float_132 res = my_sqrt(x);
      cpp_bin_float_132 ref = sqrt(x);
      if (ref != res) {
        std::cout
          << "i=" << i << std::endl
          << std::setprecision(150)
          << "x=" << x << std::endl
          << "res=" << res << std::endl
          << "ref=" << ref << std::endl;
        return 1;
      }
    }
    y = nextafter(y, cpp_bin_float_132(0));
  }
  std::cout << "[3] o.k\n";

  std::mt19937_64 rndGen;
  std::uniform_int_distribution<uint64_t> rndDistr(0, uint64_t(-1));
  auto rndFunc = std::bind ( rndDistr, std::ref(rndGen) );

  for (int i=0; i < 100000;i++) {
    cpp_bin_float_132 x(1);
    for (int wi = 0; wi < 7; ++wi)
      reinterpret_cast<uint64_t*>(&x.backend().bits())[wi] = rndFunc();
    reinterpret_cast<uint64_t*>(&x.backend().bits())[6] |= uint64_t(1) << 63;
    x.backend().exponent() = int(rndFunc() & 0x7FFF) - 0x4000;

    cpp_bin_float_132 res = my_sqrt(x);
    cpp_bin_float_132 ref = sqrt(x);
    if (ref != res) {
      std::cout
        << "i=" << i << std::endl
        << std::setprecision(150)
        << "x=" << x << std::endl
        << x.backend().exponent() << std::endl
        << int(x.backend().exponent() & -2) << std::endl
        << int(x.backend().exponent() & -2)/2 << std::endl
        << "res=" << res << std::endl
        << "ref=" << ref << std::endl;
      printf("%016llx %5d\n", reinterpret_cast<uint64_t*>(&x.backend().bits())[6], x.backend().exponent());
      printf("%016llx %5d\n", reinterpret_cast<uint64_t*>(&res.backend().bits())[6], res.backend().exponent());
      printf("%016llx %5d\n", reinterpret_cast<uint64_t*>(&ref.backend().bits())[6], ref.backend().exponent());
      return 1;
    }
  }
  std::cout << "[4] o.k\n";

  cpp_bin_float_132 xx[2];
  for (int i=0; i < 2500;i++) {
    cpp_bin_float_132 x(1);
    for (int wi = 0; wi < 7; ++wi)
      reinterpret_cast<uint64_t*>(&x.backend().bits())[wi] = rndFunc();
    reinterpret_cast<uint64_t*>(&x.backend().bits())[6] |= uint64_t(1) << 63;
    x.backend().exponent() = int(rndFunc() & 0x7FFF) - 0x4000;

    find_tip_points(x, xx);
    for (int k = 0; k < 2; ++k) {
      cpp_bin_float_132 res = my_sqrt(xx[k]);
      cpp_bin_float_132 ref = sqrt(xx[k]);
      if (ref != res) {
        std::cout
          << "i=" << i << " k=" << k << std::endl
          << std::setprecision(150)
          << "x="   << xx[k] << std::endl
          << "res=" << res << std::endl
          << "ref=" << ref << std::endl;
        return 1;
      }
    }
  }
  std::cout << "[5] o.k\n";

  return 0;
}

typedef boost::multiprecision::number<boost::multiprecision::backends::cpp_bin_float<64*16, boost::multiprecision::backends::digit_base_2, void, boost::int32_t> > float_x16_t;
static __inline float_x16_t fract(const float_x16_t& x)
{
  float_x16_t x0 = x;
  memset(
   reinterpret_cast<uint64_t*>(&x0.backend().bits()),
   0,
   sizeof(uint64_t)*9);
  return x - x0;
}


// produce two numbers whose square roots are very close to the middle between two
// representable cpp_bin_float_132 numbers, one that rounds up and another that rounds down
// both results are close (from below) to src*src
static void find_tip_points(const cpp_bin_float_132& src, cpp_bin_float_132 dst[2])
{
  typedef boost::multiprecision::number<boost::multiprecision::backends::cpp_bin_float<64*8, boost::multiprecision::backends::digit_base_2, void, boost::int32_t> > float_x8_t;
  float_x8_t x8 = src;
  reinterpret_cast<uint64_t*>(&x8.backend().bits())[0] = uint64_t(1) << 63;
  float_x16_t xx = float_x16_t(x8)*x8;

  int ulp_x_e = x8.backend().exponent()-447;
  cpp_bin_float_132 ulp = ldexp(cpp_bin_float_132(1), ulp_x_e);
  int ulp_xx_e = xx.backend().exponent()-447;
  float_x16_t ulp_xx = ldexp(float_x16_t(1), ulp_xx_e);
  float_x16_t thr    = ldexp(float_x16_t(1), ulp_xx_e-63);

  float_x16_t fract_xx = fract(xx);

  float_x16_t d1 = ulp_xx;
  float_x16_t d2 = x8*2*ulp;
  while (d2 >= d1)
    d2 -= d1;
  float_x16_t m1 = 0;
  float_x16_t m2 = 1;

  while (fract_xx > thr) {
    float_x16_t md = ceil(d1/d2);
    float_x16_t d3 = d2*md - d1;
    float_x16_t m3 = m2*md - m1;
    if (fract_xx >= d3) {
      float_x16_t mm = floor(fract_xx/d3)*m3;
      x8 -= float_x8_t(mm*ulp);
      xx  = float_x16_t(x8)*x8;
      fract_xx = fract(xx);
    }
    m1 = m2;
    m2 = m3;
    d1 = d2;
    d2 = d3;
  }
  dst[0] = cpp_bin_float_132(xx); // above tip point

  x8 -= float_x8_t(m2*ulp);
  xx = float_x16_t(x8)*x8;

  dst[1] = cpp_bin_float_132(xx); // below tip point
}
