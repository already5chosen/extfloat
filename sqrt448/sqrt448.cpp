#include <cstdint>
#include <cstring>

#include "sqrt448.h"

typedef unsigned __int128 uintx_t;

// return three middle words of square of 2N-word source number
// i.e. result = (a*a / 2**(64*(N-1)) % 2**192
// for N > 1 the result is inexact. It can bit up to 4*N-4 smaller than exact value
static inline
#ifdef __clang__
__attribute__((always_inline))
#endif
void sqrm3(uint64_t* __restrict dst, const uint64_t a[], int N)
{
  // . b b A
  // . c b b
  // . . c .
  // . . . .

  // process non-diagonal
  uintx_t x = uintx_t(a[0])*a[N*2-1];
  uintx_t sum10 = 0;
  uintx_t sum21 = x;
  for (int i = 0; i < N-1; ++i) {
    x = uintx_t(a[i])*a[2*N-3-i];   // i+j=2*N-3
    sum10 += uint64_t(x >> 64);
    x = uintx_t(a[i])*a[2*N-2-i];   // i+j=2*N-2
    sum21 += uint64_t(x >> 64);
    sum10 += uint64_t(x);
    x = uintx_t(a[i+1])*a[2*N-2-i]; // i+j=2*N-1
    sum21 += x;
    x = a[i+1]*a[2*N-1-i];          // i+j=2*N
    sum21 += x << 64;
  }
  sum10 += sum10;
  sum21 += sum21;

  // process diagonal
  x = uintx_t(a[N-1])*a[N-1];   // i+j=2*N-2
  sum21 += uint64_t(x >> 64);
  sum10 += uint64_t(x);
  x = a[N]*a[N];                // i+j=2*N
  sum21 += x << 64;

  // merge sums and store
  sum21 += uint64_t(sum10 >> 64);
  dst[0] = uint64_t(sum10);
  dst[1] = uint64_t(sum21);
  dst[2] = uint64_t(sum21>>64);
}

// Calculate y=sqr(x), where x= a*2**64 + 2**63, return y[8]
// The function is called rarely, so optimized for size, with little regard to speed
static
uint64_t sqr_tip_point(const uint64_t a[7])
{
  uintx_t acc = a[0];
  for (int K = 1; K < 7; ++K) {
    uintx_t accx = a[K];
    for (int i = 0, j = K-1; j >= 0; ++i, --j) {
      uintx_t x = uintx_t(a[i])*a[j];
      acc  += uint64_t(x);
      accx += uint64_t(x>>64);
    }
    acc = accx + uint64_t(acc>>64);
  }
  uint64_t r8 = uint64_t(acc);
  for (int i = 0, j = 6; j >= 0; ++i, --j)
    r8 += a[i]*a[j];
  return r8;
}

static uint64_t mulh(const uint64_t a, uint64_t b) {
  return uint64_t((uintx_t(a)*b) >> 64);
}

static uint64_t imulh(const int64_t a, int64_t b) {
  return int64_t((__int128(a)*b) >> 64);
}

// mulsuh - multiply signed 128-bit number by unsigned 64-bit number.
// return upper 128 bits of the product
// the 2nd input be should be >= 2**63
static __int128 mulsuh(const __int128 a, uint64_t b) {
  uint64_t a0 = uint64_t(a);
  int64_t  a1 = int64_t(a>>64);
  uintx_t ret = __int128(a1)*int64_t(b);
  ret += uintx_t(uint64_t(a1)) << 64;
  ret += uint64_t((uintx_t(a0)*b) >> 64);
  return __int128(ret);
}


// mul2sX2u - multiply signed 128-bit number by unsigned 128-bit number.
// return upper 128 bits of the product
static uintx_t mul2sX2u(const uint64_t a[2], const uint64_t b[2])
{ // start with unsigned multiplication
  uint64_t b0 = b[0], b1 = b[1];
  uintx_t sum;
  sum  = uint64_t(uintx_t(a[0])*b1 >> 64);
  sum += uint64_t(uintx_t(a[1])*b0 >> 64);
  sum += uintx_t(a[1])*b1;
  // unsigned multiplication done

  uint64_t s = int64_t(a[1]) < 0 ? uint64_t(-1) : 0;
  sum -= (uintx_t(b1 & s) << 64) | (b0 & s);
  // result converted to signed
  return sum;
}

// muluuh - multiply unsigned 128-bit number by unsigned 64-bit number.
// return upper 128 bits of the product
static uintx_t muluuh(uint64_t a0, uint64_t a1, uint64_t b) {
  return uintx_t(a1)*b + uint64_t((uintx_t(a0)*b) >> 64);
}

static
void incdecx(uint64_t* acc, int64_t incdec, int N)
{
  if (incdec < 0) {
    // decrement
    for (int i = 0; i < N; ++i) {
      uint64_t oldval = acc[i];
      acc[i] = oldval - 1;
      if (oldval != 0)
        break;
    }
  } else {
    // increment
    for (int i = 0; i < N; ++i) {
      uint64_t newval = acc[i] + 1;
      acc[i] = newval;
      if (newval != 0)
        break;
    }
  }
}

static
void addw(uint64_t* __restrict dst, const uint64_t* src, uint64_t w, int N)
{
  uintx_t sum = (uintx_t)src[0] + w;
  dst[0] = uint64_t(sum);
  for (int i = 1; i < N; ++i) {
    sum = uintx_t(src[i]) + uint64_t(sum >> 64);
    dst[i] = uint64_t(sum);
  }
}

// rsqrt64 - calculate round(ldexp(1/sqrt(ldexp(quad_t(u), -62)), 64))
// x - input in range [2**62:2**64-1]
// result is calculated with precision slightly better than +/-2
// rsqrt64(2**62) => 2**64-1
static uint64_t rsqrt64(uint64_t x)
{
  static const uint8_t rsqrt_tab[48] = {
    252, 245, 238, 232, 226, 221, 216, 211,
    207, 203, 199, 195, 192, 189, 185, 182,
    180, 177, 174, 172, 170, 167, 165, 163,
    161, 159, 157, 155, 154, 152, 150, 149,
    147, 146, 144, 143, 141, 140, 139, 137,
    136, 135, 134, 133, 132, 131, 130, 129,
  };
  int idx = (x >> 58) & 63;
  uint64_t y = rsqrt_tab[idx-16];  // y scaled by 2^8

  // 1st NR step - y = y*(3+eps - y*y*x)/2
  uint64_t yy  = y*y;                           // yy=y*y    scaled by 2^16
  uint64_t yyx = yy*(x >> (64-16));             // yyx=y*y*x scaled by 2^(62+16+16-64)=2^30
  y = (y * ((3u << 30)+(3u<<17)-yyx))>>(6+1);   // y         scaled by 2^(8+30-6)=2^32

  // 2nd NR step
  yy  = (y*y)>>32;                              // yy=y*y    scaled by 2^32
  yyx = (yy*(x >> 32)) >> 32;                   // yyx=y*y*x scaled by 2^(32+62-32-32)=2^30
  y   = (y * ((3u << 30)+(3u<<4)-yyx))>>(30+1); // y         scaled by 2^(32+30-30)=2^32

  // 3rd NR step
  // Use 2nd order polynomial: y =  y - y*(err/2 - 3/8*err**2) = y - y/2*(err - 3/4*err**2)
  yy = y*y;                                              // y*y           scaled by 2**64
  int64_t err = int64_t(uint64_t((uintx_t(yy)*x)>>62));  // err = y*y*x-1 scaled by 2**64
  uint64_t err2 = imulh(err, err);                       // err2 = err*err scaled by 2**64
  int64_t m2 = err2*3 - err*4;                           // m2 = 3*err**2 - err*4
  int64_t adj = imulh(m2, int64_t(uint64_t(y)<<29));     // adj = m2*(y/8)
  uint64_t ret = (y<<32) + adj;

  return ret;
}


static
void sub3w(uint64_t* __restrict dst, const uint64_t a[3],  const uint64_t b[3])
{
  uint64_t d0 = a[0] - b[0];
  uintx_t a21 = (uintx_t(a[2])<<64) | a[1];
  uintx_t b21 = (uintx_t(b[2])<<64) | b[1];
  uintx_t d21 = a21 - b21 - (a[0] < b[0]);
  dst[0] = d0;
  dst[1] = uint64_t(d21);
  dst[2] = uint64_t(d21 >> 64);
}

static inline
#ifdef __clang__
__attribute__((always_inline))
#endif
void NR_step(uint64_t* __restrict Sqrt, const uint64_t src[], const uint64_t Rsqrt[2], int N, int sh)
{ // improve precision of Sqrt from  128*N-eps bits to 128*(N+1)-eps bits
  uint64_t sqr[3];
  sqrm3(sqr, &Sqrt[2], N);      // sqr(sqrt) scaled by 2**190
  uint64_t err[3];
  sub3w(err, src, sqr);
  err[0] = (err[0] >> sh) | (err[1] << (64-sh));
  err[1] = (err[1] >> sh) | (err[2] << (64-sh));
  uintx_t adj = mul2sX2u(err, Rsqrt);
  uint64_t adj0 = uint64_t(uint64_t(adj) << sh);
  uint64_t adj1 = uint64_t(adj >> (64-sh));
  int64_t  adj2 = int64_t(uint64_t(adj>>64)) >> (64-sh);
  __int128 res2 = __int128(adj2) + Sqrt[2];
  Sqrt[0] = adj0;
  Sqrt[1] = adj1;
  Sqrt[2] = uint64_t(uintx_t(res2));

  int64_t incdec = int64_t(res2 >> 64);
  if (incdec != 0)
    incdecx(&Sqrt[3], incdec, N*2-1);
}

static void sqrt448(uint64_t* __restrict dst, const uint64_t* __restrict src, int exp1)
{
  uint64_t ssrc[8]; // source scaled by 2**(64*8-2)
  if (exp1) {
    ssrc[0] = 0;
    memcpy(&ssrc[1], src, sizeof(src[0])*7);
  } else {
    ssrc[0] =  src[0] << 63;
    for (int i = 1; i < 7; ++i)
      ssrc[i] = (src[i] << 63) | (src[i-1] >> 1);
    ssrc[7] = src[6] >> 1;
  }

  uint64_t Rsqrt1 = rsqrt64(ssrc[7]); // precision - 62 bits
  // Calculate sqrt() as src*invSqrt
  uint64_t Sqrt7 = mulh(Rsqrt1, ssrc[7])*2 + 2; // scaled by 2**63

  // improve precision of Sqrt to 64*2-eps bits
  uintx_t sqr = uintx_t(Sqrt7)*Sqrt7;
  __int128 err = ((uintx_t(ssrc[7]) << 64) | ssrc[6]) - sqr;
  __int128 adj = mulsuh(err, Rsqrt1);
   // Sqrt += err*invS/2
  uint64_t Sqrt6 = uint64_t(adj);
  Sqrt7 += uint64_t(uintx_t(adj)>>64);

  // improve precision of rsqrt to 64*2-eps bits
  err = muluuh(Sqrt6, Sqrt7, Rsqrt1) << 1; // err = sqrt*rsqrt-1, scaled by 2**128
  adj = mulsuh(err, Rsqrt1);               // adj = err*invS, scaled by 2**127
  uint64_t Rsqrt0 = ~uint64_t(adj);
  Rsqrt1 += ~uint64_t(uintx_t(adj)>>64);

  uint64_t Sqrt[8];
  Sqrt[6] = Sqrt6;
  Sqrt[7] = Sqrt7;
  uint64_t Rsqrt[2];
  Rsqrt[0] = Rsqrt0;
  Rsqrt[1] = Rsqrt1;
  NR_step(&Sqrt[4], &ssrc[4], Rsqrt, 1,  6); // improve precision to 64*4-eps bits
  NR_step(&Sqrt[2], &ssrc[2], Rsqrt, 2, 10); // improve precision to 64*6-eps bits
  NR_step(&Sqrt[0], &ssrc[0], Rsqrt, 3, 14); // improve precision to 64*8-eps bits
  // Sqrt scaled by 2**511

  uint64_t lsw = Sqrt[0];
  const uint64_t UINT64_MID = uint64_t(1) << 63;
  const uint64_t MAX_ERR    = 1 << 17; // probably less, but it does not cost much to be on the safe side
  if (lsw - (UINT64_MID-MAX_ERR) < MAX_ERR*2) {
    // We are very close to tip point
    // square (res*2**64+2**63), take sqr[8]
    lsw = sqr_tip_point(&Sqrt[1]) << 1;
  }
  // rounding
  addw(dst, &Sqrt[1], lsw >> 63, 7);
}

cpp_bin_float_132 my_sqrt(const cpp_bin_float_132& x) {
  if (x <= 0) {
    return 0; // do whatever you feel correct for negatives
  }
  int e = x.backend().exponent();
  cpp_bin_float_132 res(x);
  res.backend().exponent() = int(e & -2) / 2;
  sqrt448((uint64_t*)&res.backend().bits(), (uint64_t*)&x.backend().bits(), e & 1);

  return res;
}
