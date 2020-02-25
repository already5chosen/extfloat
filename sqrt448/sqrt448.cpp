#include <cstdint>
#include <cstring>

#include "sqrt448.h"

typedef unsigned __int128 uintx_t;

template<int N, int Nr = N>
static
void sqrx(
  uint64_t* r,       // [Nr*2]
  const uint64_t* a  // [N]
) {
  if (Nr < 1 || Nr > N)
    return;

  // first sum up non-diag elements
  uintx_t x01 = uintx_t(a[0])*a[1];
  r[1] = uint64_t(x01);
  uintx_t acc = uint64_t(x01>>64);

  const int N1 = Nr*2 < N ? Nr*2 : N;
  for (int k = 2; k < N1; ++k) {
    uintx_t x = uintx_t(a[0])*a[k];
    acc += uint64_t(x);
    uintx_t accx = uint64_t(x>>64);
    for (int i = 1, j = k-1; i < j; ++i, --j) {
      x = uintx_t(a[i])*a[j];
      acc  += uint64_t(x);
      accx += uint64_t(x>>64);
    }
    r[k] = uint64_t(acc);
    acc = accx + uint64_t(acc>>64);
  }
  const int N2 = Nr*2 < N*2-2 ? Nr*2 : N*2-2;
  for (int k = N1; k < N2; ++k) {
    uintx_t x = uintx_t(a[k+1-N])*a[N-1];
    acc += uint64_t(x);
    uintx_t accx = uint64_t(x>>64);
    for (int i = k+2-N, j = N-2; i < j; ++i, --j) {
      x = uintx_t(a[i])*a[j];
      acc  += uint64_t(x);
      accx += uint64_t(x>>64);
    }
    r[k] = uint64_t(acc);
    acc = accx + uint64_t(acc>>64);
  }
  if (N*2-2 < Nr*2)
    r[N*2-2] = uint64_t(acc);
  // at that point acc < 2**64

  // add diag elements
  uintx_t x00 = uintx_t(a[0])*a[0];
  r[0] = uint64_t(x00);

  x00 = uint64_t(x00>>64);
  x00 += r[1]; x00 += r[1];
  r[1] = uint64_t(x00);
  uint64_t c = uint64_t(x00>>64);

  const int N3 = Nr < N-1 ? Nr : N-1;
  for (int kh = 1; kh < N3; ++kh) {
    uintx_t x = uintx_t(a[kh])*a[kh];
    uint64_t nd0 = r[kh*2+0],  nd1 = r[kh*2+1];
    x += nd0; x += nd0;
    uintx_t sum0 = uintx_t(uint64_t(x)) + c;
    r[kh*2+0] = uint64_t(sum0);
    uintx_t sum1 = uintx_t(uint64_t(x>>64)) + uint64_t(sum0>>64);
    sum1 += nd1; sum1 += nd1;
    r[kh*2+1] = uint64_t(sum1);
    c  = uint64_t(sum1 >> 64);
  }

  if (N-1 < Nr) {
    uintx_t xhh = uintx_t(a[N-1])*a[N-1] + c;
    xhh += r[N*2-2]; xhh += r[N*2-2];
    r[N*2-2] = uint64_t(xhh);
    r[N*2-1] = uint64_t(xhh >> 64);
  }
}

static
void sqrx(uint64_t* r, uint64_t a)
{
  uintx_t x = uintx_t(a)*a;
  r[0] = uint64_t(x);
  r[1] = uint64_t(x>>64);
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

template<int N>
static
void lshift_0or1(uint64_t* __restrict buf, int sh)
{
  sh &= 1;
  for (int i = N-1; i > 0; --i)
    buf[i] = (buf[i] << sh) | ((buf[i-1] >> 63) & sh) ;
  buf[0] = buf[0] << sh;
}

template<int N>
static
void rshift_0or1(uint64_t* __restrict buf, int sh)
{
  sh &= 1;
  for (int i = 0; i < N-1; ++i)
    buf[i] = (buf[i] >> sh) | ((buf[i+1] & sh) << 63);
  buf[N-1] = buf[N-1] >> sh;
}

template<int N>
static
void xorx(uint64_t* __restrict buf, uint64_t w)
{
  for (int i = 0; i < N; ++i)
    buf[i] ^= w;
}

template<int N>
static
void addx(uint64_t* __restrict acc, const uint64_t* x, uint64_t msw)
{
  memcpy(acc, x, sizeof(uint64_t)*N);
  uintx_t sum = (uintx_t)acc[N] + x[N];
  acc[N] = uint64_t(sum);
  for (int i = N+1; i < N*2; ++i) {
    sum = (uintx_t(acc[i]) + msw) + uint64_t(sum >> 64);
    acc[i] = uint64_t(sum);
  }
}

template<int N>
static
void addw(uint64_t* __restrict dst, const uint64_t* src, uint64_t w)
{
  uintx_t sum = (uintx_t)src[0] + w;
  dst[0] = uint64_t(sum);
  for (int i = 1; i < N; ++i) {
    sum = uintx_t(src[i]) + uint64_t(sum >> 64);
    dst[i] = uint64_t(sum);
  }
}

template<int Na, int Nb, int K, int resEnd>
static
uint64_t mulh_pass(uint64_t * __restrict r, const uint64_t* a, const uint64_t* b, uint64_t carryH)
{
  typedef unsigned __int128 uintx_t;
  const int I0 = K < Na ? 0   : K+1-Na;
  const int I1 = K < Nb ? K+1 : Nb;
  uintx_t accL = r[0];
  uintx_t accH = carryH;
  uintx_t x = uintx_t(a[K-I0])*b[I0];
  accL += uint64_t(x);
  accH += uint64_t(x >> 64);
  for (int i = I0+1; i < I1; ++i) {
    x = uintx_t(a[K-i])*b[i];
    accL += uint64_t(x);
    accH += uint64_t(x >> 64);
  }
  accH += uint64_t(accL >> 64);
  r[0] = uint64_t(accL);
  r[1] = uint64_t(accH);
  const bool last = K+2 >= resEnd;
  return mulh_pass<last?0:Na, last?0:Nb, last? 0:K+1, last?0:resEnd>(r+1, a, b, uint64_t(accH >> 64));
}

template<>
uint64_t mulh_pass<0, 0, 0, 0>(uint64_t * __restrict r, const uint64_t* a, const uint64_t* b, uint64_t carryH)
{
  return carryH;
}

template<int Na, int Nb, int resBeg, int resEnd>
static
void mulh(uint64_t * __restrict r, const uint64_t* a, const uint64_t* b)
{
  typedef unsigned __int128 uintx_t;
  uintx_t acc = 0;
  if (resBeg > 0) {
    const int K  = resBeg-1;
    const int I0 = K < Na ? 0   : K+1-Na;
    const int I1 = K < Nb ? K+1 : Nb;
    acc = uint64_t((uintx_t(a[K-I0])*b[I0]) >> 64);
    for (int i = I0+1; i < I1; ++i)
      acc += uint64_t((uintx_t(a[K-i])*b[i]) >> 64);
  }
  r[0] = uint64_t(acc);
  mulh_pass<Na, Nb, resBeg, resEnd>(r, a, b, uint64_t(acc>>64));
  if (resEnd < Na+Nb) {
    const int K  = resEnd-1;
    const int I0 = K < Na ? 0   : K+1-Na;
    const int I1 = K < Nb ? K+1 : Nb;
    uint64_t accL = a[K-I0]*b[I0];
    for (int i = I0+1; i < I1; ++i)
      accL += a[K-i]*b[i];
    r[resEnd-resBeg-1] += accL;
  }
}

template<int N>
static
void subx(uint64_t* __restrict buf, const uint64_t* x)
{
  uintx_t diff = uintx_t(buf[0]) - x[0];
  buf[0] = uint64_t(diff);
  for (int i = 1; i < N; ++i) {
    diff = (uintx_t(buf[i]) - x[i]) - (uint64_t(diff >> 64) & 1);
    buf[i] = uint64_t(diff);
  }
}

// rsqrt64 - calculate round(ldexp(1/sqrt(ldexp(quad_t(u), e-63)), 64))
// x - input in range [2**63:2**64-1]
// e - exponent in range [0:1]
// result is calculated with precision better than +/-2
// rsqrt64(2**63,0) => 2**64-1
static uint64_t rsqrt64(uint64_t x, int e)
{
  static const uint8_t rsqrt_tab[16] = {
    252,  245,  238,  232,  226,  221,  216,  211,
    207,  203,  199,  195,  192,  189,  185,  182,
  };
  int idx = (x >> 59) & 15;
  uint64_t y = rsqrt_tab[idx];  // y scaled by 2^8

  // 1st NR step - y = y*(3+eps - y*y*x)/2
  uint64_t yy  = y*y;                           // yy=y*y    scaled by 2^16
  uint64_t yyx = yy*(x >> (64-15));             // yyx=y*y*x scaled by 2^(63+16+15-64)=2^30
  y = (y * ((3u << 30)+(3u<<17)-yyx))>>(6+1);   // y         scaled by 2^(8+30-6)=2^32

  // 2nd NR step
  yy  = (y*y)>>32;                              // yy=y*y    scaled by 2^32
  yyx = (yy*(x >> 32)) >> 33;                   // yyx=y*y*x scaled by 2^(32+63-32-33)=2^30
  y   = (y * ((3u << 30)+(3u<<4)-yyx))>>(30+1); // y         scaled by 2^(32+30-30)=2^32

  // handle exponent
  static const uint32_t exp_adj_tabb[2] = { uint32_t(-1), 3037000500u }; // 0xb504f333f9de6484
  y = (y * exp_adj_tabb[e]) >> 32;

  // 3rd NR step
  // Use 2nd order polynomial: y =  y - y*(err/2 - 3/8*err**2) = y - y/2*(err - 3/4*err**2)
  x = x*2;                                               // x= x - 1      scaled by 2**64
  yy = (y * y) << e;
  int64_t err = int64_t(mulh(yy, x)+yy);                 // err = y*y*x-1 scaled by 2**64
  uint64_t err2 = imulh(err, err);                       // err2 = err*err scaled by 2**64
  int64_t m2 = err2*3 - err*4;                           // m2 = 3*err**2 - err*4
  int64_t adj = imulh(m2, int64_t(uint64_t(y)<<29));     // adj = m2*(y/8)
  uint64_t ret = (y<<32) + adj;

  return ret;
}

static void sqrt448(uint64_t* __restrict dst, const uint64_t* __restrict src, int exp1)
{
  uint64_t invSx1 = rsqrt64(src[6], exp1); // precision - 62 bits

  uint64_t tmpbuf[24];
  uint64_t* invS     = &tmpbuf[8*0];
  uint64_t* invS_sqr = &tmpbuf[8*1];
  uint64_t* prod3    = &tmpbuf[8*2];
  uint64_t* prod4    = invS_sqr;
  uint64_t* Sqrt     = &tmpbuf[8*0];
  uint64_t* Sqrt_sqr = &tmpbuf[8*1];
  uint64_t* Sqrt_adj = &tmpbuf[8*2];

  // improve precision to 127-eps bits
  invS[3] = invSx1;
  sqrx(invS_sqr, invSx1);
  mulh<2, 2, 2, 4>(prod3, invS_sqr, &src[5]);
  lshift_0or1<2>(prod3, exp1);
  prod3[1] ^= uint64_t(1) << 63;
  uint64_t s = prod3[1] & (uint64_t(1) << 63) ? uint64_t(-1) : 0;
  xorx<2>(prod3, s);
  mulh<2, 1, 1, 3>(prod4, prod3, &invS[3]);
  xorx<2>(prod4, ~s);
  addx<1>(&invS[2], prod4, ~s);

  // improve precision to 255-eps bits
  sqrx<2>(invS_sqr, &invS[2]);
  mulh<4, 4, 4, 7>(prod3, invS_sqr, &src[3]);
  s = prod3[2] & (uint64_t(1) << 63) ? uint64_t(-1) : 0;
  xorx<3>(prod3, s);
  mulh<3, 2, 2, 5>(prod4, prod3, &invS[2]);
  xorx<3>(prod4, ~s);
  lshift_0or1<3>(prod4, exp1);
  addx<2>(&invS[0], prod4, ~s);

  // Calculate sqrt() as src*invSqrt, precision: 255-eps bits
  mulh<4, 4, 4, 8>(&Sqrt[4], &invS[0], &src[3]); // scaled by 2**254 or 2**255
  lshift_0or1<4>(&Sqrt[4], exp1);                // scaled by 2**255

  // improve precision to 511-eps bits by Newton step
  sqrx<4, 3>(Sqrt_sqr, &Sqrt[4]);
  lshift_0or1<5>(Sqrt_sqr, exp1 ^ 1); // sqr(sqrt(src)) scaled by 2**511
  subx<4>(&Sqrt_sqr[1], src);         // Sqrt_sqr = err = Sqrt**2 - src
  s = Sqrt_sqr[4] & (uint64_t(1) << 63) ? uint64_t(-1) : 0; // s = sign(err)
  xorx<5>(Sqrt_sqr, s);               // Sqrt_sqr = abs(err)
  mulh<5, 4, 4, 9>(Sqrt_adj, Sqrt_sqr, &invS[0]); // Sqrt_adj = abs(err)*invS/2
  rshift_0or1<5>(Sqrt_adj, exp1 ^ 1); // proper scaling
  xorx<5>(Sqrt_adj, ~s);              // Sqrt_adj = -err*invS/2
  addx<4>(Sqrt, Sqrt_adj, ~s);        // Sqrt -= err*invS/2
  // Sqrt scaled by 2**511

  uint64_t lsw = Sqrt[0];
  const uint64_t UINT64_MID = uint64_t(1) << 63;
  const uint64_t MAX_ERR    = 1 << 20; // probably less, but it does not cost much to be on the safe side
  if (lsw - (UINT64_MID-MAX_ERR) < MAX_ERR*2) {
    // We are very close to tip point
    // square (res*2**64+2**63), take sqr[8]
    lsw = sqr_tip_point(&Sqrt[1]) << 1;
  }
  // rounding
  addw<7>(dst, &Sqrt[1], lsw >> 63);
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
