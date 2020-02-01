#include <cstdint>
#include <cstring>

#include "sqrt448.h"

typedef unsigned __int128 uintx_t;

template<int N, int Nr = N>
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

void sqrx(uint64_t* r, uint64_t a)
{
  uintx_t x = uintx_t(a)*a;
  r[0] = uint64_t(x);
  r[1] = uint64_t(x>>64);
}

static uint64_t mulh(const uint64_t a, uint64_t b) {
  return uint64_t((uintx_t(a)*b) >> 64);
}

template<int N>
void lshift1(uint64_t* __restrict buf)
{
  for (int i = N-1; i > 0; --i)
    buf[i] = (buf[i] << 1) | (buf[i-1] >> 63) ;
  buf[0] = buf[0] << 1;
}

template<int N>
void lshift_0or1(uint64_t* __restrict buf, int sh)
{
  sh &= 1;
  for (int i = N-1; i > 0; --i)
    buf[i] = (buf[i] << sh) | ((buf[i-1] >> 63) & sh) ;
  buf[0] = buf[0] << sh;
}

template<int N>
void xorx(uint64_t* __restrict buf, uint64_t w)
{
  for (int i = 0; i < N; ++i)
    buf[i] ^= w;
}

template<int N>
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
void addw(uint64_t* __restrict dst, const uint64_t* src, uint64_t w)
{
  uintx_t sum = (uintx_t)src[0] + w;
  dst[0] = uint64_t(sum);
  for (int i = 1; i < N; ++i) {
    sum = uintx_t(src[i]) + uint64_t(sum >> 64);
    dst[i] = uint64_t(sum);
  }
}

static __inline uint64_t mulhx2(const uint64_t a, uint64_t b) {
  return uint64_t((uintx_t(a)*b) >> 63);
}

template<int Na, int Nb, int K, int resEnd>
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

// rsqrt64 - calculate round(ldexp(1/sqrt(ldexp(quad_t(u), e-63)), 64))
// x - input in range [2**63:2**64-1]
// e - exponent in range [0:1]
// result is calculated with precision slightly better than +/-2
// rsqrt64(2**63,0) => 2**64-1
static uint64_t rsqrt64(uint64_t x, int e)
{
  uint64_t y64 = uint64_t(-1);
  if ((x << 1) > 2) {
    static const uint8_t rsqrt_tab[17] = {
    255,  240,  226,  213,  201,  190,  180,  170,
    161,  153,  145,  137,  130,  123,  117,  111,
    105};
    int idx = (x >> 59) & 15;
    uint32_t t0 = rsqrt_tab[idx];
    uint32_t t1 = rsqrt_tab[idx+1];
    // interpolate
    uint32_t fr7 = (x >> (59-7)) & 127;
    uint32_t dt = (t0-t1)*fr7;
    uint32_t y = ((t0+257)<<7)-1-dt;
    // 1st NR step
    uint32_t yyx = ((y*y)*(x >> 32)) >> 33;
    y = (uint64_t(y) * ((3u << 30)-yyx)) >> 15;
    // 2nd NR step
    y = ((mulh((uint64_t(3)<<62)-(mulh(uint64_t(y)*y, x)>>1), uint64_t(y)<<32) >> 30)+1) >> 1;
    // 3rd NR step
    uintx_t pr3 = (uintx_t(3)<<95) - ((uintx_t(uint64_t(y)*y) * x) >> 32);
    y64 = uint64_t((pr3 >> 64)*y) + mulh(uint64_t(pr3),y);
  }
  // handle exponent
  static const uint64_t exp_adj_tab[2] = { 0, 0x95f619980c4336f7ull };
  y64 -= (mulh(y64, exp_adj_tab[e]) + 1)>>1;

  return y64;
}

static void sqrt448(uint64_t* __restrict dst, const uint64_t* __restrict src, int exp1)
{
  uint64_t invSx1 = rsqrt64(src[6], exp1); // precision - 62 bits

  uint64_t tmpbuf[24];
  uint64_t* invS     = &tmpbuf[8*0];
  uint64_t* invS_sqr = &tmpbuf[8*1];
  uint64_t* prod3    = &tmpbuf[8*2];
  uint64_t* prod4    = invS_sqr;

  // improve precision to 127-eps bits
  invS[7] = invSx1;
  sqrx(invS_sqr, invSx1);
  mulh<2, 2, 2, 4>(prod3, invS_sqr, &src[5]);
  lshift_0or1<2>(prod3, exp1);
  prod3[1] ^= uint64_t(1) << 63;
  uint64_t s = prod3[1] & (uint64_t(1) << 63) ? uint64_t(-1) : 0;
  xorx<2>(prod3, s);
  mulh<2, 1, 1, 3>(prod4, prod3, &invS[7]);
  xorx<2>(prod4, ~s);
  addx<1>(&invS[6], prod4, ~s);

  // improve precision to 255-eps bits
  sqrx<2>(invS_sqr, &invS[6]);
  mulh<4, 4, 4, 7>(prod3, invS_sqr, &src[3]);
  s = prod3[2] & (uint64_t(1) << 63) ? uint64_t(-1) : 0;
  xorx<3>(prod3, s);
  mulh<3, 2, 2, 5>(prod4, prod3, &invS[6]);
  xorx<3>(prod4, ~s);
  lshift_0or1<3>(prod4, exp1);
  addx<2>(&invS[4], prod4, ~s);

  // improve precision to 511-eps bits
  sqrx<4>(invS_sqr, &invS[4]);
  mulh<8, 7, 7, 12>(prod3, invS_sqr, &src[0]);
  lshift_0or1<5>(prod3, exp1);
  s = prod3[4] & (uint64_t(1) << 63) ? uint64_t(-1) : 0;
  xorx<5>(prod3, s);
  mulh<5, 4, 4, 9>(prod4, prod3, &invS[4]);
  xorx<5>(prod4, ~s);
  addx<4>(&invS[0], prod4, ~s);

  // Calculate sqrt() as src*invSqrt
  mulh<8, 7, 7, 15>(prod3, invS, &src[0]); // scaled by 2**510 or 2**511
  if (exp1)
    lshift1<8>(prod3); // scaled by 2**511
  // lshift_0or1<8>(prod3, exp1); // scaled by 2**511
  uint64_t lsw = prod3[0];
  const uint64_t UINT64_MID = uint64_t(1) << 63;
  const uint64_t MAX_ERR    = 1 << 20; // probably less, but it does not cost much to be on the safe side
  if (lsw - (UINT64_MID-MAX_ERR) < MAX_ERR*2) {
    // We are very close to tip point
    uint64_t* TP     = prod3;
    uint64_t* TP_sqr = &tmpbuf[8*0];
    TP[0] = UINT64_MID;
    // TP[7:0]  - Tip point TP = (res+0.5) scaled by 2**511
    sqrx<8, 5>(TP_sqr, TP);  // TP_sqr scaled by 2**1022
    lsw = TP_sqr[8] << 1;
  }
  // rounding
  addw<7>(dst, &prod3[1], lsw >> 63);
}

cpp_bin_float_132 my_sqrt(const cpp_bin_float_132& x) {
  if (x <= 0) {
    return 0; // do whatever you feel correct for negatives
  }
  int e = x.backend().exponent();
  cpp_bin_float_132 res(1);
  res.backend().exponent() = int(e & -2) / 2;
  sqrt448((uint64_t*)&res.backend().bits(), (uint64_t*)&x.backend().bits(), e & 1);

  return res;
}
