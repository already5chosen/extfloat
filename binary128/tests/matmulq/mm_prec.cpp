#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cstdint>
#include <vector>
#include <random>
#include <functional>           // for std::bind
#include <algorithm>
#include <chrono>
#include <thread>

#include <quadmath.h>
#include <mpfr.h>
#include "qmx2mpfr.h"

extern "C" {
#include "matmulq.h"
}

static __float128 make_f128(double h, double l) {
  __float128 x = h;
  return x*l*0x1.0p-53 + x;
}

static void make_mpfr(mpfr_t dst, double h, double l) {
  mpfr_set_d(dst, h, GMP_RNDN);
  mpfr_mul_d(dst, dst, l, GMP_RNDN);
  mpfr_mul_d(dst, dst, 0x1.0p-53, GMP_RNDN);
  mpfr_add_d(dst, dst, h, GMP_RNDN);
}

static
long long
matmulq_test(
 int M, int N, int K,
 const mpfr_t     A[], // [K*M]
 const mpfr_t     B[], // [N*K]
 const __float128 C[]  // [N*M]
);

int main(int argz, char** argv)
{
  int M = 401;
  int N = 401;
  int K = 667;
  int nIter_check = 1;

  for (int arg_i = 1; arg_i < argz; ++arg_i) {
    char* arg = argv[arg_i];
    static const char* prefTab[] = { "M", "N", "K", "cn" };
    const int prefTabLen = sizeof(prefTab)/sizeof(prefTab[0]);
    for (int pref_i = 0; pref_i < prefTabLen; ++pref_i) {
      const char* pref = prefTab[pref_i];
      size_t preflen = strlen(pref);
      if ( strncasecmp(pref, arg, preflen)==0 && arg[preflen]=='=') {
        // integer arguments
        char* endp;
        long val = strtol(&arg[preflen+1], &endp, 0);
        if (endp==&arg[preflen+1] || val <= 0) {
          fprintf(stderr, "Bad parameter '%s'. '%s' is not a positive number.\n", arg, &arg[preflen+1]);
          return 1;
        }
        switch (pref_i) {
          case 0: M = val; break;
          case 1: N = val; break;
          case 2: K = val; break;
          case 3: nIter_check = val; break;
          default:break;
        }
        break;
      }
    }
  }

  if (nIter_check < 1)
    nIter_check = 1;

  printf("Running matmulq with M=%d, N=%d, K=%d. %d iterations.\n", M, N, K, nIter_check);
  size_t A_SZ = K*M, B_SZ = N*K, C_SZ = N*M;
  std::vector<__float128> A(A_SZ);
  std::vector<__float128> B(B_SZ);
  std::vector<__float128> C(C_SZ);

  std::vector<mpfr_t> refA(A_SZ);
  std::vector<mpfr_t> refB(B_SZ);
  for (size_t i = 0; i < A_SZ; ++i)
    mpfr_init2(refA.data()[i], 113);
  for (size_t i = 0; i < B_SZ; ++i)
    mpfr_init2(refB.data()[i], 113);


  std::mt19937_64 rndGen;
  std::uniform_real_distribution<double> rndDistr(-1.0, 1.0);
  auto rndFunc = std::bind ( rndDistr, std::ref(rndGen) );
  long long nMism = 0;
  for (int it = 0; it < nIter_check; ++it) {
    for (size_t i = 0; i < A_SZ; ++i) {
      double h= rndFunc(), l=rndFunc();
      A[i] = make_f128(h, l);
      make_mpfr(refA[i], h, l);
    }
    for (size_t i = 0; i < B_SZ; ++i) {
      double h= rndFunc(), l=rndFunc();
      B[i] = make_f128(h, l);
      make_mpfr(refB[i], h, l);
    }
    for (size_t i = 0; i < C_SZ; ++i)
      C[i] = i;

    matmulq(M, N, K, A.data(), B.data(), C.data());
    nMism += matmulq_test(M, N, K, refA.data(), refB.data(), C.data());
  }

  if (nMism == 0)
    printf("o.k.\n");
  else
    printf("fail. %lld mismatches.\n", nMism);
  return 0;
}

static void
matmulq_test_step(
 int M, int N, int K,
 const mpfr_t     A[], // [K*M]
 const mpfr_t     B[], // [N*K]
 const __float128 C[], // [N*M]
 int              r,
 long long*       pRes)
{
  long long nMism = 0;
  MPFR_DECL_INIT(acc, 113);
  MPFR_DECL_INIT(mx,  113);
  for (int c = 0; c < N; ++c) {
    mpfr_set_zero(acc, 0);
    for (int k = 0; k < K; ++k) {
      mpfr_mul(mx, A[r*K+k], B[k*N+c], GMP_RNDN);
      mpfr_add(acc, acc, mx, GMP_RNDN);
    }
    __float128 ref;
    mpfr_to_float128(&ref, acc);
    if (ref != C[r*N+c]) {
      char resbuf[256];
      char refbuf[256];
      quadmath_snprintf(resbuf, sizeof(resbuf), "%+-45.28Qa", C[r*N+c]);
      quadmath_snprintf(refbuf, sizeof(refbuf), "%+-45.28Qa", ref);
      printf("%3d %3d %s %s\n", r, c, resbuf, refbuf);
      ++nMism;
    }
  }
  *pRes = nMism;
}

static
long long
matmulq_test(
 int M, int N, int K,
 const mpfr_t     A[], // [K*M]
 const mpfr_t     B[], // [N*K]
 const __float128 C[]  // [N*M]
)
{
  std::vector<long long> nMismAr(M);
  std::vector<std::thread> tid(M);
  for (int r = 0; r < M; ++r)
    tid[r] = std::thread(matmulq_test_step, M, N, K, A, B, C, r, &nMismAr.data()[r]);
  for (int r = 0; r < M; ++r)
    tid[r].join();
  long long nMism = 0;
  for (int r = 0; r < M; ++r)
    nMism += nMismAr[r];
  return nMism;
}
