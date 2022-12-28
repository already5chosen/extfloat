#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cstdint>
#include <vector>
#include <random>
#include <functional>           // for std::bind
#include <algorithm>
#include <chrono>

extern "C" {
#include "matmulq.h"
}

int main(int argz, char** argv)
{
  int M = 401;
  int N = 401;
  int K = 667;
  int nIter_check = 7;

  for (int arg_i = 1; arg_i < argz; ++arg_i) {
    char* arg = argv[arg_i];
      //"alpha", "beta",
      //"lda", "ldb", "ldc",
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

  int nIter_meas = nIter_check < 3 ? 3 : nIter_check;
  const int MIN_WORKING_SET_SZ = 32 * 1000 * 1000;
  int sz = (M*N + M*K + N*K)*sizeof(__float128);
  if (sz * nIter_meas < MIN_WORKING_SET_SZ)
    nIter_meas = MIN_WORKING_SET_SZ / (sz * 2) * 2 + 1;

  printf("Running matmulq with M=%d, N=%d, K=%d. %d iterations.\n", M, N, K, nIter_meas);

  size_t A_SZ = K*M, B_SZ = N*K, C_SZ = N*M;
  std::vector<__float128> A(nIter_meas*A_SZ);
  std::vector<__float128> B(nIter_meas*B_SZ);
  std::vector<__float128> C(nIter_meas*C_SZ);
  // std::vector<__float128> srcC(nIter_meas*M*N);

  std::mt19937_64 rndGen;
  std::uniform_real_distribution<double> rndDistr(-1.0, 1.0);
  auto rndFunc = std::bind ( rndDistr, std::ref(rndGen) );
  for (size_t i = 0; i < nIter_meas*A_SZ; ++i)
    A[i] = rndFunc();
  for (size_t i = 0; i < nIter_meas*B_SZ; ++i)
    B[i] = rndFunc();
  for (size_t i = 0; i < nIter_meas*C_SZ; ++i)
    C[i] = i;
  // for (int i = 0; i < nIter_meas*M*ldc; ++i)
    // srcC[i] = rndFunc();

  // std::chrono::duration<int64_t, std::nano> dt[N_REPS];
  // std::vector<std::complex<double>> vecX1(n);
  // std::this_thread::sleep_for(std::chrono::milliseconds(20));
  // ::Sleep(25);
  std::chrono::steady_clock::now();
  std::vector<long long>  dt(nIter_meas);
  // std::vector<long long>  dc(nIter_meas);
  for (int it = 0; it < nIter_meas; ++it) {
    auto t0 = std::chrono::steady_clock::now();
    // long long c0 = __rdtsc();
    // std::this_thread::sleep_for(std::chrono::milliseconds(27));
    // ::Sleep(25);
    // matmulq(M, N, K, &A.data()[A_SZ*it], &B.data()[B_SZ*0], &C.data()[C_SZ*it]);
    matmulq(M, N, K, &A.data()[A_SZ*it], &B.data()[B_SZ*it], &C.data()[C_SZ*it]);
    // long long c1 = __rdtsc();
    auto t1 = std::chrono::steady_clock::now();
    std::chrono::duration<long long, std::nano> dur = t1 - t0;
    // dc[it] = c1-c0;
    dt[it] = dur.count();
    // fprintf(stderr, "*** %d\n", it); fflush(stderr);
  }

  // for (int it = 0; it < nIter_meas; ++it)
    // printf("%2d %lld %lld %.0f\n", it, dt[it], dc[it], dc[it]/3.4);
  std::nth_element(dt.data(), dt.data() + nIter_meas/2, dt.data() + nIter_meas); // median
  long long dtMed = dt[nIter_meas/2];
  printf("%.6f msec. %.3f MFLOP/s\n", dtMed*1e-6, size_t(M)*N*K*2e3/dtMed);

  return 0;
}
