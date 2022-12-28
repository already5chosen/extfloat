#include "matmulq.h"

void matmulq(
 int M, int N, int K,
 const __float128 A[], // [K*M]
 const __float128 B[], // [N*K]
 __float128       C[]  // [N*M]
)
{
  for (int dst_r = 0; dst_r < M; ++dst_r) {
    __float128* dst        = &C[dst_r*N];
    const __float128* aSrc = &A[dst_r*K];
    const __float128* b = B;
    __float128 aVal = aSrc[0];
    for (int col = 0; col < N; ++col)
      dst[col] = b[col]*aVal;
    for (int r = 1; r < K; ++r) {
      aVal = aSrc[r];
      b += N;
      int N2 = N & (-2);
      for (int col = 0; col < N2; col+=2) {
        __float128 ab0 = b[col+0]*aVal;
        __float128 ab1 = b[col+1]*aVal;
        dst[col+0] += ab0;
        dst[col+1] += ab1;
      }
      if (N != N2)
        dst[N2] += b[N2]*aVal;
    }
  }
}
