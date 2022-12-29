#include "matmulq.h"
#include "mulsq.h"
#include "addsubvq.h"

void matmulq(
 int M, int N, int K,
 const __float128 A[], // [K*M]
 const __float128 B[], // [N*K]
 __float128       C[]  // [N*M]
)
{
  enum { MULS_N = 128 };
  for (int dst_r = 0; dst_r < M; ++dst_r) {
    __float128* dst        = &C[dst_r*N];
    const __float128* aSrc = &A[dst_r*K];
    const __float128* b = B;
    __float128 aVal = aSrc[0];
    mulSq(dst, b, aVal, N);
    for (int r = 1; r < K; ++r) {
      aVal = aSrc[r];
      b += N;
      for (int col = 0; col < N; col += MULS_N) {
        int nCol = N - col < MULS_N ? N - col : MULS_N;
        __float128 prod[MULS_N];
        mulSq(prod, &b[col], aVal, nCol);
        addvq(&dst[col], &dst[col], prod, nCol);
      }
    }
  }
}
