void matmulq(
 int M, int N, int K,
 const __float128 A[], // [K*M]
 const __float128 B[], // [N*K]
 __float128       C[]  // [N*M]
);
