#pragma once

#include <mpfr.h>

#ifdef __cplusplus
extern "C" {
#endif

void float128_to_mpfr   (mpfr_t dst, const __float128*    src);
void mpfr_to_float128   (__float128*    dst, mpfr_t src);

#ifdef __cplusplus
}
#endif
