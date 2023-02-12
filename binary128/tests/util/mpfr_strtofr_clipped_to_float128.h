
#pragma once

#include <mpfr.h>

#ifdef __cplusplus
extern "C" {
#endif

int mpfr_strtofr_clipped_to_float128(mpfr_t rop, const char *nptr, char **endptr);

#ifdef __cplusplus
}
#endif
