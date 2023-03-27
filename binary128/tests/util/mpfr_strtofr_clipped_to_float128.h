
#pragma once

#include <mpfr.h>

#ifdef __cplusplus
extern "C" {
#endif

int mpfr_strtofr_clipped_to_float128(mpfr_t rop, const char *nptr, char **endptr);
void manual_input_parser(mpfr_t rop, _Float128* rop128, const char *nptr);

#ifdef __cplusplus
}
#endif
