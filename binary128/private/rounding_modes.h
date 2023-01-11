#ifdef __amd64

#include <x86intrin.h>

static __inline
int fast_fegetround(void)
{
  return _MM_GET_ROUNDING_MODE();
}

enum {
  fast_FE_TONEAREST   = _MM_ROUND_NEAREST,
  fast_FE_DOWNWARD    = _MM_ROUND_DOWN,
  fast_FE_UPWARD      = _MM_ROUND_UP,
  fast_FE_TOWARD_ZERO = _MM_ROUND_TOWARD_ZERO,
};

#else

// Standard C function
// It is quite slow

#include <fenv.h>

static __inline
int fast_fegetround(void)
{
  return fegetround();
}

enum {
  fast_FE_TONEAREST   = FE_TONEAREST,
  fast_FE_DOWNWARD    = FE_DOWNWARD,
  fast_FE_UPWARD      = FE_UPWARD,
  fast_FE_TOWARD_ZERO = FE_TOWARD_ZERO,
};

#endif