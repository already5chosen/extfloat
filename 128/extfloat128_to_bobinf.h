#pragma once

#include "extfloat128.h"

template<typename T>
void convert_to_boost_bin_float(T* pDst, const extfloat128_t& src)
{
  if (src.m_exponent >= src.min_biased_exponent && src.m_exponent <= src.max_biased_exponent) {
    #if 0
    boost::multiprecision::uint128_t significand;
    import_bits(significand,
      reinterpret_cast<const uint8_t*>(src.m_significand),
      reinterpret_cast<const uint8_t*>(src.m_significand)+sizeof(src.m_significand),
      0, false);
    *pDst = T(boost::multiprecision::cpp_int(significand));
    pDst->backend().exponent() -= 127;
    pDst->backend().sign()    = src.m_sign;
    *pDst = ldexp(*pDst, src._get_exponent());
    #else
    boost::multiprecision::number<boost::multiprecision::cpp_bin_float<128, boost::multiprecision::backends::digit_base_2> > b(src.m_significand[1]);
    reinterpret_cast<uint64_t*>(b.backend().bits().limbs())[0] = src.m_significand[0];
    b.backend().exponent() -= 63;
    b.backend().sign()      = src.m_sign;
    *pDst = ldexp(T(b), src._get_exponent());
    #endif
  } else {
    if (src.m_exponent==src.zero_biased_exponent) {
      // zero
      *pDst = 0;
      pDst->backend().sign() = src.m_sign;
    } else if (src.m_significand[0]==0 && src.m_significand[1]==0) {
      // inf
      *pDst = std::numeric_limits<T>::infinity().backend();
      pDst->backend().sign() = src.m_sign;
    } else {
      // NaN
      *pDst = std::numeric_limits<T>::quiet_NaN();
    }
  }
}
