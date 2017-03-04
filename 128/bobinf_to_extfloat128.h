template<typename T>
void convert_from_boost_bin_float(extfloat128_t* pDst, const T& src)
{
  typedef boost::multiprecision::number<boost::multiprecision::cpp_bin_float<128, boost::multiprecision::backends::digit_base_2, void, boost::int64_t, INT_MIN, INT_MAX> > boost_quadfloat_t;
  boost_quadfloat_t qsrc(src);
  switch (fpclassify(qsrc)) {
    case FP_INFINITE:
      *pDst = extfloat128_t::inf();
      return;
    case FP_NAN:
      *pDst = extfloat128_t::nan();
      return;
    case FP_ZERO:
      *pDst = extfloat128_t::zero();
      pDst->m_sign = qsrc.backend().sign();
      return;
    default:
      break;
  }
  int64_t e = qsrc.backend().exponent();
  if (e < extfloat128_t::min_exponent_val) {
    *pDst = extfloat128_t::zero();
  } else if (e < extfloat128_t::min_exponent_val) {
    *pDst = extfloat128_t::inf();
  } else {
    pDst->_set_exponent(static_cast<int32_t>(e));
    pDst->m_significand[0] = reinterpret_cast<uint64_t*>(qsrc.backend().bits().limbs())[0];
    pDst->m_significand[1] = reinterpret_cast<uint64_t*>(qsrc.backend().bits().limbs())[1];
  }
  pDst->m_sign = qsrc.backend().sign();
}

template<typename T>
extfloat128_t from_boost_bin_float(const T& src)
{
  extfloat128_t x;
  convert_from_boost_bin_float(&x, src);
  return x;
}
