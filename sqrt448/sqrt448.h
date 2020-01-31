#include <boost/multiprecision/cpp_bin_float.hpp>

typedef boost::multiprecision::number  
  < boost::multiprecision::backends::cpp_bin_float<
   64*7,  
   boost::multiprecision::backends::digit_base_2, 
   void, boost::int32_t> 
  > cpp_bin_float_132;
  
cpp_bin_float_132 my_sqrt(const cpp_bin_float_132& x);
  
