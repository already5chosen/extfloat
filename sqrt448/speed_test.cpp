#include <iostream>
#include <cstdlib>
#include "sqrt448.h"

int main(int argz, char** argv)
{
  int Niter = 1000000;
  if (argz > 1) {
    char* endp;
    int n = strtol(argv[1], &endp, 0);
    if (endp != argv[1] && n > 0)
      Niter = n;
  }
  cpp_bin_float_132 sum = 0;
  for (int i=1; i <= Niter;i++) {
    sum += my_sqrt(cpp_bin_float_132(i));
  }
  std::cout << "sum(x=1:" <<  Niter  << "){my_sqrt(x)} = "  << sum << std::endl;
  return 0;
}
