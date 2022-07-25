#include <iostream>

#include "B.hpp"
#include "A.hpp"

int main()
{
  std::cout
    << "Called c_util: B "
    << WithSubpackages::depsB() << " A "
    << WithSubpackages::depsA()
    << "\n";  
  return 0;
}
