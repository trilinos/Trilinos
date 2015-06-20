#include <iostream>

#include "B.hpp"
#include "A.hpp"

int main()
{
  std::cout
    << "Called c_util: "
    << WithSubpackages::depsB() << " "
    << WithSubpackages::depsA()
    << "\n";  
  return 0;
}
