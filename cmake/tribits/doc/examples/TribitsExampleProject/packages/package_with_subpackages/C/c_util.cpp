#include <iostream>

#include "B.hpp"
#include "A.hpp"

int main()
{
  std::cout
    << "Called c_util: "
    << PackageWithSubpackages::depsB() << " "
    << PackageWithSubpackages::depsA()
    << "\n";  
  return 0;
}
