#include <iostream>
#include <string>

#include "A.hpp"

using namespace PackageWithSubpackages;

int main() {
  std::string label_A = getA();
  std::string deps_A  = depsA();
  std::cout << "A label is: " << label_A << std::endl;
  std::cout << "A deps are: " << deps_A << std::endl;

  return 0;
}
