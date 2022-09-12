#include <iostream>
#include <string>

#include "A.hpp"

int main() {
  std::string label_A = WithSubpackages::getA();
  std::string deps_A  = WithSubpackages::depsA();
  std::cout << "A label is: " << label_A << std::endl;
  std::cout << "A deps are: " << deps_A << std::endl;
  std::cout << "A special value: " << WithSubpackages::specialValue() << std::endl;
  return 0;
}
