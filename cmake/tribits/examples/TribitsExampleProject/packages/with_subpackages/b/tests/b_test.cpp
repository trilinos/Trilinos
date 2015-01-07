#include <iostream>
#include <string>

#include "B.hpp"

int main() {
  using namespace WithSubpackages;
  std::string label_B = getB();
  std::string deps_B  = depsB();
  std::cout << "B label is: " << label_B << std::endl;
  std::cout << "B deps are: " << deps_B << std::endl;

  return 0;
}
