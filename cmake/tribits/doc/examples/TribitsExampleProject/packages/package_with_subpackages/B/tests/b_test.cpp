#include <iostream>
#include <string>

// I can see this, because they are in subpackage B and so am I!
#include "A.hpp"
#include "B.hpp"

using namespace PackageWithSubpackages;

int main() {
  std::string label_A = getA();
  std::string label_B = getB();
  std::string deps_B  = depsB();
  std::cout << "A label is: " << label_A << std::endl;
  std::cout << "B label is: " << label_B << std::endl;
  std::cout << "B deps are: " << deps_B << std::endl;

  return 0;
}
