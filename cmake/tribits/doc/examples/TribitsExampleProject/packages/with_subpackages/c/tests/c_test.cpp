#include <iostream>
#include <string>

// I can see this, because they are in subpackage C and so am I!
#include "A.hpp"
#include "B.hpp"
#include "C.hpp"

using namespace WithSubpackages;

int main() {
  std::string label_A = getA();
  std::string label_B = getB();
  std::string label_C = getC();
  std::string deps_C  = depsC();
  std::cout << "A label is: " << label_A << std::endl;
  std::cout << "B label is: " << label_B << std::endl;
  std::cout << "C label is: " << label_C << std::endl;
  std::cout << "C deps are: " << deps_C << std::endl;

  return 0;
}
