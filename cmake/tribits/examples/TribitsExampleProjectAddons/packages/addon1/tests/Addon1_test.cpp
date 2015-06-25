#include <iostream>
#include <string>

#include "Addon1.hpp"

int main() {
  std::string label_Addon1 = Addon1::getAddon1();
  std::string deps_Addon1  = Addon1::depsAddon1();
  std::cout << "Addon1 label is: " << label_Addon1 << "\n";
  std::cout << "Addon1 deps are: " << deps_Addon1 << "\n";
  return 0;
}
