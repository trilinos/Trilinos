#include <iostream>
#include <string>

#include "ExternalPkg.hpp"

int main() {
  std::cout << "ExteranlPkg deps are: " << ExternalPkg::deps() << "\n";
  return 0;
}
