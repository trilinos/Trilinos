#include <iostream>
#include <string>

#include "InsertedPkg.hpp"

int main() {
  std::cout << "InsertedPkg deps are: " << InsertedPkg::deps() << "\n";
  return 0;
}
