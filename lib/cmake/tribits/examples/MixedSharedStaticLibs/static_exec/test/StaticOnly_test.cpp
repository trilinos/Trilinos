#include <iostream>
#include <string>

#include "StaticOnly.hpp"

int main() {
  std::cout << "StaticOnly_test returns "
            << MixedSharedStaticLibs::staticPassThrough("static") << "\n";
  return 0;
}
