#include <iostream>
#include <string>

#include "SharedOnly.hpp"

int main() {
  std::cout << "SharedOnly_test returns "
            << MixedSharedStaticLibs::sharedPassThrough("shared") << "\n";
  return 0;
}
