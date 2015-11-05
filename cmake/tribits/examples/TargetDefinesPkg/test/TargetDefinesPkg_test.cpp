#ifndef  TARGETDEFINESPKG_RETURN
# define TARGETDEFINESPKG_RETURN  "DEFAULT"
#endif

#include <iostream>
#include <string>

#include "TargetDefinesPkg.hpp"

int main() {
  std::cout << "TargetDefinesPkg_test returns "
            << TargetDefinesPkg::passThrough(TARGETDEFINESPKG_RETURN) << "\n";
  return 0;
}
