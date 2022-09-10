#include "B.hpp"

#include <iostream>
#include <string>

int main(int argc, char *argv[]) {
  std::string depsStr = WithSubpackages::getB()+" "+WithSubpackages::depsB();
  std::cout << "Util Deps: " << depsStr << "\n";
  return 0;
}
