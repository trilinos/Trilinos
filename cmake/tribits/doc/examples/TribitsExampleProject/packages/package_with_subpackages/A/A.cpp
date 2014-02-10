// include the declaration for getA()
#include "A.hpp"

// define getA()
std::string PackageWithSubpackages::getA() {
  return std::string("A");
}

// define depsA()
std::string PackageWithSubpackages::depsA() {
  return std::string("");
}
