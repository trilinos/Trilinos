// include the declaration for getB()
#include "A.hpp"
#include "B.hpp"

// define getB()
std::string PackageWithSubpackages::getB() {
  return std::string("B");
}

// define depsB()
std::string PackageWithSubpackages::depsB() {
  return std::string("A ") + depsA();
}
