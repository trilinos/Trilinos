// include the declaration for getB()
#include "C.hpp"
#include "A.hpp"
#include "B.hpp"

// define getC()
std::string PackageWithSubpackages::getC() {
  return std::string("C");
}

// define depsC()
std::string PackageWithSubpackages::depsC() {
  return std::string("B ") + depsB();
}
