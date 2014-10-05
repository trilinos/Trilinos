#include "C.hpp"
#include "A.hpp"
#include "B.hpp"

std::string PackageWithSubpackages::getC()
{
  return std::string("C");
}

std::string PackageWithSubpackages::depsC()
{
  return std::string("B ") + depsB();
}
