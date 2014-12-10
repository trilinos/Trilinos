#include "C.hpp"
#include "A.hpp"
#include "B.hpp"

std::string WithSubpackages::getC()
{
  return std::string("C");
}

std::string WithSubpackages::depsC()
{
  return std::string("B ") + depsB();
}
