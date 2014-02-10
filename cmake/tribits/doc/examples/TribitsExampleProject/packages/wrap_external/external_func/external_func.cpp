#include "external_func.hpp"

#include "SimpleCxx_HelloWorld.hpp"
#include "C.hpp"

std::string ExternalProj::external_func()
{
  return std::string("external_func ")
    + PackageWithSubpackages::getC()
    + " "
    + PackageWithSubpackages::depsC();
}
