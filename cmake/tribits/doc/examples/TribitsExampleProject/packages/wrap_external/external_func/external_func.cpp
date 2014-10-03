#include "external_func.hpp"

#include "SimpleCxx_HelloWorld.hpp"
#include "A.hpp"

std::string ExternalProj::external_func()
{
  return std::string("external_func ")
    + PackageWithSubpackages::getA()
    + " "
    + PackageWithSubpackages::depsA();
}
