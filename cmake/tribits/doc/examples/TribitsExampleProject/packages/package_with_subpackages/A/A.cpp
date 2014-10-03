#include "A.hpp"

#include "SimpleCxx_HelloWorld.hpp"

std::string PackageWithSubpackages::getA() {
  return std::string("A");
}

std::string PackageWithSubpackages::depsA() {
  return SimpleCxx::deps();
}
