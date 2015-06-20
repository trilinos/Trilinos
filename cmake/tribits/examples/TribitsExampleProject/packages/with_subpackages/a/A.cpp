#include "A.hpp"

#include "SimpleCxx_HelloWorld.hpp"

std::string WithSubpackages::getA() {
  return std::string("A");
}

std::string WithSubpackages::depsA() {
  return SimpleCxx::deps();
}
