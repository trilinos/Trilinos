#include "A.hpp"
#include "WithSubpackagesA_config.h"

#include "SimpleCxx_HelloWorld.hpp"

std::string WithSubpackages::getA() {
  return std::string("A");
}

std::string WithSubpackages::depsA() {
  return "SimpleCxx "+SimpleCxx::deps();
}

int WithSubpackages::specialValue() {
  return WITHSUBPACKAGESA_SPECIAL_VALUE;
}
