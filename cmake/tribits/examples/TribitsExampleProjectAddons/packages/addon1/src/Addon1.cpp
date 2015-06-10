#include "Addon1.hpp"

#include "SimpleCxx_HelloWorld.hpp"

std::string Addon1::getAddon1() {
  return std::string("Addon1");
}

std::string Addon1::depsAddon1() {
  return SimpleCxx::deps();
}
