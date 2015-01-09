#include "ExternalPkg.hpp"

#include "SimpleCxx_HelloWorld.hpp"

std::string ExternalPkg::deps() {
  return SimpleCxx::deps();
}
