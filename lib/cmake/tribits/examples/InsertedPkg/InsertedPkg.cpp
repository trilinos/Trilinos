#include "InsertedPkg.hpp"

#include "SimpleCxx_HelloWorld.hpp"

std::string InsertedPkg::deps() {
  return (std::string("SimpleCxx ") + SimpleCxx::deps());
}
