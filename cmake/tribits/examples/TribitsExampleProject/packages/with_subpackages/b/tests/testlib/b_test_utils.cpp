#include "b_test_utils.hpp"
#include "SimpleCxx_HelloWorld.hpp"
#include <sstream>

std::string WithSubpackages::b_test_utils()
{
  std::ostringstream oss;
  SimpleCxx::HelloWorld().printHelloWorld(oss);
  return oss.str();
}
