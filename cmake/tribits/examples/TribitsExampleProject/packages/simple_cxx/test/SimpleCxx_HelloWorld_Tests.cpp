
#include <iostream>
#include <sstream>
#include <cstdlib>
#include "SimpleCxx_HelloWorld.hpp"


#define TEST_FIND_SUBSTR_IN_STR(SUBSTR, STR) \
  { \
    const bool foundSubStr = ((STR).find(SUBSTR) != std::string::npos); \
    std::cout << "Found \"" SUBSTR "\" ? " << foundSubStr << "\n"; \
    if (!foundSubStr) success=false; \
  } \
  (void)(success)


int main() {

  bool success = true;
  std::cout << std::boolalpha;

  SimpleCxx::HelloWorld helloWorld;
  std::ostringstream oss;
  helloWorld.printHelloWorld(oss);
  std::cout << oss.str();

  TEST_FIND_SUBSTR_IN_STR("Hello World", oss.str());

#ifdef HAVE_SIMPLECXX___INT64
  TEST_FIND_SUBSTR_IN_STR("We have __int64", oss.str());
#endif

#ifdef HAVE_SIMPLECXX_DEBUG
  TEST_FIND_SUBSTR_IN_STR("Debug is enabled", oss.str());
#else
  TEST_FIND_SUBSTR_IN_STR("Release is enabled", oss.str());
#endif

  TEST_FIND_SUBSTR_IN_STR("Sqr(3) = 9", oss.str());

  if (success) {
    std::cout << "End Result: TEST PASSED\n";
  }
  else {
    std::cout << "End Result: TEST FAILED\n";
  }
  return (success ? EXIT_SUCCESS : EXIT_FAILURE);

}
