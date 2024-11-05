#include "SimpleCxx_HelloWorld_Tests.hpp"


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

#ifdef HAVE_SIMPLECXX_SIMPLETPL
  TEST_FIND_SUBSTR_IN_STR("Cube(3) = 27", oss.str());
#endif

  srand(time(NULL));
  std::cout
    << "<DartMeasurement type=\"numeric/double\" name=\"sqr_rand\">"
    << HeaderOnlyTpl::sqr(rand() % 10)  // Random number between 1 and 100
    << "</DartMeasurement>\n";
  // NOTE: The above produces a numeric test measurement that change each time
  // it calls and produce an interesting numeric plot on CDash.

  if (success) {
    std::cout << "End Result: TEST PASSED\n";
  }
  else {
    std::cout << "End Result: TEST FAILED\n";
  }
  return (success ? EXIT_SUCCESS : EXIT_FAILURE);

}
