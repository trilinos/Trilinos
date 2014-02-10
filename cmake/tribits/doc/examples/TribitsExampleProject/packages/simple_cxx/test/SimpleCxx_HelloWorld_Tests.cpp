
#include <iostream>
#include <sstream>
#include <cstdlib>
#include "SimpleCxx_HelloWorld.hpp"


int main() {

  bool success = true;
  std::cout << std::boolalpha;

  SimpleCxx::HelloWorld helloWorld;
  std::ostringstream oss;
  helloWorld.printHelloWorld(oss);
  std::cout << oss.str();
  const bool foundHelloWorld = (oss.str().find("Hello World") != std::string::npos);
  std::cout << "Found \"Hello World\" ? " << foundHelloWorld << "\n";
  if (!foundHelloWorld) success=false;


  if (success) {
    std::cout << "End Result: TEST PASSED\n";
  }
  else {
    std::cout << "End Result: TEST FAILED\n";
  }
  return (success ? EXIT_SUCCESS : EXIT_FAILURE);

}
