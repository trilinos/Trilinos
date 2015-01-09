#include <iostream>
#include "hello_world_lib.hpp"

int main() {

  bool success = true;

  const std::string rtn = HelloWorld::getHelloWorld();
  std::cout << "HelloWorld::getHelloWorld() = '"<<rtn<<"' == 'Hello World'? ";
  if (rtn == "Hello World!") {
    std::cout << "passed\n";
  }
  else {
    std::cout << "FAILED\n";
    success = false;
  }

  if (success) {
    std::cout << "All unit tests passed :-)\n";
  }
  else {
    std::cout << "At least one unit test failed :-(\n";
  }

}
