#include <iostream>

#include "SimpleCxx_HelloWorld.hpp"

int main()
{
  const SimpleCxx::HelloWorld helloWorld;
  helloWorld.printHelloWorld(std::cout);
  return 0;
}
