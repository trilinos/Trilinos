#include <iostream>

#include "Package3.hpp"

int main()
{
  std::cout << "Package3 Deps: " << Package3::deps() << "\n";
  return 0;
}
