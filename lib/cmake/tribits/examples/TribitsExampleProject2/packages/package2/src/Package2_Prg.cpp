#include <iostream>

#include "Package2.hpp"

int main()
{
  std::cout << "Package2 Deps: " << Package2::deps() << "\n";
  return 0;
}
