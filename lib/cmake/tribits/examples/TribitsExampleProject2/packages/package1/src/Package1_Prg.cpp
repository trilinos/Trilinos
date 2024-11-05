#include <iostream>
#include <string>

#include "Package1.hpp"

int main(int argc, char* argv[])
{
  std::cout << "Package1 Deps: " << Package1::deps() << "\n";
  for (int arg_i = 0; arg_i < argc; ++arg_i) {
    std::cout << argv[arg_i+1] << "\n";
  }
  return 0;
}
