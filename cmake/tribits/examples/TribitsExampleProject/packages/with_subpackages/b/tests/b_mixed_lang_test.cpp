#include "b_mixed_lang.hpp"
#include "b_test_utils.hpp"

#include <iostream>

int main() 
{
  std::cout
    << WithSubpackages::b_mixed_lang()
    << " " << WithSubpackages::b_test_utils()
    << "\n";
}
