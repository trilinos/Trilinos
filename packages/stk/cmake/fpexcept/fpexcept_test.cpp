
#include <cfenv>
#include <cmath>
#include <iostream>

int main(int argc, char** argv)
{
  bool haveFpExcept = false;
  if (math_errhandling & MATH_ERREXCEPT) {
    [[maybe_unused]] auto result = std::log(0.0);
    if (std::fetestexcept(FE_DIVBYZERO)) {
      haveFpExcept = true;
      std::cout<<"ON"; //no newline, this output will set a cmake variable
    }
  }

  if (!haveFpExcept) {
    std::cout<<"OFF";
  }

  return 0;
}

