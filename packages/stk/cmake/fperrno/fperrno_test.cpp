
#include <cfenv>
#include <cmath>
#include <iostream>

int main(int argc, char** argv)
{
  bool haveFpErrno = false;
  if (math_errhandling & MATH_ERRNO) {
    [[maybe_unused]] auto result = std::log(0.0);
    if (errno == ERANGE) {
      haveFpErrno = true;
      std::cout<<"ON"; //no newline, this output will set a cmake variable
    }
  }

  if (!haveFpErrno) {
    std::cout<<"OFF";
  }

  return 0;
}

