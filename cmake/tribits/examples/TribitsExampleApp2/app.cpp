#ifdef TRIBITSEXAPP2_HAVE_PACKAGE1
#  include "Package1.hpp"
#endif


#include <iostream>
#include <string>


void appendDepsStr(std::string &depsStr, const std::string &str)
{
  if (depsStr.length()) {
    depsStr += "; "+str;
  }
  else {
    depsStr = str;
  }
}


int main(int argc, char *argv[]) {
  // Get deps down the deps graph
  std::string depsStr;
#ifdef TRIBITSEXAPP2_HAVE_PACKAGE1
  appendDepsStr(depsStr, "Package1: "+Package1::deps());
#endif
  // NOTE: The above all call functions from the libraries and requires that
  // both the header files be found at compile time and the libraries be found
  // at link time and runtime for this these function calls to work.

  std::cout << "Full Deps: " << depsStr << "\n";

  return 0;
}
