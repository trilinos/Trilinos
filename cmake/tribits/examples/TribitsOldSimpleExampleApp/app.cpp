#include "SimpleCxx_HelloWorld.hpp"
#include "MixedLang.hpp"
#include "wsp_c/C.hpp"

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
  appendDepsStr(depsStr, "WithSubpackages:"+WithSubpackages::depsC());
  appendDepsStr(depsStr, "MixedLang:"+tribits_mixed::mixedLang());
  appendDepsStr(depsStr, "SimpleCxx:"+SimpleCxx::deps());
  // NOTE: The above all call functions from the libraries and requires that
  // both the header files be found at compile time and the libraries be found
  // at link time and runtime for this these function calls to work.

  std::cout << "Full Deps: " << depsStr << "\n";

  return 0;
}
