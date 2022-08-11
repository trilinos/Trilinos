#ifdef TRIBITSEXAPP_HAVE_SIMPLECXX
#  include "SimpleCxx_HelloWorld.hpp"
#endif
#ifdef TRIBITSEXAPP_HAVE_MIXEDLANG
#  include "MixedLang.hpp"
#endif
#ifdef TRIBITSEXAPP_HAVE_WITHSUBPACKAGESA
#  include "A.hpp"
#endif
#ifdef TRIBITSEXAPP_HAVE_WITHSUBPACKAGESB
#  include "B.hpp"
#endif
#ifdef TRIBITSEXAPP_HAVE_WITHSUBPACKAGESC
#  include "wsp_c/C.hpp"
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
#ifdef TRIBITSEXAPP_HAVE_WITHSUBPACKAGESC
  appendDepsStr(depsStr, "WithSubpackagesC:"+WithSubpackages::depsC());
#endif
#ifdef TRIBITSEXAPP_HAVE_WITHSUBPACKAGESB
  appendDepsStr(depsStr, "WithSubpackagesB:"+WithSubpackages::depsB());
#endif
#ifdef TRIBITSEXAPP_HAVE_WITHSUBPACKAGESA
  appendDepsStr(depsStr, "WithSubpackagesA:"+WithSubpackages::depsA());
#endif
#ifdef TRIBITSEXAPP_HAVE_MIXEDLANG
  appendDepsStr(depsStr, "MixedLang:"+tribits_mixed::mixedLang()); // Just print something
#endif
#ifdef TRIBITSEXAPP_HAVE_SIMPLECXX
  appendDepsStr(depsStr, "SimpleCxx:"+SimpleCxx::deps());
#endif
  // NOTE: The above all call functions from the libraries and requires that
  // both the header files be found at compile time and the libraries be found
  // at link time and runtime for this these function calls to work.

  std::cout << "Full Deps: " << depsStr << "\n";

  return 0;
}
