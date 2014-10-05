#include "B.hpp"

#ifdef HAVE_PACKAGEWITHSUBPACKAGESSUBPACKAGEB_PACKAGEWITHSUBPACKAGESSUBPACKAGEA
#  include "A.hpp"
#endif

#include "SimpleCxx_HelloWorld.hpp"


std::string PackageWithSubpackages::getB() {
  return std::string("B");
}


std::string PackageWithSubpackages::depsB() {
  std::string B_deps;
#ifdef HAVE_PACKAGEWITHSUBPACKAGESSUBPACKAGEB_PACKAGEWITHSUBPACKAGESSUBPACKAGEA
  B_deps += (std::string("A ") + depsA() + std::string(" "));
#endif
  B_deps += SimpleCxx::deps();
  return B_deps;
}
