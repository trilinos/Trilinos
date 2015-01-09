#include "B.hpp"

#ifdef HAVE_WITHSUBPACKAGESB_WITHSUBPACKAGESA
#  include "A.hpp"
#endif

#ifdef HAVE_WITHSUBPACKAGESB_EXTERNALPKG
#  include "ExternalPkg.hpp"
#endif

#include "SimpleCxx_HelloWorld.hpp"


std::string WithSubpackages::getB() {
  return std::string("B");
}


std::string WithSubpackages::depsB() {
  std::string B_deps;
#ifdef HAVE_WITHSUBPACKAGESB_WITHSUBPACKAGESA
  B_deps += (std::string("A ") + depsA() + std::string(" "));
#endif
#ifdef HAVE_WITHSUBPACKAGESB_EXTERNALPKG
  B_deps += (std::string("ExternalPkg ") + ExternalPkg::deps() + std::string(" "));
#endif
  B_deps += SimpleCxx::deps();
  return B_deps;
}
