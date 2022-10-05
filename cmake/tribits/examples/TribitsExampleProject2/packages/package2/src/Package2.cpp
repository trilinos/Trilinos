#include "Package2.hpp"

#include <sstream>

#ifdef HAVE_PACKAGE2_TPL3
#  include "Tpl3.hpp"
#endif
#include "Package1.hpp"

std::string  Package2::itsme()
{
  return "Package2";
}

std::string Package2::deps()
{
  std::ostringstream oss_deps;
  oss_deps
    << Package1::itsme() << "{" << Package1::deps() << "}"
#ifdef HAVE_PACKAGE2_TPL3
    << ", "
    << Tpl3::itsme() << "{" << Tpl3::deps() << "}"
#endif
    ;
  return oss_deps.str();
}
