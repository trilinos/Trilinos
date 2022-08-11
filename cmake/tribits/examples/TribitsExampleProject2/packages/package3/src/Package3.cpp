#include "Package3.hpp"

#include <sstream>

#include "Package1.hpp"
#ifdef HAVE_PACKAGE3_PACKAGE2
#  include "Package2.hpp"
#endif

#include "Tpl2a.hpp"
#include "Tpl2b.hpp"
#ifdef HAVE_PACKAGE3_TPL4
#  include "Tpl4.hpp"
#endif

std::string  Package3::itsme()
{
  return "Package3";
}

std::string Package3::deps()
{
  std::ostringstream oss_deps;
  oss_deps
#ifdef HAVE_PACKAGE3_PACKAGE2
    << Package2::itsme() << "{" << Package2::deps() << "}"
    << ", "
#endif
    << Package1::itsme() << "{" << Package1::deps() << "}"
    << ", "
#ifdef HAVE_PACKAGE3_TPL4
    << Tpl4::itsme() << "{" << Tpl4::deps() << "}"
    << ", "
#endif
    << Tpl2::a_itsme() << "{" << Tpl2::a_deps() << "}"
    << ", "
    << Tpl2::b_itsme() << "{" << Tpl2::b_deps() << "}"
    ;
  return oss_deps.str();
}
