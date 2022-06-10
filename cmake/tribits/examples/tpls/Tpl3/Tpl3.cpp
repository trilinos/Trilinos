#include <sstream>

#include "Tpl2a.hpp"
#include "Tpl2b.hpp"
#include "Tpl3.hpp"

std::string Tpl3::itsme()
{
  return "Tpl3";
}

std::string Tpl3::deps()
{
  std::ostringstream oss;
  oss
    << Tpl2::a_itsme() << "{" << Tpl2::a_deps() << "}"
    << ", "
    << Tpl2::b_itsme() << "{" << Tpl2::b_deps() << "}";
  return oss.str();
}
