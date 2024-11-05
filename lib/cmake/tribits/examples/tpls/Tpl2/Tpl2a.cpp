#include "Tpl2a.hpp"

#include "Tpl1.hpp"

std::string Tpl2::a_itsme()
{
  return "Tpl2a";
}

std::string Tpl2::a_deps()
{
  return Tpl1::itsme();
}
