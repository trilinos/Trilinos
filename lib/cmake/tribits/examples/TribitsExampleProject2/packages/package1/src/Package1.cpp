#include "Package1.hpp"
#include "Tpl1.hpp"

std::string Package1::itsme()
{
  return "Package1";
}

std::string Package1::deps()
{
  std::string deps;
  deps = Tpl1::itsme();
  return deps;
}
