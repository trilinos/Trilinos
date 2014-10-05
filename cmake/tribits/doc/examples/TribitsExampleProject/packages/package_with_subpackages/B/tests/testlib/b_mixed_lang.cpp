#include "b_mixed_lang.hpp"

#include "B.hpp"
#include "MixedLang.hpp"

std::string PackageWithSubpackages::b_mixed_lang()
{
  return getB()+" "+tribits_mixed::mixedLang(); 
}
