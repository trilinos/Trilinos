#include "c_b_mixed_lang.hpp"

#include "wsp_c/C.hpp"
#include "b_mixed_lang.hpp"

std::string WithSubpackages::c_b_mixed_lang()
{
  return (depsC()+" "+b_mixed_lang());
}
