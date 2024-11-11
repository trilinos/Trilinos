#include "FPExceptions.hpp"

namespace stk {
namespace util {

namespace {
void append_string(std::string& all_exceptions_string, const std::string& new_string)
{
  if (all_exceptions_string.size() == 0)
  {
    all_exceptions_string = new_string;
  } else
  {
    all_exceptions_string = all_exceptions_string + ", " + new_string;
  }
}
}

std::string get_fe_except_string(int fe_except_bitmask)
{
  std::string all_exceptions_string;
  if (fe_except_bitmask & FE_DIVBYZERO)
  {
    append_string(all_exceptions_string, "FE_DIVBYZERO");
  }

  if (fe_except_bitmask & FE_INEXACT)
  {
    append_string(all_exceptions_string, "FE_INEXACT");
  }

  if ( fe_except_bitmask & FE_INVALID)
  {
    append_string(all_exceptions_string, "FE_INVALID");
  }

  if (fe_except_bitmask & FE_OVERFLOW)
  {
    append_string(all_exceptions_string, "FE_OVERFLOW");
  }

  if (fe_except_bitmask & FE_UNDERFLOW)
  {
    append_string(all_exceptions_string, "FE_UNDERFLOW");
  }

  return all_exceptions_string;
}

}
}