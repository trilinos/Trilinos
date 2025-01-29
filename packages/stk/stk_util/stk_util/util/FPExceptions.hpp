#ifndef STK_UTIL_FPEXCEPTIONS
#define STK_UTIL_FPEXCEPTIONS

#include <fenv.h>
#include <cfenv>
#include <cmath>
#include <string>
#include <cstring>
#include <cerrno>
#include <iostream>
#include <stk_util/stk_config.h>
#include <stk_util/util/ReportHandler.hpp>

namespace stk {
namespace util {

constexpr bool have_errno()
{
#ifdef STK_HAVE_FP_ERRNO
  return math_errhandling & MATH_ERRNO;
#else
  return false;
#endif
}

constexpr bool have_errexcept()
{
#ifdef STK_HAVE_FP_EXCEPT
  return math_errhandling & MATH_ERREXCEPT;
#else
  return false;
#endif
}

constexpr int FE_EXCEPT_CHECKS = FE_ALL_EXCEPT & ~FE_INEXACT;

std::string get_fe_except_string(int fe_except_bitmask);

inline void clear_fp_errors()
{
  if constexpr (have_errexcept())
  {
    // experimental results show calling std::feclearexcept is *very*
    // expensive, so dont call it unless needed.
    if (std::fetestexcept(FE_EXCEPT_CHECKS) > 0)
    {
      std::feclearexcept(FE_EXCEPT_CHECKS);
    }
  } else if constexpr (have_errno())
  {
    errno = 0;
  }
}

inline void throw_or_warn_on_fp_error(const char* fname = nullptr, bool warn=false, std::ostream& os = std::cerr)
{
  if constexpr (have_errexcept())
  {
    int fe_except_bitmask = std::fetestexcept(FE_EXCEPT_CHECKS);
    if (fe_except_bitmask != 0)
    {
      std::string msg = std::string(fname ? fname : "") + " raised floating point error(s): " + get_fe_except_string(fe_except_bitmask);
      clear_fp_errors();
      if (warn)
      {
        os << msg << std::endl;
      } else {
        STK_ThrowRequireMsg(fe_except_bitmask == 0, msg);
      }
    }
  } else if constexpr (have_errno())
  {
    if (errno != 0)
    {
      std::string msg = std::string(fname ? fname : "") + " raised floating point error(s) " + std::strerror(errno);
      clear_fp_errors();
      if (warn)
      {
        os << msg << std::endl;
      } else
      {
        STK_ThrowRequireMsg(errno == 0, msg);
      }
    }
  }

}

inline void warn_on_fp_error(const char* fname = nullptr, std::ostream& os = std::cerr)
{
  throw_or_warn_on_fp_error(fname, true, os);
}

inline void throw_on_fp_error(const char* fname = nullptr)
{
  throw_or_warn_on_fp_error(fname, false);
}


}
}

#endif
