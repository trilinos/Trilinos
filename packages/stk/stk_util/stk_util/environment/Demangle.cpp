/*------------------------------------------------------------------------*/
/*                 Copyright 2010 - 2011 Sandia Corporation.              */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_util/environment/Demangle.hpp>
#include <stdlib.h>

#if __GNUC__ == 3 || __GNUC__ == 4
#include <cxxabi.h>
#endif

// #if defined __xlC__
// #include <demangle.h>
// #endif

namespace stk {

#ifdef STK_USE_PLATFORM_DEMANGLER

#if defined(__GNUC__)

#if (__GNUC__ == 3)
std::string
demangle(
  const char *	symbol)
{
#ifdef PURIFY_BUILD
  return symbol;
#else
  std::string   s;
  int		status = 0;

  char *demangled_symbol = abi::__cxa_demangle(symbol, 0, 0, &status);

  if (demangled_symbol) {
    s = std::string(demangled_symbol);
    free(demangled_symbol);
  }

  if (status != 0)
    s = std::string(symbol);

  return s;
#endif
}

#elif (__GNUC__ == 4)
std::string
demangle(
  const char *	symbol)
{
#ifdef PURIFY_BUILD
  return symbol;
#else
  std::string   s;

  int		status;

  char *demangled_symbol = __cxxabiv1::__cxa_demangle(symbol, 0, 0, &status);

  if (demangled_symbol) {
    s = std::string(demangled_symbol);
    free(demangled_symbol);
  }

  if (status != 0)
    s = std::string(symbol);

  return s;
#endif
}

#endif // (__GNUC__ == 3)

#elif defined __xlC__
std::string
demangle(
  const char *	symbol)
{
  return symbol;
// #ifdef PURIFY_BUILD
//   return symbol;
// #else
//   char *rest;

//   Name *name = Demangle(symbol, rest) ;

//   std::string s(name ? name->Text() : symbol);

//   delete name;

//   return s;
// #endif
}

#endif // defined __GNUC__

#else
const char *demangle(const char *symbol) {
  return symbol;
}
#endif // STK_USE_PLATFORM_DEMANGLER

} // namespace stk
