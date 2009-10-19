#include <stk_util/environment/Demangle.hpp>
#include <stdlib.h>

#if __GNUC__ == 3 || __GNUC__ == 4
#include <cxxabi.h>
#endif

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

#endif // defined(__GNUC__)

#else
const char *demangle(const char *symbol) {
  return symbol;
}
#endif // STK_USE_PLATFORM_DEMANGLER

} // namespace stk
