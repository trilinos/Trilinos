#ifndef stk_util_environment_Demangle_hpp
#define stk_util_environment_Demangle_hpp

#include <string>

#if __GNUC__ == 3 || __GNUC__ == 4
#define STK_USE_PLATFORM_DEMANGLER
#endif

namespace stk {

/**
 * @brief Function <b>demangle</b> returns the demangled C++ symbol from the mangled
 * C++ symbol.  The mangled named is obtained from the <b>type_info</b>
 * <b>name()</b> function.  From some compilers, the name is already demangled.
 *
 * @param symbol	a <b>char</b> const pointer to the symbol.
 *
 * @return		a <b>std::string</b> value of the demangled name.
 */
#ifdef STK_USE_PLATFORM_DEMANGLER
std::string demangle(const char *symbol);
#else
const char *demangle(const char *symbol);
#endif

} // namespace stk

#endif // stk_util_environment_Demangle_hpp
