dnl @synopsis AC_CXX_ABI_DEMANGLE
dnl
dnl If the gcc function abi::__cxa_demangle(...) exists,
dnl then define HAVE_GCC_ABI_DEMANGLE
dnl
dnl @version
dnl @author
dnl
AC_DEFUN([AC_CXX_ABI_DEMANGLE],
[AC_CACHE_CHECK(whether the compiler supports abi::__cxa_demangle(...),
ac_cv_cxx_abi_demangle,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <cxxabi.h>
#include <string>
namespace MyNamespace {
class MyClass {};
}
],[
const std::string
  mangledName = typeid(MyNamespace::MyClass).name();
int status;
char
  *_demangledName = abi::__cxa_demangle(mangledName.c_str(),0,0,&status);
const std::string demangledName(_demangledName);
free(_demangledName);
return ( demangledName == "MyNamespace::MyClass" ? 0 : 1 );
],
ac_cv_cxx_abi_demangle=yes, ac_cv_cxx_abi_demangle=no)
AC_LANG_RESTORE
])
if test "$ac_cv_cxx_abi_demangle" = yes; then
  AC_DEFINE(HAVE_GCC_ABI_DEMANGLE,,
            [define if the compiler supports abi::__cxa_demangle(...)])
fi
])
