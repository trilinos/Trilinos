dnl @synopsis AC_CXX_HAVE_STD
dnl
dnl If the compiler supports ISO C++ standard library (i.e., can include the
dnl files iostream, map, iomanip and cmath}), define HAVE_STD.
dnl
dnl @version $Id$
dnl @author Luc Maisonobe
dnl
AC_DEFUN([AC_CXX_HAVE_STD],
[AC_CACHE_CHECK(whether the compiler supports ISO C++ standard library,
ac_cv_cxx_have_std,
[AC_REQUIRE([AC_CXX_NAMESPACES])
 AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <iostream>
#include <map>
#include <iomanip>
#include <cmath>
#ifdef HAVE_NAMESPACES
using namespace std;
#endif],[return 0;],
 ac_cv_cxx_have_std=yes, ac_cv_cxx_have_std=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_have_std" = yes; then
  AC_DEFINE(HAVE_STD,,[define if the compiler supports ISO C++ standard library])
fi
])
