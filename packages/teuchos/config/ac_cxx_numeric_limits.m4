dnl @synopsis AC_CXX_NUMERIC_LIMITS
dnl
dnl If the compiler has <limits> and accepts std::numeric_limits<>, define HAVE_NUMERIC_LIMITS.
dnl
dnl @version
dnl @author
dnl
AC_DEFUN([AC_CXX_NUMERIC_LIMITS],
[AC_CACHE_CHECK(whether the compiler supports std::numeric_limits<>,
ac_cv_cxx_numeric_limits,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <limits>
],[
double eps = std::numeric_limits<double>::epsilon();],
 ac_cv_cxx_numeric_limits=yes, ac_cv_cxx_numeric_limits=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_numeric_limits" = yes; then
  AC_DEFINE(HAVE_NUMERIC_LIMITS,,
            [define if the compiler supports numeric_limits<>])
fi
])
