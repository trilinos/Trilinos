dnl
dnl Check to make sure nan and inf detection is supported.
dnl
dnl Author: Roger Pawlowski
dnl
dnl
AC_DEFUN([TAC_ARG_FINITE_NUMBER_CHECK],
[AC_CACHE_CHECK(whether the compiler supports isnan and isinf,
ac_cv_xyce_have_finite_number,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 ac_save_LIBS="$LIBS"
 LIBS="$LIBS -lm"
 AC_TRY_LINK([
#ifndef _ALL_SOURCE
 #define _ALL_SOURCE
#endif
#ifndef _XOPEN_SOURCE
 #define _XOPEN_SOURCE
#endif
#ifndef _XOPEN_SOURCE_EXTENDED
 #define _XOPEN_SOURCE_EXTENDED 1
#endif
#include <math.h>],[double x = 1.0; 
isnan(x); isinf(x);
return 0;],
 ac_cv_xyce_have_finite_number=yes, ac_cv_xyce_have_finite_number=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_xyce_have_finite_number" = yes; then
  AC_DEFINE(HAVE_NAN_INF_SUPPORT,,[define if the compiler supports isnan and isinf checks])
else
  echo "****************************************************"
  echo "** Warning: Your compiler doesn't support isnan() and "
  echo "** isinf().  We will supply a default checker but it "
  echo "** is *NOT* guaranteed to work on your platform!"
  echo "****************************************************"
fi

])

