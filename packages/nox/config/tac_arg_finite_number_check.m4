dnl
dnl Check to make sure nan and inf detection is supported.
dnl
dnl Author: Roger Pawlowski
dnl
dnl
AC_DEFUN([TAC_ARG_FINITE_NUMBER_CHECK],
[

AC_CACHE_CHECK(whether the compiler supports isnan(),
ac_cv_xyce_have_nan_support,
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
#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif],[double x = 1.0; 
isnan(x);
return 0;],
 ac_cv_xyce_have_nan_support=yes, ac_cv_xyce_have_nan_support=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_xyce_have_nan_support" = yes; then
  AC_DEFINE(HAVE_NAN_SUPPORT,,[define if the compiler supports the isnan() function])
else
  echo "****************************************************"
  echo "** Warning: Your compiler doesn't support isnan()."
  echo "** We will supply a default checker but it is "
  echo "** *NOT* guaranteed to work on your platform"
  echo "** unless your machine is IEEE 748/754 compliant."
  echo "****************************************************"
fi

AC_CACHE_CHECK(whether the compiler supports isinf(),
ac_cv_xyce_have_inf_support,
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
#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif],[double x = 1.0; 
isinf(x);
return 0;],
 ac_cv_xyce_have_inf_support=yes, ac_cv_xyce_have_inf_support=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_xyce_have_inf_support" = yes; then
  AC_DEFINE(HAVE_INF_SUPPORT,,[define if the compiler supports the isinf() function])
else
  echo "****************************************************"
  echo "** Warning: Your compiler doesn't support isinf()."
  echo "** We will supply a default checker but it is "
  echo "** *NOT* guaranteed to work on your platform"
  echo "** unless your machine is IEEE 748/754 compliant."
  echo "****************************************************"
fi

])

