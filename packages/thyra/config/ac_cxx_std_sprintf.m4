dnl @synopsis AC_CXX_STD_SPRINTF
dnl
dnl If the compiler recognizes std::sprintf as a function for IO,
dnl define HAVE_STD_SPRINTF.  If this test fails, use sprintf with no std prefix
dnl Note that we try to compile two versions of this routine, one using cstdio and
dnl another using stdio.h.  This approach is used to eliminate the need to test which
dnl of the two header files is present.  If one or both is usable the test will return true.
dnl

AC_DEFUN([AC_CXX_STD_SPRINTF],
[AC_CACHE_CHECK([[whether the compiler recognizes std::sprintf as supported IO function]],
ac_cv_cxx_std_sprintf,
[ AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([
#include <cstdio>
#include <string>
],
[
     int x = 100;
     char s[100];
     std::sprintf(s, "%d", x);
],
 ac_cv_cxx_std_sprintf1=yes, ac_cv_cxx_std_sprintf1=no)

AC_TRY_COMPILE([
#include <stdio.h>
#include <string>
],
[
     int x = 100;
     char s[100];
     std::sprintf(s, "%d", x);
],
 ac_cv_cxx_std_sprintf2=yes, ac_cv_cxx_std_sprintf2=no)

if (test "$ac_cv_cxx_std_sprintf1" = yes || test "$ac_cv_cxx_std_sprintf2" = yes); then
 ac_cv_cxx_std_sprintf=yes
else
 ac_cv_cxx_std_sprintf=no
fi
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_std_sprintf" = yes; then
  AC_DEFINE(HAVE_STD_SPRINTF,,[define if std::sprintf is supported])
fi
])
