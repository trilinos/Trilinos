dnl @synopsis TAC_ARG_WITH_PACKAGE(FEATURE_NAME, FEATURE_DESCRIPTION, HAVE_NAME, DEFAULT_VAL)
dnl
dnl Test for --with-${FEATURE_NAME} and set to DEFAULT_VAL value if feature not specified.
dnl Also calls AC_DEFINE to define HAVE_FEI_${HAVE_NAME} if value is not equal to "no"
dnl 
dnl Use this macro to help defining whether or not interfaces for optional 
dnl package should compiled.  For example:
dnl
dnl TAC_ARG_WITH_PACKAGE(zoltan, [Enable Zoltan interface support], ZOLTAN, no)
dnl 
dnl will test for --with-zoltan when configure is run.  If it is defined 
dnl (and not set to "no") then HAVE_ZOLTAN will be defined, 
dnl Otherwise HAVE_ZOLTAN will not be defined.
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([TAC_ARG_WITH_PACKAGE],
[
AC_ARG_WITH([$1],
AC_HELP_STRING([--with-$1],[$2 (default is [$4])]),
ac_cv_use_$1=$withval, ac_cv_use_$1=$4)

AC_MSG_CHECKING(whether to use [$1])

if test "X$ac_cv_use_$1" != "Xno"; then
  AC_MSG_RESULT(yes)  
  AC_DEFINE([HAVE_FEI_$3],,[Define if want to build with $1 enabled])
else
  AC_MSG_RESULT(no)
fi
])

