dnl @synopsis TAC_ARG_ENABLE_OPTION(FEATURE_NAME, FEATURE_DESCRIPTION, HAVE_NAME, DEFAULT_VAL)
dnl
dnl Test for --enable-${FEATURE_NAME} and set to DEFAULT_VAL value if feature not specified.
dnl Also calls AC_DEFINE to define HAVE_${HAVE_NAME} if value is not equal to "no"
dnl
dnl Use this macro to facilitate definition of options in a package.  For example:
dnl 
dnl TAC_ARG_ENABLE_OPTION(threads, [enable shared memory threads], THREADS, no)
dnl 
dnl will test for --enable-threads when configure is run.  If it is defined (and not set to "no")
dnl then HAVE_THREADS will be defined, Otherwise HAVE_THREADS will not be defined.
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([TAC_ARG_ENABLE_OPTION],
[
AC_ARG_ENABLE([$1],
AC_HELP_STRING([--enable-$1],[$2 (default is [$4])]),
ac_cv_use_$1=$enableval, ac_cv_use_$1=$4)

AC_MSG_CHECKING(whether to use [$1])

if test "X$ac_cv_use_$1" != "Xno"; then
  AC_MSG_RESULT(yes)
  AC_DEFINE([HAVE_$3],1,[Define if want to build with $1 enabled])
else
  AC_MSG_RESULT(no)
fi
])

