dnl @synopsis TAC_ARG_ENABLE_CAN_USE_PACKAGE(PACKAGE_NAME, FEATURE_DESCRIPTION, HAVE_NAME)
dnl
dnl Test for --enable-${PACKAGE_NAME} and set to no if not specified.
dnl Also calls AC_DEFINE to define HAVE_${HAVE_NAME} if value is not equal 
dnl to "no".  Note that if a package is enabled by default or because
dnl another package that is being built depends on that package, an explicit 
dnl --enable-package command is passed to the package configure scripts.
dnl We do not specify a default value here because that is determined at the
dnl Trilinos level.  Also, this way there is no maintenance required at the
dnl individual package level if a package becomes or ceases to be a package
dnl that is built by default.
dnl
dnl Use this macro to facilitate definition of options related to other 
dnl Trilinos packages that a package "can use" but does not depend on. 
dnl For example:
dnl 
dnl TAC_ARG_ENABLE_CAN_USE_PACKAGE(ifpack, [enable optional features that dependon Ifpack], IFAPCK)
dnl 
dnl will test for --enable-ifpack when configure is run.  If it is defined (and not set to "no")
dnl then HAVE_IFPACK will be defined, Otherwise HAVE_IFPACK will not be defined.
dnl
dnl @author James <jmwille@sandia.gov>
dnl
AC_DEFUN([TAC_ARG_ENABLE_CAN_USE_PACKAGE],
[
AC_ARG_ENABLE([$1],
AC_HELP_STRING([--enable-$1],[$2]),
ac_cv_use_$1=$enableval, ac_cv_use_$1=no)

AC_MSG_CHECKING(whether to build optional [$1] dependent code)

if test "X$ac_cv_use_$1" != "Xno"; then
  AC_MSG_RESULT(yes)
  AC_DEFINE([HAVE_$3],1,[Define if want to build with $1 enabled])
else
  AC_MSG_RESULT(no)
fi
])

