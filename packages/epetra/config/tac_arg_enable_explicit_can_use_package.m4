dnl @synopsis TAC_ARG_ENABLE_EXPLICIT_CAN_USE_PACKAGE(PACKAGE_NAME, FEATURE_DESCRIPTION, HAVE_NAME)
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
dnl For example, to build epetra support in to nox:
dnl 
dnl TAC_ARG_ENABLE_EXPLICIT_CAN_USE_PACKAGE(nox, epetra, 
dnl [builds epetra support into nox], EPETRA, no)
dnl 
dnl will test for --enable-nox-epetra when configure is run.  If it is 
dnl defined (and not set to "$5")
dnl then HAVE_EPETRA will be defined, Otherwise HAVE_EPETRA will not 
dnl be defined.
dnl
dnl @author Roger Pawlowski <rppawlo@sandia.gov>
dnl
AC_DEFUN([TAC_ARG_ENABLE_EXPLICIT_CAN_USE_PACKAGE],
[
AC_ARG_ENABLE([$1-$2],
AC_HELP_STRING([--enable-$1-$2],[$3]),
ac_cv_use_$2=$enableval, 
ac_cv_use_$2=$6)

AC_MSG_CHECKING(whether to build optional [$2] dependent code in [$1])

if test "X$ac_cv_use_$2" != "Xno"; then
  AC_MSG_RESULT(yes)
  AC_DEFINE([HAVE_$4],1,[Define if want to build with $1 enabled])
else
  AC_MSG_RESULT(no)
fi

AM_CONDITIONAL($5, test "X$ac_cv_use_$2" = "Xyes")

])
