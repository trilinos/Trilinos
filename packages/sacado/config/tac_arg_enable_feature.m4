dnl @synopsis TAC_ARG_ENABLE_FEATURE(FEATURE_NAME, FEATURE_DESCRIPTION, HAVE_NAME, DEFAULT_VAL)
dnl
dnl Test for --enable-${FEATURE_NAME} and set to DEFAULT_VAL value if feature not specified.
dnl Also calls AC_DEFINE to define HAVE_${HAVE_NAME} if value is not equal to "no"
dnl 
dnl Use this macro to help defining whether or not optional 
dnl features* should compiled.  For example:
dnl
dnl TAC_ARG_ENABLE_FEATURE(epetra, [Configure and build epetra], EPETRA, yes)
dnl 
dnl will test for --enable-epetra when configure is run.  If it is defined 
dnl and not set to "no" or not defined (default is "yes") then HAVE_EPETRA will
dnl be defined, if --enable-epetra is defined to be "no", HAVE_EPETRA will not
dnl be defined.
dnl
dnl *NOTE: epetra, aztecoo, komplex, ifpack, and other software found in
dnl subdirectories of Trilinos/packages are "packages" in their own right.
dnl However, these packages are also "features" of the larger package
dnl "Trilinos".  Therefore, when configuring from the Trilinos directory,
dnl it is appropriate to refer to these software packages as "features".
dnl
dnl This file was based on tac_arg_with_package.m4 by Mike Heroux
dnl @author James Willenbring <jmwille@sandia.gov>
dnl
AC_DEFUN([TAC_ARG_ENABLE_FEATURE],
[
AC_ARG_ENABLE([$1],
AC_HELP_STRING([--enable-$1],[$2 (default is [$4])]),
ac_cv_use_$1=$enableval, ac_cv_use_$1=$4)

AC_MSG_CHECKING(whether to use [$1])

if test "X$ac_cv_use_$1" != "Xno"; then
  AC_MSG_RESULT(yes)  
  AC_DEFINE([HAVE_$3],,[Define if want to build $1])
else
  AC_MSG_RESULT(no)
fi
])

