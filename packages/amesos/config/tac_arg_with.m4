dnl @synopsis TAC_ARG_WITH( VALUE_NAME, VALUE_DESCRIPTION, DEFAULT_VAL)
dnl
dnl Test for --with-${VALUE_NAME} and set to DEFAULT_VAL value if value not specified.
dnl 
dnl Use this macro to set variables, such as library names and include paths.
dnl
dnl This file was based on tac_arg_enable_feature.m4 by James Willenbring
dnl @author Ken Stanley <ksstanl@sandia.gov>
dnl
AC_DEFUN([TAC_ARG_WITH],
[
AC_ARG_WITH([$1],
AC_HELP_STRING([--with-$1],[$2 (default is [$3])]),
tac_with_$1=$withval, tac_with_$1=$3)
])

