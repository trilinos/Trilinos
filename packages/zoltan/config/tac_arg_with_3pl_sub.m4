dnl @synopsis TAC_ARG_WITH_3PL_SUB( VALUE_NAME, VALUE_SUB_NAME, VALUE_DESCRIPTION)
dnl
dnl Test for --with-${VALUE_NAME}-${VALUE_SUB_NAME} and set to no if value not specified.
dnl 
dnl Use this macro to set variables, such as library names and include paths, which
dnl include an underscore.
dnl
dnl This file was based on tac_arg_with_sub.m4 by Ken Stanley
dnl @author Jim Willenbring <jmwille@sandia.gov>
dnl
AC_DEFUN([TAC_ARG_WITH_3PL_SUB],
[
AC_ARG_WITH([$1-$2],
AC_HELP_STRING([--with-$1-$2],[$3]),
tac_with_$1_$2=$withval, tac_with_$1_$2=no)
])

