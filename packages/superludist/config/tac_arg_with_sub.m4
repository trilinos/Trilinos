dnl @synopsis TAC_ARG_WITH_SUB( VALUE_NAME, VALUE_SUB_NAME, VALUE_DESCRIPTION, DEFAULT_VAL)
dnl
dnl Test for --with-${VALUE_NAME}-${VALUE_SUB_NAME} and set to DEFAULT_VAL value if value not specified.
dnl 
dnl Use this macro to set variables, such as library names and include paths, which
dnl include an underscore.
dnl
dnl This file was based on tac_arg_enable_feature.m4 by James Willenbring
dnl @author Ken Stanley <ksstanl@sandia.gov>
dnl
AC_DEFUN([TAC_ARG_WITH_SUB],
[
AC_ARG_WITH([$1-$2],
AC_HELP_STRING([--with-$1-$2],[$3 (default is [$4])]),
tac_with_$1_$2=$withval, tac_with_$1_$2=$4)
])

