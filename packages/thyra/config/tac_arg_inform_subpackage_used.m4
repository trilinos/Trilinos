# @synopsis TAC_ARG_INFORM_SUBPACKAGE_USED(PACKAGE_NAME, OPTIONAL_DEPENDENCY_NAME, AC_DEFINE, AM_CONDITIONAL, ALLOW_IMPLICIT_ENABLE, HELP_STRING, IMPLICIT_HELP_STRING)
# 

# This macro determines if a sub-feature in another dependent package was
# enabled or not.  For example, we would like to know if support for package
# pb was enabled in package pa or not.  For this, we will assume that this is
# supported if we see --enable-pa and --enable-pb (and implicit enabling is
# allows) or --enable-pa --enable-pa-pb but not if we see --disable-pa-pb.

# Usage:
#  TAC_ARG_INFORM_SUBPACKAGE_USED(PACKAGE_NAME, 
#                                 OPTIONAL_PACKAGE_DEPENDENCY_NAME, 
#                                 AC_DEFINE, 
#                                 AM_CONDITIONAL, 
#                                 ALLOW_IMPLICIT_ENABLE, 
#                                 HELP_STRING, 
#                                 IMPLICIT_HELP_STRING)
#
# Where:
# PACKAGE_NAME - Name of package currently being enabled.
# OPTIONAL_PACKAGE_DEPENDENCY_NAME - Name of optinal package support that 
#                                    will be built in.
# AC_DEFINE - Name of define that will be put in the 
#             <PackageName>_ConfigDef.h file.  Note that a "HAVE_" 
#             will be automatically prepended to the name.
# AM_CONDITIONAL - Variable that will be set in the makefiles by an 
#                  AM_CONDITIONAL call.
# ALLOW_IMPLICIT_ENABLE - This can be used to turn off implicit enable 
#                         support.  Takes a "yes/no" argument.
# HELP_STRING - Help string for configure option:
#               --enable-<PACKAGE_NAME>-<OPTIONAL_PACKAGE_DEPENDENCY_NAME>
# IMPLICIT_HELP_STRING - Help string for implicit configure option:
#               --enable-<OPTIONAL_PACKAGE_DEPENDENCY_NAME>
#
# Results of calling this file:
#  1. An AM_CONDITIONAL will be set for the AM_CONDITIONAL argument defined
#     above.
#  2. An AC_DEFINE of the form: HAVE_<AC_DEFINE> where AC_DEFINE is the 
#     argument defined above.
#
# @author Roger Pawlowski <rppawlo@sandia.gov>
# Based on original verison by Jim Willenbring.
#
AC_DEFUN([TAC_ARG_INFORM_SUBPACKAGE_USED],
[

dnl Check for implicit enabling of base package  
AC_ARG_ENABLE([$1],
AC_HELP_STRING([--enable-$1],[$7]),
ac_cv_implicit_use_$1=$enableval, 
ac_cv_implicit_use_$1=no)

dnl Check for implicit enabling of subpackage  
AC_ARG_ENABLE([$2],
AC_HELP_STRING([--enable-$2],[$7]),
ac_cv_implicit_use_$2=$enableval, 
ac_cv_implicit_use_$2=no)

dnl If implicit enabling is used, set that as the default
if test "X$5" != "Xno"; then
  if test "X$ac_cv_implicit_use_$1" != "Xno"; then
    ac_cv_$1_using_$2_default=$ac_cv_implicit_use_$2
  else
    ac_cv_$1_using_$2_default=no
  fi
else
  ac_cv_$1_using_$2_default=no
fi

AC_ARG_ENABLE([$1-$2],
AC_HELP_STRING([--enable-$1-$2],[$6]),
ac_cv_use_$2=$enableval, 
ac_cv_use_$2=$ac_cv_$1_using_$2_default)

AC_MSG_CHECKING(whether to build optional [$2] dependent code in [$1])

if test "X$ac_cv_use_$2" != "Xno"; then
  AC_MSG_RESULT(yes)
  AC_DEFINE([HAVE_$3],1,[Define if want to build with $1 enabled])
else
  AC_MSG_RESULT(no)
fi

AM_CONDITIONAL($4, test "X$ac_cv_use_$2" = "Xyes")

])
