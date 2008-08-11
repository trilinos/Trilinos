# @synopsis TAC_ARG_ENABLE_CAN_USE_PACKAGE(PACKAGE_NAME, OPTIONAL_DEPENDENCY_NAME, AC_DEFINE, AM_CONDITIONAL, DEFAULT_BUILD, ALLOW_IMPLICIT_ENABLE, HELP_STRING, IMPLICIT_HELP_STRING)
# 

# Use this macro to facilitate definition of options related to other
# Trilinos packages that a package "can use" but does not depend on.

# This macro supports both explicit and implicit configure options.
# For example to build epetra support into nox (so that nox can call
# epetra) we can explicitly use the flag --enable-nox-epetra.  This
# requires the user to enable both packages and the support between
# them: --enable-nox --enable-epetra --enable-nox-epetra.	

# Some packages, to simplify build requirements, implicitly assume
# that epetra support in nox should be built if both --enable-nox and
# --enable-epetra are supplied.  Users can override this by using the
# explicit command --enable-nox-epetra.

# Usage:
#  TAC_ARG_ENABLE_CAN_USE_PACKAGE(PACKAGE_NAME, 
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
# For example, to force explicit configuration of epetra support in nox:
# 
#  TAC_ARG_ENABLE_EXPLICIT_CAN_USE_PACKAGE(nox, epetra, NOX_EPETRA, 
#                                          NOX_USING_EPETRA, no, 
#                                          [Builds epetra support into nox.], 
#                                          [DOES NOTHING!])
# 
# To allow both implicit and explicit configuration of epetra support in nox:
#
# TAC_ARG_ENABLE_EXPLICIT_CAN_USE_PACKAGE(nox, epetra, NOX_EPETRA, 
# NOX_USING_EPETRA, no, yes, 
# [Builds epetra support in nox.], 
# [Builds epetra support in nox.  Can be overridden with --enable-nox-epetra.])
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
AC_DEFUN([TAC_ARG_ENABLE_CAN_USE_PACKAGE],
[

dnl Check for implicit enabling of optional package  
AC_ARG_ENABLE([$2],
AC_HELP_STRING([--enable-$2],[$7]),
ac_cv_implicit_use_$2=$enableval, 
ac_cv_implicit_use_$2=no)

dnl If implicit enabling is used, set that as teh default
if test "X$5" != "Xno"; then
  ac_cv_$1_using_$2_default=$ac_cv_implicit_use_$2
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
  AC_DEFINE([HAVE_$3],1,[Define if want to build with $2 enabled])
else
  AC_MSG_RESULT(no)
fi

AM_CONDITIONAL($4, test "X$ac_cv_use_$2" = "Xyes")

])
