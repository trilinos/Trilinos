# Use this macro to facilitate circular dependencies.

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
#                                 DEFAULT_BUILD, 
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
# DEFAULT_BUILD - Default for this support package.  Takes a "yes/no" 
#                 argument.
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
#                                          NOX_USING_EPETRA, no, no, 
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
# TAC_ARG_ENABLE_CIRCULAR_DEPENDENCIES(ml, nox, AC_DEFINE, AM_CONDITIONAL, default, help_string)
#
AC_DEFUN([TAC_ARG_ENABLE_CIRCULAR_DEPENDENCIES],
[

AC_ARG_ENABLE([$1-circular-dependency-$2],
AC_HELP_STRING([--enable-$1-circular-dependency-$2],[$6]),
ac_cv_circular_dependency_$2=$enableval, 
ac_cv_circular_dependency_$2=$5)

AC_MSG_CHECKING(whether circular dependencies are enabled between [$2] and [$1])

if test "X$ac_cv_circular_dependency_$2" != "Xno"; then
  AC_MSG_RESULT(yes)
  AC_DEFINE([HAVE_CIRCULAR_DEPENDENCY_$3],1,[Define for building circular dependencies into code.])
else
  AC_MSG_RESULT(no)
fi

AM_CONDITIONAL($4, test "X$ac_cv_circular_dependency_$2" = "Xyes")

])
