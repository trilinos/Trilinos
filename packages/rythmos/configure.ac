# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

#@HEADER
# ***********************************************************************
#
#                           Rythmos Package
#                 Copyright (2006) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
#
# ***********************************************************************
#@HEADER

# ------------------------------------------------------------------------
# Initialization 
# ------------------------------------------------------------------------



# This must be the first line in configure.ac.
# Optional 3rd argument is email address for bugs.

#  #rythmos# - package name, version number, and e-mail address below
AC_INIT(rythmos, 1.0, tscoffe@sandia.gov)

# Hello World!
#  #rythmos# - package name below
echo "----------------------------------------"
echo "Running Rythmos Configure Script"
echo "----------------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.  Any file in that directory will do,
# preferably one that is unlikely to be removed or renamed.

# #np# - file name below
AC_CONFIG_SRCDIR(/src/Rythmos_Version.cpp)

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.

#
# #np# This a a subdirectory of Trilinos/packages/new_package that contains
# #np# a number of useful m4 files.  Any package that uses these m4 files
# #np# will have to add them to the repository.  The m4 files perform various
# #np# common checks and set common variable values to minimize redundant work.
#
AC_CONFIG_AUX_DIR(config)
#  #np# - Change file names in next line 
# Configure should create src/New_Package_config.h from src/New_Package_config.h.in

AM_CONFIG_HEADER(src/Rythmos_config.h:src/Rythmos_config.h.in)

# Allow users to specify their own "install" command.  If none is specified,
# the default is install-sh found in the config subdirectory.
                                                                                
AC_ARG_WITH(install,
 [AC_HELP_STRING([--with-install=INSTALL_PROGRAM],
 [Use the installation program INSTALL_PROGRAM rather the default that is provided.  For example --with-install="/path/install -p"])],
 [
   INSTALL=$withval
   INSTALL_PROGRAM=$withval
   INSTALL_SCRIPT=$withval
   INSTALL_DATA="$withval -m 644"
 ],)
                                                                                
# AM_MAINTAINER_MODE turns off maintainer-only makefile targets by
# default, and changes configure to understand a
# --enable-maintainer-mode option. --enable-maintainer-mode turns the
# maintainer-only targets back on. The maintainer-only makefile
# targets permit end users to clean automatically-generated files such
# as configure, which means they have to have autoconf and automake
# installed to repair the damage. AM_MAINTAINER_MODE makes it a bit
# harder for users to shoot themselves in the foot.

AM_MAINTAINER_MODE

# Define $build, $host, $target, etc

AC_CANONICAL_TARGET

# Use automake

#  #np# - Required version of automake.
AM_INIT_AUTOMAKE(1.8.2)

# Specify required version of autoconf.

AC_PREREQ(2.59)

# ------------------------------------------------------------------------
# Check to see if MPI enabled and if any special configuration done
# ------------------------------------------------------------------------

TAC_ARG_CONFIG_MPI

#  #np# - can eliminate compiler checks below if your package does not use the
#         language corresponding to the check.  Please note that if you use
#	  F77_FUNC to determine Fortran name mangling, you should not remove
#	  the Fortran compiler check or the check for Fortran flags.  Doing
#	  so will prevent the detection of the proper name mangling in some
#	  cases.
# ------------------------------------------------------------------------
# Checks for programs
# ------------------------------------------------------------------------

AC_PROG_CC(cc gcc)
AC_PROG_CXX(CC g++ c++ cxx)
AC_PROG_F77(f77 g77 gfortran f90 xlf90 f95)
AC_PROG_RANLIB

# ------------------------------------------------------------------------
# Checks for extra flags
# #np# - You can pull flag checks for languages not used in package.  Please 
#	  note that if you use F77_FUNC to determine Fortran name mangling, 
#	  you should not remove the Fortran compiler check or the check for 
#	  Fortran flags.  Doing so will prevent the detection of the proper 
#	  name mangling in some cases.
# ------------------------------------------------------------------------
# Check if --with-flags present, prepend any specs to FLAGS

TAC_ARG_WITH_FLAGS(ccflags, CCFLAGS)
TAC_ARG_WITH_FLAGS(cxxflags, CXXFLAGS)
TAC_ARG_WITH_FLAGS(cflags, CFLAGS)
TAC_ARG_WITH_FLAGS(fflags, FFLAGS)
TAC_ARG_WITH_LIBS
TAC_ARG_WITH_FLAGS(ldflags, LDFLAGS)

# ------------------------------------------------------------------------
# Alternate archiver
# ------------------------------------------------------------------------

TAC_ARG_WITH_AR

# ------------------------------------------------------------------------
# MPI link check
# ------------------------------------------------------------------------
TAC_ARG_CHECK_MPI

# ------------------------------------------------------------------------
# Checks for Makefile.export related systems
# ------------------------------------------------------------------------
TAC_ARG_ENABLE_EXPORT_MAKEFILES(yes)

#  #np# - Customize the checks for headers needed in package being added.
#NOTE: This section still needs to be updated for new_package, it is still 
# code from Epetra
# ------------------------------------------------------------------------
# Checks for header files - used in Epetra_ConfigDefs.h 
# ------------------------------------------------------------------------

AC_CHECK_HEADERS([cstdlib stdlib.h], [break],)
AC_CHECK_HEADERS([cstdio stdio.h], [break],)
AC_CHECK_HEADERS([cassert assert.h], [break],)
AC_CHECK_HEADERS([cmath math.h], [break],)
AC_CHECK_HEADERS([iomanip iomanip.h], [break],)
AC_CHECK_HEADERS([iostream iostream.h], [break],)
AC_CHECK_HEADERS([string string.h], [break],)

if test "X$ac_cv_header_cstdlib" != "Xyes" && test "X$ac_cv_header_stdlib_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdlib or stdlib.h.  Rythmos requires at least one or the other.]])
fi

if test "X$ac_cv_header_cstdio" != "Xyes" && test "X$ac_cv_header_stdio_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdio or stdio.h.  Rythmos requires at least one or the other.]])
fi

if test "X$ac_cv_header_cassert" != "Xyes" && test "X$ac_cv_header_assert_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cassert or assert.h.  Rythmos requires at least one or the other.]])
fi

if test "X$ac_cv_header_string" != "Xyes" && test "X$ac_cv_header_string_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find string or string.h.  Rythmos requires at least one or the other.]])
fi

if test "X$ac_cv_header_iostream" != "Xyes" && test "X$ac_cv_header_iostream_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find iostream or iostream.h.  Rythmos requires at least one or the other.]])
fi

if test "X${host}" != "Xtflop"; then
  if test "X$ac_cv_header_cmath" != "Xyes" && test "X$ac_cv_header_math_h" != "Xyes"; then
    AC_MSG_ERROR([[Cannot find cmath or math.h.  Rythmos requires at least one or the other.]])
  fi
else
  # TFLOP
  if test "X$ac_cv_header_iomanip" != "Xyes" && test "X$ac_cv_header_iomanip_h" != "Xyes"; then
    AC_MSG_ERROR([[Cannot find iomanip or iomanip.h.  Rythmos requires at least one or the other.]])
  fi
fi  

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE

# ------------------------------------------------------------------------
# Checks for compiler characteristics - Epetra_ConfigDefs.h uses these, 
# so it makes sense to leave them in - see PORTABILITY in this directory
# ------------------------------------------------------------------------
AC_CXX_NAMESPACES
AC_CXX_HAVE_STL
AC_CXX_BOOL
AC_CXX_MUTABLE

# ------------------------------------------------------------------------
# Checks for Makefile.export related systems
# ------------------------------------------------------------------------
TAC_ARG_ENABLE_EXPORT_MAKEFILES

# ---------------------------------------------------------------------
# These are pretty obscure
# ---------------------------------------------------------------------
AC_CXX_NEW_FOR_SCOPING
#AC_CXX_STD_IOS_BASE_FMTFLAGS
AC_CXX_STD_SPRINTF

# ------------------------------------------------------------------------
# Checks for special package flags
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_CAN_USE_PACKAGE(
  rythmos, nox, RYTHMOS_NOX, ENABLE_RYTHMOS_NOX, yes,
  [Build in NOX support in Rythmos.],
  [Build in NOX support in Rythmos.  Can be overridden with --disable-rythmos-nox.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(
  rythmos, stratimikos, RYTHMOS_STRATIMIKOS, ENABLE_RYTHMOS_STRATIMIKOS, yes
  ,[Build the Stratimikos adapater support in Rythmos.]
  ,[Build the Stratimikos adapater support in Rythmos.  Can be overridden with --disable-rythmos-stratimikos.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(
  rythmos, amesos, RYTHMOS_AMESOS, ENABLE_RYTHMOS_AMESOS, yes
  ,[Build the Amesos adapater support in Rythmos.]
  ,[Build the Amesos adapater support in Rythmos.  Can be overridden with --disable-rythmos-amesos.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(
  rythmos, aztecoo, RYTHMOS_AZTECOO, ENABLE_RYTHMOS_AZTECOO, yes
  ,[Build the AztecOO adapater support in Rythmos.]
  ,[Build the AztecOO adapater support in Rythmos.  Can be overridden with --disable-rythmos-aztecoo.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(
  rythmos, epetraext, RYTHMOS_EPETRAEXT, ENABLE_RYTHMOS_EPETRAEXT, yes
  ,[Build the EpetraExt adapater support in Rythmos.]
  ,[Build the EpetraExt adapater support in Rythmos.  Can be overridden with --disable-rythmos-epetraext.])

#TAC_ARG_ENABLE_CAN_USE_PACKAGE(
#  rythmos, belos, RYTHMOS_BELOS, ENABLE_RYTHMOS_BELOS, yes
#  ,[Build the Belos adapater support in Rythmos.]
#  ,[Build the Belos adapater support in Rythmos.  Can be overridden with --disable-rythmos-belos.])

TAC_ARG_INFORM_SUBPACKAGE_USED(
  epetraext, thyra, EPETRAEXT_THYRA, ENABLE_EPETRAEXT_THYRA, yes
  ,[See epetraext/configure --help]
  ,[See epetraext/configure --help])

TAC_ARG_INFORM_SUBPACKAGE_USED(
  nox, epetra, NOX_EPETRA, ENABLE_NOX_EPETRA, no
  ,[See nox/configure --help]
  ,[See nox/configure --help])

# ------------------------------------------------------------------------
# RYTHMOS interface libraries (uses config/tac_arg_enable_rythmos_interface.m4)
# ------------------------------------------------------------------------
#RYTHMOS_INTERFACE(nox, LIBRYTHMOSNOX, rythmos)
#RYTHMOS_INTERFACE(stratimikos, LIBRYTHMOSSTRATIMIKOS, rythmos)
#RYTHMOS_INTERFACE(amesos, LIBRYTHMOSAMESOS, rythmos)
#RYTHMOS_INTERFACE(aztecoo, LIBRYTHMOSAZTECOO, rythmos)
#RYTHMOS_INTERFACE(epetraext, LIBRYTHMOSEPETRAEXT, rythmos)
#RYTHMOS_INTERFACE(belos, LIBRYTHMOSBELOS, rythmos)

# 09/14/06 tscoffe:  I need nox-epetra to get the examples to work.

#AM_CONDITIONAL(ENABLE_RYTHMOS_STRATIMIKOS, test "X${ADDON_stratimikos}" != "Xno")
#AM_CONDITIONAL(HAVE_RYTHMOS_STRATIMIKOS, test "X${ADDON_stratimikos}" != "Xno")
#AM_CONDITIONAL(ENABLE_RYTHMOS_NOX, test "X${ADDON_nox}" != "Xno")

# ------------------------------------------------------------------------
# #np# newp_swahili is an example of a package feature which requires both
# a cpp #define command and one or more source or header files  - ie 
# conditional compilation which requires a change to the makefile.
# If configure is called with "--enable-newp_swahili", it will place 
# a #define HAVE_NEWP_SWAHILI in New_Package_config.h and add source 
# and header files as specified in Makefile.am 
#
# To add an option to your package, you will need to modify this file
# (i.e. configure.ac), and src/Makefile.am.  A bootstrap is also required.
# Look for SWAHILI in both files listed above.
#
# TAC_ARG_ENABLE_OPTION performs several tasks, including:
#   1)  Setting a help string
#   2)  Printing a message to standard out
#   3)  Defining HAVE_NEWP_SWAHILI (For use in New_Package_config.h)
#
# AM_CONDITIONAL makes HAVE_NEWP_SWAHILI available for use in Makefile.am
# ------------------------------------------------------------------------
#
#TAC_ARG_ENABLE_OPTION(newp_swahili, [enable swahili support], NEWP_SWAHILI, no)
#AM_CONDITIONAL(HAVE_NEWP_SWAHILI, [test "X$ac_cv_use_newp_swahili" != "Xno"])

# ------------------------------------------------------------------------
# Checks for special package flags
# ------------------------------------------------------------------------

#
#  #np# - These options are specific to Epetra and are left only as examples.
#         If uncommented, the following lines would define HAVE_FATAL_MESSAGES
#         in the (package)_config.h file.  Developers should use a 
#         more specific name such as 'new_package-fatal_messages', rather
#         than 'fatal_messages'.
#TAC_ARG_ENABLE_OPTION(fatal_messages, [enable printing of fatal error messages], 
#                   FATAL_MESSAGES, yes)
#TAC_ARG_ENABLE_OPTION(warning_messages, [enable printing of warning error messages], 
#                   WARNING_MESSAGES, no)
#TAC_ARG_ENABLE_OPTION(format_io, [Use package-defined IO formatting for output], 
#                   FORMAT_IO, yes)

#  #np# - These are specific to Epetra and are left only as examples.
#         If uncommented, the following line would define HAVE_ZOLTAN in
#         the (package)_config.h file.  Developers should use a more
#         specific name such as 'new_package-zoltan', rather than 'zoltan'.
# TAC_ARG_WITH_PACKAGE(zoltan, [Enable Zoltan interface support], ZOLTAN, no)
# TAC_ARG_WITH_PACKAGE(rtop, [Enable RTOp interface support], RTOP, no)

# ------------------------------------------------------------------------
# Checks if tests and examples should be built
# ------------------------------------------------------------------------

#  #np# - These options can disable the tests and examples of a package.
#  #np# - Packages that do not have tests or examples should #-out the 
#  #np# - option(s) that does (do) not apply.

# 09/19/06 tscoffe:  --enable-tests is checked by TAC_ARG_ENABLE_FEATURE_SUB_CHECK
#TAC_ARG_ENABLE_FEATURE(tests, [Build tests for all Trilinos packages (not all packages are sensitive to this option)], TESTS, yes)

# 09/20/06 tscoffe:  The examples are handled the same as tests below.
AC_ARG_ENABLE([examples],, ac_cv_use_examples=$enableval, ac_cv_use_examples=yes)
AC_ARG_ENABLE([rythmos-examples],
    AC_HELP_STRING([--enable-rythmos-examples],[Build Rythmos Examples (default is yes if --disable-examples is not specified), note:  the examples will not be built unless the following packages are also enabled:  epetraext, amesos, aztecoo, stratimikos, and epetraext-thyra]),ac_cv_use_rythmos_examples=$enableval, ac_cv_use_rythmos_examples=${ac_cv_use_examples})
AC_MSG_CHECKING(whether to build [rythmos-examples])

if test "X$ac_cv_use_rythmos_examples" != "Xno" &&
   test "X$ac_cv_use_epetraext" != "Xno" &&
   test "X$ac_cv_use_amesos" != "Xno" &&
   test "X$ac_cv_use_aztecoo" != "Xno" &&
   test "X$ac_cv_use_stratimikos" != "Xno" &&
   test "X$ac_cv_use_epetraext_thyra" != "Xno" ; then
  AC_MSG_RESULT(yes)
  AC_DEFINE([HAVE_RYTHMOS_EXAMPLES],,[Define if want to build rythmos-examples])
  AM_CONDITIONAL(BUILD_EXAMPLES, true)
else
  AC_MSG_RESULT(no)
  AM_CONDITIONAL(BUILD_EXAMPLES, false)
fi

# 09/20/06 tscoffe:  This is basically what TAC_ARG_ENABLE_FEATURE does but it
# checks all of my dependencies at the same time which is important for
# printing the correct message at configure time.
AC_ARG_ENABLE([tests],, ac_cv_use_tests=$enableval, ac_cv_use_tests=yes)
AC_ARG_ENABLE([rythmos-tests],
    AC_HELP_STRING([--enable-rythmos-tests],[Build Rythmos Tests (default is yes if --disable-tests is not specified), note:  all of the tests come from the examples, so if tests are enabled, the examples will be built (assuming all the dependencies are enabled)]),ac_cv_use_rythmos_tests=$enableval, ac_cv_use_rythmos_tests=${ac_cv_use_tests})
AC_MSG_CHECKING(whether to build [rythmos-tests])

if test "X$ac_cv_use_rythmos_tests" != "Xno" &&
   test "X$ac_cv_use_epetraext" != "Xno" &&
   test "X$ac_cv_use_amesos" != "Xno" &&
   test "X$ac_cv_use_aztecoo" != "Xno" &&
   test "X$ac_cv_use_stratimikos" != "Xno" &&
   test "X$ac_cv_use_epetraext_thyra" != "Xno" ; then
  AC_MSG_RESULT(yes)
  AC_DEFINE([HAVE_RYTHMOS_TESTS],,[Define if want to build rythmos-tests])
  AM_CONDITIONAL(BUILD_TESTS, true)
# 09/20/06 tscoffe:  The tests all come from the examples:
  AM_CONDITIONAL(BUILD_EXAMPLES, true)
else
  AC_MSG_RESULT(no)
  AM_CONDITIONAL(BUILD_TESTS, false)
fi

TAC_ARG_ENABLE_FEATURE(libcheck, [Check for some third-party libraries including BLAS and LAPACK.  (Cannot be disabled unless tests and examples are also disabled.)], LIBCHECK, yes)

# ------------------------------------------------------------------------
# Enable optional python wrappers
# ------------------------------------------------------------------------
TAC_ENABLE_PYTHON

# ------------------------------------------------------------------------
# Specify other directories 
# ------------------------------------------------------------------------

# 02/21/06 tscoffe:  LIBDIRS is OBSOLETE, use LDFLAGS instead.
# enable use of --with-libdirs="-Llibdir1 -Llibdir2 ..." to prepend to LDFLAGS
#TAC_ARG_WITH_LIBDIRS

# enable use of --with-incdirs="-Lincdir1 -Lincdir2 ..." to prepend to CPPFLAGS
TAC_ARG_WITH_INCDIRS

# #np# - Yet another opportunity to remove code if you aren't
# using Fortran
# Define F77_FUNC that will be used to link with Fortran subroutines, e.g., BLAS/LAPACK - trash WORKGXX 
AC_F77_WRAPPERS

# ------------------------------------------------------------------------
# Checks for libraries
# ------------------------------------------------------------------------

# If tests, examples and libcheck are disabled, we don't have to check
# for these libraries.

# #np# -
# If a package does not have tests or examples, the corresponding check(s)
# should be pulled out of the "if" statement below.
#if test "X$ac_cv_use_new_package_examples" != "Xno" || test "X$ac_cv_use_libcheck" != "Xno"; then
if test "X$ac_cv_use_rythmos_tests" != "Xno" || test "X$ac_cv_use_rythmos_examples" != "Xno" || test "X$ac_cv_use_libcheck" != "Xno"; then

# #np# - 
# To add any libraries that you want to check for, replace pow below 
# with the name of an entry point in that library.
# checks for libraries now
dnl Replace `main' with a function in -lm:
#  #np# These lines can be removed if the math lib is not needed.
AC_SEARCH_LIBS(pow,[m],,AC_MSG_ERROR(Cannot find math library))
AC_SEARCH_LIBS(sqrt,[m],,AC_MSG_ERROR(Cannot find math library))

# #np# If Lapack is also needed, see the  
# configure.ac file for epetra (Trilinos/packages/epetra)

dnl Check for user-specified names for BLAS lib
TAC_ARG_WITH_BLASLIB

  ACX_BLAS([], [AC_MSG_ERROR([[Cannot find blas library, specify a path using --with-ldflags=-I<DIR> (ex. --with-ldflags=-I/usr/path/lib) or a specific library using --with-blas=DIR/LIB (for example --with-blas=/usr/path/lib/libcxml.a)]])])

fi
# end of the list of libraries that don't need to be checked for if
# tests and examples are disabled.

# ------------------------------------------------------------------------
# Checks for compiler characteristics
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for linker characteristics
# ------------------------------------------------------------------------

# Determine libraries needed for linking with Fortran
AC_F77_LIBRARY_LDFLAGS


#  ----------------------------------------------------
#  Needs this because Dec handles templates differently
#  ----------------------------------------------------
#  #np# - Only uncomment the following lines if the package
#         creates a subdirectory called cxx_repository for
#         template object files on DEC.  (HINT: try to build on 
#         DEC/CPLANT and see if template-related errors occur)

# -- DEC CXX --
#if test "X$TRILINOS_TEMPLATE_OBJS" = "X"; then
#  case $target in
#    alpha*)
#      if test $CXX = cxx; then
#        TRILINOS_TEMPLATE_OBJS="cxx_repository/*.o"
#      fi
#    ;;
#  esac
#fi
                                                                  
#AC_MSG_CHECKING([for template-specific object files to be linked])if test "X$TRILINOS_TEMPLATE_OBJS" = "X"; then
# AC_MSG_RESULT(none)
#else
# AC_MSG_RESULT($TRILINOS_TEMPLATE_OBJS)
#fi
                                                                  
#AC_SUBST(TRILINOS_TEMPLATE_OBJS)

# ------------------------------------------------------------------------
# Checks for library functions
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for system services
# ------------------------------------------------------------------------


# ------------------------------------------------------------------------
# Perform substitutions in output files
# ------------------------------------------------------------------------

AC_SUBST(ac_aux_dir)

# ------------------------------------------------------------------------
# Output files
# ------------------------------------------------------------------------
##
#  #np# - Most packages should have an example and a full suite 
#  of tests.  New package comes only with src and example.
#  You will need to change AC_CONFIG_FILES below and Makefile.am
#  to add a new directory.
AC_CONFIG_FILES([
		Makefile 
    Makefile.export.rythmos
		src/Makefile
		example/Makefile
		example/basicExample/Makefile
		example/epetra/Makefile
		example/epetra/1Dfem/Makefile
		test/Makefile
		])
#    test/unitTests/Makefile
#    test/unitTests/InArgs/Makefile
#    test/unitTests/OutArgs/Makefile

AC_OUTPUT()

# Bye World!
#  #np# change package name
echo "----------------------------------------"
echo "Finished Running Rythmos Configure Script"
echo "----------------------------------------"
