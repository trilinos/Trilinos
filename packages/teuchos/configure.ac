# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# @HEADER
# ************************************************************************
# 
#                    Teuchos: Common Tools Package
#                 Copyright (2004) Sandia Corporation
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
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************
# @HEADER

# ------------------------------------------------------------------------
# Initialization 
# ------------------------------------------------------------------------



# This must be the first line in configure.ac.
# Optional 3rd argument is email address for bugs.

#  #np# - package name, version number, and e-mail address below
AC_INIT(teuchos, 1.0, hkthorn@sandia.gov)

# Hello World!
#  #np# - package name below
echo "----------------------------------------"
echo "Running Teuchos Configure Script"
echo "----------------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.  Any file in that directory will do,
# preferably one that is unlikely to be removed or renamed.

# #np# - file name below
AC_CONFIG_SRCDIR(/src/Teuchos_LAPACK.hpp)

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
# Configure should create src/Teuchos_config.h from src/Teuchos_config.h.in

AM_CONFIG_HEADER(src/Teuchos_config.h:src/Teuchos_config.h.in)

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

# Specify required version of automake
AM_INIT_AUTOMAKE(1.8.2)

# Specify required version of autoconf.

AC_PREREQ(2.59)

#  #np# - can eliminate compiler checks below if your package does not use the
#         language corresponding to the check
# ------------------------------------------------------------------------
# Checks for programs
# ------------------------------------------------------------------------

AC_PROG_CC(cc gcc)
AC_PROG_CXX(CC g++ gcc c++ cxx)
AC_PROG_F77(f77 g77 f90 xlf90)
AC_PROG_RANLIB

# ------------------------------------------------------------------------
# Check to see if MPI enabled and if any special configuration done
# ------------------------------------------------------------------------

TAC_ARG_CONFIG_MPI

# ------------------------------------------------------------------------
# Checks for extra flags
# #np# - You can pull flag checks for languages not used in package 
# ------------------------------------------------------------------------
# Check if --with-flags present, prepend any specs to FLAGS

TAC_ARG_WITH_FLAGS(ccflags, CCFLAGS)
TAC_ARG_WITH_FLAGS(cxxflags, CXXFLAGS)
TAC_ARG_WITH_FLAGS(cflags, CFLAGS)
TAC_ARG_WITH_FLAGS(fflags, FFLAGS)
TAC_ARG_WITH_LIBS(libs, LIBS)
TAC_ARG_WITH_FLAGS(ldflags, LDFLAGS)

# ------------------------------------------------------------------------
# Alternate archiver
# ------------------------------------------------------------------------

TAC_ARG_WITH_AR

# ------------------------------------------------------------------------
# MPI link check
# ------------------------------------------------------------------------
TAC_ARG_CHECK_MPI

#  #np# - Customize the checks for headers needed in package being added.
#NOTE: This section still needs to be updated for new_package, it is still 
# code from Epetra
# ------------------------------------------------------------------------
# Checks for header files - used in Epetra_ConfigDefs.h 
# ------------------------------------------------------------------------

AC_CHECK_HEADERS([cstdlib stdlib.h cstdio stdio.h cassert assert.h],,,)
AC_CHECK_HEADERS([cmath math.h iomanip iomanip.h iostream iostream.h],,,)
AC_CHECK_HEADERS([string string.h complex complex.h algorithm algorithm.h],,,)
AC_CHECK_HEADERS([map map.h list list.h vector vector.h],,,)
AC_CHECK_HEADERS([sstream strstream sstream.h strstream.h],,,)
AC_CHECK_HEADERS([stdexcept stdexcept.h cstdarg stdarg.h],,,)
AC_CHECK_HEADERS([typeinfo],,,)

if test "X$ac_cv_header_cstdlib" != "Xyes" && test "X$ac_cv_header_stdlib_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdlib or stdlib.h.  Teuchos requires at least one or the other.]])
fi

if test "X$ac_cv_header_cstdio" != "Xyes" && test "X$ac_cv_header_stdio_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdio or stdio.h.  Teuchos requires at least one or the other.]])
fi

if test "X$ac_cv_header_cstdarg" != "Xyes" && test "X$ac_cv_header_stdarg_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdarg or stdarg.h.  Teuchos requires at least one or the other.]])
fi

if test "X$ac_cv_header_cassert" != "Xyes" && test "X$ac_cv_header_assert_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cassert or assert.h.  Teuchos requires at least one or the other.]])
fi

if test "X$ac_cv_header_iostream" != "Xyes" && test "X$ac_cv_header_iostream_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find iostream or iostream.h.  Teuchos requires at least one or the other.]])
fi

if test "X$ac_cv_header_string" != "Xyes" && test "X$ac_cv_header_string_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find string or string.h.  Teuchos requires at least one or the other.]])
fi

if test "X$ac_cv_header_complex" != "Xyes" && test "X$ac_cv_header_complex_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find complex or complex.h.  Teuchos requires at least one or the other.]])
fi

if test "X$ac_cv_header_algorithm" != "Xyes" && test "X$ac_cv_header_algorithm_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find algorithm or algorithm.h.  Teuchos requires at least one or the other.]])
fi

if test "X$ac_cv_header_map" != "Xyes" && test "X$ac_cv_header_map_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find map or map.h.  Teuchos requires at least one or the other.]])
fi

if test "X$ac_cv_header_list" != "Xyes" && test "X$ac_cv_header_list_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find list or list.h.  Teuchos requires at least one or the other.]])
fi

if test "X$ac_cv_header_vector" != "Xyes" && test "X$ac_cv_header_vector_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find vector or vector.h.  Teuchos requires at least one or the other.]])
fi

if test "X${host}" != "Xtflop"; then
  if test "X$ac_cv_header_cmath" != "Xyes" && test "X$ac_cv_header_math_h" != "Xyes"; then
    AC_MSG_ERROR([[Cannot find cmath or math.h.  Teuchos requires at least one or the other.]])
  fi
else
  # TFLOP
  if test "X$ac_cv_header_iomanip" != "Xyes" && test "X$ac_cv_header_iomanip_h" != "Xyes"; then
    AC_MSG_ERROR([[Cannot find iomanip or iomanip.h.  Teuchos requires at least one or the other.]])
  fi
fi  

if test "X$ac_cv_header_sstream" != "Xyes" && test "X$ac_cv_header_strstream" != "Xyes" && test "X$ac_cv_header_sstream" != "Xyes" && test "X$ac_cv_header_strstream" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find sstream, strstream, sstream.h, or strstream.h.  Teuchos requires at least one of these.]])
fi

if test "X$ac_cv_header_vector" != "Xyes" && test "X$ac_cv_header_vector_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find typeinfo.  Teuchos requires this header.]])
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
# ---------------------------------------------------------------------
# These are pretty obscure
# ---------------------------------------------------------------------
AC_CXX_NEW_FOR_SCOPING
AC_CXX_STD_IOS_BASE_FMTFLAGS
AC_CXX_TEMPLATE_QUALIFIER
AC_CXX_STD_SPRINTF

# ------------------------------------------------------------------------
# Checks if tests should be built
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_FEATURE(tests, [Build tests for all Trilinos packages (not all packages are sensitive to this option)], TESTS, yes)
AM_CONDITIONAL(BUILD_TESTS, test "X$ac_cv_use_tests" != "Xno")

# ------------------------------------------------------------------------
# Checks if examples should be built
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_FEATURE(examples, [Build examples for all Trilinos packages (not all packages are sensitive to this option)], EXAMPLES, yes)
AM_CONDITIONAL(BUILD_EXAMPLES, test "X$ac_cv_use_examples" != "Xno")

# ------------------------------------------------------------------------
# Specify other directories 
# ------------------------------------------------------------------------

# enable use of --with-libdirs="-Llibdir1 -Llibdir2 ..." to prepend to LDFLAGS
TAC_ARG_WITH_LIBDIRS
# enable use of --with-incdirs="-Lincdir1 -Lincdir2 ..." to prepend to CPPFLAGS
TAC_ARG_WITH_INCDIRS

# ------------------------------------------------------------------------
# Checks for special package flags
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# #np# newp_swahili is an example of a package feature which requires both
# a cpp #define command and one or more source or header files  - ie 
# conditional compilation which requires a change to the makefile.
# If configure is called with "--enable-newp_swahili", it will place 
# a #define HAVE_NEWP_SWAHILI in Teuchos_config.h and add source 
# and header files as specified in Makefile.am 
#
# To add an option to your package, you will need to modify this file
# (i.e. configure.ac), and src/Makefile.am.  A bootstrap is also required.
# Look for SWAHILI in both files listed above.
#
# TAC_ARG_ENABLE_OPTION performs several tasks, including:
#   1)  Setting a help string
#   2)  Printing a message to standard out
#   3)  Defining HAVE_NEWP_SWAHILI (For use in Teuchos_config.h)
#
# AM_CONDITIONAL makes HAVE_NEWP_SWAHILI available for use in Makefile.am
# ------------------------------------------------------------------------
#
#TAC_ARG_ENABLE_OPTION(newp_swahili, [enable swahili support], NEWP_SWAHILI, no)
#AM_CONDITIONAL(HAVE_NEWP_SWAHILI, [test "X$ac_cv_use_newp_swahili" != "Xno"])
#
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
# TAC_ARG_WITH_PACKAGE(esi, [Enable ESI interface support], ESI, yes)
# TAC_ARG_WITH_PACKAGE(superlu, [Enable SuperLU interface support], SUPERLU, no)

# ------------------------------------------------------------------------
# Check to see if array boundschecking has been enabled
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_FEATURE_SUB( teuchos, abc, [Enable array bounds checking], 
		TEUCHOS_ARRAY_BOUNDSCHECK, no)


# ------------------------------------------------------------------------
# Checks for addons 
# ------------------------------------------------------------------------

#
#  --enable-teuchos-complex
#

TAC_ARG_ENABLE_FEATURE_SUB( teuchos, complex, [Enable Teuchos complex code], TEUCHOS_COMPLEX, no)
AM_CONDITIONAL(HAVE_TEUCHOS_COMPLEX, [test "X$ac_cv_use_teuchos_complex" != "Xno"])

#
#  --enable-teuchos-extended
#

TAC_ARG_ENABLE_FEATURE_SUB( teuchos, extended, [Enable Teuchos extended code], TEUCHOS_EXTENDED, no)
AM_CONDITIONAL(HAVE_TEUCHOS_EXTENDED, [test "X$ac_cv_use_teuchos_extended" != "Xno"])

#
#  --enable-teuchos-arprec 
# 
 
TAC_ARG_ENABLE_FEATURE_SUB( teuchos, arprec, [Enable Teuchos interface support for the ARbitrary PRECision Library], 
		TEUCHOS_ARPREC, no)  

if test "X$ac_cv_use_teuchos_arprec" != "Xno"; then

#  TAC_ARG_WITH( trilinos3pldir, [Set Trilinos3PL directory], trilinos3pldir_is_not_specified )
#  if test "X${tac_with_teuchos_trilinos3pldir}" != "Xtrilinos3pldir_is_not_specified"; then

     AC_MSG_CHECKING(for mp/mpreal.h)
     AC_TRY_CPP([#include "mp/mpreal.h"],
	[AC_MSG_RESULT(yes)],  
	[
	 AC_MSG_RESULT(no)
	 echo "-----"
     	 echo "Cannot find ARPREC library."
     	 echo "Use --with-ldflags to specify ARPREC library directory."
     	 echo "Use --with-incdirs to specify ARPREC include directory."
     	 echo "Use --with-libs to specify ARPREC library (-lmp)."
     	 echo "-----"
     	 AC_MSG_ERROR(ARPREC cannot link)
     	])
fi

# ------------------------------------------------------------------------
# Checks for libraries
# ------------------------------------------------------------------------

# #np# - 
# To add any libraries that you want to check for, replace pow below 
# with the name of an entry point in that library.
# checks for libraries now
dnl Replace `main' with a function in -lm:
#  #np# These lines can be removed if the math lib is not needed.
AC_SEARCH_LIBS(pow,[m],,AC_MSG_ERROR(Cannot find math library))
AC_CHECK_FUNCS([pow sqrt])

# #np# - Yet another opportunity to remove code if you aren't
# using Fortran
# Define F77_FUNC that will be used to link with Fortran subroutines, e.g., BLAS/LAPACK - trash WORKGXX 
AC_F77_WRAPPERS

### Include Blas/Lapack example here. ###  For now, developers should see
### The configure.ac file for epetra (Trilinos/packages/epetra) if blas and
### or Lapack are needed.

ACX_BLAS([], [AC_MSG_ERROR([[Cannot find blas library, specify a path using --with-blas=DIR/LIB (for example --with-blas=/usr/path/lib/libcxml.a)]])])
ACX_LAPACK([], [AC_MSG_ERROR([[Cannot find lapack library, specify a path using --with-lapack=DIR/LIB (ex. --with-lapack=/usr/path/lib/libcxml.a)]])])


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
if test "X$TRILINOS_TEMPLATE_OBJS" = "X"; then
  case $target in
    alpha*)
      if test $CXX = cxx; then
        TRILINOS_TEMPLATE_OBJS="cxx_repository/*.o"
      fi
    ;;
  esac
fi

AC_MSG_CHECKING([for template-specific object files to be linked])
if test "X$TRILINOS_TEMPLATE_OBJS" = "X"; then
 AC_MSG_RESULT(none)
else
 AC_MSG_RESULT($TRILINOS_TEMPLATE_OBJS)
fi

AC_SUBST(TRILINOS_TEMPLATE_OBJS)

# ------------------------------------------------------------------------
# Checks for library functions
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for system services
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Check to see if expat enabled and if any special configuration done
# ------------------------------------------------------------------------

TAC_ARG_CONFIG_EXPAT

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
		src/Makefile
		test/Makefile
		test/BLAS/Makefile
		test/LAPACK/Makefile
		test/DenseMatrix/Makefile
		test/dyn_cast/Makefile
		test/ParameterList/Makefile
		test/Containers/Makefile
		test/Time/Makefile
		test/RefCountPtr/Makefile
		test/XML/Makefile
		test/scripts/Makefile
		test/scripts/daily/Makefile
		test/scripts/daily/serial/Makefile
		example/Makefile
		example/BLAS/Makefile
		example/CommandLineProcessor/Makefile
		example/DenseMatrix/Makefile
		example/hilbert/Makefile
		example/LAPACK/Makefile
		example/ParameterList/Makefile
		example/RefCountPtr/Makefile
		example/TimeMonitor/Makefile
		])

AC_OUTPUT()

# Bye World!
#  #np# change package name
echo "----------------------------------------"
echo "Finished Running Teuchos Configure Script"
echo "----------------------------------------"
