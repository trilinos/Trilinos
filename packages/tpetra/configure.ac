# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# @HEADER
# ************************************************************************
# 
#          Tpetra: Templated Linear Algebra Services Package
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

AC_INIT(tpetra, 1.2, maherou@sandia.gov)

# Hello World!
echo "----------------------------------------"
echo "Running Tpetra Configure Script"
echo "----------------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.  Any file in that directory will do,
# preferably one that is unlikely to be removed or renamed.

AC_CONFIG_SRCDIR([src/Tpetra_Version.cpp])

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.

AC_CONFIG_AUX_DIR(config)

# Configure should create src/Tpetra_config.h from src/Tpetra_config.h.in

AM_CONFIG_HEADER(src/Tpetra_config.h:src/Tpetra_config.h.in)

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

# Specify required version of automake.
AM_INIT_AUTOMAKE(1.10 no-define tar-ustar)

# Specify required version of autoconf.
AC_PREREQ(2.61)

# ------------------------------------------------------------------------
# Check to see if MPI enabled and if any special configuration done
# ------------------------------------------------------------------------

TAC_ARG_CONFIG_MPI

# ------------------------------------------------------------------------
# Checks for programs
# ------------------------------------------------------------------------

AC_PROG_CC(cc gcc)
AC_PROG_CXX(CC g++ c++ cxx)
#AC_PROG_F77(f77 g77 gfortran f90 xlf90 f95)
AC_PROG_RANLIB

# ------------------------------------------------------------------------
# Check if --with-flags present, prepend any specs to FLAGS

TAC_ARG_WITH_FLAGS(ccflags, CCFLAGS)
TAC_ARG_WITH_FLAGS(cxxflags, CXXFLAGS)
TAC_ARG_WITH_FLAGS(cflags, CFLAGS)
#TAC_ARG_WITH_FLAGS(fflags, FFLAGS)
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

#NOTE: This section still needs to be updated for Tpetra, it is still 
# code from Epetra (via new_package)
# ------------------------------------------------------------------------
# Checks for header files - used in Tpetra_ConfigDefs.hpp
# ------------------------------------------------------------------------

AC_CHECK_HEADERS([cstdlib stdlib.h], [break],)
AC_CHECK_HEADERS([cstdio stdio.h], [break],)
AC_CHECK_HEADERS([cassert assert.h], [break],)
AC_CHECK_HEADERS([cmath math.h], [break],)
AC_CHECK_HEADERS([iomanip iomanip.h], [break],)
AC_CHECK_HEADERS([iostream iostream.h], [break],)
AC_CHECK_HEADERS([string string.h], [break],)
AC_CHECK_HEADERS([map map.h], [break],)
AC_CHECK_HEADERS([vector vector.h], [break],)
AC_CHECK_HEADERS([algorithm algo.h algorithm.h], [break],)
AC_CHECK_HEADERS([numeric algo.h], [break],)

if test "X$ac_cv_header_cstdlib" != "Xyes" && test "X$ac_cv_header_stdlib_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdlib or stdlib.h.  Tpetra requires at least one or the other.]])
fi

if test "X$ac_cv_header_cstdio" != "Xyes" && test "X$ac_cv_header_stdio_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdio or stdio.h.  Tpetra requires at least one or the other.]])
fi

if test "X$ac_cv_header_cassert" != "Xyes" && test "X$ac_cv_header_assert_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cassert or assert.h.  Tpetra requires at least one or the other.]])
fi

if test "X$ac_cv_header_cmath" != "Xyes" && test "X$ac_cv_header_math_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cmath or math.h.  Tpetra requires at least one or the other.]])
fi

if test "X$ac_cv_header_iomanip" != "Xyes" && test "X$ac_cv_header_iomanip_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find iomanip or iomanip.h.  Tpetra requires at least one or the other.]])
fi

if test "X$ac_cv_header_iostream" != "Xyes" && test "X$ac_cv_header_iostream_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find iostream or iostream.h.  Tpetra requires at least one or the other.]])
fi

if test "X$ac_cv_header_string" != "Xyes" && test "X$ac_cv_header_string_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find string or string.h.  Tpetra requires at least one or the other.]])
fi

if test "X$ac_cv_header_map" != "Xyes" && test "X$ac_cv_header_map_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find map or map.h.  Tpetra requires at least one or the other.]])
fi

if test "X$ac_cv_header_vector" != "Xyes" && test "X$ac_cv_header_vector_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find vector or vector.h.  Tpetra requires at least one or the other.]])
fi

if test "X$ac_cv_header_algorithm" != "Xyes" && test "X$ac_cv_header_algo_h" != "Xyes" && test "X$ac_cv_header_algorithm_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find algorithm, algo.h, or algorithm.h.  Tpetra requires at least one of these.]])
fi

if test "X$ac_cv_header_numeric" != "Xyes" && test "X$ac_cv_header_algo_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find numeric or algo.h.  Tpetra requires at least one or the other.]])
fi

# Checks for typedefs, structures, and compiler characteristics.
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
AC_CXX_STD_SPRINTF
AC_CXX_STD_NEW_COUNT_SYNTAX

# ------------------------------------------------------------------------
# Checks for Makefile.export related systems
# ------------------------------------------------------------------------
TAC_ARG_ENABLE_EXPORT_MAKEFILES(yes)

# ------------------------------------------------------------------------
# Checks for special package flags
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_FEATURE(tbb, [Enable use of Intel's Thread Building Blocks library], TPETRA_TBB, no)

# ------------------------------------------------------------------------
# Checks if tests should be built
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_FEATURE(tests, [Make tests for all Trilinos packages buildable with 'make tests'], TESTS, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK( tpetra, tests, [Make Tpetra tests buildable with 'make tests'], TPETRA_TESTS)
AM_CONDITIONAL(BUILD_TESTS, test "X$ac_cv_use_tpetra_tests" != "Xno")

# ------------------------------------------------------------------------
# Checks if examples should be built
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_FEATURE(examples, [Make examples for all Trilinos packages buildable with 'make examples'], EXAMPLES, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK( tpetra, examples, [Make Tpetra examples buildable with 'make examples'], TPETRA_EXAMPLES)
AM_CONDITIONAL(BUILD_EXAMPLES, test "X$ac_cv_use_tpetra_examples" != "Xno")

# Check to see if we should build the EpetraTpetraProfile example. This should only be done if
# both Epetra and TriUtils are enabled.
#TAC_ARG_ENABLE_CAN_USE_PACKAGE(tpetra, epetra, EPETRA, HAVE_EPETRA, yes, 
#[Enable examples that require Epetra.], 
#[Enable examples that require Epetra.  Can be overridden with --disable-tpetra-epetra.])
TAC_ARG_ENABLE_CAN_USE_PACKAGE(tpetra, triutils, TPETRA_TRIUTILS, ENABLE_TPETRA_TRIUTILS, yes, 
[Enable examples that require Triutils.],
[Enable examples that require Triutils.  Can be overridden with --disable-tpetra-triutils.])
AM_CONDITIONAL(BUILD_PROFILE_EXAMPLE, test "X$ac_cv_use_epetra" != "Xno" && test "X$ac_cv_use_triutils" != "Xno")

#We now build tests and examples through separate make targets, rather than
#during "make".  We still need to conditionally include the test and example
#in SUBDIRS, even though SUB_TEST and SUB_EXAMPLE will never be
#defined, so that the tests and examples are included in the distribution
#tarball.
AM_CONDITIONAL(SUB_TEST, test "X$ac_cv_use_sub_test" = "Xyes")
AM_CONDITIONAL(SUB_EXAMPLE, test "X$ac_cv_use_sub_example" = "Xyes")

TAC_ARG_ENABLE_FEATURE(libcheck, [Check for some third-party libraries.  (Cannot be disabled unless tests and examples are also disabled.)], LIBCHECK, yes)

# ------------------------------------------------------------------------
# Specify other directories 
# ------------------------------------------------------------------------

# enable use of --with-libdirs="-Llibdir1 -Llibdir2 ..." to prepend to LDFLAGS
TAC_ARG_WITH_LIBDIRS
# enable use of --with-incdirs="-Lincdir1 -Lincdir2 ..." to prepend to CPPFLAGS
TAC_ARG_WITH_INCDIRS

# ------------------------------------------------------------------------
# Checks for libraries
# ------------------------------------------------------------------------

# If tests, examples and libcheck are disabled, we don't have to check
# for these libraries.
if test "X$ac_cv_use_tpetra_tests" != "Xno" || test "X$ac_cv_use_tpetra_examples" != "Xno" || test "X$ac_cv_use_libcheck" != "Xno"; then

# To add any libraries that you want to check for, replace pow below 
# with the name of an entry point in that library.
# checks for libraries now
dnl Replace `main' with a function in -lm:
# These lines can be removed if the math lib is not needed.
AC_SEARCH_LIBS(pow,[m],,AC_MSG_ERROR(Cannot find math library))
AC_SEARCH_LIBS(sqrt,[m],,AC_MSG_ERROR(Cannot find math library))

# Define F77_FUNC that will be used to link with Fortran subroutines.
#AC_F77_WRAPPERS

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
#AC_F77_LIBRARY_LDFLAGS

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
#  You will need to change AC_CONFIG_FILES below and Makefile.am
#  to add a new directory.
AC_CONFIG_FILES([
        Makefile 
        Makefile.export.tpetra
        src/Makefile
        example/Makefile
        test/Makefile
        test/Platform/Makefile
        test/Distributor/Makefile
        test/Map/Makefile
        test/Directory/Makefile
        test/ImportExport/Makefile
        test/MultiVector/Makefile
        test/CrsMatrix/Makefile
        ])
# example/verySimple/Makefile
# example/EpetraTpetraProfile/Makefile
# example/TpetraOPs/Makefile
# test/BasicPerfTest/Makefile
# test/DistObject/Makefile
# test/CrsMatrix/Makefile
# test/Vector/Makefile

AC_OUTPUT()

# Bye World!
echo "----------------------------------------"
echo "Finished Running Tpetra Configure Script"
echo "----------------------------------------"
