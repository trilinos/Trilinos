# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# @HEADER
# ************************************************************************
# 
#        Isorropia: A Partitioning and Load Balancing Package
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
# Questions?  Contact Alan Williams (william@sandia.gov)
#                  or Erik Boman    (egboman@sandia.gov)
# 
# ************************************************************************
# @HEADER

# ------------------------------------------------------------------------
# Initialization 
# ------------------------------------------------------------------------

# This must be the first line in configure.ac.
# Optional 3rd argument is email address for bugs.

AC_INIT(isorropia, 1.0d, william@sandia.gov)

echo "----------------------------------------"
echo "Running Isorropia Configure Script"
echo "----------------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.

AC_CONFIG_SRCDIR(/test/utest/utest.cpp)

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.

AC_CONFIG_AUX_DIR(config)

# Configure should create src/Isorropia_autoheader.h from
# src/Isorropia_autoheader.h.in

AM_CONFIG_HEADER(src/Isorropia_autoheader.h:src/Isorropia_autoheader.h.in)

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

AM_INIT_AUTOMAKE(1.8.2)

# Specify required version of autoconf.

AC_PREREQ(2.59)

# Specify revision number of configure.ac that was used to generate
# configure.

dnl Commented out this line for now because it causes headaches in
dnl committing.
dnl AC_REVISION($Revision$)

# ------------------------------------------------------------------------
# Check to see if MPI enabled and if any special configuration done
# ------------------------------------------------------------------------

TAC_ARG_CONFIG_MPI

# ------------------------------------------------------------------------
# Checks for programs.
# ------------------------------------------------------------------------

AC_PROG_CC(cc gcc)
AC_PROG_CXX(CC g++ c++ cxx)
BASE_CXX=${CXX}
AC_PROG_F77(f77 g77 gfortran f90 xlf90 f95)
AC_PROG_RANLIB

# ------------------------------------------------------------------------
# Checks for --enable-debug, --enable-optimize, --enable-profile, --enable-purify
# Changes flags and compiler names accordingly
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_COMPILER_OPTIONS

# ------------------------------------------------------------------------
# Checks for extra flags
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

# ------------------------------------------------------------------------
# Checks for header files
# ------------------------------------------------------------------------

AC_HEADER_STDC
AC_CHECK_HEADERS([iostream fstream exception vector set time.h],
                 [],[], [#include <stdlib.h>],)

if test "X$ac_cv_header_iostream" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find iostream. Isorropia requires iostream.]])
fi

if test "X$ac_cv_header_exception" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find exception. Isorropia requires exception.]])
fi

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE

# ------------------------------------------------------------------------
# Checks for compiler characteristics
# ------------------------------------------------------------------------
AC_CXX_NAMESPACES
AC_CXX_HAVE_STL
AC_CXX_BOOL
AC_CXX_MUTABLE
AC_CXX_NEW_FOR_SCOPING
#AC_CXX_STD_IOS_BASE_FMTFLAGS
AC_CXX_STD_SPRINTF

# ------------------------------------------------------------------------
# Checks for special package flags
# ------------------------------------------------------------------------


# ------------------------------------------------------------------------
# Optional package dependency
# ------------------------------------------------------------------------
# The following section is required to decide if Epetra and/or Zoltan support
# should be built.

TAC_ARG_ENABLE_CAN_USE_PACKAGE(isorropia, epetra, EPETRA, HAVE_EPETRA, yes, [Build Epetra support in Isorropia.], [Build Epetra support in Isorropia.  Can be overridden with --disable-isorropia-epetra.])

TAC_ARG_ENABLE_FEATURE_SUB( isorropia, zoltan, [Enable Isorropia interface support for Zoltan], ISORROPIA_ZOLTAN, no)
AM_CONDITIONAL(BUILD_ZOLTAN, [test "X$ac_cv_use_isorropia_zoltan" != "Xno"])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(isorropia, epetraext, EPETRAEXT, HAVE_EPETRAEXT, yes, [Use EpetraExt in Isorropia.], [Use EpetraExt in Isorropia.  Can be overridden with --disable-isorropia-epetraext.])

# ------------------------------------------------------------------------
# Checks if tests and examples should be built
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_FEATURE(tests, [Build tests for all Trilinos packages (not all packages are sensitive to this option)], TESTS, yes)
TAC_ARG_ENABLE_FEATURE_SUB( isorropia, tests, [Build Isorropia tests], ISORROPIA_TESTS, yes)
AM_CONDITIONAL(BUILD_TESTS, test "X$ac_cv_use_isorropia_tests" != "Xno")

TAC_ARG_ENABLE_FEATURE(examples, [Build examples for all Trilinos packages (not all packages are sensitive to this option)], EXAMPLES, yes)
TAC_ARG_ENABLE_FEATURE_SUB( isorropia, examples, [Build Isorropia examples], ISORROPIA_EXAMPLES, yes)
AM_CONDITIONAL(BUILD_EXAMPLES, test "X$ac_cv_use_isorropia_examples" != "Xno")

TAC_ARG_ENABLE_FEATURE(libcheck, [Check for some third-party libraries including BLAS and LAPACK.  (Cannot be disabled unless tests and examples are also disabled.)], LIBCHECK, yes)

# ------------------------------------------------------------------------
# Enable optional python wrappers
# ------------------------------------------------------------------------
TAC_ENABLE_PYTHON

# ------------------------------------------------------------------------
# Interfaces to other Trilinos packages
# ------------------------------------------------------------------------

#TAC_ARG_ENABLE_CAN_USE_PACKAGE(
#  aztecoo, thyra, AZTECOO_THYRA, ENABLE_AZTECOO_THYRA, yes,
#  [Build the Thyra operator/vector adapater support in AztecOO (also rquires epetraext).],
#  [Build the Thyra operator/vector adapater support in AztecOO (also rquires epetraext).  Can be overridden with --disable-aztecoo-thyra.])

# ------------------------------------------------------------------------
# Specify other directories 
# ------------------------------------------------------------------------

# enable use of --with-libdirs="-Llibdir1 -Llibdir2 ..." to prepend to LDFLAGS
TAC_ARG_WITH_LIBDIRS
# enable use of --with-incdirs="-Lincdir1 -Lincdir2 ..." to prepend to CPPFLAGS
TAC_ARG_WITH_INCDIRS

# Define F77_FUNC that will be used to link with Fortran subroutines, e.g., BLAS/LAPACK
AC_F77_WRAPPERS

# If tests, examples and libcheck are disabled, we don't have to check
# for these libraries.
if test "X$ac_cv_use_isorropia_tests" != "Xno" || test "X$ac_cv_use_isorropia_examples" != "Xno" || test "X$ac_cv_use_libcheck" != "Xno"; then

# Check for user-specified names for LAPACK and BLAS
TAC_ARG_WITH_LAPACKLIB
TAC_ARG_WITH_BLASLIB

# find or choose to supply valid blas and lapack libraries
  ACX_BLAS([], [AC_MSG_ERROR([[Cannot find blas library, specify a path using --with-ldflags=-I<DIR> (ex. --with-ldflags=-I/usr/path/lib) or a specific library using --with-blas=DIR/LIB (for example --with-blas=/usr/path/lib/libcxml.a)]])])

  ACX_LAPACK([], [AC_MSG_ERROR([[Cannot find lapack library, specify a path using --with-ldflags=-I<DIR> (ex. --with-dlflags=-I/usr/path/lib) or a specific library using --with-lapack=DIR/LIB (ex. --with-lapack=/usr/path/lib/libcxml.a)]])])
fi
# end of the list of libraries that don't need to be checked for if
# tests and examples are disabled.

# Determine libraries needed for linking with Fortran
AC_F77_LIBRARY_LDFLAGS

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
# Perform substitutions in output files
# ------------------------------------------------------------------------

AC_SUBST(ac_aux_dir)

#-------------------------------------------------------------------------
# Add any flags that are required for external libraries that are to be
# built during the build process.
#-------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Output files
# ------------------------------------------------------------------------
#
AC_CONFIG_FILES([
    Makefile
    src/Makefile
    src/Makefile.include
    example/Makefile
    example/Makefile.include
    test/Makefile
    test/rebalance_1d_default/Makefile
    test/rebalance_1d_default/Makefile.include
    test/zoltan/Makefile
    test/zoltan/Makefile.include
    test/utest/Makefile
    test/utest/Makefile.include
    test/utils/Makefile
    test/utils/Makefile.include
    ])

AC_OUTPUT

echo "----------------------------------------"
echo "Finished Running Isorropia Configure Script"
echo "----------------------------------------"
