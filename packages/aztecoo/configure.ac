# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# @HEADER
# ************************************************************************
# 
#        AztecOO: An Object-Oriented Aztec Linear Solver Package 
#                 Copyright (2002) Sandia Corporation
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

AC_INIT(aztecoo, 3.3d, maherou@sandia.gov)

# Hello World!
echo "----------------------------------------"
echo "Running AztecOO Configure Script"
echo "----------------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.

AC_CONFIG_SRCDIR(/src/az_ilut.c)

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.

AC_CONFIG_AUX_DIR(config)

# Configure should create src/AztecOO_config.h from src/AztecOO_config.h.in

AM_CONFIG_HEADER(src/AztecOO_config.h:src/AztecOO_config.h.in)

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
AC_PROG_CXX(CC g++ gcc c++ cxx)
BASE_CXX=${CXX}
AC_PROG_F77(f77 g77 f90 xlf90)
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
# Checks for header files
# ------------------------------------------------------------------------

AC_HEADER_STDC
AC_CHECK_HEADERS([cstdlib stdlib.h], [break],)
AC_CHECK_HEADERS([cstdio stdio.h], [break],)
AC_CHECK_HEADERS([cassert assert.h], [break],)
AC_CHECK_HEADERS([cmath math.h], [break],)
AC_CHECK_HEADERS([string string.h], [break],)
AC_CHECK_HEADERS([cfloat float.h], [break],)
AC_CHECK_HEADERS([ctime sys/time.h], [break],)
AC_CHECK_HEADERS([malloc.h],,)

if test "X$ac_cv_header_cstdlib" != "Xyes" && test "X$ac_cv_header_stdlib_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdlib or stdlib.h.  AztecOO requires at least one or the other.]])
fi

if test "X$ac_cv_header_cstdio" != "Xyes" && test "X$ac_cv_header_stdio_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdio or stdio.h.  AztecOO requires at least one or the other.]])
fi

if test "X$ac_cv_header_cassert" != "Xyes" && test "X$ac_cv_header_assert_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cassert or assert.h.  AztecOO requires at least one or the other.]])
fi

if test "X$ac_cv_header_string" != "Xyes" && test "X$ac_cv_header_string_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find string or string.h.  AztecOO requires at least one or the other.]])
fi

if test "X$ac_cv_header_cmath" != "Xyes" && test "X$ac_cv_header_math_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cmath or math.h.  AztecOO requires at least one or the other.]])
fi

if test "X$ac_cv_header_cfloat" != "Xyes" && test "X$ac_cv_header_float_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cfloat or float.h.  AztecOO requires at least one or the other.]])
fi

if test "X$ac_cv_header_ctime" != "Xyes" && test "X$ac_cv_header_sys_time_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find ctime or sys/time.h.  AztecOO requires at least one or the other.]])
fi

if test "X$ac_cv_header_malloc_h" != "Xyes"; then
 AC_FUNC_MALLOC
 #if test "X$ac_cv_func_malloc_works" != "Xyes"; then
 #  AC_MSG_ERROR([malloc(...) not working, required by AztecOO])
 #fi
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

#This option is obsolete
#TAC_ARG_ENABLE_FEATURE_SUB( aztecoo, ifpack, [enable AztecOO's Ifpack dependent code], AZTECOO_IFPACK, no)

# The following section is required to decide if the AZ_lu preconditioner
# should be built.  If not, y12m is not necessary.
AC_ARG_ENABLE(aztecoo-azlu,
 [AC_HELP_STRING([--enable-aztecoo-azlu],[Enable az-lu preconditioner.  Default is no.  Requires y12m.])],
 ac_cv_use_aztecoo_azlu=$enableval, ac_cv_use_aztecoo_azlu=no)

AC_MSG_CHECKING(whether to use aztecoo-azlu)

if test "X$ac_cv_use_aztecoo_azlu" != "Xno"; then
  AC_MSG_RESULT(yes)
  AC_DEFINE([HAVE_AZLU],,[Define if you want to build AZ_lu])
else
  AC_MSG_RESULT(no)
fi

# ------------------------------------------------------------------------
# Checks if tests and examples should be built
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_FEATURE(tests, [Build tests for all Trilinos packages (not all packages are sensitive to this option)], TESTS, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK( aztecoo, tests, [Build AztecOO tests], AZTECOO_TESTS)
AM_CONDITIONAL(BUILD_TESTS, test "X$ac_cv_use_aztecoo_tests" != "Xno")

TAC_ARG_ENABLE_FEATURE(examples, [Build examples for all Trilinos packages (not all packages are sensitive to this option)], EXAMPLES, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK( aztecoo, examples, [Build AztecOO examples], AZTECOO_EXAMPLES)
AM_CONDITIONAL(BUILD_EXAMPLES, test "X$ac_cv_use_aztecoo_examples" != "Xno")

TAC_ARG_ENABLE_FEATURE(libcheck, [Check for some third-party libraries including BLAS and LAPACK.  (Cannot be disabled unless tests and examples are also disabled.)], LIBCHECK, yes)

# ------------------------------------------------------------------------
# Specify other directories 
# ------------------------------------------------------------------------

# enable use of --with-libdirs="-Llibdir1 -Llibdir2 ..." to prepend to LDFLAGS
TAC_ARG_WITH_LIBDIRS
# enable use of --with-incdirs="-Lincdir1 -Lincdir2 ..." to prepend to CPPFLAGS
TAC_ARG_WITH_INCDIRS

# Define F77_FUNC that will be used to link with Fortran subroutines, e.g., BLAS/LAPACK
AC_F77_WRAPPERS

# For now we will consider AztecOO to be dependent on Triutils,
# even though it is really only the examples.
#TAC_ARG_ENABLE_FEATURE(triutils, [Build triutils], TRIUTILS, no)
#AM_CONDITIONAL(BUILD_TRIUTILS, test "X$ac_cv_use_triutils" != "Xno")

# If tests, examples and libcheck are disabled, we don't have to check
# for these libraries.
if test "X$ac_cv_use_aztecoo_tests" != "Xno" || test "X$ac_cv_use_aztecoo_examples" != "Xno" || test "X$ac_cv_use_libcheck" != "Xno"; then

# y12m is needed only when AZ_lu is enabled
if test "X$ac_cv_use_aztecoo_azlu" != "Xno"; then
# y12m library needed
    # Find a valid y12m
    AC_LANG_SAVE
    AC_LANG_FORTRAN77
    AC_SEARCH_LIBS(y12mae,[y12m],,
    AC_MSG_ERROR([[Cannot find y12m library, specify a location using --with-ldflags-L<DIR> (ex. --with-ldflags=-L/usr/path/lib).]]))
    AC_LANG_RESTORE
fi
# Checks for library functions.
AC_SEARCH_LIBS(pow,[m],,AC_MSG_ERROR(Cannot find math library))
#Removed the following line because it caused problems on janus.
#Also, other libraries check only one function, a standard math
#lib should have sqrt provided that it has pow.
#AC_SEARCH_LIBS(sqrt,[m],,AC_MSG_ERROR(Cannot find math library))

AC_CHECK_FUNCS([floor pow sqrt strchr uname])

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
# should have this dependent on if AztecOO is built, or 
# just core Aztec - makefiles in example/
# even more critical

#
#                example/IfpackAztecOO/Makefile
#                example/azoo_iterate_hb/Makefile
#                example/aztec_app/Makefile
#                example/aztec_hb/Makefile
#
AC_CONFIG_FILES([
		Makefile
                src/Makefile
		example/Makefile
                example/UG_Ex1/Makefile
                example/simple_epetra_aztecoo/Makefile
                example/AztecOO_UserOpUserMat/Makefile
                example/AztecOO/Makefile
		test/Makefile
		test/AztecOO_SetParameters/Makefile
		test/AztecOO_UserOpUserMat/Makefile
		test/scripts/Makefile
		test/scripts/daily/Makefile
		test/scripts/daily/serial/Makefile])
AC_OUTPUT

# Bye World!
echo "----------------------------------------"
echo "Finished Running AztecOO Configure Script"
echo "----------------------------------------"
