# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# @HEADER
# ************************************************************************
# 
#       Ifpack: Object-Oriented Algebraic Preconditioner Package
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

AC_INIT(ifpack, 2.2d, maherou@sandia.gov)

# Hello World!
echo "----------------------------------------"
echo "Running Ifpack Configure Script"
echo "----------------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.

AC_CONFIG_SRCDIR(/src/Ifpack_IlukGraph.cpp)

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.

AC_CONFIG_AUX_DIR(config)

# Configure should create src/Ifpack_config.h from src/Ifpack_config.h.in

AM_CONFIG_HEADER(src/Ifpack_config.h:src/Ifpack_config.h.in)

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
# Checks for programs
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
# Special check for case of running on TFLOP
# ------------------------------------------------------------------------
AC_MSG_CHECKING([if host is TFLOP (and therefore needs special TFLOP define statement)])
if test "X${host}" = "Xtflop"; then
  AC_MSG_RESULT(yes)
  AC_DEFINE(TFLOP,,[define if compiling for ASCI Red, aka TFLOP])
else
  AC_MSG_RESULT(no)
fi

# ------------------------------------------------------------------------
# Checks for header files
# ------------------------------------------------------------------------

AC_CHECK_HEADERS([cstdlib stdlib.h], [break],)
AC_CHECK_HEADERS([cstdio stdio.h], [break],)
AC_CHECK_HEADERS([cassert assert.h], [break],)
AC_CHECK_HEADERS([cmath math.h], [break],)
AC_CHECK_HEADERS([iostream iostream.h], [break],)
AC_CHECK_HEADERS([string string.h], [break],)

if test "X$ac_cv_header_cstdlib" != "Xyes" && test "X$ac_cv_header_stdlib_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdlib or stdlib.h.  Ifpack requires at least one or the other.]])
fi

if test "X$ac_cv_header_cstdio" != "Xyes" && test "X$ac_cv_header_stdio_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdio or stdio.h.  Ifpack requires at least one or the other.]])
fi

if test "X$ac_cv_header_cassert" != "Xyes" && test "X$ac_cv_header_assert_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cassert or assert.h.  Ifpack requires at least one or the other.]])
fi

if test "X$ac_cv_header_string" != "Xyes" && test "X$ac_cv_header_string_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find string or string.h.  Ifpack requires at least one or the other.]])
fi

if test "X$ac_cv_header_iostream" != "Xyes" && test "X$ac_cv_header_iostream_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find iostream or iostream.h.  Ifpack requires at least one or the other.]])
fi

if test "X$ac_cv_header_cmath" != "Xyes" && test "X$ac_cv_header_math_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cmath or math.h.  Ifpack requires at least one or the other.]])
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
# Checks for other Trilinos packages and external libraries
# ------------------------------------------------------------------------
# Trilinos packages currently tested:
# - Epetra is required. The configure stops if Epetra is not enabled.
# - AztecOO is used in some examples.
# - Amesos is used by Ifpack_Amesos
# - Triutils is used in some examples, to generate the linear system.
# - Teuchos is used for the parameters list in several classes.
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_FEATURE_SUB( default, packages, [Build all default packages unless explicitly disabled.  When this option is disabled, only packages explicitly enabled and packages needed to satisfy the dependencies of enabled packages will be built], DEFAULT_PACKAGES, yes)

TAC_ARG_ENABLE_DEFAULT_FEATURE(epetra, [Enable EPETRA support. IFPACK depens on Epetra. Hence, if --disable-default-packages is used, the user need to specify --enable-epetra to compile IFPACK.], IFPACK_EPETRA, yes)
AM_CONDITIONAL(HAVE_IFPACK_AZTECOO, test "X$ac_cv_use_epetra" = "Xyes")
if test "X$ac_cv_use_epetra" = "Xno"; then
  echo "IFPACK depends on EPETRA. If you disable the default packages"
  echo "(--disable-default-packages), then you must enable epetra,"
  echo "using --enable-epetra."
  exit 1
fi
#

TAC_ARG_ENABLE_DEFAULT_FEATURE(aztecoo, [Enable AztecOO support. AztecOO, a Trilinos package, is used in some of the examples for solve the linear system.], IFPACK_AZTECOO, yes)
AM_CONDITIONAL(HAVE_IFPACK_AZTECOO, test "X$ac_cv_use_aztecoo" = "Xyes")
#
TAC_ARG_ENABLE_OPTION(triutils, [Enable TRIUTILS support. Triutils, a Trilinos package, is used in some IFPACK examples to construct the linear system.], IFPACK_TRIUTILS, no)
AM_CONDITIONAL(HAVE_IFPACK_TRIUTILS, test "X$ac_cv_use_triutils" = "Xyes")
#
TAC_ARG_ENABLE_OPTION(teuchos, [Enable TEUCHOS support. The parameters list of Teuchos, a Trilinos package, can be used in several IFPACK classes.], IFPACK_TEUCHOS, no)
AM_CONDITIONAL(HAVE_IFPACK_TEUCHOS, test "X$ac_cv_use_teuchos" = "Xyes")
#
TAC_ARG_ENABLE_OPTION(amesos, [Enable AMESOS support. Amesos, a Trilinos package, can be adopted as local solver for additive Schwarz preconditioners.  ], IFPACK_AMESOS, no)
AM_CONDITIONAL(HAVE_IFPACK_AMESOS, test "X$ac_cv_use_amesos" = "Xyes")

# ------------------------------------------------------------------------
# Checks for special package flags
# ------------------------------------------------------------------------
# Packages currently available:
# - METIS (used by Ifpack_METISPartitioner only)
# ------------------------------------------------------------------------

AC_ARG_ENABLE(ifpack-metis, [AC_HELP_STRING([--enable-ifpack-metis],[Enable
			     Ifpack support for METIS.  METIS can be used to
			     decompose the local graph into subgraphs. Note 1:
			     METIS is a third party library.  As such,
			     configure assumes that it has already been built.
			     Note 2: You must also specify the location of
			     headers and libraries like, for instance:
			     configure --with-libs="-lmetis"
			     --with-ldflags="-L<METIS-library-path>" and
			     --with-incdirs="<METIS-include-path>".])], ac_cv_use_ifpack_metis=$enableval, ac_cv_use_ifpack_metis=no) 
#
AC_MSG_CHECKING(whether to use ifpack-metis)
#
if test "X$ac_cv_use_ifpack_metis" != "Xno"; then
  AC_MSG_RESULT(yes)
  AC_DEFINE([HAVE_IFPACK_METIS],,[Define if you want to build Ifpack's ParameterList interface])
else
 AC_MSG_RESULT(no)
fi
#
AM_CONDITIONAL(HAVE_IFPACK_METIS, [test "X$ac_cv_use_ifpack_metis" != "Xno"])

### fill in ifdefs here - EPETRA_MPI, DEBUG, CRAY, COMPLEX_SUPPORT

# The following section is required to decide if the Ifpack ParameterList
# interface should be built.
#AC_ARG_ENABLE(ifpack-teuchos,
# [AC_HELP_STRING([--enable-ifpack-teuchos],[Enable Ifpack Teuchos::ParameterList interface. Default is no. Requires Teuchos package.])],
# ac_cv_use_ifpack_teuchos=$enableval, ac_cv_use_ifpack_teuchos=no)
#
#AC_MSG_CHECKING(whether to use ifpack-teuchos)
#
#if test "X$ac_cv_use_ifpack_teuchos" != "Xno"; then
#  AC_MSG_RESULT(yes)
#  AC_DEFINE([HAVE_IFPACK_TEUCHOS],,[Define if you want to build Ifpack's ParameterList interface])
#else
# AC_MSG_RESULT(no)
#fi
#
#AM_CONDITIONAL(HAVE_IFPACK_TEUCHOS, [test "X$ac_cv_use_ifpack_teuchos" != "Xno"])

# ------------------------------------------------------------------------
# Checks if tests should be built
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_FEATURE(tests, [Build tests for all Trilinos packages (not all packages are sensitive to this option)], TESTS, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK( ifpack, tests, [Build Ifpack tests], IFPACK_TESTS)
AM_CONDITIONAL(BUILD_TESTS, test "X$ac_cv_use_ifpack_tests" != "Xno")

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

# ------------------------------------------------------------------------
# Checks for libraries
# ------------------------------------------------------------------------

# If tests, examples and libcheck are disabled, we don't have to check
# for these libraries.
#if test "X$ac_cv_use_ifpack_tests" != "Xno" || test "X$ac_cv_use_ifpack_examples" != "Xno" || test "X$ac_cv_use_libcheck" != "Xno"; then
if test "X$ac_cv_use_ifpack_tests" != "Xno" || test "X$ac_cv_use_libcheck" != "Xno"; then

# checks for libraries now
dnl Replace `main' with a function in -lm:
AC_SEARCH_LIBS(sqrt,[m],,AC_MSG_ERROR(Cannot find math library))

dnl Check for user-specified names for LAPACK and BLAS
TAC_ARG_WITH_BLASLIB
TAC_ARG_WITH_LAPACKLIB

# find or choose to supply valid blas and lapack libraries
#TAC_ARG_ENABLE_FEATURE(referenceblasandlapack, [Build blas and lapack libs included in the distribution], BLASANDLAPACK, no)
#AM_CONDITIONAL(BUILD_BLAS, test "X$ac_cv_use_referenceblasandlapack" != "Xno")
#if test "X$ac_cv_use_referenceblasandlapack" = "Xno"; then
#  ACX_BLAS([], [AC_MSG_ERROR([[Cannot find blas library, please specify additional locations using --with-ldflags, additional lib names using --with-blas or type --with-referenceblasandlapack on the configure line in the Trilinos level directory to use the blas and lapack libraries included in the distribution]])])
  ACX_BLAS([], [AC_MSG_ERROR([[Cannot find blas library, specify a path using --with-ldflags=-I<DIR> (ex. --with-ldflags=-I/usr/path/lib) or a specific library using --with-blas=DIR/LIB (for example --with-blas=/usr/path/lib/libcxml.a)]])])
#fi
#TAC_ARG_ENABLE_FEATURE(referencelapack, [Build lapack lib included in the distribution], LAPACK, no)
#AM_CONDITIONAL(BUILD_LAPACK, test "X$ac_cv_use_referencelapack" != "Xno" || test "X$ac_cv_use_referenceblasandlapack" != "Xno")
#if test "X$ac_cv_use_referencelapack" = "Xno" && test "X$ac_cv_use_referenceblasandlapack" = "Xno"; then
#  ACX_LAPACK([], [AC_MSG_ERROR([[Cannot find lapack library, please specify additional locations using --with-ldflags, additional lib names using --with-lapacklib or type --with-referencelapack on the configure line in the Trilinos level directory to use the lapack library included in the distribution]])])
  ACX_LAPACK([], [AC_MSG_ERROR([[Cannot find lapack library, specify a path using --with-ldflags=-I<DIR> (ex. --with-dlflags=-I/usr/path/lib) or a specific library using --with-lapack=DIR/LIB (ex. --with-lapack=/usr/path/lib/libcxml.a)]])])

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
#  Only uncomment the following lines if the package
#  creates a subdirectory called cxx_repository for
#  template object files on DEC.  (HINT: try to build on
#  DEC/CPLANT and see if template-related errors occur)

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
# Perform substitutions in output files
# ------------------------------------------------------------------------

AC_SUBST(ac_aux_dir)

# ------------------------------------------------------------------------
# Output files
# ------------------------------------------------------------------------
##
#		example/Makefile 
#		example/ifpack_hb/Makefile

AC_CONFIG_FILES([
		Makefile 
                Makefile.export
		src/Makefile
		src/Makefile.include
		test/Makefile
		test/CrsRiluk/Makefile
		test/CrsRiluk/Makefile.include
		test/IlukGraph/Makefile
		test/IlukGraph/Makefile.include
		test/SetParameters/Makefile
		test/SetParameters/Makefile.include
		test/PointPreconditioner/Makefile
		test/PointPreconditioner/Makefile.include
		test/AdditiveSchwarz/Makefile
		test/AdditiveSchwarz/Makefile.include
		test/CrsAdditiveSchwarz/Makefile
		test/CrsAdditiveSchwarz/Makefile.include
		test/Container/Makefile
		test/Container/Makefile.include
		test/scripts/Makefile
          test/scripts/daily/Makefile
          test/scripts/daily/mpi/Makefile
          test/scripts/daily/serial/Makefile
          test/scripts/weekly/Makefile
          test/scripts/weekly/mpi/Makefile
          test/scripts/weekly/serial/Makefile
		])

AC_OUTPUT()

# Bye World!
echo "----------------------------------------"
echo "Finished Running Ifpack Configure Script"
echo "----------------------------------------"
