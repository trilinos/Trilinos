# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# @HEADER
# ************************************************************************
# 
#     EpetraExt: Epetra Extended - Linear Algebra Services Package 
#                 Copyright (2001) Sandia Corporation
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

AC_INIT(epetraext, 3.1d, maherou@sandia.gov)

# Hello World!
echo "----------------------------------------"
echo "Running EpetraExt Configure Script"
echo "----------------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.

AC_CONFIG_SRCDIR(/src/transform/EpetraExt_View_CrsGraph.cpp)

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.

AC_CONFIG_AUX_DIR(config)

# Configure should create src/EpetraExt_config.h from src/EpetraExt_config.h.in

AM_CONFIG_HEADER(src/EpetraExt_config.h:src/EpetraExt_config.h.in)

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
# Checks for programs
# ------------------------------------------------------------------------

AC_PROG_CC(cc gcc)
AC_PROG_CXX(CC g++ gcc c++ cxx)
BASE_CXX=${CXX}
AC_PROG_F77(f77 g77 f90 xlf90)
AC_PROG_RANLIB

# ------------------------------------------------------------------------
# Check to see if MPI enabled and if any special configuration done
# ------------------------------------------------------------------------

TAC_ARG_CONFIG_MPI

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
# Build Tests?
# ------------------------------------------------------------------------

#AC_DEFUN([EPETRAEXT_SUBPKG],
#[
#  AC_ARG_ENABLE(epetraext-$1,
#  [AC_HELP_STRING([--enable-epetraext-$1],[compile $1 sub-package])],
#  [SUBPKG_$1=$enableval],
#  [SUBPKG_$1=$2]
#  )
#
#  AC_ARG_ENABLE(epetraext-$1-tests,
#  [AC_HELP_STRING([--enable-epetraext-$1-tests],[compile $1 tests])],
#  [SUBPKG_$1_TESTS=$enableval],
#  [SUBPKG_$1_TESTS=no]
#  )
#
#  AC_ARG_ENABLE(epetraext-$1-examples,
#  [AC_HELP_STRING([--enable-epetraext-$1-examples],[compile $1 examples])],
#  [SUBPKG_$1_EXAMPLES=$enableval],
#  [SUBPKG_$1_EXAMPLES=no]
#  )
#
#  #Cannot build tests if associated sub-package not built
#  if test ${SUBPKG_$1_TESTS} = yes; then
#    SUBPKG_$1=yes;
#  fi
#
#  #Cannot build examples if associated sub-package not built
#  if test ${SUBPKG_$1_EXAMPLES} = yes; then
#    SUBPKG_$1=yes;
#  fi
#
#  AC_MSG_CHECKING(whether the sub-package $1 should be built)
#  AC_MSG_RESULT(${SUBPKG_$1})
#
#  AC_MSG_CHECKING(whether the $1 tests should be built)
#  AC_MSG_RESULT(${SUBPKG_$1_TESTS})
#
#  AC_MSG_CHECKING(whether the $1 examples should be built)
#  AC_MSG_RESULT(${SUBPKG_$1_EXAMPLES})
#])
#
#EPETRAEXT_SUBPKG(transform,yes)
#AM_CONDITIONAL(BUILD_TRANSFORM, test X${SUBPKG_transform} = Xyes)
#AM_CONDITIONAL(BUILD_TRANSFORM_TESTS, test X${SUBPKG_transform_TESTS} = Xyes)
#AM_CONDITIONAL(BUILD_TRANSFORM_EXAMPLES, test X${SUBPKG_transform_EXAMPLES} = Xyes)
#
#EPETRAEXT_SUBPKG(inout,yes)
#AM_CONDITIONAL(BUILD_INOUT, test X${SUBPKG_inout} = Xyes)
#AM_CONDITIONAL(BUILD_INOUT_TESTS, test X${SUBPKG_inout_TESTS} = Xyes)
#AM_CONDITIONAL(BUILD_INOUT_EXAMPLES, test X${SUBPKG_inout_EXAMPLES} = Xyes)

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

AC_CHECK_HEADERS([cstdlib stdlib.h cstdio stdio.h cassert assert.h],,,)
AC_CHECK_HEADERS([cmath math.h iomanip iomanip.h iostream iostream.h],,,)
AC_CHECK_HEADERS([string string.h algorithm algo.h algorithm.h],,,)

if test "X$ac_cv_header_cstdlib" != "Xyes" && test "X$ac_cv_header_stdlib_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdlib or stdlib.h.  EpetraExt requires at least one or the other.]])
fi

if test "X$ac_cv_header_cstdio" != "Xyes" && test "X$ac_cv_header_stdio_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdio or stdio.h.  EpetraExt requires at least one or the other.]])
fi

if test "X$ac_cv_header_cassert" != "Xyes" && test "X$ac_cv_header_assert_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cassert or assert.h.  EpetraExt requires at least one or the other.]])
fi

if test "X$ac_cv_header_string" != "Xyes" && test "X$ac_cv_header_string_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find string or string.h.  EpetraExt requires at least one or the other.]])
fi

if test "X$ac_cv_header_iostream" != "Xyes" && test "X$ac_cv_header_iostream_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find iostream or iostream.h.  EpetraExt requires at least one or the other.]])
fi

if test "X$ac_cv_header_algorithm" != "Xyes"  && test "X$ac_cv_header_algo_h" != "Xyes" && test "X$ac_cv_header_algorithm_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find algorithm, algo.h, or algorithm.h.  EpetraExt requires at least one of these.]])
fi

if test "X${host}" != "Xtflop"; then
  if test "X$ac_cv_header_cmath" != "Xyes" && test "X$ac_cv_header_math_h" != "Xyes"; then
    AC_MSG_ERROR([[Cannot find cmath or math.h.  EpetraExt requires at least one or the other.]])
  fi
else
  # TFLOP
  if test "X$ac_cv_header_iomanip" != "Xyes" && test "X$ac_cv_header_iomanip_h" != "Xyes"; then
    AC_MSG_ERROR([[Cannot find iomanip or iomanip.h.  EpetraExt requires at least one or the other.]])
  fi
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
AC_CXX_STL_MAP_KEY_CONST

# ------------------------------------------------------------------------
# Checks for special package flags
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_OPTION(threads, [enable shared memory threads], THREADS, no)
AM_CONDITIONAL(HAVE_THREADS, [test "X$ac_cv_use_threads" != "Xno"])

TAC_ARG_ENABLE_OPTION(fatal_messages, [enable printing of fatal error messages], 
                   FATAL_MESSAGES, yes)
TAC_ARG_ENABLE_OPTION(warning_messages, [enable printing of warning error messages], 
                   WARNING_MESSAGES, no)
TAC_ARG_ENABLE_OPTION(format_io, [Use package-defined IO formatting for output], 
                   FORMAT_IO, yes)

TAC_ARG_ENABLE_FEATURE_SUB( epetraext, transform, [Enable EpetraExt Transform subpkg], TRANSFORM, yes)
AM_CONDITIONAL(BUILD_TRANSFORM, [test "X$ac_cv_use_epetraext_transform" != "Xno"])

TAC_ARG_ENABLE_FEATURE_SUB( epetraext, inout, [Enable EpetraExt InOut subpkg], INOUT, yes)
AM_CONDITIONAL(BUILD_INOUT, [test "X$ac_cv_use_epetraext_inout" != "Xno"])

TAC_ARG_ENABLE_FEATURE_SUB( epetraext, coloring, [Enable EpetraExt Coloring subpkg], INOUT, yes)
AM_CONDITIONAL(BUILD_COLORING, [test "X$ac_cv_use_epetraext_coloring" != "Xno"])

TAC_ARG_ENABLE_FEATURE_SUB( epetraext, zoltan, [Enable EpetraExt interface support for Zoltan], ZOLTAN, no)
AM_CONDITIONAL(BUILD_ZOLTAN, [test "X$ac_cv_use_epetraext_zoltan" != "Xno"])

TAC_ARG_ENABLE_FEATURE_SUB( epetraext, btf, [Enable EpetraExt BTF subpkg], BTF, no)
AM_CONDITIONAL(BUILD_BTF, [test "X$ac_cv_use_epetraext_btf" != "Xno"])

TAC_ARG_ENABLE_FEATURE_SUB( epetraext, matlab, [Enable EpetraExt Matlab subpkg], MATLAB, no)
AM_CONDITIONAL(BUILD_MATLAB, [test "X$ac_cv_use_epetraext_matlab" != "Xno"])

# ------------------------------------------------------------------------
# Checks if tests should be built
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_FEATURE(tests, [Build tests for all Trilinos packages (not all packages are sensitive to this option)], TESTS, yes)
AM_CONDITIONAL(BUILD_TESTS, test "X$ac_cv_use_tests" != "Xno")
TAC_ARG_ENABLE_FEATURE(examples, [Build examples for all Trilinos packages (not all packages are sensitive to this option)], EXAMPLES, no)
AM_CONDITIONAL(BUILD_EXAMPLES, test "X$ac_cv_use_examples" != "Xno")
TAC_ARG_ENABLE_FEATURE(experimental, [Build experimental code for all Trilinos packages (not all packages are sensitive to this option)], EXPERIMENTAL, no)
AM_CONDITIONAL(BUILD_EXPERIMENTAL, test "X$ac_cv_use_experimental" != "Xno")

TAC_ARG_ENABLE_FEATURE_SUB( epetraext, tests, [Build tests for EpetraExt package], TESTS, yes)
AM_CONDITIONAL(BUILD_TESTS, test "X$ac_cv_use_epetraext_tests" != "Xno")
TAC_ARG_ENABLE_FEATURE_SUB( epetraext, examples, [Build examples for EpetraExt package], EXAMPLES, no)
AM_CONDITIONAL(BUILD_EXAMPLES, test "X$ac_cv_use_epetraext_examples" != "Xno")
TAC_ARG_ENABLE_FEATURE_SUB( epetraext, experimental, [Build experimental code for EpetraExt package], EXPERIMENTAL, no)
AM_CONDITIONAL(BUILD_EXPERIMENTAL, test "X$ac_cv_use_epetraext_experimental" != "Xno")

# ------------------------------------------------------------------------

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

# checks for libraries now
dnl Replace `main' with a function in -lm:
AC_SEARCH_LIBS(pow,[m],,AC_MSG_ERROR(Cannot find math library))
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
#fi

if test "X$ac_cv_use_epetraext_zoltan" != "Xno"; then
  AC_CHECK_LIB(metis,METIS_PartGraphKway,[],[AC_MSG_ERROR([Cannot find libmetis.a library needed for Zoltan subpkg, specify using --with-libs=-L<DIR>])])
  AC_CHECK_LIB(parmetis,ParMETIS_PartKway,[],[AC_MSG_ERROR([Cannot find libparmetis.a library needed for Zoltan subpkg, specify using --with-libs=-L<DIR>])])
  AC_CHECK_LIB(zoltan,Zoltan_LB_Balance,[],[AC_MSG_ERROR([Cannot find libzoltan.a library needed for Zoltan subpkg, specify using --with-libs=-L<DIR>])])
fi

# ------------------------------------------------------------------------
# Checks for compiler characteristics
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for linker characteristics
# ------------------------------------------------------------------------

# Determine libraries needed for linking with Fortran
AC_F77_LIBRARY_LDFLAGS

# -- DEC CXX --
dnl if test X${target_cpu} = Xalpha && test X${CXX} = Xcxx; then
dnl    XTRALDADD="cxx_repository/*.o"
dnl fi

dnl AC_MSG_CHECKING([for extra object files to be linked])
dnl if test X$XTRALDADD = X; then
dnl AC_MSG_RESULT(none)
dnl else
dnl AC_MSG_RESULT($XTRALDADD)
dnl fi

AC_SUBST(TRILINOS_TEMPLATE_OBJS)
AC_SUBST(XTRALDADD)

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
AC_CONFIG_FILES([
		Makefile 
		src/Makefile
                test/Makefile
                test/BTF/Makefile
                test/Composite/Makefile
		test/Matlab/Makefile
                test/MapColoring/Makefile
                test/MatrixMatrix/Makefile
                test/Overlap/Makefile
                test/Permutation/Makefile
                test/SingletonFilter/Makefile
                test/Transpose/Makefile
                test/SymmRCM/Makefile
                test/View/Makefile
                test/Zoltan/Makefile
                test/inout/Makefile
                test/scripts/Makefile
                test/scripts/daily/Makefile
                test/scripts/daily/serial/Makefile
                test/scripts/daily/mpi/Makefile
                test/scripts/weekly/Makefile
                test/scripts/weekly/serial/Makefile
                test/scripts/weekly/mpi/Makefile
                example/Makefile
                example/inout/Makefile
                example/zoltan/Makefile
		])

AC_OUTPUT()

# Bye World!
echo "----------------------------------------"
echo "Finished Running EpetraExt Configure Script"
echo "----------------------------------------"
