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

AC_INIT(epetraext, 3.7, maherou@sandia.gov)

# Hello World!
echo "----------------------------------------"
echo "Running EpetraExt Configure Script"
echo "----------------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.

AC_CONFIG_SRCDIR([src/transform/EpetraExt_View_CrsGraph.cpp])

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

AM_INIT_AUTOMAKE(1.10 no-define tar-ustar)

# Specify required version of autoconf.

AC_PREREQ(2.61)

# Specify revision number of configure.ac that was used to generate
# configure.

dnl Commented out this line for now because it causes headaches in
dnl committing.
dnl AC_REVISION($Revision$)

# ------------------------------------------------------------------------# Check to see if MPI enabled and if any special configuration done
# ------------------------------------------------------------------------                                                                          
TAC_ARG_CONFIG_MPI
                                                                          
# ------------------------------------------------------------------------
# Check to see if MPI enabled and if any special configuration done
# ------------------------------------------------------------------------                                                                          
TAC_ARG_CONFIG_MPI

# ------------------------------------------------------------------------
# Check if Fortran is enabled or not
# ------------------------------------------------------------------------

#TAC_ARG_ENABLE_OPTION(fortran, [enable Fortran support], FORTRAN_SUPPORT, yes)
#TAC_ARG_ENABLE_OPTION(fortran90, [enable Fortran 90 support (automatically enables Fortran support)], FORTRAN90_SUPPORT, no)
#AM_CONDITIONAL(FORTRAN90_DISABLED, [test "X$ac_cv_use_fortran90" != "Xyes"])
#AM_CONDITIONAL(FORTRAN_DISABLED, [test "X$ac_cv_use_fortran" != "Xyes" && test "X$ac_cv_use_fortran90" != "Xyes"])
                                                                          
# ------------------------------------------------------------------------
# Checks for programs
# ------------------------------------------------------------------------

AC_PROG_CC(cc gcc)
AC_PROG_CXX(CC g++ c++ cxx)
#if test "X$ac_cv_use_fortran90" != "Xno"; then
# AC_PROG_FC(gfortran f90 xlf90 f95)
#fi
#elif test "X$ac_cv_use_fortran" != "Xno"; then
 AC_PROG_F77(f77 g77 gfortran f90 xlf90 f95)
#fi
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
if test "X$ac_cv_use_fortran90" != "Xno"; then
 TAC_ARG_WITH_FLAGS(fcflags, FCFLAGS)
fi
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

AC_CHECK_HEADERS([cstdlib stdlib.h],,)
AC_CHECK_HEADERS([cstdio stdio.h],,)
AC_CHECK_HEADERS([cassert assert.h], [break],)
AC_CHECK_HEADERS([cmath math.h],,)
AC_CHECK_HEADERS([iomanip iomanip.h], [break],)
AC_CHECK_HEADERS([iostream iostream.h], [break],)
AC_CHECK_HEADERS([string string.h], [break],)
AC_CHECK_HEADERS([algorithm algo.h algorithm.h],,)
AC_CHECK_HEADERS([map vector],,)

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

TAC_ARG_ENABLE_FEATURE_SUB( epetraext, restrict, [Enable EpetraExt Restrict subpkg], RESTRICT, yes)
AM_CONDITIONAL(BUILD_RESTRICT, [test "X$ac_cv_use_epetraext_restrict" != "Xno"])

TAC_ARG_ENABLE_FEATURE_SUB( epetraext, coloring, [Enable EpetraExt Coloring subpkg], COLORING, yes)
AM_CONDITIONAL(BUILD_COLORING, [test "X$ac_cv_use_epetraext_coloring" != "Xno"])

TAC_ARG_ENABLE_FEATURE_SUB( epetraext, zoltan, [Enable EpetraExt interface support for Zoltan], ZOLTAN, no)
AM_CONDITIONAL(BUILD_ZOLTAN, [test "X$ac_cv_use_epetraext_zoltan" != "Xno"])

dnl Check for enabling of zoltan
AC_ARG_ENABLE([zoltan],
AC_HELP_STRING([--enable-zoltan],[Build EpetraExt interface support for Zoltan using native library.]),
ac_cv_use_native_zoltan=yes,
ac_cv_use_native_zoltan=no)

dnl Also check if the ParMETIS interface is enabled, if not disable using the native Zoltan library
TAC_ARG_WITH_PACKAGE(parmetis, [Enable Parmetis support for Zoltan.], PARMETIS, no)
if test "X$ac_cv_use_parmetis" = "Xno"; then
  ac_cv_use_native_zoltan=no
fi

dnl If both Zoltan and the ParMETIS interface are enabled, use the native Zoltan
AM_CONDITIONAL(BUILD_NATIVE_ZOLTAN, [test "X$ac_cv_use_native_zoltan" != "Xno"])

TAC_ARG_ENABLE_FEATURE_SUB( epetraext, btf, [Enable EpetraExt BTF subpkg], BTF, no)
AM_CONDITIONAL(BUILD_BTF, [test "X$ac_cv_use_epetraext_btf" != "Xno"])

TAC_ARG_ENABLE_FEATURE_SUB( epetraext, matlab, [Enable EpetraExt Matlab subpkg], MATLAB, no)
AM_CONDITIONAL(BUILD_MATLAB, [test "X$ac_cv_use_epetraext_matlab" != "Xno"])

TAC_ARG_ENABLE_FEATURE_SUB( epetraext, block, [Enable EpetraExt Block subpkg], BLOCK, yes)
AM_CONDITIONAL(BUILD_BLOCK, [test "X$ac_cv_use_epetraext_block" != "Xno"])

TAC_ARG_ENABLE_FEATURE_SUB( epetraext, modelevaluator, [Enable EpetraExt Model Evaluator subpkg], MODEL_EVALUATOR, yes)
AM_CONDITIONAL(BUILD_MODEL_EVALUATOR, [test "X$ac_cv_use_epetraext_modelevaluator" != "Xno"])

TAC_ARG_WITH_PACKAGE(petsc, [Enable EpetraExt interface support for PETSc.], PETSC,no)
AM_CONDITIONAL(EPETRAEXT_USING_PETSC, [test "X$ac_cv_use_petsc" != "Xno"])

TAC_ARG_ENABLE_FEATURE_SUB(epetraext, graphreorderings, [Enable experimental RCM and AMD graph reorderings.  This also requires --enable-epetraext-experimental.], GRAPH_REORDERINGS,no)
AM_CONDITIONAL(BUILD_GRAPH_REORDERINGS, [test "X$ac_cv_use_epetraext_graphreorderings" != "Xno"])

TAC_ARG_ENABLE_FEATURE_SUB(epetraext, hdf5, 
[Enable EpetraExt interface support for HDF5. This interface requires an
 already installed HDF5 library; the include, the library path and the library
 name must be specified using --with-incdirs, --with-ldflags and --with-libs,
 respectively.], 
                           EPETRAEXT_HDF5, no)

if test "X$ac_cv_use_epetraext_hdf5" != "Xno" ; then
  AM_CONDITIONAL(HAVE_EPETRAEXT_HDF5, true ) 
else
  AM_CONDITIONAL(HAVE_EPETRAEXT_HDF5, false ) 
fi

# ------------------------------------------------------------------------
# Check for PETSc (copied from NOX!)
# ------------------------------------------------------------------------
AC_DEFUN([PETSC_INCLUDE_CHECK],
[
  AC_MSG_CHECKING(for $1 include directory)
  if test -n "${$1_INC}" && test -f ${$1_INC}/$2; then
    AC_MSG_RESULT([${$1_INC}])
  elif test -f ${PETSC_DIR}/bmake/${PETSC_ARCH}/$2; then
    $1_INC=${PETSC_DIR}/bmake/${PETSC_ARCH}
    AC_MSG_RESULT([${$1_INC}])
  elif test -f ${PETSC_DIR}/bmake/common/$2; then
    $1_INC=${PETSC_DIR}/bmake/common
    AC_MSG_RESULT([${$1_INC}])
  elif test -f ${PETSC_DIR}/include/$2; then
    $1_INC=${PETSC_DIR}/include
    AC_MSG_RESULT([${$1_INC}])
  elif test -f ${PETSC_DIR}/include/mpiuni/$2; then
    $1_INC=${PETSC_DIR}/include/mpiuni
    AC_MSG_RESULT([${$1_INC}])
  else
    AC_MSG_RESULT(no)
    AC_MSG_ERROR([Cannot find $1 include directory containing $2. Make sure PETSC_DIR and PETSC_ARCH are defined.])
  fi

])

AM_CONDITIONAL(PETSCLIB_IS_DEFINED, [test "X${PETSC_LIB}" != "X"])
if test X$ac_cv_use_petsc = Xyes; then

  PETSC_INCLUDE_CHECK(petsc,petsc.h)
  PETSC_INCLUDE_CHECK(petsccommon,base)
  PETSC_INCLUDE_CHECK(petscbmake,petscconf.h)
  if test "X${HAVE_PKG_MPI}" = "Xno"; then
    PETSC_INCLUDE_CHECK(mpiuni,mpi.h)
  fi

  EPETRAEXT_PETSC_INCLUDES="-I${PETSC_DIR} -I${petsc_INC} -I${petsccommon_INC} -I${petscbmake_INC}"
  if test "X${HAVE_PKG_MPI}" = "Xno"; then
    EPETRAEXT_PETSC_INCLUDES=" ${EPETRAEXT_PETSC_INCLUDES} -I${mpiuni_INC}"
  fi

  AC_SUBST(EPETRAEXT_PETSC_INCLUDES)
  #JJH 8/26/2008
  # ${PETSC_DIR} and ${PETSC_ARCH} can't be used directly in
  # Makefile.export.epetraext because not all packages know about them.
  # EpetraExt does know about them, so I use AC_SUBST to get their values
  # directly into the export makefile.
  EPETRAEXT_PETSC_DIR=${PETSC_DIR}
  EPETRAEXT_PETSC_ARCH=${PETSC_ARCH}
  AC_SUBST(EPETRAEXT_PETSC_DIR)
  AC_SUBST(EPETRAEXT_PETSC_ARCH)

  AC_MSG_CHECKING(whether PETSC_LIB is defined)
  if test "X${PETSC_LIB}" = "X"; then
    #JJH 8/26/2008
    #The PETSc file "variables" defines PETSC_LIB, which contains all of
    #PETSc's link dependencies.  Unfortunately, "variables" includes another
    #file "petscconf" that redefines top_builddir (argh!). But since we need
    #that other file for PETSC_LIB, I copy both and modify to my liking.  Ugh.
    AC_MSG_RESULT(no)
    cat ${PETSC_DIR}/bmake/common/variables | sed "s/^include /#include /" > petsc_variables
    echo "include ${PWD}/petscconf" >> petsc_variables
    cat ${PETSC_DIR}/bmake/${PETSC_ARCH}/petscconf | sed "s/^top_builddir/#top_builddir/" > petscconf
  else
    AC_MSG_RESULT(yes)
    USER_DEFINED_PETSCLIB=${PETSC_LIB}
    AC_SUBST(USER_DEFINED_PETSCLIB)
  fi
fi #if test X$ac_cv_use_petsc = Xyes; then


# ------------------------------------------------------------------------
# Checks for Petsc libraries
# ------------------------------------------------------------------------

AC_DEFUN([PETSC_LIB_CHECK],
[
  AC_ARG_WITH($1-libdir,
  [AC_HELP_STRING([--with-$1-libdir],[specify location of $1 library])],
  [$1_LIBDIR=$withval])

  AC_MSG_CHECKING(whether PETSC_DIR is defined)
  if test -n "${PETSC_DIR}"; then
     AC_MSG_RESULT([yes (${PETSC_DIR})])
  else
     AC_MSG_ERROR([Cannot find PETSC_DIR environment variable. Please define.])
  fi

  AC_MSG_CHECKING(for $1)
  if test -n "${$1_LIBDIR}" && test -f ${$1_LIBDIR}/lib$1.a; then
    LDFLAGS="-L${$1_LIBDIR} ${LDFLAGS}"
    AC_MSG_RESULT(yes)
  elif test -f ${PETSC_DIR}/lib/${PETSC_ARCH}/lib$1.a; then
    LDFLAGS="-L${PETSC_DIR}/lib/${PETSC_ARCH} ${LDFLAGS}"
    AC_MSG_RESULT(yes)
  else
    AC_MSG_RESULT(no)
    AC_MSG_ERROR([Cannot find $1 library])
  fi

  LIBS="-l$1 ${LIBS}"

])
  
if test X$ac_cv_use_petsc = Xyes; then
  # Check for petsc
  PETSC_LIB_CHECK(petsc)

fi #if test X$ac_cv_use_petsc = Xyes; then

# ------------------------------------------------------------------------
# Interfaces to other Trilinos packages
# ------------------------------------------------------------------------


TAC_ARG_ENABLE_CAN_USE_PACKAGE(epetraext, triutils, EPETRAEXT_TRIUTILS, BUILD_EPETRAEXT_TRIUTILS, yes, [Build triutils support in epetraext.], [Build triutils support in epetraext. Can be overridden with --disable-epetraext-triutils.])

# ------------------------------------------------------------------------
# Checks if tests should be built
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_FEATURE(tests, [Make tests for all Trilinos packages buildable with 'make tests'], TESTS, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK( epetraext, tests, [Make EpetraExt tests buildable with 'make tests'], EPETRAEXT_TESTS)
AM_CONDITIONAL(BUILD_TESTS, test "X$ac_cv_use_epetraext_tests" != "Xno")

TAC_ARG_ENABLE_FEATURE(examples, [Make examples for all Trilinos packages buildable with 'make examples'], EXAMPLES, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK( epetraext, examples, [Make EpetraExt examples buildable with 'make examples'], EPETRAEXT_EXAMPLES)
AM_CONDITIONAL(BUILD_EXAMPLES, test "X$ac_cv_use_epetraext_examples" != "Xno")

#We now build tests and examples through separate make targets, rather than
#during "make".  We still need to conditionally include the test and example
#in SUBDIRS, even though SUB_TEST and SUB_EXAMPLE will never be
#defined, so that the tests and examples are included in the distribution
#tarball.
AM_CONDITIONAL(SUB_TEST, test "X$ac_cv_use_sub_test" = "Xyes")
AM_CONDITIONAL(SUB_EXAMPLE, test "X$ac_cv_use_sub_example" = "Xyes")

TAC_ARG_ENABLE_FEATURE(libcheck, [Check for third-party libraries.  (Cannot be disabled unless tests and examples are also disabled.)], LIBCHECK, yes)

TAC_ARG_ENABLE_FEATURE(experimental, [Build experimental code for all Trilinos packages (not all packages are sensitive to this option)], EXPERIMENTAL, no)
AM_CONDITIONAL(BUILD_EXPERIMENTAL, test "X$ac_cv_use_experimental" != "Xno")
TAC_ARG_ENABLE_FEATURE_SUB( epetraext, experimental, [Build experimental code for EpetraExt package], EXPERIMENTAL, no)
AM_CONDITIONAL(BUILD_EXPERIMENTAL, test "X$ac_cv_use_epetraext_experimental" != "Xno")

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
if test "X$ac_cv_use_epetraext_tests" != "Xno" || test "X$ac_cv_use_epetraext_examples" != "Xno" || test "X$ac_cv_use_libcheck" != "Xno"; then

# checks for libraries now
dnl Replace `main' with a function in -lm:
AC_SEARCH_LIBS(pow,[m],,AC_MSG_ERROR(Cannot find math library))
AC_SEARCH_LIBS(sqrt,[m],,AC_MSG_ERROR(Cannot find math library))

fi
# end of the list of libraries that don't need to be checked for if
# tests and examples are disabled.

if test "X$ac_cv_use_epetraext_zoltan" != "Xno" -a "X$ac_cv_use_native_zoltan" != "Xyes"; then
  dnl Check if METIS and ParMETIS libraries have proper functionality
  dnl NOTE: This is not checked if we are using the native Zoltan library, we will leave it up to Zoltan to check.
  AC_CHECK_LIB(metis,METIS_PartGraphKway,[],[AC_MSG_ERROR([Cannot find libmetis.a library needed for Zoltan subpkg, specify using --with-ldflags=-L<DIR>])])
  AC_CHECK_LIB(parmetis,ParMETIS_PartKway,[],[AC_MSG_ERROR([Cannot find libparmetis.a library needed for Zoltan subpkg, specify using --with-ldflags=-L<DIR>])])
  AC_CHECK_LIB(zoltan,Zoltan_LB_Balance,[],[AC_MSG_ERROR([Cannot find libzoltan.a library needed for Zoltan subpkg, specify using --with-ldflags=-L<DIR>])])
fi

if test "X$ac_cv_use_epetraext_matlab" != "Xno"; then
     
     AC_MSG_CHECKING(for engine.h)
     AC_TRY_CPP([#include "engine.h"], [AC_MSG_RESULT(yes)],
     	[
	 AC_MSG_RESULT(no)
	 echo "-----"
     	 echo "Cannot find Matlab libeng header file engine.h."
	 echo "For informatin on how to specify the proper Matlab headers/libraries see epetraext/doc/matlab.README."
     	 echo "Use --with-incdirs=-I<DIR> to specify the include directory containing engine.h."
     	 echo "-----"
     	 AC_MSG_ERROR(Matlab header file engine.h not found)
     	])
	
	
     AC_MSG_CHECKING(for mex.h)
     AC_TRY_CPP([#include "mex.h"], [AC_MSG_RESULT(yes)],  
	[
	 AC_MSG_RESULT(no)
	 echo "-----"
     	 echo "Cannot find Matlab libmx header file mex.h."
	 echo "For informatin on how to specify the proper Matlab headers/libraries see epetraext/doc/matlab.README."
     	 echo "Use --with-incdirs=-I<DIR> to specify the include directory containing mex.h."
     	 echo "-----"
     	 AC_MSG_ERROR(Matlab header file engine.h not found)
     	])
	
     AC_CHECK_LIB(eng,engPutVariable,[],[
     	AC_CHECK_LIB(eng,engPutArray,
		    [
#AC_DEFINE(USE_ENGPUTARRAY,[],[use older function call for Matlab versions < 6.5])
CXXFLAGS="$CXXFLAGS -DUSE_ENGPUTARRAY"
                      LIBS="$LIBS -leng"],
                    [AC_MSG_ERROR([Cannot find libeng.a library needed for Matlab subpkg, specify using --with-libs=-L<DIR> For informatin on how to specify the proper Matlab headers/libraries see epetraext/doc/matlab.REAME])]
                    )
     ])
     
     AC_CHECK_LIB(mx,mxCreateDoubleMatrix,[],[AC_MSG_ERROR([Cannot find libmx.a library needed for Matlab subpkg, specify using --with-libs=-L<DIR> For informatin on how to specify the proper Matlab headers/libraries see epetraext/doc/matlab.README.])])
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
		Makefile.export.epetraext
		src/Makefile
		test/BTF/Makefile
		test/Composite/Makefile
		test/Makefile
		test/MapColoring/Makefile
		test/Matlab/Makefile
		test/MatrixMatrix/Makefile
		test/Overlap/Makefile
		test/Permutation/Makefile
		test/SingletonFilter/Makefile
		test/SymmRCM/Makefile
		test/Transpose/Makefile
		test/View/Makefile
		test/Zoltan/Makefile
		test/inout/Makefile
		example/Makefile
		example/inout/Makefile
		example/zoltan/Makefile
		example/matlab/Makefile
		example/MapColoring/Makefile
		example/petsc/Makefile
		])

AC_OUTPUT()

# Bye World!
echo "----------------------------------------"
echo "Finished Running EpetraExt Configure Script"
echo "----------------------------------------"
