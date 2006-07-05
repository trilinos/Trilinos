# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# @HEADER
# ************************************************************************
#
#            NOX: An Object-Oriented Nonlinear Solver Package
#                 Copyright (2002) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
# Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
# (rppawlo@sandia.gov).
#
# ************************************************************************
# @HEADER

# ------------------------------------------------------------------------
# Initialization
# ------------------------------------------------------------------------

# This must be the first line in configure.ac.
# Optional 3rd argument is email address for bugs.

AC_INIT([nox],[4.0-dev],[nox-developers@software.sandia.gov])

# Hello World!
echo "----------------------------------------"
echo "Running NOX Configure Script"
echo "----------------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.

AC_CONFIG_SRCDIR(src/NOX_Solver_LineSearchBased.C)

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.

AC_CONFIG_AUX_DIR(config)

# Configure should create src/NOX_Config.h from src/NOX_Config.h.in

AM_CONFIG_HEADER(src/NOX_Config.h:src/NOX_Config.h.in)

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

# Use automake, getting name and version from variables defined by
# AC_INIT command.

AM_INIT_AUTOMAKE(1.6.3)

# Specify required version of autoconf.

AC_PREREQ(2.57)

# ------------------------------------------------------------------------
# Check for MPI C++ compiler (must be done *before* AC_PROG_CXX)
# ------------------------------------------------------------------------
TAC_ARG_CONFIG_MPI

# ------------------------------------------------------------------------
# Checks for programs
# ------------------------------------------------------------------------

AC_PROG_CC(cc gcc)
AC_PROG_CXX(CC g++ c++ cxx)
AC_PROG_CXXCPP
AC_PROG_F77(f77 g77 gfortran f90 xlf90 f95)
AC_PROG_RANLIB

# ------------------------------------------------------------------------
# Checks for extra flags
# ------------------------------------------------------------------------
# Check if --with-flags present, prepend any specs to FLAGS

TAC_ARG_WITH_FLAGS(cxxflags, CXXFLAGS)
TAC_ARG_WITH_FLAGS(ccflags, CCFLAGS)
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
# Check for doxygen dot program.  Used in doxygen config files Doxyfile.in
# ------------------------------------------------------------------------
AC_CHECK_PROG(HAVE_DOXYGEN_DOT, dot, "YES", "NO")

# ------------------------------------------------------------------------
# STL bug on solaris.
# ------------------------------------------------------------------------
# Solaris requires a const argument in constructing maps. This is wrong
# according to the standard!  DEC on the other hand fails if it is const
# so we have to wrap multimap construction with:
###      ifdef MUST_CONST_STL_MAP_KEY
# for compiling on both solaris and DEC platforms.
# The following call will set this flag if needed.
AC_CXX_STL_MAP_KEY_CONST

# ------------------------------------------------------------------------
# Checks for standard library header files
# ------------------------------------------------------------------------

#cstdlib
AC_LANG([C++])
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
AC_CHECK_HEADERS([sstream sstream.h strstream strstream.h], [break],)
AC_CHECK_HEADERS([fstream fstream.h], [break],)

if test "X$ac_cv_header_cstdlib" != "Xyes" && test "X$ac_cv_header_stdlib_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdlib or stdlib.h.  NOX requires at least one or the other.]])
fi

if test "X$ac_cv_header_cstdio" != "Xyes" && test "X$ac_cv_header_stdio_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdio or stdio.h.  NOX requires at least one or the other.]])
fi

if test "X$ac_cv_header_cassert" != "Xyes" && test "X$ac_cv_header_assert_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cassert or assert.h.  NOX requires at least one or the other.]])
fi

if test "X$ac_cv_header_string" != "Xyes" && test "X$ac_cv_header_string_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find string or string.h.  NOX requires at least one or the other.]])
fi

if test "X$ac_cv_header_iostream" != "Xyes" && test "X$ac_cv_header_iostream_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find iostream or iostream.h.  NOX requires at least one or the other.]])
fi

if test "X${host}" != "Xtflop"; then
  if test "X$ac_cv_header_cmath" != "Xyes" && test "X$ac_cv_header_math_h" != "Xyes"; then
    AC_MSG_ERROR([[Cannot find cmath or math.h.  NOX requires at least one or the other.]])
  fi
else
  # TFLOP
  if test "X$ac_cv_header_iomanip" != "Xyes" && test "X$ac_cv_header_iomanip_h" != "Xyes"; then
    AC_MSG_ERROR([[Cannot find iomanip or iomanip.h.  NOX requires at least one or the other.]])
  fi
fi  

if test "X$ac_cv_header_map" != "Xyes" && test "X$ac_cv_header_map_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find map or map.h.  NOX requires at least one or the other.]])
fi

if test "X$ac_cv_header_vector" != "Xyes" && test "X$ac_cv_header_vector_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find vector or vector.h.  NOX requires at least one or the other.]])
fi

if test "X$ac_cv_header_algorithm" != "Xyes" && test "X$ac_cv_header_algo_h" != "Xyes" && test "X$ac_cv_header_algorithm_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find algorithm, algo.h, or algorithm.h.  NOX requires at least one of these.]])
fi

if test "X$ac_cv_header_sstream" != "Xyes" && test "X$ac_cv_header_strstream" != "Xyes" && test "X$ac_cv_header_sstream" != "Xyes" && test "X$ac_cv_header_strstream" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find sstream, strstream, sstream.h, or strstream.h.  NOX requires at least one of these.]])
fi

if test "X$ac_cv_header_fstream" != "Xyes" && test "X$ac_cv_header_fstream_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find fstream or fstream.h.  NOX requires at least one or the other.]])
fi

# can we say using std::ostream?
AC_LANG([C++])
AC_MSG_CHECKING(whether "using std::ostream" is okay)
AC_COMPILE_IFELSE(
[
  AC_LANG_PROGRAM(
  [[
#ifdef HAVE_IOSTREAM
#include <iostream>
#elif defined(HAVE_IOSTREAM_H)
#include <iostream.h>
#endif
  ]],
  [[using std::ostream;]])
],
[AC_MSG_RESULT(yes)],
[
  AC_MSG_RESULT(no)
  AC_DEFINE(NO_USING_STD_OSTREAM,1,[define if "using std::ostream" is not okay])
])

# ------------------------------------------------------------------------
# Checks for typedefs, structures, and compiler characteristics.
# ------------------------------------------------------------------------
AC_C_CONST
AC_C_INLINE
AC_CXX_INVALID_TEMPLATE_QUALIFIER

# ------------------------------------------------------------------------
# Should we enable library checks?
# ------------------------------------------------------------------------
#Sets the AC_DEFINE(HAVE_LIBCHECK)  and the ac variable ac_cv_use_libcheck

TAC_ARG_ENABLE_FEATURE(libcheck, [Check for third-party (non-Trilinos) libraries including BLAS. LAPACK, and math library.  If tests or examples are enabled, this flag will be overriden and libchecks enabled.], LIBCHECK, no)

# ------------------------------------------------------------------------
# Check if the compiler supports NaN and inf checking 
# ------------------------------------------------------------------------
TAC_ARG_FINITE_NUMBER_CHECK

# ------------------------------------------------------------------------
# Enable nox debugging
# ------------------------------------------------------------------------
AC_ARG_ENABLE(nox-debug,
[AC_HELP_STRING([--enable-nox-debug],[Build NOX in debug mode - FOR DEVELOPERS ONLY!])],
[NOX_DEBUG_MODE=$enableval],
[NOX_DEBUG_MODE=no]
)

AC_MSG_CHECKING(whether nox is being built in debug mode)
if test "X$NOX_DEBUG_MODE" != "Xno"; then
  AC_DEFINE([HAVE_NOX_DEBUG],1,[Define to build NOX in debug mode.])
  AC_MSG_RESULT(yes)
else
  AC_MSG_RESULT(no)
fi

# ------------------------------------------------------------------------
# Determine if nox is being built as part of Trilinos or as a stand-alone
# package.  This will affect the default settings such as building loca.
# ------------------------------------------------------------------------
AC_ARG_WITH(trilinos,
[AC_HELP_STRING([--with-trilinos],[This flag is set by the top level of Trilinos. If NOX is being built as part of a Trilinos build, this will automatically be enabled. Defaults in NOX and LOCA will be different depeneding on this flag. Default is no.])],
[INSIDE_TRILINOS=$withval],
[INSIDE_TRILINOS=no])

AC_MSG_CHECKING(whether nox is being built as part of a Trilinos build)
AC_MSG_RESULT(${INSIDE_TRILINOS})
AM_CONDITIONAL(INSIDE_TRILINOS, test X${INSIDE_TRILINOS} = Xyes)

# ------------------------------------------------------------------------
# LOCA Library Check
# ------------------------------------------------------------------------
# If inside trilinos, the default is to build loca.  If not,
# the default is not to build it.
# We also must turn off loca builds if inside trilinos and
# --disable-default-packages is set.
AC_ARG_ENABLE([default-packages],
AC_HELP_STRING([--enable-default-packages, Flag used to toggle the automatic building of a list of default packages in Trilinos. Both NOX and LOCA are enabled if default packages are being built. Default is yes.]),
ac_cv_use_default_packages=$enableval, ac_cv_use_default_packages=yes)

AC_MSG_CHECKING(whether we should enable default libs)

if test "X$ac_cv_use_default_packages" != "Xno"; then
  AC_MSG_RESULT(yes)
else
  AC_MSG_RESULT(no)
fi

if test X${ac_cv_use_default_packages} != Xyes; then
  LOCA_DEFAULT=no
else
  LOCA_DEFAULT=${INSIDE_TRILINOS}
fi

AC_ARG_ENABLE(loca,
[AC_HELP_STRING([--enable-loca],[Compile the loca libraries. Default is yes if part of a Trilinos build and no if nox is being built stand alone.])],
[BUILD_LIBLOCA=${enableval}],
[BUILD_LIBLOCA=${LOCA_DEFAULT}]
)
AC_MSG_CHECKING(whether we should build the loca library)
AC_MSG_RESULT([${BUILD_LIBLOCA}])
AM_CONDITIONAL(BUILD_LIBLOCA, test X${BUILD_LIBLOCA} = Xyes)

# ------------------------------------------------------------------------
# NOX interface libraries (uses config/tac_arg_enable_nox_interface.m4)
# ------------------------------------------------------------------------

NOX_INTERFACE(lapack, LIBNOXLAPACK, nox)
NOX_INTERFACE(epetra, LIBNOXEPETRA, nox)
NOX_INTERFACE(petsc, LIBNOXPETSC, nox)
NOX_INTERFACE(belos, LIBNOXBELOS, nox)

# ------------------------------------------------------------------------
# LOCA interface libraries
# ------------------------------------------------------------------------
# Loca support is enabled if nox support and loca are enabled
AM_CONDITIONAL(BUILD_LIBLOCALAPACK, test X${ADDON_lapack} = Xyes && 
                                    test X${BUILD_LIBLOCA} = Xyes)
AM_CONDITIONAL(BUILD_LIBLOCAEPETRA, test X${ADDON_epetra} = Xyes && 
                                    test X${BUILD_LIBLOCA} = Xyes)

# ------------------------------------------------------------------------
# Library check variable.  This can be set to yes by examples or tests
# in the next two sections.
# ------------------------------------------------------------------------
REQUIRE_LIB_CHECK=no

# ------------------------------------------------------------------------
# Checks for Regression Testing
# ------------------------------------------------------------------------

  AC_DEFUN([NOX_TESTING],
  [
    AC_MSG_CHECKING(whether we should build the nox $1 regression tests)
    if test X${ADDON_TESTS} = Xyes && test X${ADDON_$1} = Xyes; then
      AC_MSG_RESULT(yes)
      ADDON_$1_TESTS=yes

      AC_MSG_CHECKING(whether we should build the loca $1 regression tests)
      if test X${BUILD_LIBLOCA} = Xyes; then
	ADDON_loca_$1_TESTS=yes
        AC_MSG_RESULT(yes)
      else
	ADDON_loca_$1_TESTS=no
        AC_MSG_RESULT(no)
      fi

    else
      AC_MSG_RESULT(no)
      ADDON_$1_TESTS=no
      AC_MSG_CHECKING(whether we should build the loca $1 regression tests)
      AC_MSG_RESULT(no)
      ADDON_loca_$1_TESTS=no
    fi

    if test X${ADDON_$1_TESTS} = Xyes; then
      REQUIRE_LIB_CHECK=yes
    fi

  ])

  # See if testing is enabled
  AC_MSG_CHECKING(whether we should build the regression tests)
  AC_ARG_ENABLE(nox-tests,
  [AC_HELP_STRING([--enable-nox-tests],[compile regression tests, over-rides the trilinos flag "--enable-tests" so that nox tests can be built separately from other trilinos package tests])],
  [ADDON_TESTS=$enableval],
  [
     AC_ARG_ENABLE(tests,
     [AC_HELP_STRING([--enable-tests],[compile regression tests, will be ignored if "--enable-nox-tests" is specified])],
     [ADDON_TESTS=$enableval],
     [ADDON_TESTS=no]
     )
  ]
  )
  AC_MSG_RESULT([${ADDON_TESTS}])

  # Basic tests that are not dependent on a support package
  if test X${ADDON_TESTS} = Xyes; then
    REQUIRE_LIB_CHECK=yes
  fi
  AM_CONDITIONAL(BUILD_NOX_TESTS, test X${ADDON_TESTS} = Xyes)
  AM_CONDITIONAL(BUILD_LOCA_TESTS, test X${ADDON_TESTS} = Xyes && 
                                   test X${BUILD_LIBLOCA} = Xyes)

  # lapack testing
  NOX_TESTING(lapack)
  AM_CONDITIONAL(BUILD_LAPACK_TESTS, test X${ADDON_lapack_TESTS} = Xyes)
  AM_CONDITIONAL(BUILD_LOCA_LAPACK_TESTS, test X${ADDON_loca_lapack_TESTS} = Xyes)
  AC_SUBST(ADDON_lapack_TESTS)
  AC_SUBST(ADDON_loca_lapack_TESTS)

  # epetra testing
  NOX_TESTING(epetra)
  AM_CONDITIONAL(BUILD_EPETRA_TESTS, test X${ADDON_epetra_TESTS} = Xyes)
  AM_CONDITIONAL(BUILD_LOCA_EPETRA_TESTS, test X${ADDON_loca_epetra_TESTS} = Xyes)
  AC_SUBST(ADDON_epetra_TESTS)
  AC_SUBST(ADDON_loca_epetra_TESTS)

  # petsc testing
  NOX_TESTING(petsc)
  AM_CONDITIONAL(BUILD_PETSC_TESTS, test X${ADDON_petsc_TESTS} = Xyes)
  AM_CONDITIONAL(BUILD_LOCA_PETSC_TESTS, test X${ADDON_loca_petsc_TESTS} = Xyes)
  AC_SUBST(ADDON_petsc_TESTS)
  AC_SUBST(ADDON_loca_petsc_TESTS)

  # belos testing
  NOX_TESTING(belos)
  AM_CONDITIONAL(BUILD_BELOS_TESTS, test X${ADDON_belos_TESTS} = Xyes)
  AC_SUBST(ADDON_belos_TESTS)

  # Determine whether to build the new test utils library
  AM_CONDITIONAL(BUILD_TESTUTILS, test X${REQUIRE_LIB_CHECK} = Xyes)

# ------------------------------------------------------------------------
# Checks for Examples
# ------------------------------------------------------------------------

  AC_DEFUN([NOX_EXAMPLES],
  [
    AC_MSG_CHECKING(whether we should build the nox $1 examples)
    if test X${ADDON_EXAMPLES} = Xyes && test X${ADDON_$1} = Xyes; then
      AC_MSG_RESULT(yes)
      ADDON_$1_EXAMPLES=yes

      AC_MSG_CHECKING(whether we should build the loca $1 examples)
      if test X${BUILD_LIBLOCA} = Xyes; then
	ADDON_loca_$1_EXAMPLES=yes
        AC_MSG_RESULT(yes)
      else
	ADDON_loca_$1_EXAMPLES=no
        AC_MSG_RESULT(no)
      fi

    else
      AC_MSG_RESULT(no)
      ADDON_$1_EXAMPLES=no
      AC_MSG_CHECKING(whether we should build the loca $1 examples)
      AC_MSG_RESULT(no)
      ADDON_loca_$1_EXAMPLES=no
    fi

    if test X${ADDON_$1_EXAMPLES} = Xyes; then
      REQUIRE_LIB_CHECK=yes
    fi

  ])

  # See if examples are enabled
  AC_MSG_CHECKING(whether we should build the examples)
  AC_ARG_ENABLE(nox-examples,
  [AC_HELP_STRING([--enable-nox-examples],[compile examples, over-rides the trilinos flag "--enable-examples" so that nox examples can be built separately from other trilinos package examples])],
  [ADDON_EXAMPLES=$enableval],
  [
     AC_ARG_ENABLE(examples,
     [AC_HELP_STRING([--enable-examples],[compile examples, will be ignored if "--enable-nox-examples" is specified])],
     [ADDON_EXAMPLES=$enableval],
     [ADDON_EXAMPLES=no]
     )
  ]
  )
  AC_MSG_RESULT([${ADDON_EXAMPLES}])

  # lapack 
  NOX_EXAMPLES(lapack)
  AM_CONDITIONAL(BUILD_LAPACK_EXAMPLES, test X${ADDON_lapack_EXAMPLES} = Xyes)
  AM_CONDITIONAL(BUILD_LOCA_LAPACK_EXAMPLES, test X${ADDON_loca_lapack_EXAMPLES} = Xyes)
  AC_SUBST(ADDON_lapack_EXAMPLES)
  AC_SUBST(ADDON_loca_lapack_EXAMPLES)

  # epetra 
  NOX_EXAMPLES(epetra)
  AM_CONDITIONAL(BUILD_EPETRA_EXAMPLES, test X${ADDON_epetra_EXAMPLES} = Xyes)
  AM_CONDITIONAL(BUILD_LOCA_EPETRA_EXAMPLES, test X${ADDON_loca_epetra_EXAMPLES} = Xyes)
  AC_SUBST(ADDON_epetra_EXAMPLES)
  AC_SUBST(ADDON_loca_epetra_EXAMPLES)

  # petsc (no loca tests)
  NOX_EXAMPLES(petsc)
  AM_CONDITIONAL(BUILD_PETSC_EXAMPLES, test X${ADDON_petsc_EXAMPLES} = Xyes)
  AC_SUBST(ADDON_petsc_EXAMPLES)

  # belos
  NOX_EXAMPLES(belos)
  AM_CONDITIONAL(BUILD_BELOS_EXAMPLES, test X${ADDON_belos_EXAMPLES} = Xyes)
  AC_SUBST(ADDON_belos_EXAMPLES)
  
# ------------------------------------------------------------------------
# LOCA-Anasazi support
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_CAN_USE_PACKAGE(loca, anasazi, LOCA_ANASAZI, BUILD_LOCAANASAZI, yes, [Build anasazi support into the loca library.], [Implicitly build anasazi support into the loca library. Can be overridden with --<en|dis>able-loca-anasazi.])

# ------------------------------------------------------------------------
# LOCA-MF support (MF is NOT a Trilinos package - do not use the 
# tac_arg_can_use_package.m4)
# ------------------------------------------------------------------------
AC_ARG_WITH(loca-mf,
[AC_HELP_STRING([--with-loca-mf=PATH],[Enable Multifario (MF) support in LOCA, optionally providing the path to the base multifario directory (e.g., /usr/local/multifario). Default is no.])],
[
  if test "X${withval}" != "Xyes"; then
     MF_PATH=${withval}
  fi
  BUILD_LOCA_MF_SUPPORT=yes
],
[BUILD_LOCA_MF_SUPPORT=no]
)
AC_MSG_CHECKING(whether we should build MF support into loca library)
AC_MSG_RESULT([${BUILD_LOCA_MF_SUPPORT}])
AM_CONDITIONAL(BUILD_LIBLOCAMF, test X${BUILD_LOCA_MF_SUPPORT} = Xyes)

AC_ARG_WITH(mf-incdir,
[AC_HELP_STRING([--with-mf-incdir=PATH],[Path to Multifario include directory])],
[MF_INCDIR=${withval}],
[
  if test "X${MF_PATH}" != "X"; then
     MF_INCDIR="${MF_PATH}/include"
  else
     MF_INCDIR=
  fi
]
)

AC_ARG_WITH(mf-libdir,
[AC_HELP_STRING([--with-mf-libdir=PATH],[Path to Multifario library directory])],
[MF_LIBDIR=${withval}],
[
  if test "X${MF_PATH}" != "X"; then
     MF_LIBDIR="${MF_PATH}/lib"
  else
     MF_LIBDIR=
  fi
]
)

AC_ARG_WITH(mf-libs,
[AC_HELP_STRING([--with-mf-libs=LIBS],[Multifario libraries (e.g., "-lMF -lExpCmp")])],
[MF_LIBS=${withval}],
[MF_LIBS="-lMF -lExpCmp"]
)

AC_SUBST(MF_INCDIR)
AC_SUBST(MF_LIBDIR)
AC_SUBST(MF_LIBS) 

# ------------------------------------------------------------------------
# NOX-EpetraExt support
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_CAN_USE_PACKAGE(nox, epetraext, NOX_EPETRAEXT, BUILD_NOX_EPETRAEXT_SUPPORT, yes, [Build epetraext support into the noxepetra library.], [Implicitly build epetraext support into the noxepetra library. Can be overridden with --<en|dis>able-nox-epetraext.])

# ------------------------------------------------------------------------
# NOX-Epetra Multiphysics support
# ------------------------------------------------------------------------
AM_CONDITIONAL(BUILD_EPETRA_MULTIPHYSICS, false)

# ------------------------------------------------------------------------
# NOX-Epetra-ML support
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_CAN_USE_PACKAGE(nox, ml, NOX_ML_EPETRA, BUILD_NOX_EPETRA_ML_SUPPORT, yes, [Build ml support into the noxepetra library.], [Implicitly build ml support into the noxepetra library.  Can be overridden with --<en|dis>able-nox-ml.])

# ------------------------------------------------------------------------
# NOX/LOCA-Python Check.
# ------------------------------------------------------------------------
# RPP: 11/30/2005 - disabling all python in NOX.  The refactor from
# last month added smart pointers that phyton currently can't handle.
# Until Bill Spotz fixes this for Trilinos in general, we can't
# support python.  
BUILD_NOX_PYTHON=no

#TAC_ENABLE_PYTHON

# Separate builds for loca and nox
#AC_ARG_ENABLE(nox-python,
#[AC_HELP_STRING([--disable-nox-python],[Used ONLY to disable the python supportin nox.  It cannot be used to enable python support.  Use the flag --enable-python to enable support.])],
#[
#  BUILD_NOX_PYTHON=${enableval}
#  if test X${BUILD_NOX_PYTHON} != Xno; then
#    echo "*****************************"
#    echo "NOX Configure ERROR - the flag \"--enable-nox-python\" is not"
#    echo "allowed!  To enable python support in nox, use --enable-python."
#    echo "The flag \"--disable-nox-python\" is allowed to specifically"
#    echo "disable nox python support while allowing python support in"
#    echo "other Trilinos packages."
#    echo "*****************************"
#    AC_MSG_ERROR([The flag \"--enable-nox-python\" is not allowed.])
#  fi
#],
#[
#  AC_ARG_ENABLE(python,
#  [],
#  [BUILD_NOX_PYTHON=${enableval}],
#  [BUILD_NOX_PYTHON=no]
#  )
#]
#)

if test X${BUILD_NOX_PYTHON} = Xyes && test X${BUILD_LOCA} = Xyes; then
  ## RPP: For now we have completely disabled python support in loca
  ##BUILD_LOCA_PYTHON=yes
  BUILD_LOCA_PYTHON=no
else
  BUILD_LOCA_PYTHON=no
fi
AM_CONDITIONAL(BUILD_NOX_PYTHON, test X${BUILD_NOX_PYTHON} = Xyes)
AM_CONDITIONAL(BUILD_LOCA_PYTHON, test X${BUILD_LOCA_PYTHON} = Xyes)

# Determine which addon libraries to support in Python
IMPORT_LAPACK=0
if test X${BUILD_PYTHON} = Xyes && test X${ADDON_lapack} = Xyes; then
  IMPORT_LAPACK=1
fi
AC_SUBST(IMPORT_LAPACK)

IMPORT_EPETRA=0
if test X${BUILD_PYTHON} = Xyes && test X${ADDON_epetra} = Xyes; then
  IMPORT_EPETRA=1
fi
AC_SUBST(IMPORT_EPETRA)

# ------------------------------------------------------------------------
# User defined directories 
# ------------------------------------------------------------------------
# enable use of --with-incdirs="-Lincdir1 -Lincdir2 ..." to prepend to CPPFLAGS
TAC_ARG_WITH_INCDIRS

# ------------------------------------------------------------------------
# Checks for basic libraries for examples and/or lapack support
# ------------------------------------------------------------------------
TAC_ARG_WITH_BLASLIB
TAC_ARG_WITH_LAPACKLIB

# Define F77_FUNC - used for correctly munging Fortran function calls
if test X${ADDON_lapack} = Xyes || test X${REQUIRE_LIB_CHECK} = Xyes || test X${BUILD_LIBLOCA} = Xyes; then

  AC_F77_WRAPPERS

fi

if test X${REQUIRE_LIB_CHECK} = Xyes || test X${ac_cv_use_libcheck} = Xyes; then

  # Math library
  AC_SEARCH_LIBS(pow,[m],,AC_MSG_WARN(Cannot find math library))
  AC_CHECK_FUNCS([pow sqrt])

  # Fortran libraries
  AC_F77_LIBRARY_LDFLAGS
  LIBS=" ${FLIBS} ${LIBS}"

  # Find valid BLAS
  ACX_BLAS([LIBS="${BLAS_LIBS} ${LIBS}"],[AC_MSG_ERROR([Cannot find BLAS library])])

  # Find valid LAPACK
  ACX_LAPACK([LIBS=" ${LAPACK_LIBS} ${LIBS}"],[AC_MSG_ERROR([Cannot find LAPACK library])])

  # Find valid LAPACK DGGEV function
  ACX_DGGEV

fi

# ------------------------------------------------------------------------
# Checks for Trilinos header files
# No other Trilinos library does this.  We only do it to support 
# building nox separate from Trilinos.
# ------------------------------------------------------------------------

# NOX always requires Teuchos!!!!
TRILINOS_HEADER(teuchos,Teuchos_config.h,src)
TRILINOS_HEADER(teuchos,Teuchos_Flops.hpp,src)

if test X${ADDON_epetra} = Xyes; then

  # We have to check for the standard Trilinos headers as well as
  # those created by configure.
  TRILINOS_HEADER(epetra,Epetra_config.h,src)
  TRILINOS_HEADER(epetra,Epetra_DataAccess.h,src)
  TRILINOS_HEADER(aztecoo,AztecOO_config.h,src)
  TRILINOS_HEADER(aztecoo,AztecOO_StatusType.h,src)
  TRILINOS_HEADER(ifpack,Ifpack_config.h,src)
  TRILINOS_HEADER(ifpack,Ifpack_ScalingType.h,src)

  if test X${ac_cv_use_epetraext} = Xyes; then
    TRILINOS_HEADER(epetraext,EpetraExt_config.h,src)
    TRILINOS_HEADER(epetraext,EpetraExt_ConfigDefs.h,src)
  fi

  if test X${ac_cv_use_ml} = Xyes; then
    TRILINOS_HEADER(ml,ml_config.h,src)
    TRILINOS_HEADER(ml,ml_common.h,src/Include)
  fi

fi

if test X${ADDON_belos} = Xyes; then
  TRILINOS_HEADER(belos,Belos_config.h,src)
  TRILINOS_HEADER(belos,BelosTypes.hpp,src)
  AC_DEFINE(HAVE_NOX_BELOS,1,[define that Belos is activated in nox])
fi

if test X${ac_cv_use_anasazi} = Xyes; then
  TRILINOS_HEADER(anasazi,Anasazi_config.h,src)
  TRILINOS_HEADER(anasazi,AnasaziTypes.hpp,src)
fi

if test X${BUILD_PYTHON} = Xyes && test X${ADDON_epetra} = Xyes; then
# Commenting out the check for PyTrilinos_config.h for right now since
# PyTrilinos is configured after NOX, so this file may not exist
#  TRILINOS_HEADER(PyTrilinos,PyTrilinos_config.h,src)
  TRILINOS_HEADER(PyTrilinos,FileStream.h,src)
fi

# ------------------------------------------------------------------------
# Petsc include file check for building noxpetsc library
# ------------------------------------------------------------------------

AC_ARG_WITH(petsc-incdir,
[  --with-petsc-incdir     specify location of PETSc header files],
[petsc_INC=$withval])


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

if test X${ADDON_petsc} = Xyes; then

  PETSC_INCLUDE_CHECK(petsc,petsc.h)
  PETSC_INCLUDE_CHECK(petsccommon,base)
  PETSC_INCLUDE_CHECK(petscbmake,petscconf.h)
  if test "X${HAVE_PKG_MPI}" = "Xno"; then
    PETSC_INCLUDE_CHECK(mpiuni,mpi.h)
  fi

  PETSC_CXXFLAGS="-I${petsc_INC} -I${petsccommon_INC} -I${petscbmake_INC}"
  if test "X${HAVE_PKG_MPI}" = "Xno"; then
    PETSC_CXXFLAGS=" ${PETSC_CXXFLAGS} -I${mpiuni_INC}"
  fi

  AC_SUBST(PETSC_CXXFLAGS)

fi

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

if test X${ADDON_petsc_EXAMPLES} = Xyes ||
   test X${ADDON_petsc_TESTS} = Xyes; then

  # Check for petsc
  PETSC_LIB_CHECK(petsc)

fi

# ------------------------------------------------------------------------
# Checks for linker characteristics
# ------------------------------------------------------------------------

## RPP: 5/23/2003 - CXXREPOS is not standard.  All of Trilinos uses
##                  the flag TRILINOS_TEMPLATE_OBJECTS instead
##                  of CXXREPOS.  I added the check for
##                  TRILINOS_TEMPLATE_OBJECTS to be consistent
##                  with the rest of trilinos (otherwise trilinos
##                  users must set redundant flags when configuring).
AC_MSG_CHECKING([for trilinos template-specific object files to be linked])
if test "X$TRILINOS_TEMPLATE_OBJS" = "X"; then
 AC_MSG_RESULT(none)
else
 AC_MSG_RESULT($TRILINOS_TEMPLATE_OBJS)
 CXXREPOS=$TRILINOS_TEMPLATE_OBJS
fi

# -- DEC CXX --
if test "X$TRILINOS_TEMPLATE_OBJS" = "X"; then
  case $target in
    alpha*)
      if test $CXX = cxx; then
        CXXREPOS="cxx_repository/*.o"
	TRILINOS_TEMPLATE_OBJS="cxx_repository/*.o"
      fi
    ;;
  esac
fi

AC_MSG_CHECKING([for C++ repository files])
if test X$CXXREPOS = X; then
AC_MSG_RESULT(none)
else
AC_MSG_RESULT($CXXREPOS)
fi

AC_SUBST(CXXREPOS)

# ------------------------------------------------------------------------
# Should we compile 'predevelopment'  versions of the code
# ------------------------------------------------------------------------
AC_ARG_ENABLE(prerelease,
[AC_HELP_STRING([--enable-prerelease],[Enable compilation and documentation of prerelease code.  Prerelease code is code under development that we restrict from release to the public.  It is avaliable to developers with no support - use at your own risk.])],
[WITH_PRERELEASE=$enableval],
[WITH_PRERELEASE=no])

if test X${WITH_PRERELEASE} = Xyes; then
   CPPFLAGS="${CPPFLAGS} -DWITH_PRERELEASE"
fi

if test X${WITH_PRERELEASE} = Xyes; then
  WITH_PRERELEASE_SUBST=WITH_PRERELEASE
else
  WITH_PRERELEASE_SUBST=
fi

AM_CONDITIONAL(BUILD_PRERELEASE, test X${WITH_PRERELEASE} = Xyes)
AC_SUBST(WITH_PRERELEASE_SUBST)

AC_MSG_CHECKING([whether or not we should build prerelease code])
AC_MSG_RESULT($WITH_PRERELEASE)

# ------------------------------------------------------------------------
# Perform substitutions in output files
# ------------------------------------------------------------------------

AC_SUBST(ac_aux_dir)
AC_SUBST(abs_builddir)
AC_SUBST(abs_top_builddir)

AC_SUBST(MPI_LIBS)

# ------------------------------------------------------------------------
# Output files
# ------------------------------------------------------------------------

AC_CONFIG_FILES([
		Makefile
		Makefile.export.nox
		doc/Makefile
		doc/Doxyfile
		src/Makefile
                src-lapack/Makefile
		src-epetra/Makefile
                src-petsc/Makefile
		src-belos/Makefile
		src-loca/Makefile
		src-loca/doc/Makefile
		src-loca/src/Makefile
                src-loca/src-epetra/Makefile
		src-loca/src-lapack/Makefile
		src-loca/src-mf/Makefile
		examples/Makefile
		examples/lapack/Makefile
		examples/lapack/NOX_SimpleExamples/Makefile
		examples/lapack/LOCA_Chan/Makefile
		examples/lapack/LOCA_Brusselator/Makefile
		examples/lapack/LOCA_1DFDPitchfork/Makefile
		examples/epetra/Makefile
		examples/epetra/NOX_1DFEMNonlinear/Makefile
		examples/epetra/NOX_1DFEMNonlinearColoring/Makefile
		examples/epetra/LOCA_Tcubed/Makefile
		examples/epetra/LOCA_Brusselator_xyzt/Makefile
                examples/petsc/Makefile
                examples/petsc/DS6.5.1/Makefile
                examples/petsc/1DFiniteDifference/Makefile
                test/Makefile
                test/basic/Makefile
                test/utils/Makefile
                test/lapack/Makefile
		test/lapack/NOX/Makefile
                test/lapack/NOX_NewTest/Makefile
                test/lapack/NOX_Group/Makefile
		test/lapack/LOCA_TestProblems/Makefile
		test/lapack/LOCA_UnitTests/Makefile
                test/lapack/LOCA_ChanProblem/Makefile
                test/epetra/Makefile
                test/epetra/NOX_NewTest/Makefile
                test/epetra/NOX_Vector/Makefile
                test/epetra/NOX_Operators/Makefile
                test/epetra/NOX_Group/Makefile
		test/epetra/LOCA_TestProblems/Makefile
		test/epetra/LOCA_UnitTests/Makefile
		test/epetra/LOCA_Tcubed/Makefile
		test/epetra/LOCA_Pitchfork/Makefile
                test/epetra/1Dfem/Makefile
                test/epetra/MultiPhysics/Makefile
                test/epetra/DS6.5.1/Makefile
                test/epetra/Brusselator/Makefile
                test/epetra/FiniteDifferenceColoring/Makefile
                test/petsc/Makefile
		python/Makefile
		python/src/Makefile
		python/src/TRILINOS_HOME_DIR
		python/src/__init__.py
		python/src/PyTrilinos/__init__.py
		python/src/.deps/dummy
		python/test/Makefile
		python/test-epetra/Makefile
		python/test-lapack/Makefile
		python/example-epetra/Makefile
		python/example-lapack/Makefile
		python/src-loca/Makefile
		python/src-loca/src/Makefile
		python/src-loca/src/TRILINOS_HOME_DIR
		python/src-loca/src/__init__.py
		python/src-loca/src/.deps/dummy
		python/src-loca/src/PyTrilinos/__init__.py
		python/src-loca/test/Makefile
		python/src-loca/test-epetra/Makefile
		python/src-loca/test-lapack/Makefile
		python/src-loca/example-epetra/Makefile
		python/src-loca/example-lapack/Makefile
		])

AC_OUTPUT

# ------------------------------------------------------------------------
# Bye World!
# ------------------------------------------------------------------------

echo "------------------------------------------"
echo "Finished Running NOX Configure Script"
echo "------------------------------------------"
