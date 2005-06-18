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
AC_PROG_CXX(CC g++ gcc c++ cxx)
AC_PROG_CXXCPP
AC_PROG_F77(f77 g77 f90 xlf90)
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

# Libtool flags
TAC_ARG_WITH_FLAGS(libtool-compiler-flags, LIBTOOL_COMPILER_FLAGS)
TAC_ARG_WITH_FLAGS(libtool-linker-flags, LIBTOOL_LINKER_FLAGS)
TAC_ARG_WITH_FLAGS(libtool-cclinker-flags, LIBTOOL_CCLINKER_FLAGS)

# libtoolize LIBTOOL flags (i.e., add -Xcompiler, -Xlinker, -XCClinker flags
# respectively)
AC_LIBTOOLIZE_COMPILER_FLAGS(LIBTOOL_COMPILER_FLAGS)
AC_LIBTOOLIZE_LINKER_FLAGS(LIBTOOL_LINKER_FLAGS)
AC_LIBTOOLIZE_CCLINKER_FLAGS(LIBTOOL_CCLINKER_FLAGS)

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

# ------------------------------------------------------------------------
# Determine if nox is being built as part of Trilinos or as a stand-alone
# package.  This will affect the default settings such as building loca.
# ------------------------------------------------------------------------
AC_ARG_WITH(trilinos,
[AC_HELP_STRING([--with-trilinos],[building inside Trilinos directory structure])],
[INSIDE_TRILINOS=$withval],
[INSIDE_TRILINOS=no])

AC_MSG_CHECKING(whether nox is being built as part of a Trilinos build)
AC_MSG_RESULT(${INSIDE_TRILINOS})

# ------------------------------------------------------------------------
# Which Add-Ons To Build
# ------------------------------------------------------------------------

BUILD_ANY_EXAMPLES=no

AC_DEFUN([NOX_INTERFACE],
[
  AC_ARG_ENABLE(nox-$1,
  [AC_HELP_STRING([--enable-nox-$1],[compile $1 interface libraries])],
  [ADDON_$1=$enableval],
  [ADDON_$1=no]
  )

  AC_ARG_ENABLE(nox-$1-examples,
  [AC_HELP_STRING([--enable-nox-$1-examples],[compile $1 interface examples])],
  [ADDON_$1_EXAMPLES=${enableval}],
  [ADDON_$1_EXAMPLES=no]
  )

  # Cannot build the examples if we don't build the library
  if test ${ADDON_$1_EXAMPLES} = yes; then
     BUILD_ANY_EXAMPLES=yes
     ADDON_$1=yes;
  fi

  AC_MSG_CHECKING(whether we should build the $1 interface libraries)
  AC_MSG_RESULT(${ADDON_$1})

  AC_MSG_CHECKING(whether we should build the $1 interface examples)
  AC_MSG_RESULT(${ADDON_$1_EXAMPLES})

])


# ------------------------------------------------------------------------
# NOX-LAPACK library?
# ------------------------------------------------------------------------

NOX_INTERFACE(lapack)
AM_CONDITIONAL(BUILD_LIBNOXLAPACK, test X${ADDON_lapack} = Xyes)
AM_CONDITIONAL(BUILD_LAPACK_EXAMPLES, test X${ADDON_lapack_EXAMPLES} = Xyes)

# ------------------------------------------------------------------------
# NOX-Epetra library?
# ------------------------------------------------------------------------

NOX_INTERFACE(epetra)
AM_CONDITIONAL(BUILD_LIBNOXEPETRA, test X${ADDON_epetra} = Xyes)
AM_CONDITIONAL(BUILD_EPETRA_EXAMPLES, test X${ADDON_epetra_EXAMPLES} = Xyes)

# ------------------------------------------------------------------------
# NOX-PETSC library?
# ------------------------------------------------------------------------

NOX_INTERFACE(petsc)
AM_CONDITIONAL(BUILD_LIBNOXPETSC, test X${ADDON_petsc} = Xyes)
AM_CONDITIONAL(BUILD_PETSC_EXAMPLES, test X${ADDON_petsc_EXAMPLES} = Xyes)

# ------------------------------------------------------------------------
# NOX-TSF library?
# ------------------------------------------------------------------------

NOX_INTERFACE(tsf)
AM_CONDITIONAL(BUILD_LIBNOXTSF, test X${ADDON_tsf} = Xyes)
AM_CONDITIONAL(BUILD_TSF_EXAMPLES, test X${ADDON_tsf_EXAMPLES} = Xyes)

# ------------------------------------------------------------------------
# NOX-Belos library?
# ------------------------------------------------------------------------

NOX_INTERFACE(belos)
AM_CONDITIONAL(BUILD_LIBNOXBELOS, test X${ADDON_belos} = Xyes)
AM_CONDITIONAL(BUILD_BELOS_EXAMPLES, test X${ADDON_belos_EXAMPLES} = Xyes)

# ------------------------------------------------------------------------
# LOCA Library Check
# ------------------------------------------------------------------------
# If inside trilinos, the default is to build loca.  If not,
# the default is not to build it.
# We also must turn off loca builds if inside trilinos and
# --disable-default-packages is set.
AC_ARG_ENABLE([default-packages],
AC_HELP_STRING([--enable-default-packages, flag in trilinos to Turn off default packages including nox and loca]),
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
[AC_HELP_STRING([--enable-loca],[compile the loca libraries])],
[BUILD_LIBLOCA=${enableval}],
[BUILD_LIBLOCA=${LOCA_DEFAULT}]
)
AC_MSG_CHECKING(whether we should build the loca library)
AC_MSG_RESULT([${BUILD_LIBLOCA}])
AM_CONDITIONAL(BUILD_LIBLOCA, test X${BUILD_LIBLOCA} = Xyes)

# ------------------------------------------------------------------------
# NOX/LOCA-Python Check.
# ------------------------------------------------------------------------
AC_ARG_ENABLE(python,
 [AC_HELP_STRING([--enable-python],[enable Python support])],
 [ADDON_python=$enableval],
 [ADDON_python=no]
)

AC_MSG_CHECKING(whether we should build the python support)
if test X${ADDON_python} = Xyes; then
  BUILD_PYTHON=yes
  AC_MSG_RESULT(yes)
  AM_CONDITIONAL(BUILD_PYTHON, true)

  # Get python include/lib paths
  AC_PYTHON_DEVEL

  # Check that Python.h is available
  save_CPPFLAGS=$CPPFLAGS
  CPPFLAGS="$save_CPPFLAGS $PYTHON_CPPFLAGS"
  AC_LANG([C++])
  AC_CHECK_HEADER(
  [Python.h],
  break,
  AC_MSG_ERROR([You must have Python.h in order to build the Python support!!]))
  CPPFLAGS="$save_CPPFLAGS"
else
  AC_MSG_RESULT(no)
  BUILD_PYTHON=no
  AM_CONDITIONAL(BUILD_PYTHON, false)
fi

# Determine which addon libraries to support in Python
IMPORT_LAPACK=0
if test X${BUILD_PYTHON} = Xyes && test X${ADDON_lapack} = Xyes; then
  IMPORT_LAPACK=1
fi
AC_SUBST(IMPORT_LAPACK)

IMPORT_EPETRA=0
if test X${BUILD_PYTHON} = Xyes && test X${ADDON_epetra} = Xyes; then
  IMPORT_EPETRA=1

  # Check for Numeric when building the Epetra module
  AC_PYTHON_MODULE(Numeric,yes)
fi
AC_SUBST(IMPORT_EPETRA)

# ------------------------------------------------------------------------
# Configure libtool
# ------------------------------------------------------------------------
# The default is to not build LOCA as a shared library, but we use libtool
# either way.  If Python support is enabled, then we must build LOCA
# as a shared library
if test X${BUILD_LIBLOCA} = Xyes || test X${BUILD_PYTHON} = Xyes; then

  # Disable shared libraries by default
  AC_DISABLE_SHARED

  # Enable shared libraries if Python is enabled
  if test X${BUILD_PYTHON} = Xyes; then
    enable_shared=yes
  fi

  # Configure libtool
  AC_PROG_LIBTOOL

fi

# ------------------------------------------------------------------------
# Configure swig
# ------------------------------------------------------------------------

AC_ARG_WITH(swig,
	    [AC_HELP_STRING([--with-swig],[enable swig and set swig binary])],
	    [WITH_SWIG=yes
	     if test X${withval} != Xyes ; then
	     	SWIG=$withval
	     fi],
            [WITH_SWIG=no])
AC_MSG_CHECKING(whether we should enable swig support)
AC_MSG_RESULT([${WITH_SWIG}])
AM_CONDITIONAL(HAVE_SWIG,test X${WITH_SWIG} = Xyes)
AC_SUBST(GNU_HAVE_SWIG, ${WITH_SWIG})

if test X${BUILD_PYTHON} = Xyes && test X${WITH_SWIG} = Xyes; then
  AC_PROG_SWIG(1.3.23)
  SWIG_ENABLE_CXX
  SWIG_MULTI_MODULE_SUPPORT
  SWIG_PYTHON
fi

# ------------------------------------------------------------------------
# Multivector check
# ------------------------------------------------------------------------

REQUIRE_TEUCHOS=no

AC_ARG_ENABLE(nox-multivectors,
 [AC_HELP_STRING([--enable-nox-multivectors],[enable NOX multivector support])],
 [BUILD_NOX_MULTIVECS=$enableval],
 [BUILD_NOX_MULTIVECS=no]
)

# If LOCA is enabled, must enable multivector support
if test X${BUILD_LIBLOCA} = Xyes; then
   BUILD_NOX_MULTIVECS=yes
fi

# If NOX-Belos support is enabled, must enable multivector support
if test X${BUILD_LIBNOXBELOS} = Xyes; then
   BUILD_NOX_MULTIVECS=yes
fi

# If multivec support is enabled, must have teuchos
if test X${BUILD_NOX_MULTIVECS} = Xyes; then
   REQUIRE_TEUCHOS=yes
   HAVE_NOX_MULTIVECS_SUBST=HAVE_NOX_MULTIVECS
   AC_DEFINE(HAVE_NOX_MULTIVECS,1,[define that multivector support is activated in nox])
else
   HAVE_NOX_MULTIVECS_SUBST=
fi

AC_MSG_CHECKING(whether we should enable multivector support)
AC_MSG_RESULT(${BUILD_NOX_MULTIVECS})
AM_CONDITIONAL(BUILD_NOX_MULTIVECS, test X${BUILD_NOX_MULTIVECS} = Xyes)
AC_SUBST(HAVE_NOX_MULTIVECS_SUBST)

# ------------------------------------------------------------------------
# Teuchos::any check
# ------------------------------------------------------------------------

AC_ARG_ENABLE(nox-any,
 [AC_HELP_STRING([--enable-nox-any],[enable Teuchos::Any support in NOX])],
 [BUILD_NOX_ANY=$enableval],
 [BUILD_NOX_ANY=no]
)



# If LOCA is enabled, must enable any support
if test X${BUILD_LIBLOCA} = Xyes; then
   BUILD_NOX_ANY=yes
fi

# If TSF is enabled, must enable any and Teuchos support
if test X${ADDON_tsf} = Xyes; then
  BUILD_NOX_ANY=yes
   REQUIRE_TEUCHOS=yes
   HAVE_NOX_ANY_SUBST=HAVE_NOX_ANY
   AC_DEFINE(HAVE_NOX_ANY,1,[define that Techos::any support is activated in nox])
else
   HAVE_NOX_ANY_SUBST=
fi

# If any support is enabled, must have teuchos
if test X${BUILD_NOX_ANY} = Xyes; then
   REQUIRE_TEUCHOS=yes
   HAVE_NOX_ANY_SUBST=HAVE_NOX_ANY
   AC_DEFINE(HAVE_NOX_ANY,1,[define that Techos::any support is activated in nox])
else
   HAVE_NOX_ANY_SUBST=
fi

AC_MSG_CHECKING(whether we should enable Teuchos::any support)
AC_MSG_RESULT(${BUILD_NOX_ANY})
AM_CONDITIONAL(BUILD_NOX_ANY, test X${BUILD_NOX_ANY} = Xyes)
AC_SUBST(HAVE_NOX_ANY_SUBST)

# ------------------------------------------------------------------------
# LOCA-LAPACK Library Check
# ------------------------------------------------------------------------

AC_MSG_CHECKING(whether we should build the loca-lapack interface library)
if test X${BUILD_LIBLOCA} = Xyes && test X${ADDON_lapack} = Xyes; then
  AC_MSG_RESULT(yes)
  AM_CONDITIONAL(BUILD_LIBLOCALAPACK, true)
else
  AC_MSG_RESULT(no)
  AM_CONDITIONAL(BUILD_LIBLOCALAPACK, false)
fi

# ------------------------------------------------------------------------
# Checks for building loca lapack examples
# ------------------------------------------------------------------------

AC_MSG_CHECKING(whether we should build the loca-lapack examples)
if test X${BUILD_LIBLOCA} = Xyes && test X${ADDON_lapack_EXAMPLES} = Xyes; then
  AC_MSG_RESULT(yes)
  AM_CONDITIONAL(BUILD_LOCA_LAPACK_EXAMPLES, true)
else
  AC_MSG_RESULT(no)
  AM_CONDITIONAL(BUILD_LOCA_LAPACK_EXAMPLES, false)
fi

# ------------------------------------------------------------------------
# LOCA-Epetra Library Check
# ------------------------------------------------------------------------

AC_MSG_CHECKING(whether we should build the loca-epetra interface library)
if test X${BUILD_LIBLOCA} = Xyes && test X${ADDON_epetra} = Xyes; then
  AC_MSG_RESULT(yes)
  AM_CONDITIONAL(BUILD_LIBLOCAEPETRA, true)
else
  AC_MSG_RESULT(no)
  AM_CONDITIONAL(BUILD_LIBLOCAEPETRA, false)
fi

# ------------------------------------------------------------------------
# Checks for building loca epetra examples
# ------------------------------------------------------------------------

AC_MSG_CHECKING(whether we should build the loca-epetra examples)
if test X${BUILD_LIBLOCA} = Xyes && test X${ADDON_epetra_EXAMPLES} = Xyes; then
  AC_MSG_RESULT(yes)
  AM_CONDITIONAL(BUILD_LOCA_EPETRA_EXAMPLES, true)
else
  AC_MSG_RESULT(no)
  AM_CONDITIONAL(BUILD_LOCA_EPETRA_EXAMPLES, false)
fi

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
      BUILD_ANY_EXAMPLES=yes
    fi

  ])

  # See if testing is enabled
  # We
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

  # python testing
  NOX_TESTING(python)
  AM_CONDITIONAL(BUILD_PYTHON_TESTS, test X${ADDON_python_TESTS} = Xyes)
  AM_CONDITIONAL(BUILD_LOCA_PYTHON_TESTS, test X${ADDON_loca_python_TESTS} = Xyes)
  AC_SUBST(ADDON_python_TESTS)
  AC_SUBST(ADDON_loca_python_TESTS)

# ------------------------------------------------------------------------
# LOCA-Anasazi support
# ------------------------------------------------------------------------

AC_ARG_WITH(loca-anasazi,
[AC_HELP_STRING([--with-loca-anasazi],[compile anasazi support into libloca])],
[BUILD_LOCA_ANASAZI_SUPPORT=${withval}],
[BUILD_LOCA_ANASAZI_SUPPORT=no]
)
AC_MSG_CHECKING(whether we should build anasazi support into loca library)
AC_MSG_RESULT([${BUILD_LOCA_ANASAZI_SUPPORT}])
AM_CONDITIONAL(BUILD_LOCAANASAZI, test X${BUILD_LOCA_ANASAZI_SUPPORT} = Xyes)

# ------------------------------------------------------------------------
# LOCA-MF support
# ------------------------------------------------------------------------

AC_ARG_WITH(loca-mf,
[AC_HELP_STRING([--with-loca-mf],[compile MF support into libloca])],
[BUILD_LOCA_MF_SUPPORT=${withval}],
[BUILD_LOCA_MF_SUPPORT=no]
)
AC_MSG_CHECKING(whether we should build MF support into loca library)
AC_MSG_RESULT([${BUILD_LOCA_MF_SUPPORT}])
AM_CONDITIONAL(BUILD_LIBLOCAMF, test X${BUILD_LOCA_MF_SUPPORT} = Xyes)

# ------------------------------------------------------------------------
# User defined directories 
# ------------------------------------------------------------------------

# enable use of --with-libdirs="-Llibdir1 -Llibdir2 ..." to prepend to LDFLAGS
TAC_ARG_WITH_LIBDIRS
# enable use of --with-incdirs="-Lincdir1 -Lincdir2 ..." to prepend to CPPFLAGS
TAC_ARG_WITH_INCDIRS

# ------------------------------------------------------------------------
# Checks for basic libraries for examples and/or lapack support
# ------------------------------------------------------------------------

if test X${ADDON_lapack} = Xyes || test X${BUILD_ANY_EXAMPLES} = Xyes || test X${BUILD_LIBLOCA} = Xyes; then

  # Define F77_FUNC - used for correctly munging Fortran function calls
  AC_F77_WRAPPERS

fi

if test X${BUILD_ANY_EXAMPLES} = Xyes; then

  AC_SEARCH_LIBS(pow,[m],,AC_MSG_WARN(Cannot find math library))
  AC_CHECK_FUNCS([pow sqrt])

  dnl Checks for libraries.
  AC_F77_LIBRARY_LDFLAGS
  LIBS=" ${FLIBS} ${LIBS}"

  # Find valid BLAS
  ACX_BLAS([LIBS="${BLAS_LIBS} ${LIBS}"],[AC_MSG_ERROR([Cannot find BLAS library])])

  # Find valid LAPACK
  ACX_LAPACK([LIBS=" ${LAPACK_LIBS} ${LIBS}"],[AC_MSG_ERROR([Cannot find LAPACK library])])

  # Find valid LAPACK DGGEV function
  ACX_DGGEV

  # User defined libraries
  TAC_ARG_WITH_LIBS

fi

# ------------------------------------------------------------------------
# Trilinos include file check for building noxepetra library
# ------------------------------------------------------------------------

AC_ARG_WITH(trilinos-basedir,
[AC_HELP_STRING([--with-trilinos-basedir],[specify location of trilinos source base directory, e.g., ${home}/Trilinos])],
[trilinos_basedir=$withval])

AC_ARG_WITH(trilinos-builddir,
[AC_HELP_STRING([--with-trilinos-builddir],[specify location of trilinos build directory, e.g., ${home}/Trilinos/build])],
[trilinos_builddir=$withval])

AC_ARG_WITH(epetra-libdir,
[AC_HELP_STRING([--with-epetra-libdir],[specify location of epetra library])],
[epetra_LIBDIR=$withval])

AC_ARG_WITH(epetra-incdir,
[AC_HELP_STRING([--with-epetra-incdir],[specify location of all epetra header files])],
[epetra_INCDIR=$withval])

AC_ARG_WITH(epetra-libdir,
[AC_HELP_STRING([--with-epetra-libdir],[specify location of epetra library])],
[epetra_LIBDIR=$withval])

AC_ARG_WITH(nox-epetraext,
[AC_HELP_STRING([--with-nox-epetraext],[compile epetraext support into the noxepetra library])],
[BUILD_NOX_EPETRAEXT_SUPPORT=${withval}],
[BUILD_NOX_EPETRAEXT_SUPPORT=no]
)
AC_MSG_CHECKING(whether we should build epetraext support into the noxepetra library)
AC_MSG_RESULT([${BUILD_NOX_EPETRAEXT_SUPPORT}])

AC_ARG_WITH(epetraext-incdir,
[AC_HELP_STRING([--with-epetraext-incdir],[specify location of epetraext header files])],
[epetraext_INCDIR=$withval])

AC_ARG_WITH(epetraext-libdir,
[AC_HELP_STRING([--with-epetraext-libdir],[specify location of epetraext library])],
[epetraext_LIBDIR=$withval])

AC_ARG_WITH(nox-epetra-ml,
[AC_HELP_STRING([--with-nox-epetra-ml],[compile ml support into the noxepetra library])],
[BUILD_NOX_EPETRA_ML_SUPPORT=${withval}],
[BUILD_NOX_EPETRA_ML_SUPPORT=no]
)
AC_MSG_CHECKING(whether we should build ml support into the noxepetra library)
AC_MSG_RESULT([${BUILD_NOX_EPETRA_ML_SUPPORT}])

AC_ARG_WITH(ml-incdir,
[AC_HELP_STRING([--with-ml-incdir],[specify location of ml header files])],
[ml_INCDIR=$withval])

AC_ARG_WITH(ml-libdir,
[AC_HELP_STRING([--with-ml-libdir],[specify location of ml library])],
[ml_LIBDIR=$withval])

AC_ARG_WITH(aztecoo-incdir,
[AC_HELP_STRING([--with-aztecoo-incdir],[specify location of all aztecoo header files])],
[aztecoo_INCDIR=$withval])

AC_ARG_WITH(aztecoo-libdir,
[AC_HELP_STRING([--with-aztecoo-libdir],[specify location of aztecoo library])],
[aztecoo_LIBDIR=$withval])

AC_ARG_WITH(ifpack-incdir,
[AC_HELP_STRING([--with-ifpack-incdir],[specify location of all ifpack header files])],
[ifpack_INCDIR=$withval])

AC_ARG_WITH(ifpack-libdir,
[AC_HELP_STRING([--with-ifpack-libdir],[specify location of ifpack library])],
[ifpack_LIBDIR=$withval])

AC_ARG_WITH(teuchos-incdir,
[AC_HELP_STRING([--with-teuchos-incdir],[specify location of all teuchos header files])],
[teuchos_INCDIR=$withval])

AC_ARG_WITH(teuchos-libdir,
[AC_HELP_STRING([--with-teuchos-libdir],[specify location of teuchos library])],
[teuchos_LIBDIR=$withval])

AC_ARG_WITH(anasazi-incdir,
[AC_HELP_STRING([--with-anasazi-incdir],[specify location of all anasazi header files])],
[anasazi_INCDIR=$withval])

AC_ARG_WITH(anasazi-libdir,
[AC_HELP_STRING([--with-anasazi-libdir],[specify location of anasazi library])],
[anasazi_LIBDIR=$withval])

AC_ARG_WITH(belos-incdir,
[AC_HELP_STRING([--with-belos-incdir],[specify location of all belos header files])],
[belos_INCDIR=$withval])

AC_ARG_WITH(belos-libdir,
[AC_HELP_STRING([--with-belos-libdir],[specify location of belos library])],
[belos_LIBDIR=$withval])

AC_ARG_WITH(mf-incdir,
[AC_HELP_STRING([--with-mf-incdir],[specify location of all MF header files])],
[MF_INCDIR=$withval])

AC_ARG_WITH(mf-libdir,
[AC_HELP_STRING([--with-mf-libdir],[specify location of MF library])],
[MF_LIBDIR=$withval])

# Relative path to Trilinos include files, to be used in Makefiles
REL_TRILINOS_CPPFLAGS=

# Absolute path to Trilinos include files from location that configure
# is run - to be used in configuration tests
ABS_TRILINOS_CPPFLAGS=

# Determine if anything needs to be added to REL/ABS_TRILINOS_FLAGS.
# Does not modify CPPFLAGS.
# TRILINOS_HEADER now takes three arguments:
# TRILINOS_HEADER(<library name>,<file name>,<subdirectory>)
# This is needed for multiple libraries being supplied by the
# same trilinos package (anasazi).
AC_DEFUN([TRILINOS_HEADER],
[
  save_CPPFLAGS=$CPPFLAGS
  done=no

  # First try the directories we already have
  CPPFLAGS="$save_CPPFLAGS $ABS_TRILINOS_CPPFLAGS"
  AC_PREPROC_IFELSE([AC_LANG_SOURCE(
  [[
#undef PACKAGE
#undef PACKAGE_NAME
#undef PACKAGE_BUGREPORT
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION
#include "$2"
  ]])],[done=yes],[])
  CPPFLAGS=$save_CPPFLAGS

  # Next  try the specified trilinos base directory (if any)
  if test $done = no && test -n "${trilinos_basedir}"; then
    CPPFLAGS="-I${trilinos_basedir}/packages/$1/$3 $save_CPPFLAGS $ABS_TRILINOS_CPPFLAGS"
    AC_PREPROC_IFELSE([AC_LANG_SOURCE(
  [[
#undef PACKAGE
#undef PACKAGE_NAME
#undef PACKAGE_BUGREPORT
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION
#include "$2"
  ]])],[ done=yes
      REL_TRILINOS_CPPFLAGS="-I${trilinos_basedir}/packages/$1/$3 ${REL_TRILINOS_CPPFLAGS}"
      ABS_TRILINOS_CPPFLAGS="-I${trilinos_basedir}/packages/$1/$3 ${ABS_TRILINOS_CPPFLAGS}"
    ],[])
    CPPFLAGS=$save_CPPFLAGS
  fi

  # Next  try the specified trilinos build directory (if any)
  if test $done = no && test -n "${trilinos_builddir}"; then
    CPPFLAGS="-I${trilinos_builddir}/packages/$1/$3 $save_CPPFLAGS $ABS_TRILINOS_CPPFLAGS"
    AC_PREPROC_IFELSE([AC_LANG_SOURCE(
  [[
#undef PACKAGE
#undef PACKAGE_NAME
#undef PACKAGE_BUGREPORT
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION
#include "$2"
  ]])],[ done=yes
      REL_TRILINOS_CPPFLAGS="-I${trilinos_builddir}/packages/$1/$3 ${REL_TRILINOS_CPPFLAGS}"
      ABS_TRILINOS_CPPFLAGS="-I${trilinos_builddir}/packages/$1/$3 ${ABS_TRILINOS_CPPFLAGS}"
    ],[])
    CPPFLAGS=$save_CPPFLAGS
  fi

  # Next  try the specified incdir (if any)
  if test $done = no && test -n "${$1_INCDIR}"; then
    CPPFLAGS="-I${$1_INCDIR} $save_CPPFLAGS $ABS_TRILINOS_CPPFLAGS"
    AC_PREPROC_IFELSE([AC_LANG_SOURCE(
  [[
#undef PACKAGE
#undef PACKAGE_NAME
#undef PACKAGE_BUGREPORT
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION
#include "$2"
  ]])],[ done=yes
      REL_TRILINOS_CPPFLAGS="-I${$1_INCDIR} ${REL_TRILINOS_CPPFLAGS}"
      ABS_TRILINOS_CPPFLAGS="-I${$1_INCDIR} ${ABS_TRILINOS_CPPFLAGS}"
    ],[])
    CPPFLAGS=$save_CPPFLAGS
  fi

  # Next try the Trilinos build dir
  if test $done = no; then
    CPPFLAGS="-I../$1/$3 $save_CPPFLAGS $ABS_TRILINOS_CPPFLAGS"
    AC_PREPROC_IFELSE([AC_LANG_SOURCE(
  [[
#undef PACKAGE
#undef PACKAGE_NAME
#undef PACKAGE_BUGREPORT
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION
#include "$2"
  ]])],[ done=yes
      REL_TRILINOS_CPPFLAGS="-I\$(top_builddir)/../$1/$3 ${REL_TRILINOS_CPPFLAGS}"
      ABS_TRILINOS_CPPFLAGS="-I../$1/$3 ${ABS_TRILINOS_CPPFLAGS}"
    ],[])
    CPPFLAGS=$save_CPPFLAGS
  fi

  # Next try the Trilinos source dir
  if test $done = no && test -n "$srcdir"; then
    CPPFLAGS="-I$srcdir/../$1/$3 $save_CPPFLAGS $ABS_TRILINOS_CPPFLAGS"
    AC_PREPROC_IFELSE([AC_LANG_SOURCE(
  [[
#undef PACKAGE
#undef PACKAGE_NAME
#undef PACKAGE_BUGREPORT
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION
#include "$2"
  ]])],[ done=yes
      REL_TRILINOS_CPPFLAGS="-I\$(top_srcdir)/../$1/$3 ${REL_TRILINOS_CPPFLAGS}"
      ABS_TRILINOS_CPPFLAGS="-I$srcdir/../$1/$3 ${ABS_TRILINOS_CPPFLAGS}"
    ],[])
    CPPFLAGS=$save_CPPFLAGS
  fi

  if test $done = no; then
    echo "------"
    echo "Cannot preprocess the include file $2, part of the $1 package."
    echo "Try --with-$1-incdir=<dir> to specify a single directory or"
    echo "--with-cppflags=\"-I<dir1> -I<dir2>\" for multiple directories."
    echo "------"
    AC_MSG_ERROR([Cannot find $1 include directory containing $2])
  fi
])

if test X${ADDON_epetra} = Xyes; then

  # We have to check for the standard Trilinos headers as well as
  # those created by configure.
  TRILINOS_HEADER(epetra,Epetra_config.h,src)
  TRILINOS_HEADER(epetra,Epetra_DataAccess.h,src)
  TRILINOS_HEADER(aztecoo,AztecOO_config.h,src)
  TRILINOS_HEADER(aztecoo,AztecOO_StatusType.h,src)
  TRILINOS_HEADER(ifpack,Ifpack_config.h,src)
  TRILINOS_HEADER(ifpack,Ifpack_ScalingType.h,src)

  REQUIRE_TEUCHOS=yes

  if test X${BUILD_NOX_EPETRAEXT_SUPPORT} = Xyes; then
    TRILINOS_HEADER(epetraext,EpetraExt_config.h,src)
    TRILINOS_HEADER(epetraext,EpetraExt_ConfigDefs.h,src)
    TRILINOS_HEADER(epetraext,EpetraExt_Transform.h,src/transform)
    TRILINOS_HEADER(epetraext,EpetraExt_BlockUtility.h,src/block)
    TRILINOS_CXXFLAGS="${REL_TRILINOS_CPPFLAGS}"
  fi

  if test X${BUILD_NOX_EPETRAEXT_SUPPORT} = Xyes; then
    TRILINOS_HEADER(epetraext,EpetraExt_MapColoring.h,src/coloring)
  fi

  if test X${BUILD_NOX_EPETRA_ML_SUPPORT} = Xyes; then
    TRILINOS_HEADER(ml,ml_config.h,src)
    TRILINOS_HEADER(ml,ml_common.h,src/Include)
    TRILINOS_HEADER(ml,ml_memory.h,src/Utils)
    TRILINOS_HEADER(ml,ml_comm.h,src/Comm)
    TRILINOS_HEADER(ml,ml_op_utils.h,src/Operator)
    TRILINOS_HEADER(ml,ml_check.h,src/Coarsen)
    TRILINOS_HEADER(ml,ml_mapper.h,src/FEGrid)
    TRILINOS_HEADER(ml,ml_solver.h,src/Smoother)
    TRILINOS_HEADER(ml,ml_viz_stats.h,src/Main)
    TRILINOS_HEADER(ml,ml_krylov.h,src/Krylov)
    TRILINOS_HEADER(triutils,Trilinos_Util_CommandLineParser.h,src)
    TRILINOS_CXXFLAGS="${REL_TRILINOS_CPPFLAGS}"
  fi

  TRILINOS_CXXFLAGS="${REL_TRILINOS_CPPFLAGS}"

fi

if test X${REQUIRE_TEUCHOS} = Xyes; then
  TRILINOS_HEADER(teuchos,Teuchos_config.h,src)
  TRILINOS_HEADER(teuchos,Teuchos_Flops.hpp,src)
  TRILINOS_CXXFLAGS="${REL_TRILINOS_CPPFLAGS}"
fi

if test X${ADDON_belos} = Xyes; then
   TRILINOS_HEADER(belos,Belos_config.h,src)
   TRILINOS_HEADER(belos,BelosTypes.hpp,src)
   TRILINOS_CXXFLAGS="${REL_TRILINOS_CPPFLAGS}"
fi

if test X${BUILD_LOCA_ANASAZI_SUPPORT} = Xyes; then
  TRILINOS_HEADER(anasazi,Anasazi_config.h,src)
  TRILINOS_HEADER(anasazi,AnasaziReturnType.hpp,src)
  TRILINOS_CXXFLAGS="${REL_TRILINOS_CPPFLAGS}"
fi

if test X${BUILD_PYTHON} = Xyes && test X${ADDON_epetra} = Xyes; then
# Commenting out the check for PyTrilinos_config.h for right now since
# PyTrilinos is configured after NOX, so this file may not exist
#  TRILINOS_HEADER(PyTrilinos,PyTrilinos_config.h,src)
  TRILINOS_HEADER(PyTrilinos,PyTrilinos.h,src)
  TRILINOS_CXXFLAGS="${REL_TRILINOS_CPPFLAGS}"
fi

AC_SUBST(MF_INCDIR)

AC_SUBST(TRILINOS_CXXFLAGS)

# ------------------------------------------------------------------------
# Checks for Trilinos libraries
# ------------------------------------------------------------------------

# TRILINOS_LIB_CHECK now takes three arguments:
# TRILINOS_LIB_CHECK(<library name>,<package name>,<subdirectory>)
# This is needed for multiple libraries being supplied by the
# same trilinos package (anasazi).
AC_DEFUN([TRILINOS_LIB_CHECK],
[
  LIBS="-l$1 ${LIBS}"
  save_LDFLAGS=$LDFLAGS
  done=no

  # First, see if we're inside Trilinos.
  # If so, just assume the libraries will be built later.
  if test X$INSIDE_TRILINOS = Xyes; then
    PWD=`pwd`
    LDFLAGS="-L$PWD/../$2/$3 $save_LDFLAGS"
    done=auto
  fi

  # Next, see if we can link without doing anything
  if test $done = no; then
    AC_LANG([C++])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([[]], [[]])],[done=yes],[])
    LDFLAGS=$save_LDFLAGS
  fi

  # Next try the specified Trilinos base dir (if any)
   if test $done = no && test -n "${trilinos_basedir}"; then
     LDFLAGS="-L${trilinos_basedir}/packages/$1/$3 $save_LDFLAGS"
     AC_LINK_IFELSE([AC_LANG_PROGRAM([[]], [[]])],[done=yes],[LDFLAGS=$save_LDFLAGS])
   fi

  # Next try the specified Trilinos build dir (if any)
   if test $done = no && test -n "${trilinos_builddir}"; then
     LDFLAGS="-L${trilinos_builddir}/packages/$1/$3 $save_LDFLAGS"
     AC_LINK_IFELSE([AC_LANG_PROGRAM([[]], [[]])],[done=yes],[LDFLAGS=$save_LDFLAGS])
   fi

  # Next try the specified libdir (if any)
   if test $done = no && test -n "${$1_LIBDIR}"; then
     LDFLAGS="-L${$1_LIBDIR} $save_LDFLAGS"
     AC_LINK_IFELSE([AC_LANG_PROGRAM([[]], [[]])],[done=yes],[LDFLAGS=$save_LDFLAGS])
   fi

  # Next try the Trilinos build dir
   if test $done = no; then
    PWD=`pwd`
    LDFLAGS="-L$PWD/../$2/$3 $save_LDFLAGS"
    AC_LINK_IFELSE([AC_LANG_PROGRAM([[]], [[]])],[done=yes],[LDFLAGS=$save_LDFLAGS])
   fi

  AC_MSG_CHECKING([for $1 library])

  if test $done = auto; then
    AC_MSG_RESULT(no check - assuming library will be built by Trilinos)
  elif test $done = yes; then
    AC_MSG_RESULT(yes)
  else
    AC_MSG_RESULT(no)
    echo "------"
    echo "Cannot find the library $1."
    echo "Try --with-$1-libdir=<dir> to specify a single directory or"
    echo "--with-ldflags=\"-L<dir1> -L<dir2>\" for multiple directories."
    echo "------"
    AC_MSG_ERROR([Cannot find $1 library.])
  fi

])

#
# We have to be careful about the order we do these checks in, since this
# determines the order in which the libraries appear on the link line.
# TRILINOS_LIB_CHECK prepends each library to LIBS, so we have to do the
# checks in the opposite order
#

if test X${BUILD_NOX_EPETRA_ML_SUPPORT} = Xyes ||
   test X${ADDON_epetra_EXAMPLES} = Xyes ||
   test X${ADDON_epetra_TESTS} = Xyes; then
  REQUIRE_TEUCHOS=yes
fi

if test X${REQUIRE_TEUCHOS} = Xyes; then
  TRILINOS_LIB_CHECK(teuchos,teuchos,src)
  AC_DEFINE(HAVE_NOX_TEUCHOS,1,[define that Teuchos is activated in nox])
fi

if test X${ADDON_epetra_EXAMPLES} = Xyes ||
   test X${ADDON_epetra_TESTS} = Xyes; then

  TRILINOS_LIB_CHECK(epetra,epetra,src)
  TAC_ARG_ENABLE_CAN_USE_PACKAGE(amesos, [Test for enabled amesos.], NOX_AMESOS)
  if test "X$ac_cv_use_amesos" = "Xyes"; then
    TRILINOS_LIB_CHECK(amesos,amesos,src)
  fi
  TRILINOS_LIB_CHECK(ifpack,ifpack,src)
  TRILINOS_LIB_CHECK(aztecoo,aztecoo,src)

  if test X${BUILD_NOX_EPETRAEXT_SUPPORT} = Xyes; then
    TRILINOS_LIB_CHECK(epetraext,epetraext,src)
    AC_DEFINE(HAVE_NOX_EPETRAEXT,1,[define that epetraext is activated in nox])
  fi

fi

if test X${BUILD_NOX_EPETRA_ML_SUPPORT} = Xyes; then
  TRILINOS_LIB_CHECK(ml,ml,src)
  AC_DEFINE(HAVE_NOX_ML_EPETRA,1,[define that ml is activated in nox])
fi

if test X${BUILD_LOCA_ANASAZI_SUPPORT} = Xyes; then
  TRILINOS_LIB_CHECK(anasazi,anasazi,src)
  AC_DEFINE(HAVE_LOCA_ANASAZI,1,[define that anasazi is activated in loca])
fi

if test X${ADDON_belos} = Xyes; then
  TRILINOS_LIB_CHECK(belos,belos,src)
  AC_DEFINE(HAVE_NOX_BELOS,1,[define that Belos is activated in nox])
fi

AC_SUBST(MF_LIBDIR)

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
  elif test -f ${PETSC_DIR}/src/sys/src/mpiuni/$2; then
    $1_INC=${PETSC_DIR}/src/sys/src/mpiuni
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

  AC_MSG_CHECKING(whether PETSC_BOPT is defined)
  if test -n "${PETSC_BOPT}"; then
     AC_MSG_RESULT([yes (${PETSC_BOPT})])
  else
     AC_MSG_ERROR([Cannot find PETSC_BOPT environment variable. Please define.])  fi

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
  elif test -f ${PETSC_DIR}/lib/lib${PETSC_BOPT}/${PETSC_ARCH}/lib$1.a; then
    LDFLAGS="-L${PETSC_DIR}/lib/lib${PETSC_BOPT}/${PETSC_ARCH} ${LDFLAGS}"
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
# Checks for types
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for structures
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for libraries
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for compiler characteristics
# ------------------------------------------------------------------------

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
[AC_HELP_STRING([--enable-prerelease],[Enable compilation and documentation of prerelease code.])],
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
# Checks for library functions
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for system services
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Perform substitutions in output files
# ------------------------------------------------------------------------

AC_SUBST(ac_aux_dir)
AC_SUBST(abs_builddir)
AC_SUBST(abs_top_builddir)

AC_SUBST(MPI_LIBS)
AC_SUBST(LIBTOOL_COMPILER_FLAGS)
AC_SUBST(LIBTOOL_LINKER_FLAGS)
AC_SUBST(LIBTOOL_CCLINKER_FLAGS)

# ------------------------------------------------------------------------
# Output files
# ------------------------------------------------------------------------

AC_CONFIG_FILES([
		Makefile
		doc/Makefile
		doc/Doxyfile
		src/Makefile
                src-lapack/Makefile
		example-lapack/Makefile
		src-tsf/Makefile
		example-tsf/Makefile
		src-epetra/Makefile
		example-epetra/Makefile
		example-epetra/1Dfem_linear/Makefile
		example-epetra/1Dfem_nonlinear/Makefile
		example-epetra/1Dfem_nonlinearColoring/Makefile
		example-epetra/1Dtransfem_nonlinear/Makefile
		example-epetra/Brusselator_Split/Makefile
		example-epetra/Brusselator_EpetraNew/Makefile
                src-petsc/Makefile
                example-petsc/Makefile
                example-petsc/DS6.5.1/Makefile
                example-petsc/1Dfd_nonlinear/Makefile
		src-belos/Makefile
		example-belos/Makefile
		example-belos/epetra/Makefile
		example-belos/epetra/1Dfem_nonlinear/Makefile
		python/Makefile
		python/GNUmakefile
		python/__init__.py
		src-loca/Makefile
		src-loca/doc/Makefile
		src-loca/src/Makefile
                src-loca/src-epetra/Makefile
		src-loca/src-lapack/Makefile
		src-loca/src-mf/Makefile
                src-loca/example-epetra/Makefile
                src-loca/example-epetra/1Dfem/Makefile
                src-loca/example-epetra/Tcubed/Makefile
                src-loca/example-epetra/Tcubed_epetra_new/Makefile
                src-loca/example-epetra/Brusselator_xyzt/Makefile
		src-loca/example-lapack/Makefile
		src-loca/example-lapack/Chan/Makefile
		src-loca/example-lapack/1Dfd_pitchfork/Makefile
		src-loca/example-lapack/Brusselator/Makefile
		src-loca/python/Makefile
		src-loca/python/GNUmakefile
		src-loca/python/__init__.py
                test/Makefile
                test/scripts/Makefile
                test/scripts/daily/Makefile
                test/scripts/daily/serial/Makefile
		test/scripts/daily/serial/noxTestAll
                test/scripts/daily/mpi/Makefile
		test/scripts/daily/mpi/noxTestAllMpi
                test/utils/Makefile
                test/lapack/Makefile
		test/lapack/NOX/Makefile
                test/lapack/NOX_NewTest/Makefile
                test/lapack/NOX_Group/Makefile
		test/lapack/LOCA_TestProblems/Makefile
		test/lapack/LOCA_UnitTests/Makefile
                test/lapack/LOCA_ChanProblem/Makefile
		test/lapack/LOCA_python/Makefile
		test/lapack/LOCA_python/setpath.py
                test/epetra/Makefile
                test/epetra/NOX_NewTest/Makefile
                test/epetra/NOX_Vector/Makefile
                test/epetra/NOX_Group/Makefile
		test/epetra/LOCA_TestProblems/Makefile
		test/epetra/LOCA_TestProblems/Tcubed/Makefile
		test/epetra/LOCA_UnitTests/Makefile
                test/epetra/LOCA_NewTest/Makefile
                test/epetra/LOCA_Group/Makefile
		test/epetra/LOCA_Pitchfork/Makefile
                test/epetra/1Dfem/Makefile
                test/epetra/MultiPhysics/Makefile
                test/epetra/DS6.5.1/Makefile
                test/petsc/Makefile
		])

AC_OUTPUT


# ------------------------------------------------------------------------
# Bye World!
# ------------------------------------------------------------------------

echo "----------------------------------------"
echo "Finished Running NOX Configure Script"
echo "----------------------------------------"

