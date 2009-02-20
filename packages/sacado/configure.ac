## $Id$ 
## $Source$ 
# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# @HEADER
# ************************************************************************
# 
#                           Sacado Package
#                 Copyright (2006) Sandia Corporation
# 
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
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
# Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
# (etphipp@sandia.gov).
# 
# ************************************************************************
# @HEADER

# ------------------------------------------------------------------------
# Initialization 
# ------------------------------------------------------------------------

# This must be the first line in configure.ac.
# Optional 3rd argument is email address for bugs.
AC_INIT([sacado], [9.1-dev], [sacado-developers@software.sandia.gov])

# Hello World!
echo "----------------------------------------"
echo "Running Sacado Configure Script"
echo "----------------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.  Any file in that directory will do,
# preferably one that is unlikely to be removed or renamed.
AC_CONFIG_SRCDIR([src/Sacado_Version.hpp])

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.
AC_CONFIG_AUX_DIR(config)

# Configure should create src/Sacado_config.h from src/Sacado_config.h.in
AM_CONFIG_HEADER(src/Sacado_config.h:src/Sacado_config.h.in)

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

# ------------------------------------------------------------------------
# Check to see if MPI enabled and if any special configuration done
# ------------------------------------------------------------------------

TAC_ARG_CONFIG_MPI

# ------------------------------------------------------------------------
# Checks for programs
# ------------------------------------------------------------------------

AC_PROG_CC(cc gcc)
AC_PROG_CXX(CC g++ gcc c++ cxx)
#AC_PROG_F77(f77 g77 f90 xlf90)
AC_PROG_RANLIB
AC_PROG_AWK

# ------------------------------------------------------------------------
# Checks for extra flags
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

# ------------------------------------------------------------------------
# Checks for header files - used in Sacado_ConfigDefs.h 
# ------------------------------------------------------------------------

AC_LANG([C++])
AC_CHECK_HEADERS(
 [cstdlib cstdio cassert cmath cstring iomanip iostream string algorithm new vector valarray map iterator typeinfo],, 
 AC_MSG_ERROR([[Cannot find standard header (see above).]]))

# ------------------------------------------------------------------------
# Checks for Makefile.export related systems
# ------------------------------------------------------------------------
TAC_ARG_ENABLE_EXPORT_MAKEFILES(yes)

# ------------------------------------------------------------------------
# Checks for special package flags
# ------------------------------------------------------------------------

# Check for --enable-sacado-debug
AC_ARG_ENABLE([sacado-debug], 
              [AC_HELP_STRING([--enable-sacado-debug], 
                              [Enable extra debugging checks that may decrease performance (default is no)])],
              [ac_cv_use_sacado_debug=$enableval], 
              [ac_cv_use_sacado_debug=no])
AC_MSG_CHECKING(whether to use sacado-debug)
if test "X$ac_cv_use_sacado_debug" != "Xno"; then
  AC_MSG_RESULT(yes)  
  AC_DEFINE([SACADO_DEBUG],,[Define to enable extra debugging checks])
else
  AC_MSG_RESULT(no)
fi

TAC_ARG_ENABLE_CAN_USE_PACKAGE(
  sacado, teuchos, SACADO_TEUCHOS, ENABLE_SACADO_TEUCHOS, yes
  ,[Build in Teuchos support in Sacado.]
  ,[Build in Teuchos support in Sacado.  Can be overridden with --disable-sacado-teuchos.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(
  sacado, stokhos, SACADO_STOKHOS, ENABLE_SACADO_STOKHOS, yes
  ,[Build in Stokhos support in Sacado.]
  ,[Build in Stokhos support in Sacado.  Can be overridden with --disable-sacado-stokhos.])

# Check for --[en|dis]able-sacado-paramlib
case "${target_os}" in
  solaris*)
    sacado_paramlib_default=no
    ;;
  *)
    sacado_paramlib_default=yes
esac
AC_ARG_ENABLE([sacado-paramlib], 
              [AC_HELP_STRING([--enable-sacado-paramlib], 
                              [Enable building Sacado parameter library (default is yes except on Sun platform)])],
              [ac_cv_use_sacado_paramlib=$enableval], 
              [ac_cv_use_sacado_paramlib=$sacado_paramlib_default])
AC_MSG_CHECKING(whether to enable the sacado parameter library)
if test "X$ac_cv_use_sacado_paramlib" != "Xno"; then
  AC_MSG_RESULT(yes)  
else
  AC_MSG_RESULT(no)
fi
AM_CONDITIONAL(BUILD_PARAMLIB, test "X$ac_cv_use_sacado_paramlib" != "Xno")

# Check for --[en|dis]able-sacado-nestedfad
case "${target_os}" in
  solaris*)
    sacado_nested_fad_default=no
    ;;
  *)
    sacado_nested_fad_default=yes
esac
AC_ARG_ENABLE([sacado-nestedfad], 
              [AC_HELP_STRING([--enable-sacado-nestedfad], 
                              [Enable nested Fad inside of Fad types (default is yes except on Sun platform)])],
              [ac_cv_use_sacado_nested_fad=$enableval], 
              [ac_cv_use_sacado_nested_fad=$sacado_nested_fad_default])
AC_MSG_CHECKING(whether to enable nested fad types)
if test "X$ac_cv_use_sacado_nested_fad" != "Xno"; then
  AC_MSG_RESULT(yes)  
else
  AC_MSG_RESULT(no)
fi
AM_CONDITIONAL(BUILD_NESTED_FAD, test "X$ac_cv_use_sacado_nested_fad" != "Xno")

#
#  --enable-sacado-complex
#

TAC_ARG_ENABLE_FEATURE_SUB( sacado, complex, [Enable Sacado complex code], SACADO_COMPLEX, yes)
AM_CONDITIONAL(HAVE_SACADO_COMPLEX, [test "X$ac_cv_use_sacado_complex" != "Xno"])

# ------------------------------------------------------------------------
# Checks if tests and examples should be built
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_FEATURE(tests, [Make tests for all Trilinos packages buildable with 'make tests'], TESTS, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK(sacado, tests, [Make Sacado tests buildable with 'make tests'], SACADO_TESTS)
AM_CONDITIONAL(BUILD_TESTS, test "X$ac_cv_use_sacado_tests" != "Xno")

TAC_ARG_ENABLE_FEATURE(examples, [Make examples for all Trilinos packages buildable with 'make examples'], EXAMPLES, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK(sacado, examples, [Make Sacado examples buildable with 'make examples'], SACADO_EXAMPLES)
AM_CONDITIONAL(BUILD_EXAMPLES, test "X$ac_cv_use_sacado_examples" != "Xno")

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
# specify execincludedir
AC_ARG_WITH_EXECINCLUDEDIR

# enable use of --with-libdirs="-Llibdir1 -Llibdir2 ..." to prepend to LDFLAGS
TAC_ARG_WITH_LIBDIRS
# enable use of --with-incdirs="-Lincdir1 -Lincdir2 ..." to prepend to CPPFLAGS
TAC_ARG_WITH_INCDIRS

# #np# - Yet another opportunity to remove code if you aren't
# using Fortran
# Define F77_FUNC that will be used to link with Fortran subroutines.
#AC_F77_WRAPPERS

# ------------------------------------------------------------------------
# Checks for libraries
# ------------------------------------------------------------------------

# If tests, examples and libcheck are disabled, we don't have to check
# for these libraries.

if test "X$ac_cv_use_sacado_tests" != "Xno" || test "X$ac_cv_use_sacado_examples" != "Xno" || test "X$ac_cv_use_libcheck" != "Xno"; then

# checks for libraries now
dnl Replace `main' with a function in -lm:
AC_SEARCH_LIBS(pow,[m],,AC_MSG_ERROR(Cannot find math library))
AC_SEARCH_LIBS(sqrt,[m],,AC_MSG_ERROR(Cannot find math library))

# #np# If Lapack is also needed, see the  
# configure.ac file for epetra (Trilinos/packages/epetra)

fi
# end of the list of libraries that don't need to be checked for if
# tests and examples are disabled.

# Check for cppunit for tests
if test "X$ac_cv_use_sacado_tests" != "Xno"; then
  AM_PATH_CPPUNIT(1.10.2,[ac_cv_have_cppunit=yes],[ac_cv_have_cppunit=no])
  if test X${ac_cv_have_cppunit} = Xno; then
    AC_MSG_WARN([Cannot find cppunit and tests are enabled.  Disabling unit tests.  Configure with --with-cppunit-prefix=PATH to specify location of cppunit installation to enable unit tests])
  fi
fi
AM_CONDITIONAL(HAVE_CPPUNIT, test X${ac_cv_have_cppunit} = Xyes)

# Check for ADOL-C
AC_ARG_WITH(adolc, 
            [AC_HELP_STRING([--with-adolc=PATH],
                            [Enable ADOL-C support for Taylor polynomial tests, optionally providing the path to the base ADOL-C directory.  Default is no.])],
            [
              if test "X${withval}" != "Xyes"; then
                 ADOLC_PATH=${withval}
              fi
              ENABLE_ADOLC_SUPPORT=yes
            ],
            [ENABLE_ADOLC_SUPPORT=no]
            )
if test "X$ac_cv_use_sacado_tests" != "Xno"; then
  AC_MSG_CHECKING([whether we should enable ADOL-C support])
  AC_MSG_RESULT([${ENABLE_ADOLC_SUPPORT}])
fi
AM_CONDITIONAL(ENABLE_ADOLC, test X${ENABLE_ADOLC_SUPPORT} = Xyes)
if test X${ENABLE_ADOLC_SUPPORT} = Xyes; then
  AC_DEFINE([HAVE_ADOLC],,[Define if ADOL-C is enabled])
fi
  
AC_ARG_WITH(adolc-incdir,
            [AC_HELP_STRING([--with-adolc-incdir=PATH],
                            [Path to ADOL-C include directory])],
            [ADOLC_INCDIR=${withval}],
            [
              if test "X${ADOLC_PATH}" != "X"; then
                ADOLC_INCDIR="${ADOLC_PATH}/include"
              else
                ADOLC_INCDIR=
              fi
            ])

AC_ARG_WITH(adolc-libdir,
            [AC_HELP_STRING([--with-adolc-libdir=PATH],
                            [Path to ADOL-C library directory])],
            [ADOLC_LIBDIR=${withval}],
            [
              if test "X${ADOLC_PATH}" != "X"; then
                ADOLC_LIBDIR="${ADOLC_PATH}/lib"
              else
                ADOLC_LIBDIR=
              fi
            ])

AC_ARG_WITH(adolc-libs,
            [AC_HELP_STRING([--with-adolc-libs=LIBS],
            [ADOL-C libraries (e.g., "-ladaolc")])],
            [ADOLC_LIBS=${withval}],
            [ADOLC_LIBS="-ladolc"])

AC_SUBST(ADOLC_INCDIR)
AC_SUBST(ADOLC_LIBDIR)
AC_SUBST(ADOLC_LIBS)

# ------------------------------------------------------------------------
# Checks for compiler characteristics
# ------------------------------------------------------------------------

# Check for bug (e.g., with Sun's CC 5.7) in handling std:: within a namespace
AC_CXX_STDCCBUG

# ------------------------------------------------------------------------
# Checks for linker characteristics
# ------------------------------------------------------------------------

# Determine libraries needed for linking with Fortran
#AC_F77_LIBRARY_LDFLAGS


#  ----------------------------------------------------
#  Needs this because Dec handles templates differently
#  ----------------------------------------------------

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
AC_CONFIG_FILES([
		Makefile 
		Makefile.export.sacado
		src/Makefile
		Fad/Makefile
		example/Makefile
		test/Makefile
		test/utils/Makefile
		test/TestSuite/Makefile
		test/tradoptest/Makefile
		test/tradoptest/testgen.as
		test/performance/Makefile
		])

AC_OUTPUT()

# Bye World!
echo "----------------------------------------"
echo "Finished Running Sacado Configure Script"
echo "----------------------------------------"
