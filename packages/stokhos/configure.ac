## $Id$ 
## $Source$ 
# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# @HEADER
# ************************************************************************
# 
#                     Stokhos Package
#                 Copyright (2008) Sandia Corporation
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
# Questions? Contact Eric T. Phipps (etphipp@sandia.gov) 
# 
# ************************************************************************
# @HEADER

# ------------------------------------------------------------------------
# Initialization 
# ------------------------------------------------------------------------

# This must be the first line in configure.ac.
# Optional 3rd argument is email address for bugs.
AC_INIT(stokhos, 1.0d, etphipp@sandia.gov)

# Hello World!
echo "----------------------------------------"
echo "Running Stokhos Configure Script"
echo "----------------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.  Any file in that directory will do,
# preferably one that is unlikely to be removed or renamed.
AC_CONFIG_SRCDIR([src/Stokhos_Version.hpp])

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.
AC_CONFIG_AUX_DIR(config)

# Configure should create src/Stokhos_config.h from src/Stokhos_config.h.in
AM_CONFIG_HEADER(src/Stokhos_config.h:src/Stokhos_config.h.in)

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
#  - Required version of automake.
AM_INIT_AUTOMAKE(1.10 no-define tar-ustar)

# Specify required version of autoconf.
AC_PREREQ(2.61)

# ------------------------------------------------------------------------
# Check if Fortran is enabled or not
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_OPTION(fortran, [enable Fortran support], FORTRAN_SUPPORT, yes)
AM_CONDITIONAL(FORTRAN_ENABLED, [test "X$ac_cv_use_fortran" != "Xno"])
# ------------------------------------------------------------------------
# Check to see if MPI enabled and if any special configuration done
# ------------------------------------------------------------------------

TAC_ARG_CONFIG_MPI

# ------------------------------------------------------------------------
# Checks for programs
# ------------------------------------------------------------------------

AC_PROG_CC(cc gcc)
AC_PROG_CXX(CC g++ c++ cxx)
AC_PROG_F77(f77 g77 gfortran f90 xlf90 f95)
AC_PROG_RANLIB

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
# Checks for header files - used in Stokhos_ConfigDefs.h 
# ------------------------------------------------------------------------

AC_LANG([C++])
AC_CHECK_HEADERS(
 [cstdlib cstdio cassert cmath cstring iomanip iostream string vector new],, 
 AC_MSG_ERROR([[Cannot find standard header (see above).]]))

# ------------------------------------------------------------------------
# Checks for Makefile.export related systems
# ------------------------------------------------------------------------
TAC_ARG_ENABLE_EXPORT_MAKEFILES(yes)

# ------------------------------------------------------------------------
# Checks for special package flags
# ------------------------------------------------------------------------

# Check for --enable-stokhos-debug
AC_ARG_ENABLE([stokhos-debug], 
              [AC_HELP_STRING([--enable-stokhos-debug], 
                              [Enable extra debugging checks that may decrease performance (default is no)])],
              [ac_cv_use_stokhos_debug=$enableval], 
              [ac_cv_use_stokhos_debug=no])
AC_MSG_CHECKING(whether to use stokhos-debug)
if test "X$ac_cv_use_stokhos_debug" != "Xno"; then
  AC_MSG_RESULT(yes)  
  AC_DEFINE([STOKHOS_DEBUG],,[Define to enable extra debugging checks])
else
  AC_MSG_RESULT(no)
fi

# ------------------------------------------------------------------------
# Checks if tests and examples should be built
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_FEATURE(tests, [Make tests for all Trilinos packages buildable with 'make tests'], TESTS, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK( stokhos, tests, [Make Stokhos tests buildable with 'make tests'], STOKHOS_TESTS)
AM_CONDITIONAL(BUILD_TESTS, test "X$ac_cv_use_stokhos_tests" != "Xno")

TAC_ARG_ENABLE_FEATURE(examples, [Make examples for all Trilinos packages buildable with 'make examples'], EXAMPLES, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK( stokhos, examples, [Make Stokhos examples buildable with 'make examples'], STOKHOS_EXAMPLES)
AM_CONDITIONAL(BUILD_EXAMPLES, test "X$ac_cv_use_stokhos_examples" != "Xno")

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

# #np# - Yet another opportunity to remove code if you aren't
# using Fortran
# Define F77_FUNC that will be used to link with Fortran subroutines. - trash WORKGXX 
#AC_F77_WRAPPERS

# ------------------------------------------------------------------------
# Checks for libraries
# ------------------------------------------------------------------------

# If tests, examples and libcheck are disabled, we don't have to check
# for these libraries.

if test "X$ac_cv_use_stokhos_tests" != "Xno" || test "X$ac_cv_use_stokhos_examples" != "Xno" || test "X$ac_cv_use_libcheck" != "Xno"; then

# checks for libraries now
dnl Replace `main' with a function in -lm:
AC_SEARCH_LIBS(pow,[m],,AC_MSG_ERROR(Cannot find math library))
AC_SEARCH_LIBS(sqrt,[m],,AC_MSG_ERROR(Cannot find math library))

fi
# end of the list of libraries that don't need to be checked for if
# tests and examples are disabled.

#Read in the path to the Dakota installation
AC_ARG_WITH(dakota,
AC_HELP_STRING([--with-dakota=PATH],[Specify the base path of an installed version of Dakota.  .  (ex --with-dakota=/home/user/Dakota)]),
ac_cv_use_dakota=$withval, ac_cv_use_dakota=no)

AC_MSG_CHECKING(if we are using Dakota)

if test "X$ac_cv_use_dakota" != "Xno"; then
  AC_MSG_RESULT(yes: $ac_cv_use_dakota)
  dakota_installdir=$ac_cv_use_dakota
  AC_DEFINE([HAVE_STOKHOS_DAKOTA],,[Define if Dakota is enabled])
else
  AC_MSG_RESULT(no)
  DAKOTA_DIR=none
fi

AM_CONDITIONAL(HAVE_STOKHOS_DAKOTA, test "X$ac_cv_use_dakota" != "Xno")
AC_SUBST(dakota_installdir)

#Read in the path to the Fortran UQTK installation
AC_ARG_WITH(fortran-uqtk,
	    AC_HELP_STRING([--with-fortran-uqtk=PATH],
			   [Enable support for Fortran UQTK library, optionally providing the installed library path (e.g., --with-fortran-uqtk=/home/user/fluq_dist/lib/lnx86).  Default is no.]),
	    [ac_cv_use_fortran_uqtk=$withval], 
	    [ac_cv_use_fortran_uqtk=no])

AC_MSG_CHECKING([if we are using Fortran UQTK])
if test "X$ac_cv_use_fortran_uqtk" != "Xno"; then
  AC_MSG_RESULT([yes: $ac_cv_use_fortran_uqtk])
  fortan_uqtk_installdir=$ac_cv_use_fortran_uqtk
  AC_DEFINE([HAVE_STOKHOS_FORUQTK],,[Define if Fortran UQTK is enabled])
else
  AC_MSG_RESULT([no])
fi
AM_CONDITIONAL(HAVE_STOKHOS_FORUQTK, test "X$ac_cv_use_fortran_uqtk" != "Xno")
AC_SUBST(fortan_uqtk_installdir)

AC_ARG_WITH(fortran-uqtk-libs,
            [AC_HELP_STRING([--with-fortran-uqtk-libs=LIBS],
            [Fortan UQTK libraries (e.g., "-luqgen -ldvode-clean -ltiming -lcvode -lgh -lslatec -llapack -lblas ")])],
            [fortran_uqtk_libs=${withval}],
            [fortran_uqtk_libs="-luqgen -ldvode-clean -ltiming -lcvode -lgh -lslatec -llapack -lblas"])
AC_SUBST(fortran_uqtk_libs)

AM_CONDITIONAL(ENABLE_STOKHOS_QUADRATURE, test "X$ac_cv_use_dakota" != "Xno" -o "X$ac_cv_use_fortran_uqtk" != "Xno")

# ------------------------------------------------------------------------
# Checks for linker characteristics
# ------------------------------------------------------------------------

# Determine libraries needed for linking with Fortran
#AC_F77_LIBRARY_LDFLAGS


#  ----------------------------------------------------
#  Needs this because Dec handles templates differently
#  ----------------------------------------------------

# -- DEC CXX --
#if test "X$TRILINOS_TEMPLATE_OBJS" = "X"; then
#  case $target in
#    alpha*)
#      if test $CXX = cxx; then
#        TRILINOS_TEMPLATE_OBJS="cxx_repository/*.o"
#      fi
#    ;;
#  esac
#fi
                                                                  
#AC_MSG_CHECKING([for template-specific object files to be linked])if test "X$TRILINOS_TEMPLATE_OBJS" = "X"; then
# AC_MSG_RESULT(none)
#else
# AC_MSG_RESULT($TRILINOS_TEMPLATE_OBJS)
#fi
                                                                  
#AC_SUBST(TRILINOS_TEMPLATE_OBJS)

# ------------------------------------------------------------------------
# Perform substitutions in output files
# ------------------------------------------------------------------------
AC_SUBST(ac_aux_dir)

# ------------------------------------------------------------------------
# Output files
# ------------------------------------------------------------------------
AC_CONFIG_FILES([
		Makefile
		Makefile.export.stokhos
		src/Makefile
		example/Makefile
		example/generate-sample-output.pl.stub
		test/Makefile
		test/Hermite/Makefile
		])

AC_OUTPUT()

# Bye World!
echo "---------------------------------------------"
echo "Finished Running Stokhos Configure Script"
echo "---------------------------------------------"
