# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# @HEADER
# ************************************************************************
# 
#         Stratimikos: Thyra-based strategies for linear solvers
#                Copyright (2006) Sandia Corporation
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
# Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
# 
# ************************************************************************
# @HEADER

# ------------------------------------------------------------------------
# Initialization 
# ------------------------------------------------------------------------

# This must be the first line in configure.ac.
# Optional 3rd argument is email address for bugs.

AC_INIT(stratimikos, Development Version, rabartl@sandia.gov)

echo "----------------------------------------"
echo "Running stratimikos Configure Script"
echo "----------------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.  Any file in that directory will do,
# preferably one that is unlikely to be removed or renamed.

AC_CONFIG_SRCDIR([src/Thyra_DefaultRealLinearSolverBuilder.hpp])

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.

#
#   I wimped out here.  I still think that ../../config can work 
#   for any package with a frim dependency on Trilinos.
#   But, using ../../config puts some ../../config files into a list
#   of files that it looked like it might delete under certain circumstances,
#   or perhaps put into some tar ball.  Anyway, it scared me off.
#   In order to get rid of the local config directory, you have to 
#   reference the config directory in Makefile.am and in bootstrap
#   (which should use aclocal -I ../../config)
#
AC_CONFIG_AUX_DIR(config)

# Configure should create src/Stratimikos_Config.h from src/Stratimikos_Config.h.in

AM_CONFIG_HEADER(src/Stratimikos_Config.h:src/Stratimikos_Config.h.in)
AM_CONFIG_HEADER(src/Stratimikos_InternalConfig.h:src/Stratimikos_InternalConfig.h.in)

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
AC_PROG_CXX(CC g++ c++ cxx)
#AC_PROG_F77(f77 g77 gfortran f90 xlf90 f95)
AC_PROG_RANLIB

# ------------------------------------------------------------------------
# Checks for extra flags
# ------------------------------------------------------------------------
# Check if --with-flags present, prepend any specs to FLAGS

TAC_ARG_WITH_FLAGS(ccflags, CCFLAGS)
TAC_ARG_WITH_FLAGS(cxxflags, CXXFLAGS)
TAC_ARG_WITH_FLAGS(cflags, CFLAGS)
#TAC_ARG_WITH_FLAGS(fflags, FFLAGS)
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

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE

# ------------------------------------------------------------------------
# Checks for compiler characteristics - Epetra_ConfigDefs.h uses these, 
# so it makes sense to leave them in - see PORTABILITY in this directory
# ------------------------------------------------------------------------

# ---------------------------------------------------------------------
# These are pretty obscure
# ---------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for special package flags
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_CAN_USE_PACKAGE(
  stratimikos, amesos, STRATIMIKOS_AMESOS, ENABLE_STRATIMIKOS_AMESOS, yes,
  [Enable Thyra adapters for Amesos in Stratimikos.],
  [Enable Thyra adapters for Amesos in Stratimikos.  Can be overridden with --disable-stratimikos-amesos.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(
  stratimikos, aztecoo, STRATIMIKOS_AZTECOO, ENABLE_STRATIMIKOS_AZTECOO, yes,
  [Enable Thyra adapters for AztecOO in Stratimikos.],
  [Enable Thyra adapters for AztecOO in Stratimikos.  Can be overridden with --disable-stratimikos-aztecoo.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(
  stratimikos, belos, STRATIMIKOS_BELOS, ENABLE_STRATIMIKOS_BELOS, yes,
  [Enable Thyra adapters for Belos in Stratimikos.],
  [Enable Thyra adapters for Belos in Stratimikos.  Can be overridden with --disable-stratimikos-belos.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(
  stratimikos, ifpack, STRATIMIKOS_IFPACK, ENABLE_STRATIMIKOS_IFPACK, yes,
  [Enable Thyra adapters for Ifpack in Stratimikos.],
  [Enable Thyra adapters for Ifpack in Stratimikos.  Can be overridden with --disable-stratimikos-ifpack.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(
  stratimikos, ml, STRATIMIKOS_ML, ENABLE_STRATIMIKOS_ML, yes,
  [Enable Thyra adapters for ML in Stratimikos.],
  [Enable Thyra adapters for ML in Stratimikos.  Can be overridden with --disable-stratimikos-ml.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(
  stratimikos, triutils, STRATIMIKOS_TRIUTILS, ENABLE_STRATIMIKOS_TRIUTILS, yes,
  [Build tests that use Triutils. Can be overridden with --disable-stratimikos-triutils.],
  [Build tests that use Triutils.])

TAC_ARG_INFORM_SUBPACKAGE_USED(
  thyra, epetra, STRATIMIKOS_THYRA_EPETRA, ENABLE_STRATIMIKOS_THYRA_EPETRA, yes,
  [see thyra/configure --help],
  [see thyra/configure --help])

TAC_ARG_INFORM_SUBPACKAGE_USED(
  thyra, tpetra, STRATIMIKOS_THYRA_TPETRA, ENABLE_STRATIMIKOS_THYRA_TPETRA, yes,
  [see thyra/configure --help],
  [see thyra/configure --help])

TAC_ARG_INFORM_SUBPACKAGE_USED(
  thyra, epetraext, STRATIMIKOS_THYRA_EPETRAEXT, ENABLE_STRATIMIKOS_THYRA_EPETRAEXT, yes,
  [See thyra/configure --help],
  [See thyra/configure --help])

TAC_ARG_INFORM_SUBPACKAGE_USED(
  ml, nox, STRATIMIKOS_ML_NOX, ENABLE_STRATIMIKOS_ML_NOX, yes,
  [See ml/configure --help],
  [See ml/configure --help])

TAC_ARG_ENABLE_FEATURE_SUB(teuchos, complex, [Build examples/tests which use complex scalar type], TEUCHOS_COMPLEX, yes)
AM_CONDITIONAL(BUILD_TEUCHOS_COMPLEX, [test "X$ac_cv_use_teuchos_complex" != "Xno"])

# ------------------------------------------------------------------------
# Checks if examples should be built
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_FEATURE(tests, [Make tests for all Trilinos packages buildable with 'make tests'], TESTS, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK( stratimikos, tests, [Make Stratimikos tests buildable with 'make tests'], STRATIMIKOS_TESTS)
AM_CONDITIONAL(BUILD_TESTS, test "X$ac_cv_use_stratimikos_tests" != "Xno")

TAC_ARG_ENABLE_FEATURE(examples, [Make examples for all Trilinos packages buildable with 'make examples'], EXAMPLES, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK( stratimikos, examples, [Make Stratimikos examples buildable with 'make examples'], STRATIMIKOS_EXAMPLES)
AM_CONDITIONAL(BUILD_EXAMPLES, test "X$ac_cv_use_stratimikos_examples" != "Xno")

#We now build tests and examples through separate make targets, rather than
#during "make".  We still need to conditionally include the test and example
#in SUBDIRS, even though SUB_TEST and SUB_EXAMPLE will never be
#defined, so that the tests and examples are included in the distribution
#tarball.
AM_CONDITIONAL(SUB_TEST, test "X$ac_cv_use_sub_test" = "Xyes")
AM_CONDITIONAL(SUB_EXAMPLE, test "X$ac_cv_use_sub_example" = "Xyes")

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

if test "X$ac_cv_use_stratimikos_tests" != "Xno" || test "X$ac_cv_use_stratimikos_examples" != "Xno" || test "X$ac_cv_use_libcheck" != "Xno"; then

# #np# - 
# To add any libraries that you want to check for, replace pow below 
# with the name of an entry point in that library.
# checks for libraries now
dnl Replace `main' with a function in -lm:
#  #np# These lines can be removed if the math lib is not needed.
AC_SEARCH_LIBS(pow,[m],,AC_MSG_ERROR(Cannot find math library))
AC_SEARCH_LIBS(sqrt,[m],,AC_MSG_ERROR(Cannot find math library))

fi

# end of the list of libraries that don't need to be checked for if
# tests and examples are disabled.

# ------------------------------------------------------------------------
# Checks for compiler characteristics
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for Makefile.export related systems
# ------------------------------------------------------------------------
TAC_ARG_ENABLE_EXPORT_MAKEFILES(yes)

# ------------------------------------------------------------------------
# Checks for linker characteristics
# ------------------------------------------------------------------------

# Determine libraries needed for linking with Fortran
#AC_F77_LIBRARY_LDFLAGS

#  ----------------------------------------------------
#  Needs this because Dec handles templates differently
#  ----------------------------------------------------
#  #np# - Only uncomment the following lines if the package
#         creates a subdirectory called cxx_repository for
#         template object files on DEC.  (HINT: try to build on 
#         DEC/CPLANT and see if template-related errors occur)

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
    Makefile.export.stratimikos
    src/Makefile
    example/Makefile
    example/generate-sample-output.pl.stub
    test/Makefile
    adapters/Makefile
    adapters/amesos/Makefile
    adapters/amesos/src/Makefile
    adapters/amesos/test/Makefile
    adapters/aztecoo/Makefile
    adapters/aztecoo/src/Makefile
    adapters/aztecoo/test/Makefile
    adapters/belos/Makefile
    adapters/belos/src/Makefile
    adapters/belos/test/Makefile
    adapters/belos/test/MVOPTester/Makefile
    adapters/belos/test/LOWSFactory/Makefile
    adapters/belos/example/Makefile
    adapters/belos/example/LOWSFactory/Makefile
    adapters/belos/example/LOWSFactory/Epetra/Makefile
    adapters/belos/example/LOWSFactory/Tpetra/Makefile
    adapters/ifpack/Makefile
    adapters/ifpack/src/Makefile
    adapters/ml/Makefile
    adapters/ml/src/Makefile
    ])

AC_OUTPUT()

# Bye World!
echo "---------------------------------------------"
echo "Finished Running Stratimikos Configure Script"
echo "---------------------------------------------"
