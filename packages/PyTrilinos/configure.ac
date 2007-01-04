# @HEADER
# ************************************************************************
# 
#                PyTrilinos: Python Interface to Trilinos
#                   Copyright (2005) Sandia Corporation
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
# Questions? Contact Bill Spotz (wfspotz@sandia.gov) 
# 
# ************************************************************************
# @HEADER

# ------------------------------------------------------------------------
# Initialization 
# ------------------------------------------------------------------------

AC_INIT(PyTrilinos, 4.0, wfspotz@sandia.gov)

# Hello World!
echo "-----------------------------------"
echo "Running PyTrilinos Configure Script"
echo "-----------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.  Any file in that directory will do,
# preferably one that is unlikely to be removed or renamed.

AC_CONFIG_SRCDIR(/src/setup.py)

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.

AC_CONFIG_AUX_DIR(config)

# Get all of the enabled and disabled packages

TAC_ARG_ENABLE_FEATURE_SUB( default, packages, [Build all default packages unless explicitly disabled.  When this option is disabled, only packages explicitly enabled and packages needed to satisfy the dependencies of enabled packages will be built], DEFAULT_PACKAGES, yes)

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, teuchos, PYTRILINOS_TEUCHOS,
       HAVE_PYTRILINOS_TEUCHOS, yes,
       [Build python interface to teuchos.],
       [Build python interface to teuchos.  Can be overridden with --disable-pytrilinos-teuchos.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, thyra, PYTRILINOS_THYRA,
       HAVE_PYTRILINOS_THYRA, yes,
       [Build python interface to thyra.],
       [Build python interface to thyra.  Can be overridden with --disable-pytrilinos-thyra.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, epetra, PYTRILINOS_EPETRA,
       HAVE_PYTRILINOS_EPETRA, yes,
       [Build python interface to epetra.],
       [Build python interface to epetra.  Can be overridden with --disable-pytrilinos-epetra.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, triutils, PYTRILINOS_TRIUTILS,
       HAVE_PYTRILINOS_TRIUTILS, yes,
       [Build python interface to triutils.],
       [Build python interface to triutils.  Can be overridden with --disable-pytrilinos-triutils.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, epetraext, PYTRILINOS_EPETRAEXT,
       HAVE_PYTRILINOS_EPETRAEXT, yes,
       [Build python interface to epetraext.],
       [Build python interface to epetraext.  Can be overridden with --disable-pytrilinos-epetraext.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, aztecoo, PYTRILINOS_AZTECOO,
       HAVE_PYTRILINOS_AZTECOO, yes,
       [Build python interface to aztecoo.],
       [Build python interface to aztecoo.  Can be overridden with --disable-pytrilinos-aztecoo.])

TAC_ARG_ENABLE_FEATURE_SUB(aztecoo, teuchos,
       [Build python interface to aztecoo-teuchos.], 
       PYTRILINOS_AZTECOO_TEUCHOS, yes)
AM_CONDITIONAL(HAVE_PYTRILINOS_AZTECOO_TEUCHOS, test "X$ac_cv_use_aztecoo_teuchos" = "Xyes")

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, galeri, PYTRILINOS_GALERI,
       HAVE_PYTRILINOS_GALERI, yes,
       [Build python interface to galeri.],
       [Build python interface to galeri.  Can be overridden with --disable-pytrilinos-galeri.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, amesos, PYTRILINOS_AMESOS,
       HAVE_PYTRILINOS_AMESOS, yes,
       [Build python interface to amesos.],
       [Build python interface to amesos.  Can be overridden with --disable-pytrilinos-amesos.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, ifpack, PYTRILINOS_IFPACK,
       HAVE_PYTRILINOS_IFPACK, yes,
       [Build python interface to ifpack.],
       [Build python interface to ifpack.  Can be overridden with --disable-pytrilinos-ifpack.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, anasazi, PYTRILINOS_ANASAZI,
       HAVE_PYTRILINOS_ANASAZI, yes,
       [Build python interface to anasazi.],
       [Build python interface to anasazi.  Can be overridden with --disable-pytrilinos-anasazi.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, ml, PYTRILINOS_ML,
       HAVE_PYTRILINOS_ML, yes,
       [Build python interface to ml.],
       [Build python interface to ml.  Can be overridden with --disable-pytrilinos-ml.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, nox, PYTRILINOS_NOX,
       HAVE_PYTRILINOS_NOX, yes,
       [Build python interface to nox.],
       [Build python interface to nox.  Can be overridden with --disable-pytrilinos-nox.])

TAC_ARG_ENABLE_FEATURE_SUB(nox, epetra,
       [Build python interface to nox-epetra.], 
       PYTRILINOS_NOX_EPETRA, no)
AM_CONDITIONAL(HAVE_PYTRILINOS_NOX_EPETRA, test "X$ac_cv_use_nox_epetra" = "Xyes")

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, loca, PYTRILINOS_LOCA,
       HAVE_PYTRILINOS_LOCA, yes,
       [Build python interface to loca.],
       [Build python interface to loca.  Can be overridden with --disable-pytrilinos-loca.])

# Configure should create src/PyTrilinos_config.h from src/PyTrilinos_config.h.in

AM_CONFIG_HEADER(src/PyTrilinos_config.h:src/PyTrilinos_config.h.in)

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

# ------------------------------------------------------------------------
# Check to see if MPI enabled and if any special configuration done
# ------------------------------------------------------------------------

TAC_ARG_CONFIG_MPI

# ------------------------------------------------------------------------
# MPI link check
# ------------------------------------------------------------------------
TAC_ARG_CHECK_MPI

# ------------------------------------------------------------------------
# Checks if tests and examples should be built
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_FEATURE(tests, [Build tests for all Trilinos packages (not all packages are sensitive to this option)], TESTS, yes)
TAC_ARG_ENABLE_FEATURE_SUB(pytrilinos, tests, [Enable PyTrilinos tests], PYTRILINOS_TESTS)
AM_CONDITIONAL(BUILD_TESTS, [test "X$ac_cv_use_pytrilinos_tests" != "Xno" ])

TAC_ARG_ENABLE_FEATURE(examples, [Build examples for all Trilinos packages (not all packages are sensitive to this option)], EXAMPLES, yes)
TAC_ARG_ENABLE_FEATURE_SUB(pytrilinos, examples, [Enable PyTrilinos examples], PYTRILINOS_EXAMPLES)
AM_CONDITIONAL(BUILD_EXAMPLES, test "X$ac_cv_use_pytrilinos_examples" != "Xno")

#TAC_ARG_ENABLE_FEATURE(libcheck, [Check for some third-party libraries including BLAS and LAPACK.  (Cannot be disabled unless tests and examples are also disabled.)], LIBCHECK, yes)

# ------------------------------------------------------------------------
# Specify other directories 
# ------------------------------------------------------------------------

# enable use of --with-libdirs="-Llibdir1 -Llibdir2 ..." to prepend to LDFLAGS
TAC_ARG_WITH_LIBDIRS
# enable use of --with-incdirs="-Iincdir1 -Iincdir2 ..." to prepend to CPPFLAGS
TAC_ARG_WITH_INCDIRS

# ------------------------------------------------------------------------
# Configure python
# ------------------------------------------------------------------------

AZ_PYTHON_WITH( )     # Check for --with-python, --without-python
az_python_use=true
AM_CONDITIONAL(PYTHON_USE, test x"$az_python_use" = x"true")

AC_PATH_PROG([PYTHON], python)

if test -z "$PYTHON"; then
   AC_MSG_ERROR(no path to python found)
fi
AZ_PYTHON_PREFIX( )

# Ensure that we have python version 2.3 or greater (for bool type)
AZ_PYTHON_VERSION_ENSURE( [2.3] )

# Python compiler and linker flags
AZ_PYTHON_CSPEC
AZ_PYTHON_LSPEC

# Check that Python.h is available
save_CPPFLAGS=$CPPFLAGS
CPPFLAGS="$save_CPPFLAGS $PYTHON_CSPEC"
AC_LANG([C++])
AC_CHECK_HEADER(
	[Python.h],
	break,
	AC_MSG_ERROR([You must have Python.h in order to build the Python support!!]))
CPPFLAGS="$save_CPPFLAGS"

# Check for numpy python module and user_array/UserArray submodule
AC_PYTHON_MODULE(numpy,yes)
AC_PYTHON_MODULE(numpy.lib.user_array)
if test "X$HAVE_PYMOD_NUMPY_LIB_USER_ARRAY" == "Xno"; then
  AC_PYTHON_MODULE(numpy.lib.UserArray,yes)
fi

# If user specifies prefix, use it for the PYTHON_PREFIX
if test "$prefix" != "$ac_default_prefix"; then
  PYTHON_PREFIX=$prefix
fi
if test "$exec_prefix" != "$ac_default_prefix"; then
  PYTHON_EXECPREFIX=$exec_prefix
fi

# ------------------------------------------------------------------------
# If the python wrappers are to be built, then SWIG (Simple Wrapper
# Interface Generator) is required
# ------------------------------------------------------------------------

if test -n "$PYTHON"; then

  # Check for --with-swig[=path]
  AC_MSG_CHECKING(for --with-swig)
  AC_ARG_WITH(swig,
              [AC_HELP_STRING([--with-swig@<:@=SWIG@:>@],
                              [enable swig and set swig binary])],
              [AC_MSG_RESULT(yes)
               WITH_SWIG=yes
               if test X${withval} != Xyes ; then
                 SWIG=$withval
               fi],
              [AC_MSG_RESULT(no)
               AC_CHECK_PROG(WITH_SWIG,swig,yes,no)])

  # Report error if no swig found
  if test ${WITH_SWIG} = no; then
     AC_MSG_ERROR(
     [Python wrappers require swig (Simple Wrapper Interface Generator).
      See http://www.swig.org])
  fi

  # SWIG configuration
  AC_PROG_SWIG(1.3.28)
  SWIG_ENABLE_CXX
  SWIG_MULTI_MODULE_SUPPORT
  SWIG_PYTHON
fi
AM_CONDITIONAL(HAVE_SWIG,test X${WITH_SWIG} = Xyes)
AC_SUBST(GNU_HAVE_SWIG, ${WITH_SWIG})

# ------------------------------------------------------------------------
# Output files
# ------------------------------------------------------------------------

AC_CONFIG_FILES([
		Makefile 
		src/Makefile
		src/.deps/dummy
		src/PyTrilinos/.dummy
		src/PyTrilinos/NOX/.dummy
		src/PyTrilinos/NOX/Epetra/.dummy
		src/PyTrilinos/NOX/LAPACK/.dummy
		test/Makefile
		example/Makefile
		])

AC_OUTPUT()

# Bye World!
echo "--------------------------------------------"
echo "Finished Running PyTrilinos Configure Script"
echo "--------------------------------------------"
