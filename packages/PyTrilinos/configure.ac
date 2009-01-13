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

AC_INIT(PyTrilinos, 4.2d, wfspotz@sandia.gov)

# Hello World!
echo "-----------------------------------"
echo "Running PyTrilinos Configure Script"
echo "-----------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.  Any file in that directory will do,
# preferably one that is unlikely to be removed or renamed.

AC_CONFIG_SRCDIR([src/setup.py])

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.

AC_CONFIG_AUX_DIR(config)

# ------------------------------------------------------------------------
# Check for enabled packages that have python wrappers, and will have
# their libraries converted to shared
# ------------------------------------------------------------------------
TAC_ARG_ENABLE_FEATURE_SUB( default, packages, [Build all default packages unless explicitly disabled.  When this option is disabled, only packages explicitly enabled and packages needed to satisfy the dependencies of enabled packages will be built], DEFAULT_PACKAGES, yes)

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, teuchos, PYTRILINOS_TEUCHOS,
       HAVE_PYTRILINOS_TEUCHOS, yes,
       [Build python interface and shared library for teuchos.],
       [Build python interface and shared library for teuchos.  Can be overridden with --disable-pytrilinos-teuchos.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, epetra, PYTRILINOS_EPETRA,
       HAVE_PYTRILINOS_EPETRA, yes,
       [Build python interface and shared library for epetra.],
       [Build python interface and shared library for epetra.  Can be overridden with --disable-pytrilinos-epetra.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, triutils, PYTRILINOS_TRIUTILS,
       HAVE_PYTRILINOS_TRIUTILS, yes,
       [Build python interface and shared library for triutils.],
       [Build python interface and shared library for triutils.  Can be overridden with --disable-pytrilinos-triutils.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, epetraext, PYTRILINOS_EPETRAEXT,
       HAVE_PYTRILINOS_EPETRAEXT, yes,
       [Build python interface and shared library for epetraext.],
       [Build python interface and shared library for epetraext.  Can be overridden with --disable-pytrilinos-epetraext.])

TAC_ARG_ENABLE_FEATURE_SUB(epetraext, hdf5,
       [Build python interface for hdf5 file support for epetraext.], 
       PYTRILINOS_EPETRAEXT_HDF5, no)
AM_CONDITIONAL(HAVE_PYTRILINOS_EPETRAEXT_HDF5, test "X$ac_cv_use_epetraext_hdf5" = "Xyes")

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, aztecoo, PYTRILINOS_AZTECOO,
       HAVE_PYTRILINOS_AZTECOO, yes,
       [Build python interface and shared library for aztecoo.],
       [Build python interface and shared library for aztecoo.  Can be overridden with --disable-pytrilinos-aztecoo.])

TAC_ARG_ENABLE_FEATURE_SUB(aztecoo, teuchos,
       [Build python interface and shared library for aztecoo-teuchos.], 
       PYTRILINOS_AZTECOO_TEUCHOS, yes)
AM_CONDITIONAL(HAVE_PYTRILINOS_AZTECOO_TEUCHOS, test "X$ac_cv_use_aztecoo_teuchos" = "Xyes")

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, galeri, PYTRILINOS_GALERI,
       HAVE_PYTRILINOS_GALERI, yes,
       [Build python interface and shared library for galeri.],
       [Build python interface and shared library for galeri.  Can be overridden with --disable-pytrilinos-galeri.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, amesos, PYTRILINOS_AMESOS,
       HAVE_PYTRILINOS_AMESOS, yes,
       [Build python interface and shared library for amesos.],
       [Build python interface and shared library for amesos.  Can be overridden with --disable-pytrilinos-amesos.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, ifpack, PYTRILINOS_IFPACK,
       HAVE_PYTRILINOS_IFPACK, yes,
       [Build python interface and shared library for ifpack.],
       [Build python interface and shared library for ifpack.  Can be overridden with --disable-pytrilinos-ifpack.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, komplex, PYTRILINOS_KOMPLEX,
       HAVE_PYTRILINOS_KOMPLEX, yes,
       [Build python interface and shared library for komplex.],
       [Build python interface and shared library for komplex.  Can be overridden with --disable-pytrilinos-komplex.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, anasazi, PYTRILINOS_ANASAZI,
       HAVE_PYTRILINOS_ANASAZI, yes,
       [Build python interface and shared library for anasazi.],
       [Build python interface and shared library for anasazi.  Can be overridden with --disable-pytrilinos-anasazi.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, ml, PYTRILINOS_ML,
       HAVE_PYTRILINOS_ML, yes,
       [Build python interface and shared library for ml.],
       [Build python interface and shared library for ml.  Can be overridden with --disable-pytrilinos-ml.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, nox, PYTRILINOS_NOX,
       HAVE_PYTRILINOS_NOX, yes,
       [Build python interface and shared library for nox.],
       [Build python interface and shared library for nox.  Can be overridden with --disable-pytrilinos-nox.])

TAC_ARG_ENABLE_FEATURE_SUB(nox, epetra,
       [Build python interface and shared library for nox-epetra.], 
       PYTRILINOS_NOX_EPETRA, no)
AM_CONDITIONAL(HAVE_PYTRILINOS_NOX_EPETRA, test "X$ac_cv_use_nox_epetra" = "Xyes")

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, loca, PYTRILINOS_LOCA,
       HAVE_PYTRILINOS_LOCA, yes,
       [Build python interface and shared library for loca.],
       [Build python interface shared library for loca.  Can be overridden with --disable-pytrilinos-loca.])

# ------------------------------------------------------------------------
# Check for enabled packages that do not have python wrappers, but
# will have their libraries converted to shared
# ------------------------------------------------------------------------
TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, rtop, PYTRILINOS_RTOP,
       HAVE_PYTRILINOS_RTOP, yes,
       [Build shared library for rtop.],
       [Build shared library for rtop.  Can be overridden with --disable-pytrilinos-rtop.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, kokkos, PYTRILINOS_KOKKOS,
       HAVE_PYTRILINOS_KOKKOS, yes,
       [Build shared library for kokkos.],
       [Build shared library for kokkos.  Can be overridden with --disable-pytrilinos-kokkos.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, tpetra, PYTRILINOS_TPETRA,
       HAVE_PYTRILINOS_TPETRA, yes,
       [Build shared library for tpetra.],
       [Build shared library for tpetra.  Can be overridden with --disable-pytrilinos-tpetra.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, thyra, PYTRILINOS_THYRA,
       HAVE_PYTRILINOS_THYRA, yes,
       [Build shared libraries for thyra.],
       [Build shared libraries for thyra.  Can be overridden with --disable-pytrilinos-thyra.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, isorropia, PYTRILINOS_ISORROPIA,
       HAVE_PYTRILINOS_ISORROPIA, yes,
       [Build shared library for isorropia.],
       [Build shared library for isorropia.  Can be overridden with --disable-pytrilinos-isorropia.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, claps, PYTRILINOS_CLAPS,
       HAVE_PYTRILINOS_CLAPS, yes,
       [Build shared library for claps.],
       [Build shared library for claps.  Can be overridden with --disable-pytrilinos-claps.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, belos, PYTRILINOS_BELOS,
       HAVE_PYTRILINOS_BELOS, yes,
       [Build shared library for belos.],
       [Build shared library for belos.  Can be overridden with --disable-pytrilinos-belos.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, moertel, PYTRILINOS_MOERTEL,
       HAVE_PYTRILINOS_MOERTEL, yes,
       [Build shared library for moertel.],
       [Build shared library for moertel.  Can be overridden with --disable-pytrilinos-moertel.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, stratimikos, PYTRILINOS_STRATIMIKOS,
       HAVE_PYTRILINOS_STRATIMIKOS, yes,
       [Build shared libraries for stratimikos.],
       [Build shared libraries for stratimikos.  Can be overridden with --disable-pytrilinos-stratimikos.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, meros, PYTRILINOS_MEROS,
       HAVE_PYTRILINOS_MEROS, yes,
       [Build shared library for meros.],
       [Build shared library for meros.  Can be overridden with --disable-pytrilinos-meros.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, sacado, PYTRILINOS_SACADO,
       HAVE_PYTRILINOS_SACADO, yes,
       [Build shared library for sacado.],
       [Build shared library for sacado.  Can be overridden with --disable-pytrilinos-sacado.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, rythmos, PYTRILINOS_RYTHMOS,
       HAVE_PYTRILINOS_RYTHMOS, yes,
       [Build shared library for rythmos.],
       [Build shared library for rythmos.  Can be overridden with --disable-pytrilinos-rythmos.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(pytrilinos, moocho, PYTRILINOS_MOOCHO,
       HAVE_PYTRILINOS_MOOCHO, yes,
       [Build shared libraries for moocho.],
       [Build shared libraries for moocho.  Can be overridden with --disable-pytrilinos-moocho.])

# ------------------------------------------------------------------------
# Allow users to specify their own "install" command.  If none is specified,
# the default is install-sh found in the config subdirectory.
# ------------------------------------------------------------------------
AC_ARG_WITH(install,
 [AC_HELP_STRING([--with-install=INSTALL_PROGRAM],
 [Use the installation program INSTALL_PROGRAM rather the default that is provided.  For example --with-install="/path/install -p"])],
 [
   INSTALL=$withval
   INSTALL_PROGRAM=$withval
   INSTALL_SCRIPT=$withval
   INSTALL_DATA="$withval -m 644"
 ],)

# ------------------------------------------------------------------------
# AM_MAINTAINER_MODE turns off maintainer-only makefile targets by
# default, and changes configure to understand a
# --enable-maintainer-mode option. --enable-maintainer-mode turns the
# maintainer-only targets back on. The maintainer-only makefile
# targets permit end users to clean automatically-generated files such
# as configure, which means they have to have autoconf and automake
# installed to repair the damage. AM_MAINTAINER_MODE makes it a bit
# harder for users to shoot themselves in the foot.
# ------------------------------------------------------------------------
AM_MAINTAINER_MODE

# ------------------------------------------------------------------------
# Define $build, $host, $target, etc
# ------------------------------------------------------------------------
AC_CANONICAL_TARGET

# ------------------------------------------------------------------------
# Use automake
# ------------------------------------------------------------------------
AM_INIT_AUTOMAKE(1.10 no-define tar-ustar)
AC_PREREQ(2.61)

# ------------------------------------------------------------------------
# Check to see if MPI enabled and if any special configuration done
# ------------------------------------------------------------------------
TAC_ARG_CONFIG_MPI

# ------------------------------------------------------------------------
# MPI link check
# ------------------------------------------------------------------------
TAC_ARG_CHECK_MPI

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

AZ_PYTHON_ENABLE( )   # Check for --enable-python, --disable-python
AZ_PYTHON_WITH( )     # Check for --with-python, --without-python
AZ_PYTHON_PATH( )     # Look for python assuming none has been previously found

if test -z "$PYTHON"; then
   AC_MSG_ERROR(no path to python found)
fi
AZ_PYTHON_PREFIX( )

# Ensure that we have python version 2.3 or greater (for bool type)
AZ_PYTHON_VERSION_ENSURE( [2.3] )

# ------------------------------------------------------------------------
# Checks if tests and examples should be built
# ------------------------------------------------------------------------
TAC_ARG_ENABLE_FEATURE(tests, [Make tests for all Trilinos packages buildable with 'make tests'], TESTS, yes)
TAC_ARG_ENABLE_FEATURE_SUB(pytrilinos, tests, [Make PyTrilinos tests buildable with 'make tests'], PYTRILINOS_TESTS)
AM_CONDITIONAL(BUILD_TESTS, [test "X$ac_cv_use_pytrilinos_tests" != "Xno" ])

TAC_ARG_ENABLE_FEATURE(examples, [Make examples for all Trilinos packages buildable with 'make examples'], EXAMPLES, yes)
TAC_ARG_ENABLE_FEATURE_SUB(pytrilinos, examples, [Make PyTrilinos examples buildable with 'make examples'], PYTRILINOS_EXAMPLES)
AM_CONDITIONAL(BUILD_EXAMPLES, test "X$ac_cv_use_pytrilinos_examples" != "Xno")

#We now build tests and examples through separate make targets, rather than
#during "make".  We still need to conditionally include the test and example
#in SUBDIRS, even though SUB_TEST and SUB_EXAMPLE will never be
#defined, so that the tests and examples are included in the distribution
#tarball.
AM_CONDITIONAL(SUB_TEST, test "X$ac_cv_use_sub_test" = "Xyes")
AM_CONDITIONAL(SUB_EXAMPLE, test "X$ac_cv_use_sub_example" = "Xyes")

# ------------------------------------------------------------------------
# It may seem a little strange to check for whether PyTrilinos is
# enabled here, but there is a use case whereby a user specifies
# --enable-shared but does not enable PyTrilinos.  In this case, make
# will be executed in the shared directory, but not the src, test or
# example directory.
# ------------------------------------------------------------------------
TAC_ARG_ENABLE_FEATURE(pytrilinos, [Build PyTrilinos, a Python wrapper package], PYTRILINOS, no)
AM_CONDITIONAL(BUILD_PYTRILINOS, test "X$ac_cv_use_pytrilinos" != "Xno")
if test "X$ac_cv_use_pytrilinos" != "Xno"; then

# ------------------------------------------------------------------------
# Configure should create src/PyTrilinos_config.h from src/PyTrilinos_config.h.in
# ------------------------------------------------------------------------
AM_CONFIG_HEADER(src/PyTrilinos_config.h:src/PyTrilinos_config.h.in)

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
  AC_PROG_SWIG(1.3.31)
  SWIG_ENABLE_CXX
  SWIG_MULTI_MODULE_SUPPORT
  SWIG_PYTHON
fi

# This end-if denotes the end of configuration code executed only if
# PyTrilinos is enabled.  That is to say, configuration code that is
# disabled for shared-library-only builds.
fi

AM_CONDITIONAL(HAVE_SWIG,test X${WITH_SWIG} = Xyes)
AC_SUBST(GNU_HAVE_SWIG, ${WITH_SWIG})

# ------------------------------------------------------------------------
# Output files
# ------------------------------------------------------------------------

AC_CONFIG_FILES([
		Makefile
		util/MakefileVariables.py
		util/PyTrilinosExtension.py
		util/SharedUtils.py
		util/pkg_info.py
		shared/Makefile
		src/Makefile
		src/pyfragments.swg
		src/.deps/dummy
		src/PyTrilinos/.dummy
		src/PyTrilinos/NOX/.dummy
		src/PyTrilinos/NOX/Epetra/.dummy
		src/PyTrilinos/NOX/LAPACK/.dummy
		src/PyTrilinos/LOCA/.dummy
		src/PyTrilinos/LOCA/Epetra/.dummy
		test/Makefile
		example/Makefile
		])

AC_OUTPUT()

# Bye World!
echo "--------------------------------------------"
echo "Finished Running PyTrilinos Configure Script"
echo "--------------------------------------------"
