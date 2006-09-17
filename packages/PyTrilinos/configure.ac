# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# @HEADER
# ************************************************************************
# 
#                  PyTrilinos: Rapid Prototyping Package
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
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************
# @HEADER

# ------------------------------------------------------------------------
# Initialization 
# ------------------------------------------------------------------------

AC_INIT(PyTrilinos, 3.0, wfspotz@sandia.gov)

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
AM_CONDITIONAL(BUILD_TESTS, test "X$ac_cv_use_tests" != "Xno")

TAC_ARG_ENABLE_FEATURE(examples, [Build examples for all Trilinos packages (not all packages are sensitive to this option)], EXAMPLES, yes)
AM_CONDITIONAL(BUILD_EXAMPLES, test "X$ac_cv_use_examples" != "Xno")

#TAC_ARG_ENABLE_FEATURE(libcheck, [Check for some third-party libraries including BLAS and LAPACK.  (Cannot be disabled unless tests and examples are also disabled.)], LIBCHECK, yes)

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

# Check for numpy python module
AC_PYTHON_MODULE(numpy,yes)

# If user specifies prefix, use it for the PYTHON_PREFIX
if test "$prefix" != $ac_default_prefix; then
    PYTHON_PREFIX=$prefix
fi
if test "$exec_prefix" != $ac_default_prefix; then
    PYTHON_EXECPREFIX=$exec_prefix
fi

# ------------------------------------------------------------------------
# Output files
# ------------------------------------------------------------------------

AC_CONFIG_FILES([
		Makefile 
		example/Makefile
		src/Makefile
		src/TRILINOS_HOME_DIR
		test/Makefile
		])

AC_OUTPUT()

# Bye World!
echo "--------------------------------------------"
echo "Finished Running PyTrilinos Configure Script"
echo "--------------------------------------------"
