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

AC_INIT(nox, 0.1b)

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.

AC_CONFIG_SRCDIR(src/NOX_Solver_Newton.C)

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.

AC_CONFIG_AUX_DIR(config)

# Configure should create src/NOX_Config.h from src/NOX_Config.h.in

AM_CONFIG_HEADER(src/NOX_Config.h:src/NOX_Config.h.in)

# Define $build, $host, $target, etc

AC_CANONICAL_TARGET

# Use automake, getting name and version from variables defined by
# AC_INIT command.

AM_INIT_AUTOMAKE(${PACKAGE_NAME}, ${PACKAGE_VERSION})

# Specify required version of autoconf.

AC_PREREQ(2.52)

# Specify revision number of configure.ac that was used to generate
# configure.

AC_REVISION($Revision$)

# ------------------------------------------------------------------------
# Configure command line options
# ------------------------------------------------------------------------

AC_ARG_ENABLE(debug,
[  --enable-debug          Compile with debugging options],
[
case $enableval in
  yes)
    DEBUG=yes
    CXXFLAGS="${CXXFLAGS} -g"
    CFLAGS="${CFLAGS} -g"
  ;;
  opt)
    DEBUG=opt
    CXXFLAGS="${CXXFLAGS} -O2"
    CFLAGS="${CFLAGS} -O2"
  ;;
  no)
    DEBUG=no
  ;;
  *)
    AC_MSG_ERROR([Invalid value for --enable-debug ($enableval)])
  ;;
esac
],
[DEBUG=no]
)
AC_MSG_CHECKING(whether debugging is on)
AC_MSG_RESULT([${DEBUG}])

AC_ARG_ENABLE(mpi,
[  --enable-mpi            Enable MPI.],
[
case $enableval in
  yes)
    USE_MPI=yes 
    USE_MPI_COMPILER=yes
  ;;
  no)
    USE_MPI=no
  ;;
  *)
    AC_MSG_ERROR([Invalid value for --enable-mpi ($enableval)])
  ;;
esac
],
[USE_MPI=no]
)

AC_ARG_WITH(mpiinc,
[  --with-mpiinc           Specify MPI includes (enables MPI)],
[
USE_MPI=yes
USE_MPI_COMPILER=no
AC_MSG_CHECKING(user-defined MPI includes)
AC_MSG_RESULT([${withval}])
CPPFLAGS="${CPPFLAGS} ${withval}"
]
)

AC_ARG_WITH(mpilib,
[  --with-mpilib           Specify MPI libraries (enables MPI)],
[
USE_MPI=yes
USE_MPI_COMPILER=no
AC_MSG_CHECKING(user-defined MPI libraries)
AC_MSG_RESULT([${withval}])
LIBS="${LIBS} ${withval}"
]
)

AC_MSG_CHECKING(whether we are using mpi)
AC_MSG_RESULT([${USE_MPI}])

# ------------------------------------------------------------------------
# Checks for programs
# ------------------------------------------------------------------------

if test X${USE_MPI} = Xyes; then

  if test -z "$CC" && test X${USE_MPI_COMPILER} = Xyes; then
    AC_MSG_CHECKING(for mpicc)
    AC_CHECK_PROG(MPICC_EXISTS, mpicc, yes, no)
    if test X${MPICC_EXISTS} = Xyes; then
      CC=mpicc
    fi
  fi

  if test -z "$CXX" && test X${USE_MPI_COMPILER} = Xyes; then
    AC_MSG_CHECKING(for mpiCC)
    AC_CHECK_PROG(MPICXX_EXISTS, mpiCC, yes, no)
    if test X${MPICXX_EXISTS} = Xyes; then
      CXX=mpiCC
    fi
  fi

fi

AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB

# ------------------------------------------------------------------------
# MPI check
# ------------------------------------------------------------------------

if test X${USE_MPI} = Xyes; then

  AC_LANG_SAVE

  AC_LANG_C
  AC_MSG_CHECKING(whether MPI will link using C compiler)
  AC_TRY_LINK([#include <mpi.h>],
  [int c; char** v; MPI_Init(&c,&v);],
  [AC_MSG_RESULT(yes)],
  [AC_MSG_RESULT(no)
   AC_MSG_ERROR(MPI cannot link)]
  )

  AC_LANG_CPLUSPLUS 
  AC_MSG_CHECKING(whether MPI will link using C++ compiler)
  AC_TRY_LINK([#include <mpi.h>],
  [int c; char** v; MPI_Init(&c,&v);],
  [AC_MSG_RESULT(yes)], 
  [AC_MSG_RESULT(no)  
   AC_MSG_ERROR(MPI cannot link)]
  )

  AC_LANG_RESTORE

fi

# ------------------------------------------------------------------------
# Checks for libraries
# ------------------------------------------------------------------------

AC_CHECK_LIB([m],[pow])
AC_CHECK_FUNCS([pow])

# Check that all libraries will link with C++
if test ${LIBS}; then
  AC_MSG_CHECKING([[whether libraries (${LIBS}) will link]])
  AC_LANG_CPLUSPLUS 
  AC_TRY_LINK([],[],
  [AC_MSG_RESULT(yes)],
  [
  AC_MSG_RESULT(no)
  AC_MSG_ERROR(cannot link)
  ]
  )
fi

# ------------------------------------------------------------------------
# Checks for header files
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for types
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for structures
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for compiler characteristics
# ------------------------------------------------------------------------

AC_CXX_NAMESPACES
AC_CXX_MUTABLE
AC_CXX_HAVE_STL
AC_CXX_HAVE_STD

AC_MSG_CHECKING([for special compiler options])
case $target in
  rs6000-ibm-aix3.2.* | rs6000-ibm-aix4.* | powerpc-ibm-aix4.*)
    # -qnolm squashes warnings about compiler licences
    COPTIONS="-qnolm"
    # -qrtti allows for dynamic casts
    CXXOPTIONS="-qnolm -qrtti"
    AC_MSG_RESULT([rs6000 or powerpc])
  ;;
  *)
    AC_MSG_RESULT([none])
  ;;
esac

CFLAGS="${CFLAGS} ${COPTIONS}"
CXXFLAGS="${CXXFLAGS} ${CXXOPTIONS}"
CPPFLAGS="${CPPFLAGS} ${XTRAFLAGS} ${XTRAINC}"

# ------------------------------------------------------------------------
# Checks for library functions
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for system services
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for MPI/SERIAL status (only affects the examples)
# ------------------------------------------------------------------------

AC_TRILINOS

if test X${USE_TRILINOS} = Xyes; then
   AM_CONDITIONAL(BUILD_NOXEPETRA, true)
else
   AM_CONDITIONAL(BUILD_NOXEPETRA, false)
fi

# ------------------------------------------------------------------------
# Perform substitutions in output files
# ------------------------------------------------------------------------

AC_SUBST(ac_aux_dir)

# ------------------------------------------------------------------------
# Output files
# ------------------------------------------------------------------------

AC_CONFIG_FILES([
		Makefile 
		src/Makefile 
		example/Makefile 
		example/epetra_DS6.5.1/Makefile
		example/epetra_1Dfem_linear/Makefile
		example/epetra_1Dfem_nonlinear/Makefile
		])

AC_OUTPUT()


