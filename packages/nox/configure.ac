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

AC_INIT(nox, 0.2a, nox-developers@csmr.ca.sandia.gov)

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.

AC_CONFIG_SRCDIR(src/NOX_Solver_LineSearchBased.C)

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.

AC_CONFIG_AUX_DIR(config)

# Configure should create src/NOX_Config.h from src/NOX_Config.h.in

AM_CONFIG_HEADER(src/NOX_Config.h:src/NOX_Config.h.in)

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

AM_INIT_AUTOMAKE(${PACKAGE_NAME}, ${PACKAGE_VERSION})

# Specify required version of autoconf.

AC_PREREQ(2.52)

# ------------------------------------------------------------------------
# Check for user overrides (to be used later)
# ------------------------------------------------------------------------

if test -n "${CXXFLAGS}"; then
  USERDEF_CXXFLAGS=yes
fi

if test -n "${CPPFLAGS}"; then
  USERDEF_CPPFLAGS=yes
fi

# ------------------------------------------------------------------------
# Checks for extra flags, etc
# ------------------------------------------------------------------------

AC_ARG_WITH(cppflags,
[  --with-cppflags         add extra preprocessor flags (e.g., "-Dflag -Idir")],
[EXTRA_CPPFLAGS=${withval}]
)

AC_ARG_WITH(cxxflags,
[  --with-cxxflags         add extra C++ compiler flags],
[EXTRA_CXXFLAGS=${withval}]
)

# ------------------------------------------------------------------------
# Alternate archiver
# ------------------------------------------------------------------------

AC_ARG_WITH(ar,
[  --with-ar               specify archiver command (default is "ar cru")],
[
AC_MSG_CHECKING(user-defined archiver)
AC_MSG_RESULT([${withval}])
ALTERNATE_AR="${withval}"
]
)

# ------------------------------------------------------------------------
# MPI configure command line option
# ------------------------------------------------------------------------
AC_ARG_WITH(mpi-cxx,
[  --with-mpi-cxx          specify MPI C++ compiler (uses mpiCC if no argument)],
[
HAVE_PKG_MPI=yes 
if test X${withval} = Xyes; then
  MPI_CXX=mpiCC
else
  MPI_CXX=${withval}
fi
AC_MSG_CHECKING(user-defined MPI C++ compiler)
AC_MSG_RESULT([${MPI_CXX}])
]
)

AC_ARG_WITH(mpi,
[  --with-mpi              specify root directory of MPI installation],
[
HAVE_PKG_MPI=yes
MPI_DIR=${withval}
AC_MSG_CHECKING(MPI directory)
AC_MSG_RESULT([${MPI_DIR}])
]
)

AC_ARG_WITH(mpi-libs,
[  --with-mpi-libs         specify MPI libraries],
[
HAVE_PKG_MPI=yes
MPI_LIBS=${withval}
AC_MSG_CHECKING(user-defined MPI libraries)
AC_MSG_RESULT([${MPI_LIBS}])
]
)

AC_ARG_WITH(mpi-include,
[  --with-mpi-include      specify include directory for MPI],
[
HAVE_PKG_MPI=yes
MPI_INC=${withval}
AC_MSG_CHECKING(user-defined MPI includes)
AC_MSG_RESULT([${MPI_INC}])
]
)

AC_ARG_WITH(mpi-libdir,
[  --with-mpi-libdir       specify location of MPI libraries],
[
HAVE_PKG_MPI=yes
MPI_LIBDIR=${withval}
AC_MSG_CHECKING(user-defined MPI library directory)
AC_MSG_RESULT([${MPI_LIBDIR}])
]
)

AC_MSG_CHECKING(whether we are using MPI)
AC_MSG_RESULT([${HAVE_PKG_MPI}])

# ------------------------------------------------------------------------
# Check for MPI C++ compiler (must be done *before* AC_PROG_CXX)
# ------------------------------------------------------------------------

if test -n "${MPI_CXX}"; then

  if test -f ${MPI_CXX}; then   
    MPI_CXX_EXISTS=yes
  else
    AC_CHECK_PROG(MPI_CXX_EXISTS, ${MPI_CXX}, yes, no)
  fi

  if test "X${MPI_CXX_EXISTS}" = "Xyes"; then
    CXX=${MPI_CXX}
  else
    AC_MSG_ERROR([MPI C++ compiler (${MPI_CXX}) not found.])
  fi

fi

# ------------------------------------------------------------------------
# Checks for programs
# ------------------------------------------------------------------------

AC_PROG_CXX(CC c++ g++ gcc cxx)
AC_PROG_CXXCPP
AC_PROG_RANLIB

# ------------------------------------------------------------------------
# Checks for special compiler options
# ------------------------------------------------------------------------

AC_ARG_ENABLE(debug,
[  --enable-debug          compile with debugging],
[DEBUG=${enableval}],
[DEBUG="default"]
)
AC_MSG_CHECKING(whether debugging is on)
AC_MSG_RESULT([${DEBUG}])

AC_ARG_ENABLE(opt,
[  --enable-opt            compile with optimization],
[OPT=${enableval}],
[OPT=default]
)
AC_MSG_CHECKING([whether optimization is on])
AC_MSG_RESULT([${OPT}])

# options needed only for optimization
CXXFLAGS_OPT=-O
# options needed only for debugging
CXXFLAGS_DBG=-g
# options that are always needed
CXXFLAGS_MISC=

case $target in

  # IBM test
  rs6000-ibm-aix3.2.* | rs6000-ibm-aix4.* | powerpc-ibm-aix4.*)

    case $CXX in
      xlC* ) 
      CXXFLAGS_MISC="-qnolm -qrtti"
      ;;
    esac

  ;;

  # SGI/IRIX test
  mips-sgi-irix* ) 

    if test "${CXX}" = "CC"; then
     CXXFLAGS_MISC="-LANG:std -LANG:ansi-for-init-scope=ON -64 -DSGI64"
     F77=f77
     FFLAGS="-64"
     SPECIAL_AR="CC -ar -W,l -W,s -o"    
    fi

  ;;

esac

AC_MSG_CHECKING([for C++ flags])

if test "X${USERDEF_CXXFLAGS}" = "Xyes"; then

  AC_MSG_RESULT([${CXXFLAGS} (user-specified via CXXFLAGS)])

else

  if test ${DEBUG} = yes; then

  CXXFLAGS=${CXXFLAGS_DBG}

  elif test ${OPT} = yes; then

  CXXFLAGS=${CXXFLAGS_OPT}

  fi

  if test -n "${CXXFLAGS_MISC}"; then
    CXXFLAGS="${CXXFLAGS} ${CXXFLAGS_MISC}"
  fi			  

  AC_MSG_RESULT([${CXXFLAGS}])

fi

AC_MSG_CHECKING([for extra C++ compiler flags])
if test -z "${EXTRA_CPPFLAGS}"; then
  AC_MSG_RESULT(none)
else
  AC_MSG_RESULT([${EXTRA_CPPFLAGS}])
  CXXFLAGS="${CXXFLAGS} ${EXTRA_CPPFLAGS}"
fi

AC_MSG_CHECKING([for extra preprocessor flags])
if test -z "${EXTRA_CXXFLAGS}"; then
  AC_MSG_RESULT(none)
else
  AC_MSG_RESULT([${EXTRA_CXXFLAGS}])
  CPPFLAGS="${CPPFLAGS} ${EXTRA_CXXFLAGS}"
fi

# ------------------------------------------------------------------------
# Set up alternate archiver
# ------------------------------------------------------------------------

if test -n "${SPECIAL_AR}" && test -z "${ALTERNATE_AR}"; then
  ALTERNATE_AR="${SPECIAL_AR}"
fi

AC_MSG_CHECKING(for special archiver command)
if test -n "${ALTERNATE_AR}"; then
   AC_MSG_RESULT([${ALTERNATE_AR}])
   AM_CONDITIONAL(USE_ALTERNATE_AR, true)
else
   AC_MSG_RESULT([none])
   AM_CONDITIONAL(USE_ALTERNATE_AR, false)
fi
AC_SUBST(ALTERNATE_AR)



# ------------------------------------------------------------------------
# MPI cpp check
# ------------------------------------------------------------------------

if test "X${HAVE_PKG_MPI}" = "Xyes"; then

  if test -n "${MPI_DIR}" && test -z "${MPI_INC}"; then
    MPI_INC="${MPI_DIR}/include"
  fi

  if test -n "${MPI_INC}"; then
    CPPFLAGS="${CPPFLAGS} -I${MPI_INC}"
  fi

  AC_LANG_CPLUSPLUS 
  AC_MSG_CHECKING(for mpi.h)
  AC_TRY_CPP(mpi.h,
    [AC_MSG_RESULT(yes)], 
    [AC_MSG_RESULT(no)  
     AC_MSG_ERROR(cannot find mpi.h)]
  )

fi

# ------------------------------------------------------------------------
# Checks for header files
# ------------------------------------------------------------------------

#iostream
AC_LANG_CPLUSPLUS 
AC_MSG_CHECKING([for iostream])
AC_TRY_CPP(iostream,
[
  AC_MSG_RESULT(yes)
  AC_DEFINE(HAVE_IOSTREAM,,[define if we can include <iostream>])
], 
[
  AC_MSG_RESULT(no)
  AC_MSG_CHECKING([for iostream.h])
  AC_TRY_CPP(iostream.h,
  [
    AC_MSG_RESULT(yes)
    AC_DEFINE(HAVE_IOSTREAM_H,,[define if we can include <iostream.h>])
  ], 
  [
    AC_MSG_RESULT(no)
    AC_MSG_NOTICE([You do not appear to have iostream nor iostream.h!!])
  ]
  )
]
)

#iomanip
AC_LANG_CPLUSPLUS 
AC_MSG_CHECKING([for iomanip])
AC_TRY_CPP(iomanip,
[
  AC_MSG_RESULT(yes)
  AC_DEFINE(HAVE_IOMANIP,,[define if we can include <iomanip>])
], 
[
  AC_MSG_RESULT(no)
  AC_MSG_CHECKING([for iomanip.h])
  AC_TRY_CPP(iomanip.h,
  [
    AC_MSG_RESULT(yes)
    AC_DEFINE(HAVE_IOMANIP_H,,[define if we can include <iomanip.h>])
  ], 
  [
    AC_MSG_RESULT(no)
    AC_MSG_NOTICE([You do not appear to have iomanip nor iomanip.h!!])
  ]
  )
]
)

#cmath

# Note that the check for cmath had to be changed from AC_TRY_CPP to
# AC_TRY_COMPILE to accomodate the TFLOPS cross-compilation.

AC_LANG_CPLUSPLUS 
AC_MSG_CHECKING(for cmath)
AC_TRY_COMPILE(
[#include <cmath>],
[],
[
  AC_MSG_RESULT(yes)
  AC_DEFINE(HAVE_CMATH,,[define if we can include <cmath>])
], 
[
  AC_MSG_RESULT(no)
  AC_MSG_CHECKING([for math.h])
  AC_TRY_CPP(math.h,
  [
    AC_MSG_RESULT(yes)
    AC_DEFINE(HAVE_MATH_H,,[define if we can include <math.h>])
  ], 
  [
    AC_MSG_RESULT(no)
    AC_MSG_NOTICE([You do not appear to have cmath nor math.h!!])
  ]
  )
]
)

# can we say using std::ostream?
AC_LANG_CPLUSPLUS 
AC_MSG_CHECKING(whether "using std::ostream" is okay)
AC_TRY_COMPILE(
[
#ifdef HAVE_IOSTREAM
#include <iostream>
#elif defined(HAVE_IOSTREAM_H)
#include <iostream.h>
#endif
],
[using std::ostream;],
[AC_MSG_RESULT(yes)], 
[
  AC_MSG_RESULT(no)  
  AC_DEFINE(NO_USING_STD_OSTREAM,,[define if "using std::ostream" is okay])
]
)

# ------------------------------------------------------------------------
# With Trilinos
# ------------------------------------------------------------------------


AC_ARG_ENABLE(epetra,
[  --enable-epetra         compile noxepetra library],
[BUILD_LIBNOXEPETRA=$enableval],
[BUILD_LIBNOXEPETRA=no]
)

TRILINOS_DIR=${TRILINOS_HOME}

AC_ARG_WITH(trilinos,
[  --with-trilinos         specify Trilinos root directory
                          (uses TRILINOS_HOME if no argument)],
[
  BUILD_LIBNOXEPETRA=yes
  if test ${withval} = yes; then
    TRILINOS_DIR=${TRILINOS_HOME}
  else
    TRILINOS_DIR=${withval}
  fi  
]
)

AC_MSG_CHECKING(whether we should build the noxepetra library)
AC_MSG_RESULT(${BUILD_LIBNOXEPETRA})

# ------------------------------------------------------------------------
# With Petsc
# ------------------------------------------------------------------------


AC_ARG_ENABLE(petsc,
[  --enable-petsc         compile noxpetsc library],
[BUILD_LIBNOXPETSC=$enableval],
[BUILD_LIBNOXPETSC=no]
)

#PETSC_BOPT=g_c++
#PETSC_ARCH=linux

AC_MSG_CHECKING(whether we should build the noxpetsc library)
AC_MSG_RESULT(${BUILD_LIBNOXPETSC})

# ------------------------------------------------------------------------
# Trilinos include file check for building noxepetra library
# ------------------------------------------------------------------------

AC_DEFUN(TRILINOS_INCLUDE_CHECK,
[
  AC_ARG_WITH($1-include,
  [  --with-$1-include    specify location of $1 header files],
  [$1_INC=$withval])

  AC_MSG_CHECKING(for $1 include directory)
  if test -n "${$1_INC}" && test -f ${$1_INC}/$2; then
    AC_MSG_RESULT([${$1_INC}])
  elif test -f ${TRILINOS_DIR}/include/$2; then
    $1_INC=${TRILINOS_DIR}/include
    AC_MSG_RESULT([${$1_INC}])
  elif test -f ${TRILINOS_DIR}/include/$1/$2; then
    $1_INC=${TRILINOS_DIR}/include/$1
    AC_MSG_RESULT([${$1_INC}])
  elif test -f ${TRILINOS_DIR}/packages/$1/src/$2; then
    $1_INC=${TRILINOS_DIR}/packages/$1/src
    AC_MSG_RESULT([${$1_INC}])
  else
    AC_MSG_RESULT(no)
    AC_MSG_ERROR([Cannot find $1 include directory containing $2])
  fi  

])

if test X${BUILD_LIBNOXEPETRA} = Xyes; then

  TRILINOS_INCLUDE_CHECK(epetra,Epetra_Vector.h)
  TRILINOS_INCLUDE_CHECK(aztecoo,AztecOO.h)

  if test "X${HAVE_PKG_MPI}" = "Xyes"; then
    TRILINOS_COMM="MPI"
  else
    TRILINOS_COMM="SERIAL"
  fi

  TRILINOS_CXXFLAGS="-DEPETRA_${TRILINOS_COMM} -I${epetra_INC} -I${aztecoo_INC}"

fi

AC_SUBST(TRILINOS_CXXFLAGS)

# ------------------------------------------------------------------------
# Petsc include file check for building noxpetsc library
# ------------------------------------------------------------------------

AC_DEFUN(PETSC_INCLUDE_CHECK,
[
  AC_ARG_WITH($1-include,
  [  --with-$1-include    specify location of $1 header files],
  [$1_INC=$withval])

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
  else
    AC_MSG_RESULT(no)
    AC_MSG_ERROR([Cannot find $1 include directory containing $2. Make sure PETSC_DIR and PETSC_ARCH are defined.])
  fi

])

if test X${BUILD_LIBNOXPETSC} = Xyes; then

  PETSC_INCLUDE_CHECK(petsc,petsc.h)
  PETSC_INCLUDE_CHECK(petsccommon,base)
  PETSC_INCLUDE_CHECK(petscbmake,petscconf.h)

#  if test "X${HAVE_PKG_MPI}" = "Xyes"; then
#    TRILINOS_COMM="MPI"
#  else
#    TRILINOS_COMM="SERIAL"
#  fi
#
  PETSC_CXXFLAGS="-I${petsc_INC} -I${petsccommon_INC} -I${petscbmake_INC}"

  AC_SUBST(PETSC_CXXFLAGS)

fi

#AC_SUBST(PETSC_CXXFLAGS)

# ------------------------------------------------------------------------
# Checks for building simple examples
# ------------------------------------------------------------------------

AC_ARG_ENABLE(examples,
[  --enable-examples       compile simple examples],
[BUILD_EXAMPLES=${enableval}],
[BUILD_EXAMPLES=no]
)

AC_MSG_CHECKING(whether we should build the basic examples)
AC_MSG_RESULT(${BUILD_EXAMPLES})
AM_CONDITIONAL(BUILD_EXAMPLES, test X${BUILD_EXAMPLES} = Xyes)

# ------------------------------------------------------------------------
# Checks for building epetra examples
# ------------------------------------------------------------------------

AC_ARG_ENABLE(epetra-examples,
[  --enable-epetra-examples     compile epetra examples],
[BUILD_EPETRA_EXAMPLES=${enableval}],
[BUILD_EPETRA_EXAMPLES=no]
)

# Cannot build the example if we don't build the library
if test X${BUILD_LIBNOXEPETRA} != Xyes; then
   BUILD_EPETRA_EXAMPLES=no
fi

AC_MSG_CHECKING(whether we should build the noxepetra examples)
AC_MSG_RESULT(${BUILD_EPETRA_EXAMPLES})

# ------------------------------------------------------------------------
# Checks for building petsc examples
# ------------------------------------------------------------------------

AC_ARG_ENABLE(petsc-examples,
[  --enable-petsc-examples     compile petsc examples],
[BUILD_PETSC_EXAMPLES=${enableval}],
[BUILD_PETSC_EXAMPLES=no]
)

# Cannot build the example if we don't build the library
if test X${BUILD_LIBNOXPETSC} != Xyes; then
   BUILD_PETSC_EXAMPLES=no
fi

AC_MSG_CHECKING(whether we should build the noxpetsc examples)
AC_MSG_RESULT(${BUILD_PETSC_EXAMPLES})

# ------------------------------------------------------------------------
# Checks for basic libraries for examples
# ------------------------------------------------------------------------

if test ${BUILD_EXAMPLES} = yes || test ${BUILD_EPETRA_EXAMPLES} = yes; then

  AC_SEARCH_LIBS(pow,[m],,AC_MSG_ERROR(Cannot find math library))
  AC_SEARCH_LIBS(sqrt,[m],,AC_MSG_ERROR(Cannot find math library))

  dnl Checks for libraries.
  AC_F77_LIBRARY_LDFLAGS
  LIBS="${LIBS} ${FLIBS}"

  # Find valid BLAS
  AC_LANG_FORTRAN77
  AC_SEARCH_LIBS(dcopy,[cxml blas],,AC_MSG_ERROR([Cannot find BLAS library]))

  # Find valid LAPACK
  AC_LANG_FORTRAN77
  AC_SEARCH_LIBS(dgetrs,[cxml lapack complib.sgimath],,AC_MSG_ERROR([Cannot find LAPACK library]))

fi

# ------------------------------------------------------------------------
# Check for MPI library
# ------------------------------------------------------------------------

if test "X${HAVE_PKG_MPI}" = "Xyes" && test X${BUILD_EPETRA_EXAMPLES} = Xyes; then

  if test -n "${MPI_DIR}" && test -z "${MPI_LIBDIR}"; then
    MPI_LIBDIR="${MPI_DIR}/lib"
  fi

  if test -n "${MPI_LIBDIR}"; then
    LDFLAGS="${LDFLAGS} -L${MPI_LIBDIR}"
  fi

  if test -z "${MPI_LIBS}" && test -n "${MPI_LIBDIR}"; then
    MPI_LIBS="-lmpi"
  fi

  if test -n "${MPI_LIBS}"; then
    LIBS="${MPI_LIBS} ${LIBS}"
  fi

  AC_LANG_CPLUSPLUS 
  AC_MSG_CHECKING(whether MPI will link using C++ compiler)
  AC_TRY_LINK([#include <mpi.h>],
  [int c; char** v; MPI_Init(&c,&v);],
  [AC_MSG_RESULT(yes)], 
  [AC_MSG_RESULT(no)  
   AC_MSG_ERROR(MPI cannot link)]
  )

fi

# ------------------------------------------------------------------------
# Checks for Trilinos libraries
# ------------------------------------------------------------------------

AC_DEFUN(TRILINOS_LIB_CHECK,
[
  AC_ARG_WITH($1-libdir,
  [  --with-$1-libdir    specify location of $1 library],
  [$1_LIBDIR=$withval])

  AC_MSG_CHECKING(for $1)
  if test -n "${$1_LIBDIR}" && test -f ${$1_LIBDIR}/lib$1.a; then
    LDFLAGS="-L${$1_LIBDIR} ${LDFLAGS}"
    AC_MSG_RESULT(yes)
  elif test -f ${TRILINOS_DIR}/lib/lib$1.a; then
    LDFLAGS="-L${TRILINOS_DIR}/lib ${LDFLAGS}"
    AC_MSG_RESULT(yes)
  elif test -f ${TRILINOS_DIR}/lib/$1/lib$1.a; then
    LDFLAGS="-L${TRILINOS_DIR}/lib/$1 ${LDFLAGS}"
    AC_MSG_RESULT(yes)
  elif test -f ${TRILINOS_DIR}/lib/${TRILINOS_TARGET}/lib$1.a; then
    LDFLAGS="-L${TRILINOS_DIR}/lib/${TRILINOS_TARGET} ${LDFLAGS}"
    AC_MSG_RESULT(yes)
  else
    AC_MSG_RESULT(no)
    AC_MSG_ERROR([Cannot find $1 library])
  fi  

  LIBS="-l$1 ${LIBS}"

])

if test X${BUILD_EPETRA_EXAMPLES} = Xyes; then

  # Determine TRILINOS_TARGET
  AC_MSG_CHECKING(whether TRILINOS_ARCH is defined)
  if test -n "${TRILINOS_ARCH}"; then
    AC_MSG_RESULT([yes (${TRILINOS_ARCH})])
  else
    case $target in
      rs6000-ibm-aix3.2.* | rs6000-ibm-aix4.* | powerpc-ibm-aix4.*)
   	TRILINOS_ARCH="AIX"
	;;
      alpha-unknown-linux-gnu)
 	TRILINOS_ARCH="LINUXALPHA"
 	;;	
      *)
 	TRILINOS_ARCH="LINUX"
   	;;
    esac
    AC_MSG_RESULT([no (using ${TRILINOS_ARCH})])
  fi
 
  AC_MSG_CHECKING(whether TRILINOS_ID is defined)
  if test -n "${TRILINOS_ID}"; then
     AC_MSG_RESULT([yes (${TRILINOS_ID})])
  else
     AC_MSG_RESULT(no)
  fi
 
  AC_MSG_CHECKING(whether TRILINOS_TARGET is defined)
  if test -n "${TRILINOS_TARGET}"; then
     AC_MSG_RESULT([yes (${TRILINOS_TARGET})])
  else
     TRILINOS_TARGET="${TRILINOS_ARCH}.${TRILINOS_COMM}${TRILINOS_ID}"
     AC_MSG_RESULT([no (using ${TRILINOS_TARGET})])
  fi

  # Find a valid y12m
  AC_LANG_FORTRAN77

  AC_ARG_WITH(y12m-libdir,
  [  --with-y12m-libdir      specify location of y12m library],
  [y12m_LIBDIR=$withval])

  AC_MSG_CHECKING(for y12m)
  if test -n "${y12m_LIBDIR}" && test -f ${y12m_LIBDIR}/liby12m.a; then
    LDFLAGS="-L${y12m_LIBDIR} ${LDFLAGS}"
    AC_MSG_RESULT(yes)
  elif test -f ${TRILINOS_DIR}/external/y12m/lib/${TRILINOS_TARGET}/liby12m.a; then
    LDFLAGS="-L${TRILINOS_DIR}/external/y12m/lib/${TRILINOS_TARGET} ${LDFLAGS}"
    AC_MSG_RESULT(yes)
  else
    AC_MSG_RESULT(no)
    AC_MSG_ERROR([Cannot find y12m library])
  fi  

  AC_SEARCH_LIBS(y12mae,[y12m],,AC_MSG_ERROR([Cannot find y12m library]))


  # Check for ifpack
  TRILINOS_LIB_CHECK(ifpack)
  TRILINOS_LIB_CHECK(aztecoo)
  TRILINOS_LIB_CHECK(epetra)

fi

# ------------------------------------------------------------------------
# Checks for TRILINOS status 
# ------------------------------------------------------------------------

AM_CONDITIONAL(BUILD_LIBNOXEPETRA, test X${BUILD_LIBNOXEPETRA} = Xyes)
AM_CONDITIONAL(BUILD_EPETRA_EXAMPLES, test X${BUILD_EPETRA_EXAMPLES} = Xyes)

# ------------------------------------------------------------------------
# Checks for Petsc libraries
# ------------------------------------------------------------------------

AC_DEFUN(PETSC_LIB_CHECK,
[
  AC_ARG_WITH($1-libdir,
  [  --with-$1-libdir    specify location of $1 library],
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
    
if test X${BUILD_PETSC_EXAMPLES} = Xyes; then
  
  # Check for petsc
  PETSC_LIB_CHECK(petsc)

fi

# ------------------------------------------------------------------------
# Checks for PETSC status 
# ------------------------------------------------------------------------

AM_CONDITIONAL(BUILD_LIBNOXPETSC, test X${BUILD_LIBNOXPETSC} = Xyes)
AM_CONDITIONAL(BUILD_PETSC_EXAMPLES, test X${BUILD_PETSC_EXAMPLES} = Xyes)

# ------------------------------------------------------------------------
# LOCA Library Check
# ------------------------------------------------------------------------

AC_ARG_ENABLE(loca,
[  --enable-loca         compile the loca library],
[BUILD_LIBLOCA=${enableval}],
[BUILD_LIBLOCA=no]
)
AC_MSG_CHECKING(whether we should build the loca library)
AC_MSG_RESULT([${BUILD_LIBLOCA}])
AM_CONDITIONAL(BUILD_LIBLOCA, test X${BUILD_LIBLOCA} = Xyes)

# ------------------------------------------------------------------------
# LOCA-Epetra Library Check
# ------------------------------------------------------------------------

AC_ARG_ENABLE(loca-epetra,
[  --enable-loca-epetra      compile the locaepetra library],
[BUILD_LIBLOCAEPETRA=${enableval}],
[BUILD_LIBLOCAEPETRA=no]
)

AC_MSG_CHECKING(whether we should build the locaepetra library)
AC_MSG_RESULT([${BUILD_LIBLOCAEPETRA}])
AM_CONDITIONAL(BUILD_LIBLOCAEPETRA, test X${BUILD_LIBLOCAEPETRA} = Xyes)

# ------------------------------------------------------------------------
# Checks for types
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for structures
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for libraries
# ------------------------------------------------------------------------

# THIS TEST ISN'T WORKING - WILL COME BACK TO IT LATER.
# Check that all libraries will link with C++
dnl if test X${LIBS} != X; then
dnl   AC_MSG_CHECKING([[whether libraries (${LIBS}) will link]])
dnl   AC_LANG_CPLUSPLUS 
dnl   AC_TRY_LINK([],[],
dnl   [AC_MSG_RESULT(yes)],
dnl   [
dnl   AC_MSG_RESULT(no)
dnl   AC_MSG_ERROR(cannot link)
dnl   ]
dnl   )
dnl fi

# ------------------------------------------------------------------------
# Checks for compiler characteristics
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for linker characteristics
# ------------------------------------------------------------------------

# -- DEC CXX --
case $target in
  alpha*)
    if test CXX = cxx; then
      CXXREPOS="cxx_repository/*.o"
    fi
  ;;
esac

AC_MSG_CHECKING([for C++ repository files])
if test X$CXXREPOS = X; then
AC_MSG_RESULT(none)
else
AC_MSG_RESULT($CXXREPOS)
fi

AC_SUBST(CXXREPOS)

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
		doc/Makefile
		src/Makefile
		src-epetra/Makefile 
		example/Makefile
		example-epetra/Makefile 
		example-epetra/DS6.5.1/Makefile
		example-epetra/1Dfem_linear/Makefile
		example-epetra/1Dfem_nonlinear/Makefile
		example-epetra/1Dtransfem_nonlinear/Makefile
                src-petsc/Makefile
                example-petsc/Makefile
                example-petsc/DS6.5.1/Makefile
		src-loca/Makefile
                src-loca/src-epetra/Makefile
		])

AC_OUTPUT()


