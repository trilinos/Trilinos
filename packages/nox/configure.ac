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

AC_INIT(nox, 1.0a-dev, nox-developers@software.sandia.gov)

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
# Check for user overrides (some to be used later)
# ------------------------------------------------------------------------

if test -n "${CXXFLAGS}"; then
  USERDEF_CXXFLAGS=yes
fi

if test -n "${CFLAGS}"; then
  USERDEF_CFLAGS=yes
fi

if test -n "${CPPFLAGS}"; then
  USERDEF_CPPFLAGS=yes
fi

if test -n "${F77}"; then
  USERDEF_F77=yes
fi

if test -n "${FFLAGS}"; then
  USERDEF_FFLAGS=yes
fi

if test -n "${LDFLAGS}"; then
  USERDEF_LDFLAGS=yes
fi

# ------------------------------------------------------------------------
# Checks for extra flags, etc
# ------------------------------------------------------------------------

# Since this is the first "with" command, we skip a line

AC_ARG_WITH(cppflags,
[AC_HELP_STRING([--with-cppflags],[add extra preprocessor flags (e.g., "-Dflag -Idir")])],
[EXTRA_CPPFLAGS=${withval}]
)


AC_ARG_WITH(cxxflags,
[AC_HELP_STRING([--with-cxxflags],[add extra C++ compiler flags])],
[EXTRA_CXXFLAGS=${withval}]
)


AC_ARG_WITH(ldflags,
[AC_HELP_STRING([--with-ldflags],[add extra linker (typically -L) flags])],
[EXTRA_LDFLAGS=${withval}]
)

# ------------------------------------------------------------------------
# Alternate archiver
# ------------------------------------------------------------------------


AC_ARG_WITH(ar,
[AC_HELP_STRING([--with-ar],[specify archiver command @<:@"ar cru"@:>@])],
[
AC_MSG_CHECKING(user-defined archiver)
AC_MSG_RESULT([${withval}])
ALTERNATE_AR="${withval}"
]
)

# ------------------------------------------------------------------------
# MPI configure command line option
# ------------------------------------------------------------------------

AC_ARG_ENABLE(mpi,
[AC_HELP_STRING([--enable-mpi],[MPI support])],
[HAVE_PKG_MPI=$enableval],
[HAVE_PKG_MPI=no]
)

AC_ARG_WITH(mpi-cxx,
[AC_HELP_STRING([--with-mpi-cxx],[use MPI C++ compiler (enables MPI) @<:@mpiCC@:>@])],
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
[AC_HELP_STRING([--with-mpi=MPIROOT],[use MPI root directory (enables MPI)])],
[
  HAVE_PKG_MPI=yes
  MPI_DIR=${withval}
  AC_MSG_CHECKING(MPI directory)
  AC_MSG_RESULT([${MPI_DIR}])
]
)


AC_ARG_WITH(mpi-libs,
[AC_HELP_STRING([--with-mpi-libs="LIBS"],[MPI libraries @<:@"-lmpi"@:>@])],
[
  MPI_LIBS=${withval}
  AC_MSG_CHECKING(user-defined MPI libraries)
  AC_MSG_RESULT([${MPI_LIBS}])
]
)


AC_ARG_WITH(mpi-incdir,
[AC_HELP_STRING([--with-mpi-incdir=DIR],[MPI include directory @<:@MPIROOT/include@:>@])],
[
  MPI_INC=${withval}
  AC_MSG_CHECKING(user-defined MPI includes)
  AC_MSG_RESULT([${MPI_INC}])
]
)


AC_ARG_WITH(mpi-libdir,
[AC_HELP_STRING([--with-mpi-libdir=DIR],[MPI library directory @<:@MPIROOT/lib@:>@])],
[
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
    echo "-----"
    echo "Cannot find MPI C++ compiler ${MPI_CXX}."
    echo "Specify with --with-mpi-cxx."
    echo "-----"
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
[AC_HELP_STRING([--enable-debug],[compile with debugging])],
[DEBUG=${enableval}],
[DEBUG="default"]
)
AC_MSG_CHECKING(whether debugging is on)
AC_MSG_RESULT([${DEBUG}])

AC_ARG_ENABLE(opt,
[AC_HELP_STRING([--enable-opt],[compile with optimization])],
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


# Checks for special cases...
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
     F77_MISC=f77
     FFLAGS_MISC="-64"
     SPECIAL_AR="CC -ar -W,l -W,s -o"    
    fi

  ;;

  # Solaris test
  sparc-sun-solaris* ) 

    if test "${CXX}" = "CC"; then
     SPECIAL_AR="CC -xar -o"    
    fi

  ;;

esac


# Process the C++ flags
AC_MSG_CHECKING([for C++ flags])

# Test for user-defined CXXFLAGS
if test "X${USERDEF_CXXFLAGS}" = "Xyes"; then

  AC_MSG_RESULT([${CXXFLAGS} (user-specified via CXXFLAGS)])

# Otherwise, determine the flags automatically
else

  # Test for adding debug flags
  if test ${DEBUG} = yes; then

    CXXFLAGS=${CXXFLAGS_DBG}

  # Otherwise, test for adding optimization flags
  elif test ${OPT} = yes; then

    CXXFLAGS=${CXXFLAGS_OPT}

  fi

  # Add on any additional automatically determined flags
  if test -n "${CXXFLAGS_MISC}"; then
    CXXFLAGS="${CXXFLAGS} ${CXXFLAGS_MISC}"
  fi			  

  # Print the result
  AC_MSG_RESULT([${CXXFLAGS}])

fi

# Check for user-defined F77 compiler (used only for compiler checks)
if test "X${USERDEF_F77}" != "Xyes"; then
  if test -n "${F77_MISC}"; then
    F77="${F77_MISC}"
  fi			  
fi

# Test for user-defined FFLAGS flags
AC_MSG_CHECKING([for F77 flags])

if test "X${USERDEF_FFLAGS}" = "Xyes"; then

  AC_MSG_RESULT([${FFLAGS} (user-specified via FFLAGS)])

# Otherwise, determine the flags automatically
else

  if test -n "${FFLAGS_MISC}"; then
    FFLAGS="${FFLAGS} ${FFLAGS_MISC}"
  fi			  

  AC_MSG_RESULT([${FFLAGS}])

fi

# Add on any user-defined CXXFLAGS
AC_MSG_CHECKING([for extra C++ compiler flags])
if test -z "${EXTRA_CXXFLAGS}"; then
  AC_MSG_RESULT(none)
else
  AC_MSG_RESULT([${EXTRA_CXXFLAGS}])
  CXXFLAGS="${CXXFLAGS} ${EXTRA_CXXFLAGS}"
fi

# Add on any user-defined CPPFLAGS
AC_MSG_CHECKING([for extra preprocessor flags])
if test -z "${EXTRA_CPPFLAGS}"; then
  AC_MSG_RESULT(none)
else
  AC_MSG_RESULT([${EXTRA_CPPFLAGS}])
  CPPFLAGS="${CPPFLAGS} ${EXTRA_CPPFLAGS}"
fi

# Add on any user-defined LDFLAGS
AC_MSG_CHECKING([for extra linker flags])
if test -z "${EXTRA_LDFLAGS}"; then
  AC_MSG_RESULT(none)
else
  AC_MSG_RESULT([${EXTRA_LDFLAGS}])
  LDFLAGS="${LDFLAGS} ${EXTRA_LDFLAGS}"
fi

# ------------------------------------------------------------------------
# Set up alternate archiver (option -with-ar above)
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
  AC_TRY_CPP([#include "mpi.h"],
    [
     AC_MSG_RESULT(yes)
     AC_DEFINE(HAVE_MPI,,[define that mpi is being used])
    ], 
    [
     AC_MSG_RESULT(no)  
     echo "---"
     echo "Cannot find header file mpi.h."
     echo "Either compile without mpi (no --with-mpi directives)"
     echo "or specify --with-mpi-incdir=<dir>."
     echo "---"
     AC_MSG_ERROR(cannot find mpi.h)
    ]
  )

fi

# ------------------------------------------------------------------------
# Checks for header files
# ------------------------------------------------------------------------

#iostream
AC_LANG_CPLUSPLUS 
AC_MSG_CHECKING([for iostream])
AC_TRY_CPP([#include <iostream>],
[
  AC_MSG_RESULT(yes)
  AC_DEFINE(HAVE_IOSTREAM,,[define if we can include <iostream>])
], 
[
  AC_MSG_RESULT(no)
  AC_MSG_CHECKING([for iostream.h])
  AC_TRY_CPP([#include "iostream.h"],
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
AC_TRY_CPP([#include "iomanip"],
[
  AC_MSG_RESULT(yes)
  AC_DEFINE(HAVE_IOMANIP,,[define if we can include <iomanip>])
], 
[
  AC_MSG_RESULT(no)
  AC_MSG_CHECKING([for iomanip.h])
  AC_TRY_CPP([#include "iomanip.h"],
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
  AC_TRY_CPP([#include "math.h"],
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
   echo "-----"
   echo "Cannot link simple MPI program."
   echo "Try --with-mpi-cxx to specify MPI C++ compile script."
   echo "Or try --with-mpi-libs, --with-mpi-incdir, --with-mpi-libdir"
   echo "to specify all the specific MPI compile options."
   echo "-----"
   AC_MSG_ERROR(MPI cannot link)]
  )

fi

# ------------------------------------------------------------------------
# Which Add-Ons To Build
# ------------------------------------------------------------------------

BUILD_ANY_EXAMPLES=no
AC_DEFUN(NOX_INTERFACE,
[
  
  AC_ARG_ENABLE(nox-$1,
  [AC_HELP_STRING([--enable-nox-$1],[compile $1 interace libraries])],
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
# NOX-PETSC
# ------------------------------------------------------------------------

NOX_INTERFACE(petsc)
AM_CONDITIONAL(BUILD_LIBNOXPETSC, test X${ADDON_petsc} = Xyes)
AM_CONDITIONAL(BUILD_PETSC_EXAMPLES, test X${ADDON_petsc_EXAMPLES} = Xyes)

# ------------------------------------------------------------------------
# LOCA Library Check
# ------------------------------------------------------------------------

AC_ARG_ENABLE(loca,
[AC_HELP_STRING([--enable-loca],[compile the loca library])],
[BUILD_LIBLOCA=${enableval}],
[BUILD_LIBLOCA=no]
)
AC_MSG_CHECKING(whether we should build the loca library)
AC_MSG_RESULT([${BUILD_LIBLOCA}])
AM_CONDITIONAL(BUILD_LIBLOCA, test X${BUILD_LIBLOCA} = Xyes)


# ------------------------------------------------------------------------
# LOCA-LAPACK Library Check
# ------------------------------------------------------------------------

AC_MSG_CHECKING(whether we should build the localapack library)
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

AC_MSG_CHECKING(whether we should build the locaepetra library)
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
# Checks for basic libraries for examples and/or lapack support
# ------------------------------------------------------------------------

if test ${BUILD_ANY_EXAMPLES} = yes; then

  AC_SEARCH_LIBS(pow,[m],,AC_MSG_ERROR(Cannot find math library))
  AC_SEARCH_LIBS(sqrt,[m],,AC_MSG_ERROR(Cannot find math library))

  dnl Checks for libraries.
  AC_F77_LIBRARY_LDFLAGS
  LIBS=" ${FLIBS} ${LIBS}"

  # Define F77_FUNC - used for correctly munging Fortran function calls
  AC_F77_WRAPPERS

  # Find valid BLAS
  ACX_BLAS([LIBS="${BLAS_LIBS} ${LIBS}"],[AC_MSG_ERROR([Cannot find BLAS library])])

  # Find valid LAPACK
  ACX_LAPACK([LIBS=" ${LAPACK_LIBS} ${LIBS}"],[AC_MSG_ERROR([Cannot find LAPACK library])])

fi

# ------------------------------------------------------------------------
# Trilinos include file check for building noxepetra library
# ------------------------------------------------------------------------

AC_ARG_WITH(trilinos,
[AC_HELP_STRING([--with-trilinos],[building inside Trilinos directory structure])],
[INSIDE_TRILINOS=$withval])

AC_ARG_WITH(epetra-incdir,
[AC_HELP_STRING([--with-epetra-incdir],[specify location of all epetra header files])],
[epetra_INCDIR=$withval])

AC_ARG_WITH(aztecoo-incdir,
[AC_HELP_STRING([--with-aztecoo-incdir],[specify location of all aztecoo header files])],
[aztecoo_INCDIR=$withval])

AC_ARG_WITH(ifpack-incdir,
[AC_HELP_STRING([--with-ifpack-incdir],[specify location of all ifpack header files])],
[ifpack_INCDIR=$withval])

AC_ARG_WITH(epetra-libdir,
[AC_HELP_STRING([--with-epetra-libdir],[specify location of epetra library])],
[epetra_LIBDIR=$withval])

AC_ARG_WITH(aztecoo-libdir,
[AC_HELP_STRING([--with-aztecoo-libdir],[specify location of aztecoo library])],
[aztecoo_LIBDIR=$withval])

AC_ARG_WITH(ifpack-libdir,
[AC_HELP_STRING([--with-ifpack-libdir],[specify location of ifpack library])],
[ifpack_LIBDIR=$withval])

AC_ARG_WITH(y12m-libdir,
[AC_HELP_STRING([--with-y12m-libdir],[specify location of y12m library])],
[y12m_LIBDIR=$withval])


# Relative path to Trilinos include files, to be used in Makefiles
REL_TRILINOS_CPPFLAGS=

# Absolute path to Trilinos include files from location that configure
# is run - to be used in configuration tests
ABS_TRILINOS_CPPFLAGS=

# Determine if anything needs to be added to REL/ABS_TRILINOS_FLAGS.
# Does not modify CPPFLAGS.
AC_DEFUN(TRILINOS_HEADER,
[
  save_CPPFLAGS=$CPPFLAGS
  done=no

  # First try the directories we already have
  CPPFLAGS="$save_CPPFLAGS $ABS_TRILINOS_CPPFLAGS"
  AC_TRY_CPP([#include "$2"],[done=yes],)
  CPPFLAGS=$save_CPPFLAGS

  # Next  try the specified incdir (if any)
  if test $done = no && test -n ${$1_INCDIR}; then
    CPPFLAGS="-I${$1_INCDIR} $save_CPPFLAGS $ABS_TRILINOS_CPPFLAGS"
    AC_TRY_CPP([#include "$2"],
    [ done=yes
      REL_TRILINOS_CPPFLAGS="-I${$1_INCDIR} ${REL_TRILINOS_CPPFLAGS}"
      ABS_TRILINOS_CPPFLAGS="-I${$1_INCDIR} ${ABS_TRILINOS_CPPFLAGS}"
    ],[])
    CPPFLAGS=$save_CPPFLAGS
  fi   

  # Next try the Trilinos build dir
  if test $done = no; then
    CPPFLAGS="-I../$1/src $save_CPPFLAGS $ABS_TRILINOS_CPPFLAGS"
    AC_TRY_CPP([#include "$2"],
    [ done=yes
      REL_TRILINOS_CPPFLAGS="-I\$(top_builddir)/../$1/src ${REL_TRILINOS_CPPFLAGS}"
      ABS_TRILINOS_CPPFLAGS="-I../$1/src ${ABS_TRILINOS_CPPFLAGS}"
    ],[])
    CPPFLAGS=$save_CPPFLAGS
  fi   

  # Next try the Trilinos source dir
  if test $done = no && test -n $srcdir; then
    CPPFLAGS="-I$srcdir/../$1/src $save_CPPFLAGS $ABS_TRILINOS_CPPFLAGS"
    AC_TRY_CPP([#include "$2"],
    [ done=yes
      REL_TRILINOS_CPPFLAGS="-I\$(top_srcdir)/../$1/src ${REL_TRILINOS_CPPFLAGS}"
      ABS_TRILINOS_CPPFLAGS="-I$srcdir/../$1/src ${ABS_TRILINOS_CPPFLAGS}"
    ],[])					       
    CPPFLAGS=$save_CPPFLAGS
  fi   

  if test $done = no; then
    echo "------"
    echo "Cannot find the include file $2, part of the $1 package."
    echo "Try --with-$1-incdir=<dir> to specify a single directory or"
    echo "--with-cppflags=\"-I<dir1> -I<dir2>\" for multiple directories."
    echo "------"
    AC_MSG_ERROR([Cannot find $1 include directory containing $2])
  fi 
])

if test X${ADDON_epetra} = Xyes; then

  # We have to check for the standard Trilinos headers as well as
  # those created by configure.
  TRILINOS_HEADER(epetra,Epetra_config.h)
  TRILINOS_HEADER(epetra,Epetra_Vector.h)
  TRILINOS_HEADER(aztecoo,AztecOO_config.h)
  TRILINOS_HEADER(aztecoo,AztecOO.h)
  TRILINOS_HEADER(ifpack,Ifpack_config.h)
  TRILINOS_HEADER(ifpack,Ifpack_CrsRiluk.h)

  TRILINOS_CXXFLAGS="${REL_TRILINOS_CPPFLAGS}"

fi

AC_SUBST(TRILINOS_CXXFLAGS)

# ------------------------------------------------------------------------
# Checks for Trilinos libraries
# ------------------------------------------------------------------------

AC_DEFUN(TRILINOS_LIB_CHECK,
[
  LIBS="-l$1 ${LIBS}"
  save_LDFLAGS=$LDFLAGS
  done=no

  # First, see if we're inside Trilinos. 
  # If so, just assume the libraries will be built later. 
  if test X$INSIDE_TRILINOS = Xyes; then
    PWD=`pwd`
    LDFLAGS="-L$PWD/../$1/src $save_LDFLAGS"
    done=auto
  fi

  # Next, see if we can link without doing anything
  if test $done = no; then
    AC_LANG_CPLUSPLUS
    AC_TRY_LINK(,,[done=yes],)
    LDFLAGS=$save_LDFLAGS
  fi

  # Next try the specified libdir (if any)
   if test $done = no && test -n ${$1_LIBDIR}; then
     LDFLAGS="-L${$1_LIBDIR} $save_LDFLAGS"
     AC_TRY_LINK(,,[done=yes],[LDFLAGS=$save_LDFLAGS])
   fi   

  # Next try the Trilinos build dir
   if test $done = no; then
    PWD=`pwd`
    LDFLAGS="-L$PWD/../$1/src $save_LDFLAGS"
    AC_TRY_LINK(,,[done=yes],[LDFLAGS=$save_LDFLAGS])
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

if test X${ADDON_epetra_EXAMPLES} = Xyes || test X${BUILD_LOCA_EPETRA_EXAMPLES} = Xyes; then

  # Find the dreaded y12 library
  LIBS="-ly12m ${LIBS}"
  save_LDFLAGS=$LDFLAGS
  done=no

  # First, see if we can link without doing anything
  AC_LANG_CPLUSPLUS
  AC_TRY_LINK(,,[done=yes],)
  LDFLAGS=$save_LDFLAGS

  # Next try the specified libdir (if any)
   if test $done = no && test -n ${y12m_LIBDIR}; then
     LDFLAGS="-L${y12m_LIBDIR} $save_LDFLAGS"
     AC_TRY_LINK(,,[done=yes],[LDFLAGS=$save_LDFLAGS])
   fi   

  # Next try the Trilinos build dir for reference y12m
   if test $done = no; then
    PWD=`pwd`
    LDFLAGS="-L$PWD/../../external/y12m/src $save_LDFLAGS"
    AC_TRY_LINK(,,[done=yes],[LDFLAGS=$save_LDFLAGS])
   fi   

  AC_MSG_CHECKING([for y12m library])
  
  if test $done = yes; then
    AC_MSG_RESULT(yes)
  else
    AC_MSG_RESULT(no)
    echo "---"
    echo "Cannot find y12m library."
    echo "Try --with-y12m-libdir=<dir> to specify the library location."
    echo "If requiring the y12m library annoys you, vote for Bug 106 by visiting"
    echo "http://software.sandia.gov/bugzilla/show_bug.cgi?id=106."
    echo "---"
    AC_MSG_ERROR([Cannot find y12m library])
  fi  

  TRILINOS_LIB_CHECK(epetra)
  TRILINOS_LIB_CHECK(ifpack)
  TRILINOS_LIB_CHECK(aztecoo)

fi

# ------------------------------------------------------------------------
# Petsc include file check for building noxpetsc library
# ------------------------------------------------------------------------


AC_ARG_WITH(petsc-incdir,
[  --with-petsc-incdir     specify location of PETSc header files],
[petsc_INC=$withval])


AC_DEFUN(PETSC_INCLUDE_CHECK,
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
  else
    AC_MSG_RESULT(no)
    AC_MSG_ERROR([Cannot find $1 include directory containing $2. Make sure PETSC_DIR and PETSC_ARCH are defined.])
  fi

])

if test X${ADDON_petsc} = Xyes; then

  PETSC_INCLUDE_CHECK(petsc,petsc.h)
  PETSC_INCLUDE_CHECK(petsccommon,base)
  PETSC_INCLUDE_CHECK(petscbmake,petscconf.h)

  PETSC_CXXFLAGS="-I${petsc_INC} -I${petsccommon_INC} -I${petscbmake_INC}"

  AC_SUBST(PETSC_CXXFLAGS)

fi

# ------------------------------------------------------------------------
# Checks for Petsc libraries
# ------------------------------------------------------------------------

AC_DEFUN(PETSC_LIB_CHECK,
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
    
if test X${ADDON_petsc_EXAMPLES} = Xyes; then
  
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

# -- DEC CXX --
case $target in
  alpha*)
    if test $CXX = cxx; then
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
		doc/Doxyfile
		src/Makefile
                src-lapack/Makefile
		example-lapack/Makefile
		src-epetra/Makefile 
		example-epetra/Makefile 
		example-epetra/DS6.5.1/Makefile
		example-epetra/1Dfem_linear/Makefile
		example-epetra/1Dfem_nonlinear/Makefile
		example-epetra/1Dfem_nonlinearColoring/Makefile
		example-epetra/1Dtransfem_nonlinear/Makefile
                src-petsc/Makefile
                example-petsc/Makefile
                example-petsc/DS6.5.1/Makefile
		src-loca/Makefile
		src-loca/src/Makefile
                src-loca/src-epetra/Makefile
		src-loca/src-lapack/Makefile
                src-loca/example-epetra/Makefile
                src-loca/example-epetra/1Dfem/Makefile
		src-loca/example-lapack/Makefile
		src-loca/example-lapack/Chan/Makefile
		])

AC_OUTPUT()


# ------------------------------------------------------------------------
# Bye World!
# ------------------------------------------------------------------------

echo "----------------------------------------"
echo "Finished Running NOX Configure Script"
echo "----------------------------------------"

