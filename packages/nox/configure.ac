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

dnl Commented out this line for now because it causes headaches in
dnl committing.
dnl AC_REVISION($Revision$)

# ------------------------------------------------------------------------
# Check for user overrides
# ------------------------------------------------------------------------

if test -z "$CC"; then
  USERDEF_CC=yes
fi

if test -z "$CFLAGS"; then
  USERDEF_CFLAGS=yes
fi

if test -z "$CXX"; then
  USERDEF_CXX=yes
fi

if test -z "$CXXFLAGS"; then
  USERDEF_CXXFLAGS=yes
fi

# ------------------------------------------------------------------------
# Extra flags and includes
# ------------------------------------------------------------------------

AC_ARG_WITH(flags,
[  --with-flags            Specify add'l preprocessor flags ("-Dflag -Idir")],
[
XTRAFLAGS=$withval
echo Using extra flags: $withval
]
)

# ------------------------------------------------------------------------
# Debug configure command line option
# ------------------------------------------------------------------------

AC_ARG_ENABLE(debug,
[  --enable-debug          Compile with debugging options],
[
DEBUG=$enableval
case $enableval in

  yes) CXXFLAGS="${CXXFLAGS} -g"
       CFLAGS="${CFLAGS} -g"
       ;;

  no)  CXXFLAGS="${CXXFLAGS} -O2"
       CFLAGS="${CFLAGS} -O2"
       ;;
esac
],
[DEBUG=default]
)
AC_MSG_CHECKING(whether debugging is on)
AC_MSG_RESULT([${DEBUG}])

# ------------------------------------------------------------------------
# AR override command line option
# ------------------------------------------------------------------------

AC_ARG_WITH(ar,
[  --with-ar               Specify archiver (default is "ar cru")],
[
AC_MSG_CHECKING(user-defined archiver)
AC_MSG_RESULT([${withval}])
USE_ALTERNATE_AR=yes
ALTERNATE_AR="${withval}"
]
)

# ------------------------------------------------------------------------
# MPI configure command line option
# ------------------------------------------------------------------------

AC_ARG_ENABLE(mpi,
[  --enable-mpi            Enable MPI],
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

AC_MSG_CHECKING(whether we are using MPI)
AC_MSG_RESULT([${USE_MPI}])

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

# ------------------------------------------------------------------------
# Package-specific command line option
# ------------------------------------------------------------------------

AC_ARG_WITH(pkgdir,
[  --with-pkgdir           Specify location of add'l packages (i.e., epetra)],
[
USE_PKGDIR=yes
USE_TRILINOS=yes
PKGDIR=$withval
]
)

AC_ARG_ENABLE(trilinos,
[  --disable-trilinos      Disable TRILINOS],
[
case $enableval in
  yes)
    USE_TRILINOS=yes
  ;;
  no)
    USE_TRILINOS=no
  ;;
  *)
    AC_MSG_ERROR([Invalid value for --disable-trilinos ($enableval)])
  ;;
esac
],
[USE_TRILINOS=yes]
)

# ------------------------------------------------------------------------
# Checks for MPI compilers (mpicc/mpiCC)
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

# ------------------------------------------------------------------------
# Checks for programs
# ------------------------------------------------------------------------

AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB

# ------------------------------------------------------------------------
# Checks for extra flags, etc
# ------------------------------------------------------------------------


case $target in

  # IBM test
  rs6000-ibm-aix3.2.* | rs6000-ibm-aix4.* | powerpc-ibm-aix4.*)

    case $CC in
      xlc* | cc* ) COPTIONS="-qnolm";;
    esac

    case $CXX in
      xlC* ) CXXOPTIONS="-qnolm -qrtti";;
    esac

  ;;

  # SGI/IRIX test
  mips-sgi-irix* ) 

    if $CXX = CC; then
     COPTIONS="-LANG:std -LANG:ansi-for-init-scope=ON -64 -DSGI64"
     
     if X$USE_ALTERNATE_AR != Xyes; then
       USE_ALTERNATE_AR=yes
       ALTERNATE_AR="CC -ar -W,l -W,s -o"
     fi

    fi

  ;;

esac

AC_MSG_CHECKING([for extra preprocessor flags])
if test X$XTRAFLAGS = X; then
AC_MSG_RESULT(none)
else
AC_MSG_RESULT([$XTRAFLAGS])
fi

AC_MSG_CHECKING([for extra C flags])
if test X$COPTIONS = X; then
AC_MSG_RESULT(none)
else
AC_MSG_RESULT($COPTIONS)
fi

AC_MSG_CHECKING([for extra C++ flags])
if test X"$CXXOPTIONS" = X; then
AC_MSG_RESULT(none)
else
AC_MSG_RESULT($CXXOPTIONS)
fi

CFLAGS="${CFLAGS} ${COPTIONS}"
CXXFLAGS="${CXXFLAGS} ${CXXOPTIONS}"
CPPFLAGS="${CPPFLAGS} ${XTRAFLAGS}"

# ------------------------------------------------------------------------
# MPI link check
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
# Checks for header files
# ------------------------------------------------------------------------

#iostream
AC_LANG_CPLUSPLUS 
AC_MSG_CHECKING([for iostream])
AC_TRY_LINK([#include <iostream>],[],
[
  AC_MSG_RESULT(yes)
  AC_DEFINE(HAVE_IOSTREAM,,[define if we can include <iostream>])
], 
[
  AC_MSG_RESULT(no)
  AC_MSG_CHECKING([for iostream.h])
  AC_TRY_LINK([#include <iostream.h>],[],
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
AC_TRY_LINK([#include <iomanip>],[],
[
  AC_MSG_RESULT(yes)
  AC_DEFINE(HAVE_IOMANIP,,[define if we can include <iomanip>])
], 
[
  AC_MSG_RESULT(no)
  AC_MSG_CHECKING([for iomanip.h])
  AC_TRY_LINK([#include <iomanip.h>],[],
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
AC_LANG_CPLUSPLUS 
AC_MSG_CHECKING([for cmath])
AC_TRY_LINK([#include <cmath>],[],
[
  AC_MSG_RESULT(yes)
  AC_DEFINE(HAVE_CMATH,,[define if we can include <cmath>])
], 
[
  AC_MSG_RESULT(no)
  AC_MSG_CHECKING([for math.h])
  AC_TRY_LINK([#include <math.h>],[],
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
AC_LANG_SAVE

AC_LANG_CPLUSPLUS 
AC_MSG_CHECKING(whether "using std::ostream" is okay)
AC_TRY_LINK(
[
#ifdef HAVE_IOSTREAM
#include <iostream>
#elif defined(HAVE_IOSTREAM_H)
#include <iostream.h>
#else
// Try to include it anyway
#include <iostream>
#endif
],
[
using std::ostream;
],
[AC_MSG_RESULT(yes)], 
[
  AC_MSG_RESULT(no)  
  AC_DEFINE(NO_USING_STD_OSTREAM,,[define is "using std::ostream" is okay])
]
)

AC_LANG_RESTORE


# ------------------------------------------------------------------------
# Checks for libraries
# ------------------------------------------------------------------------

AC_SEARCH_LIBS(pow,[m])
AC_SEARCH_LIBS(sqrt,[m])
AC_CHECK_FUNCS([pow sqrt])


# ------------------------------------------------------------------------
# Trilinos
# ------------------------------------------------------------------------

# Case I: Using standard Trilinos directory structure
# (a) Check if TRILINOS_HOME is valid
if test ${USE_TRILINOS} = yes && test X${USE_PKGDIR} != Xyes; then

   AC_MSG_CHECKING(whether TRILINOS_HOME is defined)
   if test $TRILINOS_HOME; then
      AC_MSG_RESULT([yes (${TRILINOS_HOME})])
   else
      AC_MSG_RESULT([no (using ${HOME}/Trilinos)])
      TRILINOS_HOME=${HOME}/Trilinos
   fi

   AC_MSG_CHECKING(whether TRILINOS_HOME is valid)
   if test -d ${TRILINOS_HOME}; then
      AC_MSG_RESULT(yes)
   else
      AC_MSG_RESULT([no (WARNING: DISABLING TRILINOS!!!!)])
      USE_TRILINOS=no
   fi
fi

# Case I: Using standard Trilinos directory structure
# (b): Set Trilinos variables
if test ${USE_TRILINOS} = yes && test X${USE_PKGDIR} != Xyes; then

   AC_MSG_CHECKING(whether TRILINOS_ARCH is defined)
   if test ${TRILINOS_ARCH}; then
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

   AC_MSG_CHECKING(whether TRILINOS_COMM is defined)
   if test ${USE_MPI} = yes; then
      TRILINOS_COMM="MPI"
   else
      TRILINOS_COMM="SERIAL"
   fi
   AC_MSG_RESULT([overriding with ${TRILINOS_COMM}])

   AC_MSG_CHECKING(whether TRILINOS_ID is defined)
   if test ${TRILINOS_ID}; then
      AC_MSG_RESULT([yes (${TRILINOS_ID})])
   else
      AC_MSG_RESULT(no)
   fi

   AC_MSG_CHECKING(whether TRILINOS_TARGET is defined)
   if test ${TRILINOS_TARGET}; then
      AC_MSG_RESULT([yes (${TRILINOS_TARGET})])
   else
      TRILINOS_TARGET="${TRILINOS_ARCH}.${TRILINOS_COMM}${TRILINOS_ID}"
      AC_MSG_RESULT([no (using ${TRILINOS_TARGET})])
   fi

   TRILINOS_EPETRA_INC=${TRILINOS_HOME}/packages/epetra/src
   TRILINOS_AZTECOO_INC=${TRILINOS_HOME}/packages/aztecoo/src
   TRILINOS_LIB_DIR=${TRILINOS_HOME}/lib/${TRILINOS_TARGET}
fi


# Case II: Using Trilinos ARCH directory (ala Xyce)
# (a) Check is trilinos-arch-dir is valid
if test ${USE_TRILINOS} = yes && test X${USE_PKGDIR} = Xyes; then

   AC_MSG_CHECKING(whether PKGDIR is valid)
   if test -d ${PKGDIR}; then
      AC_MSG_RESULT(yes)
   else
      AC_MSG_NOTICE(no)
      AC_MSG_NOTICE([DISABLING TRILINOS!!])
      USE_TRILINOS=no
   fi

fi

# Case II: Using Trilinos ARCH directory (ala Xyce)
# (b) Set Trilinos variables
if test ${USE_TRILINOS} = yes && test X${USE_PKGDIR} = Xyes; then

   if test ${USE_MPI} = yes; then
      TRILINOS_COMM="MPI"
   else
      TRILINOS_COMM="SERIAL"
   fi

   TRILINOS_EPETRA_INC=${PKGDIR}/include/epetra/
   TRILINOS_AZTECOO_INC=${PKGDIR}/include/aztecoo/
   TRILINOS_LIB_DIR=${PKGDIR}/lib/

fi

# Set up flags for test compiles and links below
if test ${USE_TRILINOS} = yes; then
   LDFLAGS="-L${TRILINOS_LIB_DIR} $LDFLAGS"
fi

# Find a valid y12m
if test ${USE_TRILINOS} = yes; then
  AC_LANG_FORTRAN77
  AC_SEARCH_LIBS(y12mae,[y12m],,DONT_HAVE_LIBY12M="no")
  if test -e "$DONT_HAVE_LIBY12M"; then
    echo "------------------------------------------------------------------------"
    echo "Could not link to the Y12M library. This library is distributed with"
    echo "Trilinos as an external package. Specify LDFLAGS on the configure line"
    echo "to specify the name and/or path to the library."
    echo "------------------------------------------------------------------------"
    AC_MSG_NOTICE([DISABLING TRILINOS!!])
    USE_TRILINOS=no
  fi
fi

# Find g2c if it exists (don't seem to need this if it's not there)
if test ${USE_TRILINOS} = yes; then
  AC_LANG_C
  AC_SEARCH_LIBS(s_copy,[g2c])
fi

# Find valid BLAS
if test ${USE_TRILINOS} = yes; then
  AC_LANG_FORTRAN77
  AC_SEARCH_LIBS(dcopy,[cxml blas],,DONT_HAVE_BLAS="no")
  if test -e "$DONT_HAVE_LIBBLAS"; then
    echo "------------------------------------------------------------------------"
    echo "Could not link to the BLAS library. It can be obtained at"
    echo "http://www.netlib.org/lapack. Specify LDFLAGS on the configure line to"
    echo "specify the name and/or path to the library."
    echo "------------------------------------------------------------------------"
    AC_MSG_NOTICE([DISABLING TRILINOS!!])
    USE_TRILINOS=no
  fi
fi

# Find valid LAPACK
if test ${USE_TRILINOS} = yes; then
  AC_LANG_FORTRAN77
  AC_SEARCH_LIBS(dgetrs,[cxml lapack],,DONT_HAVE_LIBLAPACK="no")
  if test -e "$DONT_HAVE_LIBLAPACK"; then
    echo "------------------------------------------------------------------------"
    echo "Could not link to the LAPACK library. It can be obtained at"
    echo "http://www.netlib.org/lapack. Specify LDFLAGS on the configure line to"
    echo "specify the name and/or path to the library."
    echo "------------------------------------------------------------------------"
    AC_MSG_NOTICE([DISABLING TRILINOS!!])
    USE_TRILINOS=no
  fi

fi

# Check for ifpack
if test ${USE_TRILINOS} = yes; then
  AC_MSG_CHECKING(for ifpack)
  if test -f ${TRILINOS_LIB_DIR}/libifpack.a; then
    AC_MSG_RESULT(yes)
  else
    AC_MSG_RESULT(no)
    AC_MSG_NOTICE([DISABLING TRILINOS!!])
    USE_TRILINOS=no
  fi  
fi

# Check for aztecoo
if test ${USE_TRILINOS} = yes; then
  AC_MSG_CHECKING(for aztecoo)
  if test -f ${TRILINOS_LIB_DIR}/libaztecoo.a; then
    AC_MSG_RESULT(yes)
  else
    AC_MSG_RESULT(no)
    AC_MSG_NOTICE([DISABLING TRILINOS!!])
    USE_TRILINOS=no
  fi  
fi

# Check for epetra
if test ${USE_TRILINOS} = yes; then
  AC_MSG_CHECKING(for epetra)
  if test -f ${TRILINOS_LIB_DIR}/libepetra.a; then
    AC_MSG_RESULT(yes)
  else
    AC_MSG_RESULT(no)
    AC_MSG_NOTICE([DISABLING TRILINOS!!])
    USE_TRILINOS=no
  fi  
fi

# If we make it through all the tests, define the correct libs and includes!!
if test ${USE_TRILINOS} = yes; then
  LIBS="-lepetra -laztecoo -lifpack ${LIBS}"
  TRILINOS_CXXFLAGS="-DEPETRA_${TRILINOS_COMM} -I${TRILINOS_EPETRA_INC} -I${TRILINOS_AZTECOO_INC}"
fi

AC_SUBST(TRILINOS_CXXFLAGS)

# ------------------------------------------------------------------------
# Checks for TRILINOS status 
# ------------------------------------------------------------------------

if test X${USE_TRILINOS} = Xyes; then
   AM_CONDITIONAL(BUILD_NOXEPETRA, true)
else
   AM_CONDITIONAL(BUILD_NOXEPETRA, false)
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

AC_CXX_NAMESPACES
AC_CXX_MUTABLE
AC_CXX_HAVE_STL
AC_CXX_HAVE_STD

# ------------------------------------------------------------------------
# Checks for linker characteristics
# ------------------------------------------------------------------------

# -- DEC CXX --
if test X${target_cpu} = Xalpha && test X${CXX} = Xcxx; then
   XTRALDADD="cxx_repository/*.o"
fi

AC_MSG_CHECKING([for extra object files to be linked])
if test X$XTRALDADD = X; then
AC_MSG_RESULT(none)
else
AC_MSG_RESULT($XTRALDADD)
fi

AC_SUBST(XTRALDADD)

# ------------------------------------------------------------------------
# Checks for library functions
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for system services
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Alternate archiver
# ------------------------------------------------------------------------

AC_MSG_CHECKING(for special archiver command)
if test X${USE_ALTERNATE_AR} = Xyes; then
   AC_MSG_RESULT([yes (${ALTERNATE_AR})])
   AM_CONDITIONAL(USE_ALTERNATE_AR, true)
else
   AC_MSG_RESULT([none])
   AM_CONDITIONAL(USE_ALTERNATE_AR, false)
fi
AC_SUBST(ALTERNATE_AR)

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


